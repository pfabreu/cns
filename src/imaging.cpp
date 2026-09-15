#include "celestial/imaging.hpp"

#include <algorithm>
#include <cmath>
#include <random>

#include "celestial/star_catalog.hpp"

namespace celestial {
namespace {

/// A precomputed table of standard normal deviates. Generating 4.7M
/// std::normal_distribution samples per 2.35 Mpx frame dominates the whole
/// simulator; a table plus a per-frame offset is ~20x faster and is
/// statistically fine for a sensor noise model.
const std::vector<float>& normalTable() {
  static const std::vector<float> kTable = [] {
    std::vector<float> t(1 << 20);
    std::mt19937 rng(987654321u);
    std::normal_distribution<float> g(0.0f, 1.0f);
    for (float& v : t) v = g(rng);
    return t;
  }();
  return kTable;
}

constexpr double kArcsecPerRad = 180.0 * 3600.0 / M_PI;

/// Slerp-free small-angle attitude interpolation. Over one frame period the
/// change is a fraction of a degree, so linear interpolation of Euler angles
/// is well inside the error budget and much cheaper.
FrameTruth lerpPose(const FrameTruth& a, const FrameTruth& b, double s) {
  FrameTruth f = a;
  auto wrap = [](double d) { return std::atan2(std::sin(d), std::cos(d)); };
  f.roll = a.roll + s * wrap(b.roll - a.roll);
  f.pitch = a.pitch + s * wrap(b.pitch - a.pitch);
  f.yaw = a.yaw + s * wrap(b.yaw - a.yaw);
  f.pos.lat = a.pos.lat + s * (b.pos.lat - a.pos.lat);
  f.pos.lon = a.pos.lon + s * (b.pos.lon - a.pos.lon);
  f.pos.alt = a.pos.alt;
  return f;
}

/// Gauss-Newton fit of a circular Gaussian plus a flat background.
/// Params: A, x0, y0, sigma, B. Returns false if it fails to converge.
bool fitGaussian(const Image& img, double u_init, double v_init, int hw,
                 double sigma_init, double& u_out, double& v_out) {
  const int cx = static_cast<int>(std::round(u_init - 0.5));
  const int cy = static_cast<int>(std::round(v_init - 0.5));
  const int x0 = std::max(0, cx - hw), x1 = std::min(img.width - 1, cx + hw);
  const int y0 = std::max(0, cy - hw), y1 = std::min(img.height - 1, cy + hw);
  const int n = (x1 - x0 + 1) * (y1 - y0 + 1);
  if (n < 12) return false;

  // Seed background from the window border, amplitude from the peak.
  double border = 0.0; int nb = 0, peak = 0;
  for (int y = y0; y <= y1; ++y) {
    for (int x = x0; x <= x1; ++x) {
      const int val = img.at(x, y);
      if (val > peak) peak = val;
      if (x == x0 || x == x1 || y == y0 || y == y1) { border += val; ++nb; }
    }
  }
  if (nb == 0) return false;
  double p[5] = {std::max(1.0, peak - border / nb), u_init, v_init,
                 std::max(0.4, sigma_init), border / nb};

  for (int iter = 0; iter < 25; ++iter) {
    Eigen::Matrix<double, 5, 5> H = Eigen::Matrix<double, 5, 5>::Zero();
    Eigen::Matrix<double, 5, 1> g = Eigen::Matrix<double, 5, 1>::Zero();
    const double s2 = p[3] * p[3];
    for (int y = y0; y <= y1; ++y) {
      for (int x = x0; x <= x1; ++x) {
        const double dx = (x + 0.5) - p[1], dy = (y + 0.5) - p[2];
        const double r2 = dx * dx + dy * dy;
        const double E = std::exp(-r2 / (2.0 * s2));
        const double model = p[0] * E + p[4];
        const double res = img.at(x, y) - model;
        Eigen::Matrix<double, 5, 1> J;
        J(0) = E;
        J(1) = p[0] * E * dx / s2;
        J(2) = p[0] * E * dy / s2;
        J(3) = p[0] * E * r2 / (s2 * p[3]);
        J(4) = 1.0;
        H += J * J.transpose();
        g += J * res;
      }
    }
    H.diagonal() *= 1.0 + 1e-3;  // light damping
    const Eigen::Matrix<double, 5, 1> d = H.ldlt().solve(g);
    if (!d.allFinite()) return false;
    for (int k = 0; k < 5; ++k) p[k] += d(k);
    p[3] = std::clamp(p[3], 0.3, 20.0);
    if (std::abs(d(1)) < 1e-4 && std::abs(d(2)) < 1e-4) break;
  }
  if (std::abs(p[1] - u_init) > hw || std::abs(p[2] - v_init) > hw) return false;
  u_out = p[1];
  v_out = p[2];
  return true;
}

}  // namespace

SensorModel SensorModel::ideal() {
  SensorModel s;
  s.exposure_s = 0.0;  // no blur
  s.blur_substeps = 1;
  s.read_noise_e = 0.0;
  s.dark_e_per_s = 0.0;
  s.sky_mag_per_arcsec2 = 99.0;  // effectively no sky
  s.shot_noise = false;
  s.prnu_frac = 0.0;
  s.hot_pixel_rate = 0.0;
  return s;
}

Eigen::Vector3d pixelToRay(double u, double v, const Camera& cam) {
  const double f = cam.focalPx();
  return Eigen::Vector3d((u - cam.width / 2.0) / f, (v - cam.height / 2.0) / f,
                         1.0)
      .normalized();
}

bool rayToPixel(const Eigen::Vector3d& v_cam, const Camera& cam, double& u,
                double& v) {
  if (v_cam.z() <= 1e-9) return false;
  const double f = cam.focalPx();
  u = f * v_cam.x() / v_cam.z() + cam.width / 2.0;
  v = f * v_cam.y() / v_cam.z() + cam.height / 2.0;
  return true;
}

// ---------------------------------------------------------------------------
// Render
// ---------------------------------------------------------------------------

RenderedFrame renderFrame(const FrameTruth& a, const FrameTruth& b,
                          double frame_period_s, const StarField& field,
                          const std::vector<double>& vmags, const Camera& cam,
                          const Eigen::Matrix3d& C_b_c_true,
                          const SensorModel& sensor, unsigned frame_seed) {
  RenderedFrame out;
  out.image.width = cam.width;
  out.image.height = cam.height;
  const int npix = cam.width * cam.height;

  std::vector<double> accum(npix, 0.0);  // electrons

  const int steps = std::max(1, sensor.blur_substeps);
  const double expose_frac =
      (frame_period_s > 0.0)
          ? std::clamp(sensor.exposure_s / frame_period_s, 0.0, 1.0)
          : 0.0;
  const double smear_exposure =
      (sensor.exposure_s > 0.0) ? sensor.exposure_s : 0.0;
  // Photon collection time. Normally equal to the exposure; overridden only by
  // the blur-isolation diagnostic.
  const double exposure = (sensor.photometric_exposure_s > 0.0)
                              ? sensor.photometric_exposure_s
                              : (smear_exposure > 0.0 ? smear_exposure : 0.05);

  std::vector<Eigen::Vector3d> ecef;

  // Accumulators for the per-star ground truth.
  struct Acc {
    double wsum = 0, uw = 0, vw = 0;
    double u0 = 0, v0 = 0, u1 = 0, v1 = 0;
    bool started = false;
    bool truncated = false;
  };
  std::vector<Acc> acc(field.size());
  std::vector<double> star_e(field.size(), 0.0);

  for (int s = 0; s < steps; ++s) {
    const double frac = (steps == 1) ? 0.0
                                     : expose_frac * static_cast<double>(s) /
                                           static_cast<double>(steps - 1);
    const FrameTruth p = lerpPose(a, b, frac);
    field.ecefDirections(p.epoch, ecef);

    const Eigen::Matrix3d C_ecef_ned = nedBasisEcef(p.pos);
    const Eigen::Matrix3d cam_from_ned =
        (eulerToDcm(p.roll, p.pitch, p.yaw) * C_b_c_true).transpose();

    for (size_t i = 0; i < ecef.size(); ++i) {
      const Eigen::Vector3d v_cam =
          cam_from_ned * (C_ecef_ned.transpose() * ecef[i]);
      double u, v;
      if (!rayToPixel(v_cam, cam, u, v)) continue;

      // Total electrons from this star over the whole exposure, split evenly
      // across sub-steps.
      const double e_total =
          sensor.zeropoint_e_per_s * std::pow(10.0, -0.4 * vmags[i]) * exposure;
      const double e_step = e_total / steps;

      Acc& A = acc[i];
      const bool on_sensor =
          u >= -3 && u < cam.width + 3 && v >= -3 && v < cam.height + 3;
      if (on_sensor) {
        if (!A.started) { A.u0 = u; A.v0 = v; A.started = true; }
        A.u1 = u; A.v1 = v;
        A.wsum += e_step; A.uw += e_step * u; A.vw += e_step * v;
        star_e[i] += e_step;
      } else if (A.started) {
        A.truncated = true;
      }
      if (!on_sensor) continue;

      // Deposit a Gaussian PSF over a local window.
      const double sig = std::max(sensor.psf_sigma_px, 0.3);
      const int rad = static_cast<int>(std::ceil(4.0 * sig));
      const int u0i = std::max(0, static_cast<int>(std::floor(u)) - rad);
      const int u1i = std::min(cam.width - 1, static_cast<int>(std::ceil(u)) + rad);
      const int v0i = std::max(0, static_cast<int>(std::floor(v)) - rad);
      const int v1i = std::min(cam.height - 1, static_cast<int>(std::ceil(v)) + rad);
      const double inv2s2 = 1.0 / (2.0 * sig * sig);
      const double norm = e_step / (2.0 * M_PI * sig * sig);
      for (int y = v0i; y <= v1i; ++y) {
        const double dy = y + 0.5 - v;
        for (int x = u0i; x <= u1i; ++x) {
          const double dx = x + 0.5 - u;
          accum[y * cam.width + x] += norm * std::exp(-(dx * dx + dy * dy) * inv2s2);
        }
      }
    }
  }

  // Labels
  for (size_t i = 0; i < field.size(); ++i) {
    if (!acc[i].started || acc[i].wsum <= 0.0) continue;
    StarLabel L;
    L.id = field.id(i);
    L.vmag = vmags[i];
    L.u = acc[i].uw / acc[i].wsum;
    L.v = acc[i].vw / acc[i].wsum;
    L.u0 = acc[i].u0; L.v0 = acc[i].v0;
    L.u1 = acc[i].u1; L.v1 = acc[i].v1;
    L.total_e = star_e[i];
    L.truncated = acc[i].truncated;
    out.labels.push_back(L);
  }

  // --- background, fixed pattern, noise, digitisation ----------------------
  const std::vector<float>& ntab = normalTable();
  const size_t nmask = ntab.size() - 1;
  size_t ni = (static_cast<size_t>(sensor.seed) * 2654435761u +
               static_cast<size_t>(frame_seed) * 40503u) & nmask;
  auto gauss = [&]() { ni = (ni + 1) & nmask; return double(ntab[ni]); };

  // Sky: convert surface brightness to electrons per pixel.
  const double px_arcsec = cam.pixelIfov() * kArcsecPerRad;
  const double sky_mag_px =
      sensor.sky_mag_per_arcsec2 - 2.5 * std::log10(px_arcsec * px_arcsec);
  const double sky_e =
      sensor.zeropoint_e_per_s * std::pow(10.0, -0.4 * sky_mag_px) * exposure;
  const double dark_e = sensor.dark_e_per_s * exposure;

  // Hot pixels are a FIXED pattern: same pixels every frame. They pass a
  // k-sigma threshold reliably and track perfectly, which is exactly the kind
  // of false detection RANSAC exists to reject.
  std::mt19937 fixed_rng(sensor.seed);
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  std::vector<float> prnu(npix, 1.0f);
  std::vector<float> hot(npix, 0.0f);
  if (sensor.prnu_frac > 0.0 || sensor.hot_pixel_rate > 0.0) {
    std::normal_distribution<double> pg(1.0, sensor.prnu_frac);
    for (int i = 0; i < npix; ++i) {
      if (sensor.prnu_frac > 0.0) prnu[i] = static_cast<float>(pg(fixed_rng));
      if (uni(fixed_rng) < sensor.hot_pixel_rate) {
        hot[i] = static_cast<float>(sensor.hot_pixel_e_per_s * exposure);
      }
    }
  }

  // --- weather ------------------------------------------------------------
  // Cloud: a few octaves of value noise, smoothly varying, ATTENUATING the
  // stars and LIFTING the background. Flare: broad radial gradients anchored
  // off-axis and fixed in the camera frame.
  std::vector<float> transmit, glow;
  if (sensor.cloud_amount > 0 || sensor.flare_amount > 0) {
    transmit.assign(npix, 1.f);
    glow.assign(npix, 0.f);
    std::mt19937 wr(sensor.seed ^ 0x51ed270bu);
    std::uniform_real_distribution<double> ph(0.0, 2 * M_PI);
    // Value noise as a small sum of sinusoids: crude, cheap, and smooth --
    // which is all that matters for standing in as cloud structure.
    struct Wave { double kx, ky, px, py, amp; };
    std::vector<Wave> waves;
    for (int o = 0; o < 4; ++o) {
      const double k = 2 * M_PI / (sensor.cloud_scale_px / (1 << o));
      std::uniform_real_distribution<double> dir(0.0, M_PI);
      const double a = dir(wr);
      waves.push_back({k * std::cos(a), k * std::sin(a), ph(wr), ph(wr),
                       1.0 / (1 << o)});
    }
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    const double fx = u01(wr) * out.image.width * 1.6 - 0.3 * out.image.width;
    const double fy = u01(wr) * out.image.height * 1.6 - 0.3 * out.image.height;
    for (int y = 0; y < out.image.height; ++y) {
      for (int x = 0; x < out.image.width; ++x) {
        const size_t i = size_t(y) * out.image.width + x;
        if (sensor.cloud_amount > 0) {
          double v = 0, norm = 0;
          for (const auto& w : waves) {
            v += w.amp * std::sin(w.kx * x + w.px) * std::sin(w.ky * y + w.py);
            norm += w.amp;
          }
          const double c = std::clamp(0.5 + 0.5 * v / norm, 0.0, 1.0);
          const double opacity = sensor.cloud_amount * c;
          if (sensor.cloud_attenuates) transmit[i] *= float(1.0 - opacity);
          glow[i] += float(opacity * sensor.cloud_glow_e_per_s * exposure);
        }
        if (sensor.flare_amount > 0) {
          const double r = std::hypot(x - fx, y - fy) /
                           double(out.image.width);
          glow[i] += float(sensor.flare_amount * 900.0 * exposure *
                           std::exp(-r * r / 0.18));
        }
      }
    }
  }

  const double max_adu = std::pow(2.0, sensor.bit_depth) - 1.0;
  out.image.data.resize(npix);
  for (int i = 0; i < npix; ++i) {
    double star_e = accum[i];
    double extra = 0.0;
    if (!transmit.empty()) { star_e *= transmit[i]; extra = glow[i]; }
    double e = (star_e + sky_e + extra) * prnu[i] + dark_e + hot[i];
    if (sensor.shot_noise && e > 0.0) {
      if (e > 50.0) {
        e += std::sqrt(e) * gauss();  // normal approximation
      } else {
        // Small lambda: table-free Knuth Poisson is fine here, it is rare.
        double L = std::exp(-e), p = 1.0; int k = 0;
        do { ++k; p *= (double(ni = (ni * 1103515245u + 12345u) & 0x7fffffff) /
                        double(0x7fffffff)); } while (p > L && k < 200);
        e = k - 1;
      }
    }
    if (sensor.read_noise_e > 0.0) e += sensor.read_noise_e * gauss();
    // Saturate on the high side only. Do NOT clamp at zero here: a real
    // sensor's BIAS OFFSET exists precisely so that negative read-noise
    // excursions remain representable. Clipping first collapses the noise
    // distribution to a half-Gaussian with a spike at zero, which drives the
    // MAD to zero and destroys any k-sigma threshold downstream.
    e = std::min(e, sensor.full_well_e);
    const double adu = sensor.bias_adu + e / sensor.gain_e_per_adu;
    out.image.data[i] =
        static_cast<uint16_t>(std::clamp(adu, 0.0, max_adu));
  }
  return out;
}

// ---------------------------------------------------------------------------
// Detect
// ---------------------------------------------------------------------------

namespace {

/// POLAR MATCHED FILTER: TRIED, REMOVED. Kept as a note so it is not retried.
///
/// To first order the smear field is a rotation about the principal point plus
/// a uniform translation, and a translation can be absorbed into the centre:
/// the whole field is a pure rotation about the instantaneous centre
/// `c = (w_x, w_y) * f / w_z`. About THAT point every streak subtends the same
/// angle, so in polar coordinates the kernel is a fixed width and a running sum
/// makes it O(1) regardless of streak length. Elegant, and it does not work
/// here.
///
/// The centre is not near the image. In a coordinated turn the body rates are
/// roughly `q = W sin(bank)`, `r = W cos(bank)`, so with 23 deg of bank the
/// translation terms are ~40% of the rotation term and the centre lands ~810 px
/// OUTSIDE the frame. The polar grid must then span r up to ~1950 px, and to
/// resolve one pixel of arc at that radius needs ~12000 angular bins. The
/// warped image comes out LARGER than the original, and peak detection over it
/// costs more than the thing it replaced: measured, 1585 ms against 247 ms for
/// the naive filter.
///
/// Polar is right when the rotation centre is inside or near the frame. For a
/// banked orbit it is not. Decimation is what worked -- see below.

/// Detection on the matched-filter surface.
///
/// Connected components are WRONG here and it is worth saying why: the filter
/// integrates along the streak, so its output is spatially smoothed and its
/// noise is correlated. A median+MAD threshold on a smoothed map floods --
/// measured, 800-2600 false components per frame against 0 for the plain
/// detector. The signal in a matched-filter output is a PEAK, so detect peaks:
/// a local maximum above threshold, one detection each.
///
/// Centroiding is done on the ORIGINAL image over a window covering the local
/// streak. The filtered map is a detection surface, not a photometric one.
/// SExtractor-style mesh background, subtracted with a positive offset so the
/// result stays representable in the unsigned image type.
Image meshSubtract(const Image& img, int cell) {
  const int W = img.width, H = img.height;
  const int nx = (W + cell - 1) / cell, ny = (H + cell - 1) / cell;
  std::vector<float> node(size_t(nx) * ny, 0.f);
  std::vector<uint16_t> buf;

  for (int cy = 0; cy < ny; ++cy) {
    for (int cx = 0; cx < nx; ++cx) {
      buf.clear();
      for (int y = cy * cell; y < std::min(H, (cy + 1) * cell); y += 2)
        for (int x = cx * cell; x < std::min(W, (cx + 1) * cell); x += 2)
          buf.push_back(img.at(x, y));
      if (buf.empty()) continue;
      // Sigma-clip so that sources in the cell do not drag the estimate up.
      double med = 0;
      for (int it = 0; it < 3; ++it) {
        std::nth_element(buf.begin(), buf.begin() + buf.size() / 2, buf.end());
        med = buf[buf.size() / 2];
        std::vector<double> ad;
        ad.reserve(buf.size());
        for (auto v : buf) ad.push_back(std::abs(double(v) - med));
        std::nth_element(ad.begin(), ad.begin() + ad.size() / 2, ad.end());
        const double sig = 1.4826 * ad[ad.size() / 2] + 1e-6;
        std::vector<uint16_t> keep;
        keep.reserve(buf.size());
        for (auto v : buf)
          if (std::abs(double(v) - med) < 3 * sig) keep.push_back(v);
        if (keep.size() < 8) break;
        buf.swap(keep);
      }
      node[size_t(cy) * nx + cx] = float(med);
    }
  }

  // Median filter across nodes: one cell containing something bright must not
  // pull down its whole neighbourhood when interpolated.
  std::vector<float> sm = node;
  for (int cy = 0; cy < ny; ++cy) {
    for (int cx = 0; cx < nx; ++cx) {
      std::vector<float> w;
      for (int j = -1; j <= 1; ++j)
        for (int i = -1; i <= 1; ++i) {
          const int y = cy + j, x = cx + i;
          if (x < 0 || y < 0 || x >= nx || y >= ny) continue;
          w.push_back(node[size_t(y) * nx + x]);
        }
      std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
      sm[size_t(cy) * nx + cx] = w[w.size() / 2];
    }
  }

  Image out;
  out.width = W;
  out.height = H;
  out.data.assign(size_t(W) * H, 0);
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const double gx =
          std::min(double(nx - 1), std::max(0.0, (x + 0.5) / cell - 0.5));
      const double gy =
          std::min(double(ny - 1), std::max(0.0, (y + 0.5) / cell - 0.5));
      const int x0 = int(gx), y0 = int(gy);
      const int x1 = std::min(nx - 1, x0 + 1), y1 = std::min(ny - 1, y0 + 1);
      const double fx = gx - x0, fy = gy - y0;
      const double b = (1 - fx) * (1 - fy) * sm[size_t(y0) * nx + x0] +
                       fx * (1 - fy) * sm[size_t(y0) * nx + x1] +
                       (1 - fx) * fy * sm[size_t(y1) * nx + x0] +
                       fx * fy * sm[size_t(y1) * nx + x1];
      // Offset keeps read-noise excursions below background representable,
      // for the same reason renderFrame carries a bias offset.
      out.data[size_t(y) * W + x] =
          uint16_t(std::clamp(double(img.at(x, y)) - b + 1000.0, 0.0, 65535.0));
    }
  }
  return out;
}


}  // namespace

std::vector<Detection> detectStars(const Image& raw,
                                   const DetectorConfig& cfg) {
  // Mesh background first, if enabled: everything downstream -- threshold,
  // centroid -- then sees a flat field.
  const Image bg_removed =
      cfg.bg_mesh_px > 0 ? meshSubtract(raw, cfg.bg_mesh_px) : Image{};
  const Image& img = cfg.bg_mesh_px > 0 ? bg_removed : raw;
  std::vector<Detection> dets;
  const int n = img.width * img.height;
  if (n == 0) return dets;

  double bg = 0.0, sigma = 1.0;
  if (cfg.robust_background) {
    // Median + MAD. Robust to the stars themselves and to a vignetted corner,
    // unlike the mean/stddev the paper uses.
    //
    // Estimated on a subsample: the background is smooth and stars occupy a
    // tiny fraction of the frame, so every 17th pixel gives the same answer
    // ~17x faster. Full-frame nth_element was the simulator's bottleneck.
    const int stride = 17;
    std::vector<uint16_t> s;
    s.reserve(n / stride + 1);
    for (int i = 0; i < n; i += stride)
      s.push_back(img.data[i]);
    std::nth_element(s.begin(), s.begin() + s.size() / 2, s.end());
    bg = s[s.size() / 2];
    for (size_t i = 0; i < s.size(); ++i)
      s[i] = static_cast<uint16_t>(std::abs(int(s[i]) - int(bg)));
    std::nth_element(s.begin(), s.begin() + s.size() / 2, s.end());
    sigma = std::max(1.0, 1.4826 * s[s.size() / 2]);
  } else {
    double sum = 0, sum2 = 0;
    for (int i = 0; i < n; ++i) { sum += img.data[i]; sum2 += double(img.data[i]) * img.data[i]; }
    bg = sum / n;
    sigma = std::max(1.0, std::sqrt(std::max(0.0, sum2 / n - bg * bg)));
  }
  const double thresh = bg + cfg.threshold_k * sigma;

  // Spatially varying floor from a detection prior, if one was supplied. Never
  // below threshold_k_low, so there is always a real significance requirement
  // and a confident prediction over empty sky still yields nothing.
  const bool use_prior =
      cfg.prior.size() == size_t(img.width) * size_t(img.height) &&
      cfg.threshold_k_low > 0 && cfg.threshold_k_low < cfg.threshold_k;
  auto floorAt = [&](int idx) -> double {
    if (!use_prior) return thresh;
    const double p = std::clamp(double(cfg.prior[size_t(idx)]), 0.0, 1.0);
    return bg + (cfg.threshold_k -
                 (cfg.threshold_k - cfg.threshold_k_low) * p) * sigma;
  };

  // 8-connected component labelling, iterative flood fill.
  std::vector<uint8_t> seen(n, 0);
  std::vector<int> stack;
  std::vector<int> comp;
  for (int start = 0; start < n; ++start) {
    if (seen[start] || img.data[start] <= floorAt(start)) continue;
    stack.clear(); comp.clear();
    stack.push_back(start); seen[start] = 1;
    while (!stack.empty()) {
      const int p = stack.back(); stack.pop_back();
      comp.push_back(p);
      const int x = p % img.width, y = p / img.width;
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          if (!dx && !dy) continue;
          const int nx = x + dx, ny = y + dy;
          if (nx < 0 || ny < 0 || nx >= img.width || ny >= img.height) continue;
          const int q = ny * img.width + nx;
          if (seen[q] || img.data[q] <= floorAt(q)) continue;
          seen[q] = 1;
          stack.push_back(q);
        }
      }
    }
    if (static_cast<int>(comp.size()) < cfg.min_pixels ||
        static_cast<int>(comp.size()) > cfg.max_pixels) continue;

    // Background-subtracted weighted centroid over the WHOLE component.
    // (The paper uses a fixed 3x3 window, which does not contain a streak.)
    double w = 0, uw = 0, vw = 0;
    for (int p : comp) {
      const double val = img.data[p] - bg;
      const double x = p % img.width + 0.5, y = p / img.width + 0.5;
      w += val; uw += val * x; vw += val * y;
    }
    if (w <= 0) continue;
    Detection d;
    d.u = uw / w; d.v = vw / w; d.flux = w;
    d.n_pixels = static_cast<int>(comp.size());

    if (cfg.centroid == CentroidMethod::GaussianFit) {
      double gu, gv;
      // Seed from the CoG; fall back to it if the fit fails.
      if (fitGaussian(img, d.u, d.v, cfg.fit_halfwidth,
                      std::max(1.0, std::sqrt(double(comp.size()) / M_PI) / 2.0),
                      gu, gv)) {
        d.u = gu; d.v = gv;
      }
    }

    // Second moments -> elongation, a free streak/point discriminator.
    double sxx = 0, syy = 0, sxy = 0;
    for (int p : comp) {
      const double val = img.data[p] - bg;
      const double dx = (p % img.width + 0.5) - d.u;
      const double dy = (p / img.width + 0.5) - d.v;
      sxx += val * dx * dx; syy += val * dy * dy; sxy += val * dx * dy;
    }
    sxx /= w; syy /= w; sxy /= w;
    const double tr = sxx + syy;
    const double det = sxx * syy - sxy * sxy;
    const double disc = std::sqrt(std::max(0.0, tr * tr / 4.0 - det));
    const double l1 = tr / 2.0 + disc, l2 = tr / 2.0 - disc;
    d.elongation = (l2 > 1e-6) ? std::sqrt(l1 / l2) : 99.0;
    dets.push_back(d);
  }

  return dets;
}

// ---------------------------------------------------------------------------
// Match
// ---------------------------------------------------------------------------

int matchDetections(const std::vector<Detection>& dets, const Epoch& epoch,
                    const Eigen::Matrix3d& C_l_b_est,
                    const Eigen::Matrix3d& C_b_c_assumed, const Camera& cam,
                    const StarField& field, const Geodetic& assumed_pos,
                    const MatcherConfig& mcfg, FrameData& out) {
  out.id.clear();
  out.v_cam.clear();
  out.epoch = epoch;
  out.C_l_b_est = C_l_b_est;
  if (dets.empty()) return 0;

  // Project the catalogue with the ESTIMATED attitude and mounting, at the
  // assumed position. No ground truth is used.
  std::vector<Eigen::Vector3d> ecef;
  field.ecefDirections(epoch, ecef);
  const Eigen::Matrix3d C_ecef_ned = nedBasisEcef(assumed_pos);
  const Eigen::Matrix3d cam_from_ned = (C_l_b_est * C_b_c_assumed).transpose();

  std::vector<double> pu, pv;
  std::vector<int> pid;
  for (size_t i = 0; i < ecef.size(); ++i) {
    const Eigen::Vector3d v_cam =
        cam_from_ned * (C_ecef_ned.transpose() * ecef[i]);
    double u, v;
    if (!rayToPixel(v_cam, cam, u, v)) continue;
    const double m = mcfg.radius_px;
    if (u < -m || v < -m || u > cam.width + m || v > cam.height + m) continue;
    pu.push_back(u); pv.push_back(v); pid.push_back(field.id(i));
  }

  // MUTUAL nearest neighbour. The previous rule -- nearest prediction, then
  // reject if a second lies within `ambiguity_ratio` times as far -- collapses
  // when the catalogue is deep and the radius is wide: with ~130 predictions
  // in frame at a 120 px radius there is almost always a second candidate, so
  // most genuine matches are discarded (observed: 12 detected, 6 matched).
  //
  // Requiring the pair to be each other's nearest is scale-free: it stays
  // correct as the catalogue deepens or the search widens, and it enforces
  // one-to-one for free.
  std::vector<int> best_for_pred(pu.size(), -1);
  std::vector<double> best_d2(pu.size(), 1e18);
  for (size_t di = 0; di < dets.size(); ++di) {
    for (size_t k = 0; k < pu.size(); ++k) {
      const double dx = pu[k] - dets[di].u, dy = pv[k] - dets[di].v;
      const double r2 = dx * dx + dy * dy;
      if (r2 < best_d2[k]) { best_d2[k] = r2; best_for_pred[k] = int(di); }
    }
  }

  int matched = 0;
  for (const Detection& d : dets) {
    int best = -1, second = -1;
    double bd = 1e18, sd = 1e18;
    for (size_t k = 0; k < pu.size(); ++k) {
      const double dx = pu[k] - d.u, dy = pv[k] - d.v;
      const double r2 = dx * dx + dy * dy;
      if (r2 < bd) { sd = bd; second = best; bd = r2; best = static_cast<int>(k); }
      else if (r2 < sd) { sd = r2; second = static_cast<int>(k); }
    }
    if (best < 0 || bd > mcfg.radius_px * mcfg.radius_px) continue;
    // Mutual: the prediction's nearest detection must be this one.
    const size_t di = size_t(&d - dets.data());
    if (best_for_pred[best] != int(di)) continue;
    (void)second; (void)sd;
    out.id.push_back(pid[best]);
    out.v_cam.push_back(pixelToRay(d.u, d.v, cam));
    ++matched;
  }
  return matched;
}

int buildFrameData(const Image& img, const Epoch& epoch,
                   const Eigen::Matrix3d& C_l_b_est,
                   const Eigen::Matrix3d& C_b_c_assumed, const Camera& cam,
                   const StarField& field, const Geodetic& assumed_pos,
                   const DetectorConfig& dcfg, const MatcherConfig& mcfg,
                   FrameData& out) {
  return matchDetections(detectStars(img, dcfg), epoch, C_l_b_est,
                         C_b_c_assumed, cam, field, assumed_pos, mcfg, out);
}

}  // namespace celestial
