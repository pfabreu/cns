#include "celestial/horizon.hpp"

#include <Eigen/Geometry>
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>

namespace celestial {
namespace {

/// Moving average of the tilt measurement, held per HorizonState.
void averageTilt(TiltMeasurement& z, HorizonState& st, double window_s,
                 double dt) {
  if (window_s <= 0 || dt <= 0) return;
  const int n = std::max(1, int(window_s / dt + 0.5));
  st.tilt_hist.push_back(z.tilt);
  while (int(st.tilt_hist.size()) > n) st.tilt_hist.pop_front();
  Eigen::Vector2d m = Eigen::Vector2d::Zero();
  for (const auto& v : st.tilt_hist) m += v;
  const double k = double(st.tilt_hist.size());
  z.tilt = m / k;
  z.R /= k;   // the mean of k samples is that much more certain
}

}  // namespace
namespace {

constexpr double kEarthR = 6371008.8;

std::mt19937& rng(unsigned seed) {
  static std::mt19937 g(seed);
  return g;
}

}  // namespace

double horizonDip(double altitude_m, double refraction_k) {
  if (altitude_m <= 0) return 0.0;
  // Geometric dip sqrt(2h/R), reduced by refraction bending the ray downward.
  return std::sqrt(2.0 * altitude_m / kEarthR) * std::sqrt(1.0 - 2.0 * refraction_k);
}

TiltMeasurement measureHorizon(const Eigen::Matrix3d& C_true,
                               const Eigen::Matrix3d& C_est, double altitude_m,
                               double t, const HorizonConfig& cfg,
                               HorizonState& state) {
  const double dt_est =
      (state.last_meas_t >= 0 && t > state.last_meas_t) ? t - state.last_meas_t
                                                        : 0.1;
  state.last_meas_t = t;
  TiltMeasurement out;
  if (cfg.cameras.empty() || altitude_m <= 0) return out;

  auto& g = rng(cfg.seed);
  std::normal_distribution<double> nd(0.0, 1.0);

  // --- propagate the shared atmosphere state -----------------------------
  if (state.last_t < 0) {
    state.anom_dip = cfg.anom_dip_sigma * nd(g);
    state.anom_grad = cfg.anom_grad_sigma * nd(g);
  } else {
    const double dt = std::max(0.0, t - state.last_t);
    const double a = std::exp(-dt / std::max(1e-6, cfg.anom_dip_tau));
    const double q = std::sqrt(std::max(0.0, 1.0 - a * a));
    state.anom_dip = a * state.anom_dip + q * cfg.anom_dip_sigma * nd(g);
    state.anom_grad = a * state.anom_grad + q * cfg.anom_grad_sigma * nd(g);
  }
  state.last_t = t;

  // --- the quantity being measured ---------------------------------------
  // delta rotates the estimated body frame onto the true one, so its x and y
  // components ARE the IMU tilt error about the forward and lateral axes.
  const Eigen::AngleAxisd aa(C_est.transpose() * C_true);
  const Eigen::Vector3d delta = aa.axis() * aa.angle();
  const double tau_fwd = delta.x();
  const double tau_lat = delta.y();

  // --- per-camera measurement --------------------------------------------
  double sum_lat = 0, sum_fwd = 0;
  double var_lat = 0, var_fwd = 0;
  int n = 0;

  for (const HorizonCamera& c : cfg.cameras) {
    const double s = c.aft ? -1.0 : 1.0;  // aft view flips both axes

    // Per-camera anomalous component: the fore and aft tangent points are
    // ~200 km apart, so the anomaly is only partly shared. This is what stops
    // the fore/aft pair from cancelling perfectly, and it is the term that
    // sets the floor on that configuration.
    const double anom_own =
        cfg.anom_decorrelated_frac * cfg.anom_dip_sigma * nd(g);

    // Line-fit statistics over `width` columns.
    // intercept -> dip;  slope -> roll.
    const double sig_dip_edge = c.edge_sigma_px / std::sqrt(double(c.width)) * c.ifov();
    const double sig_slope_edge =
        c.edge_sigma_px * std::sqrt(12.0 / c.width) / c.width;
    // Across-field gradient of anomalous dip appears as a line rotation.
    const double sig_slope_anom = cfg.anom_grad_sigma / c.hfov;

    // DIP AXIS: absorbs the FULL common-mode anomalous dip.
    const double meas_lat = s * (tau_lat * s + state.anom_dip + anom_own) +
                            c.mount_lat + sig_dip_edge * nd(g);
    // SLOPE AXIS: immune to the common-mode term; sees only its gradient.
    const double meas_fwd = s * (tau_fwd * s) + c.mount_fwd +
                            (state.anom_grad / c.hfov) + sig_slope_edge * nd(g);

    sum_lat += meas_lat;
    sum_fwd += meas_fwd;
    var_lat += sig_dip_edge * sig_dip_edge;
    var_fwd += sig_slope_edge * sig_slope_edge + sig_slope_anom * sig_slope_anom;
    ++n;
  }

  out.tilt(0) = sum_fwd / n;
  out.tilt(1) = sum_lat / n;

  out.R.setZero();
  out.R(0, 0) = var_fwd / (n * n);
  out.R(1, 1) = var_lat / (n * n);

  // The common-mode anomalous dip cancels with a fore/aft pair but survives
  // with a single camera, so it enters R only in the single-camera case.
  // With a pair, only the DIFFERENTIAL anomaly remains: modelled as a fixed
  // fraction of the common-mode magnitude, since the two tangent points are
  // ~200 km apart and only weakly correlated at that separation.
  bool has_fore = false, has_aft = false;
  for (const HorizonCamera& c : cfg.cameras) (c.aft ? has_aft : has_fore) = true;
  if (has_fore && has_aft) {
    const double resid = 0.33 * cfg.anom_dip_sigma;
    out.R(1, 1) += resid * resid;
  } else {
    out.R(1, 1) += cfg.anom_dip_sigma * cfg.anom_dip_sigma;
  }

  out.ok = true;

  // Short moving average: a Lepton-class sensor is noise-limited, so this is
  // worth ~37% on the fix. See HorizonConfig::tilt_average_s -- and note the
  // window must stay SHORT, because a lag during an orbit is a
  // heading-correlated error the orbit average cannot remove.
  averageTilt(out, state, cfg.tilt_average_s, dt_est);
  return out;
}

Eigen::Matrix3d fuseTilt(const Eigen::Matrix3d& C_est, const TiltMeasurement& z,
                         double imu_tilt_sigma, Eigen::Matrix2d* posterior) {
  if (!z.ok) {
    if (posterior) *posterior = Eigen::Matrix2d::Identity() *
                                (imu_tilt_sigma * imu_tilt_sigma);
    return C_est;
  }
  // Prior: the IMU tilt error is zero-mean with this sigma. The horizon
  // measures it. Standard linear-Gaussian update on a 2-vector.
  Eigen::Matrix2d P = Eigen::Matrix2d::Identity() * (imu_tilt_sigma * imu_tilt_sigma);
  const Eigen::Matrix2d S = P + z.R;
  const Eigen::Matrix2d K = P * S.inverse();
  const Eigen::Vector2d dtheta = K * z.tilt;
  if (posterior) *posterior = (Eigen::Matrix2d::Identity() - K) * P;

  // delta was defined as C_est^T * C_true, so the correction multiplies on
  // the right.
  //
  // ZEROING THE BODY-Z COMPONENT IS NOT ENOUGH, and the original comment here
  // claimed "yaw is untouched", which is false whenever the aircraft is banked.
  // Yaw is rotation about the LOCAL VERTICAL, not about body z, and the two
  // differ by the bank angle -- so a correction about a horizontal BODY axis
  // carries a component along the nav vertical and moves the heading.
  //
  // Measured with a PERFECT sensor and NO atmosphere, one orbit at 14.3 deg
  // mean bank: -0.176 deg of yaw change per frame, consistently signed. That
  // drove the celestial compass residual from 0.089 to 0.310 deg and made the
  // horizon a net LOSS on the transit despite improving every fix.
  //
  // This is the mirror of the trap `mountTiltOnly` exists to avoid, which
  // documents it from the other side: "Level it is body z; banked it tilts by
  // the bank angle." Project the correction orthogonal to the local vertical
  // expressed in body axes, and the compass is left alone.
  Eigen::Vector3d rv(dtheta(0), dtheta(1), 0.0);
  const Eigen::Vector3d vert_body = C_est.transpose() * Eigen::Vector3d::UnitZ();
  if (vert_body.norm() > 1e-9) {
    const Eigen::Vector3d v = vert_body.normalized();
    rv -= rv.dot(v) * v;
  }
  const double ang = rv.norm();
  const Eigen::Matrix3d dC =
      (ang < 1e-12) ? Eigen::Matrix3d::Identity()
                    : Eigen::AngleAxisd(ang, rv / ang).toRotationMatrix();
  return C_est * dC;
}

Image renderHorizonView(const Eigen::Matrix3d& C_true, double altitude_m,
                        const HorizonCamera& cam, const HorizonState& state,
                        unsigned seed) {
  Image img;
  img.width = cam.width;
  img.height = cam.width * 4 / 5;  // 640x512, Boson aspect
  img.data.assign(size_t(img.width) * img.height, 0);

  const double dip = horizonDip(altitude_m) + state.anom_dip;
  // TRUE focal length; the extractor will assume the nominal one.
  const double f = (cam.width / 2.0) / std::tan(cam.hfov / 2.0) *
                   (1.0 + cam.focal_err);

  // Camera axes in body: boresight along +x (fore) or -x (aft), image x to
  // starboard, image y down.
  Eigen::Matrix3d C_b_h;
  const double s = cam.aft ? -1.0 : 1.0;
  C_b_h.col(0) = Eigen::Vector3d(0, s, 0);
  C_b_h.col(1) = Eigen::Vector3d(0, 0, 1);
  C_b_h.col(2) = Eigen::Vector3d(s, 0, 0);
  const Eigen::Matrix3d C_l_h = C_true * C_b_h;

  std::mt19937 g(seed ? seed : 1u);
  std::normal_distribution<double> nd(0.0, 1.0);

  // Column FPN: fixed per-column offset, the dominant microbolometer mode.
  std::vector<double> fpn(img.width, 0.0);
  if (cam.fpn_counts > 0) {
    std::mt19937 fg(seed ? seed * 7919u : 7919u);
    std::normal_distribution<double> fn(0.0, cam.fpn_counts);
    for (int x = 0; x < img.width; ++x) fpn[x] = fn(fg);
  }

  for (int y = 0; y < img.height; ++y) {
    for (int x = 0; x < img.width; ++x) {
      double nx = (x + 0.5 - img.width / 2.0) / f;
      double ny = (y + 0.5 - img.height / 2.0) / f;
      // Radial distortion of the true optics (first-order Brown, inverted
      // approximately -- adequate for the magnitudes involved).
      if (cam.k1 != 0.0) {
        const double r2 = nx * nx + ny * ny;
        const double sc = 1.0 - cam.k1 * r2;
        nx *= sc; ny *= sc;
      }
      const Eigen::Vector3d v_cam(nx, ny, 1.0);
      const Eigen::Vector3d v_ned = C_l_h * v_cam.normalized();
      const double el = nedToElevation(v_ned);
      // Sea is warm and textured; sky is cold and smooth.
      // Soft limb: the sea/sky transition is spread by sea state and blur,
      // not a step. blur_ang converts the pixel blur into an angle.
      const double blur_ang = std::max(1e-9, 1.5 / f);
      const double mix = 0.5 * (1.0 - std::tanh((el + dip) / blur_ang));
      const double base = 1200.0 + 1400.0 * mix;
      const double tex = 25.0 + 65.0 * mix;
      const double val = base + tex * nd(g) + fpn[x] -
                         300.0 * std::min(0.0, el * 8.0);
      img.data[size_t(y) * img.width + x] =
          uint16_t(std::clamp(val, 0.0, 4095.0));
    }
  }
  return img;
}


HorizonLine detectHorizon(const Image& img, const HorizonCamera& /*cam*/) {
  HorizonLine out;
  if (img.width < 16 || img.height < 8) return out;

  // 1. Per-column vertical gradient peak, refined to subpixel by a parabola
  //    through the peak and its neighbours.
  std::vector<double> xs, ys;
  xs.reserve(img.width);
  ys.reserve(img.width);
  for (int x = 0; x < img.width; ++x) {
    int best_y = -1;
    double best_g = 0.0;
    for (int y = 1; y < img.height - 1; ++y) {
      const double gr = double(img.data[size_t(y + 1) * img.width + x]) -
                        double(img.data[size_t(y - 1) * img.width + x]);
      if (std::abs(gr) > std::abs(best_g)) { best_g = gr; best_y = y; }
    }
    if (best_y < 2 || best_y > img.height - 3) continue;
    if (std::abs(best_g) < 50.0) continue;  // no usable edge in this column

    auto G = [&](int y) {
      return std::abs(double(img.data[size_t(y + 1) * img.width + x]) -
                      double(img.data[size_t(y - 1) * img.width + x]));
    };
    const double a = G(best_y - 1), b = G(best_y), c = G(best_y + 1);
    const double den = a - 2 * b + c;
    const double sub = (std::abs(den) > 1e-9) ? 0.5 * (a - c) / den : 0.0;
    xs.push_back(x + 0.5 - img.width / 2.0);
    ys.push_back(best_y + 0.5 + std::clamp(sub, -1.0, 1.0));
  }
  if (xs.size() < 16) return out;

  // 2. THEIL-SEN: median of pairwise slopes. Robust to cloud tops, ships and
  //    whitecaps, which produce outlier columns that would drag a least
  //    squares fit. Subsampled stride keeps this cheap.
  std::vector<double> slopes;
  const size_t n = xs.size();
  const size_t stride = std::max<size_t>(1, n / 300);
  slopes.reserve((n / stride) * (n / stride) / 2 + 8);
  for (size_t i = 0; i < n; i += stride) {
    for (size_t j = i + stride; j < n; j += stride) {
      const double dx = xs[j] - xs[i];
      if (std::abs(dx) < 1e-6) continue;
      slopes.push_back((ys[j] - ys[i]) / dx);
    }
  }
  if (slopes.size() < 8) return out;
  std::nth_element(slopes.begin(), slopes.begin() + slopes.size() / 2, slopes.end());
  const double m = slopes[slopes.size() / 2];

  std::vector<double> icpt(n);
  for (size_t i = 0; i < n; ++i) icpt[i] = ys[i] - m * xs[i];
  std::nth_element(icpt.begin(), icpt.begin() + n / 2, icpt.end());
  const double b0 = icpt[n / 2];

  // 3. Covariance from the RESIDUALS, so a ragged limb widens sigma instead
  //    of silently biasing the answer.
  std::vector<double> res(n);
  for (size_t i = 0; i < n; ++i) res[i] = std::abs(ys[i] - (b0 + m * xs[i]));
  std::nth_element(res.begin(), res.begin() + n / 2, res.end());
  const double mad = 1.4826 * res[n / 2];

  double sxx = 0;
  for (double x : xs) sxx += x * x;
  out.ok = true;
  out.v_centre = b0;
  out.slope = m;
  out.residual_mad = mad;
  out.n_inliers = int(n);
  out.sigma_v = mad / std::sqrt(double(n));
  out.sigma_slope = (sxx > 0) ? mad / std::sqrt(sxx) : mad;
  return out;
}

HorizonLine predictHorizonLine(const Eigen::Matrix3d& C, double altitude_m,
                               const HorizonCamera& cam) {
  HorizonLine out;
  const double dip = horizonDip(altitude_m);
  const int height = cam.width * 4 / 5;
  const double f = (cam.width / 2.0) / std::tan(cam.hfov / 2.0);
  const double s = cam.aft ? -1.0 : 1.0;

  Eigen::Matrix3d C_b_h;
  C_b_h.col(0) = Eigen::Vector3d(0, s, 0);
  C_b_h.col(1) = Eigen::Vector3d(0, 0, 1);
  C_b_h.col(2) = Eigen::Vector3d(s, 0, 0);
  const Eigen::Matrix3d C_l_h = C * C_b_h;

  // Sample two columns and solve for the row where elevation == -dip.
  auto rowAt = [&](double xn) {
    double lo = -height, hi = height;
    for (int it = 0; it < 40; ++it) {
      const double mid = 0.5 * (lo + hi);
      const Eigen::Vector3d v(xn, mid / f, 1.0);
      const double el = nedToElevation(C_l_h * v.normalized());
      (el > -dip ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi) + height / 2.0;
  };
  const double xl = -cam.width / 4.0, xr = cam.width / 4.0;
  const double yl = rowAt(xl / f), yr = rowAt(xr / f);
  out.ok = true;
  out.slope = (yr - yl) / (xr - xl);
  out.v_centre = 0.5 * (yl + yr);
  return out;
}

TiltMeasurement measureHorizonCV(const Eigen::Matrix3d& C_true,
                                 const Eigen::Matrix3d& C_est,
                                 double altitude_m, double t,
                                 const HorizonConfig& cfg, HorizonState& state) {
  TiltMeasurement out;
  if (cfg.cameras.empty()) return out;

  auto& g = rng(cfg.seed);
  std::normal_distribution<double> nd(0.0, 1.0);
  if (state.last_t < 0) {
    state.anom_dip = cfg.anom_dip_sigma * nd(g);
    state.anom_grad = cfg.anom_grad_sigma * nd(g);
  } else {
    const double dt = std::max(0.0, t - state.last_t);
    const double a = std::exp(-dt / std::max(1e-6, cfg.anom_dip_tau));
    const double q = std::sqrt(std::max(0.0, 1.0 - a * a));
    state.anom_dip = a * state.anom_dip + q * cfg.anom_dip_sigma * nd(g);
    state.anom_grad = a * state.anom_grad + q * cfg.anom_grad_sigma * nd(g);
  }
  state.last_t = t;

  double sum_lat = 0, sum_fwd = 0, var_lat = 0, var_fwd = 0;
  int n = 0;
  for (const HorizonCamera& c : cfg.cameras) {
    HorizonCamera cc = c;
    cc.fpn_counts = c.fpn_counts * cfg.fpn_drift;

    // Per-camera anomalous component. The fore and aft tangent points are
    // ~200 km apart, so the anomaly is only PARTLY shared. Without this the
    // pair cancels perfectly and the model flatters itself -- the same trap
    // as in the analytic path.
    HorizonState local = state;
    local.anom_dip += cfg.anom_decorrelated_frac * cfg.anom_dip_sigma * nd(g);

    const Image im = renderHorizonView(C_true, altitude_m, cc, local,
                                       unsigned(cfg.seed + n * 977));
    const HorizonLine meas = detectHorizon(im, cc);
    if (!meas.ok) continue;
    const HorizonLine pred = predictHorizonLine(C_est, altitude_m, cc);

    const double f = (c.width / 2.0) / std::tan(c.hfov / 2.0);
    const double s = c.aft ? -1.0 : 1.0;
    // Line lower in the frame than predicted => nose pitched up relative to
    // the estimate. The aft camera sees both axes reversed.
    //
    // SIGN on the slope: image x is to starboard and y is down, so a positive
    // (right-wing-down) roll lifts the right of the horizon toward SMALLER
    // row indices. Positive roll therefore gives NEGATIVE slope.
    const double d_lat = s * (meas.v_centre - pred.v_centre) / f + c.mount_lat;
    const double d_fwd = -s * (meas.slope - pred.slope) + c.mount_fwd;

    sum_lat += d_lat;
    sum_fwd += d_fwd;
    var_lat += (meas.sigma_v / f) * (meas.sigma_v / f);
    var_fwd += meas.sigma_slope * meas.sigma_slope;
    ++n;
  }
  if (!n) return out;

  out.tilt(0) = sum_fwd / n;
  out.tilt(1) = sum_lat / n;
  out.R.setZero();
  out.R(0, 0) = var_fwd / (n * n) +
                (cfg.anom_grad_sigma / cfg.cameras[0].hfov) *
                    (cfg.anom_grad_sigma / cfg.cameras[0].hfov);
  out.R(1, 1) = var_lat / (n * n);
  bool fore = false, aft = false;
  for (const HorizonCamera& c : cfg.cameras) (c.aft ? aft : fore) = true;
  const double resid = (fore && aft) ? cfg.anom_decorrelated_frac * cfg.anom_dip_sigma
                                     : cfg.anom_dip_sigma;
  out.R(1, 1) += resid * resid;
  out.ok = true;
  return out;
}

}  // namespace celestial
