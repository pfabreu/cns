// Training-data generator for the learned detector.
//
//   ./build/gen_dataset out/ --frames 2000 --seed 1
//
// Writes, per frame, a 16-bit PGM image and an 8-bit PGM mask. The LABELS
// ALREADY EXIST: `StarLabel` carries the true centroid, both streak endpoints
// and a truncation flag, so a segmentation mask is the (u0,v0)->(u1,v1) segment
// rasterised and dilated by the PSF. Nothing here is new science; it is
// plumbing that reuses the renderer.
//
// SPLIT BY TRAJECTORY, NOT BY FRAME. Consecutive frames are nearly identical,
// so a random frame-level split leaks the validation set into training and
// reports a number that means nothing. Each --seed produces one trajectory;
// generate several and hold whole seeds out.
//
// See "Learned detector" in NOTES-private.md, and in particular the
// warning that a network trained on cloud we invented will beat a classical
// filter on cloud we invented. This generator is for building the pipeline,
// not for proving the network works.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "celestial/imaging.hpp"
#include "celestial/orbit.hpp"
#include "celestial/star_catalog.hpp"

using namespace celestial;

namespace {

bool writePgm16(const std::string& path, const Image& img) {
  FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) return false;
  std::fprintf(f, "P5\n%d %d\n65535\n", img.width, img.height);
  std::vector<uint8_t> row(size_t(img.width) * 2);
  for (int y = 0; y < img.height; ++y) {
    for (int x = 0; x < img.width; ++x) {
      const uint16_t v = img.at(x, y);
      row[2 * x] = uint8_t(v >> 8);      // PGM is big-endian
      row[2 * x + 1] = uint8_t(v & 0xff);
    }
    std::fwrite(row.data(), 1, row.size(), f);
  }
  std::fclose(f);
  return true;
}

bool writePgm8(const std::string& path, const std::vector<uint8_t>& m, int w,
               int h) {
  FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) return false;
  std::fprintf(f, "P5\n%d %d\n255\n", w, h);
  std::fwrite(m.data(), 1, m.size(), f);
  std::fclose(f);
  return true;
}

/// Mask: the streak segment, dilated by the PSF. Truncated streaks are left
/// OUT -- a partial star is not something the network should learn to call a
/// star, and the detector discards them too.
std::vector<uint8_t> maskFromLabels(const RenderedFrame& rf, int w, int h,
                                    double sigma) {
  std::vector<uint8_t> m(size_t(w) * h, 0);
  const int rad = std::max(1, int(std::ceil(2.0 * sigma)));
  for (const StarLabel& L : rf.labels) {
    if (L.truncated) continue;
    const double dx = L.u1 - L.u0, dy = L.v1 - L.v0;
    const int n = std::max(2, int(std::ceil(std::hypot(dx, dy))) + 1);
    for (int i = 0; i < n; ++i) {
      const double t = double(i) / (n - 1);
      const int cx = int(std::lround(L.u0 + t * dx));
      const int cy = int(std::lround(L.v0 + t * dy));
      for (int j = -rad; j <= rad; ++j)
        for (int k = -rad; k <= rad; ++k) {
          const int x = cx + k, y = cy + j;
          if (x < 0 || y < 0 || x >= w || y >= h) continue;
          if (k * k + j * j <= rad * rad) m[size_t(y) * w + x] = 255;
        }
    }
  }
  return m;
}

}  // namespace

int main(int argc, char** argv) {
  std::string out = "dataset";
  int frames = 200;
  unsigned seed = 1;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--frames" && i + 1 < argc) frames = std::atoi(argv[++i]);
    else if (a == "--seed" && i + 1 < argc) seed = unsigned(std::atoi(argv[++i]));
    else if (a[0] != '-') out = a;
  }

  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u01(0.0, 1.0);

  // One trajectory per seed, randomised so the set spans the conditions the
  // detector will actually meet.
  OrbitConfig cfg;
  cfg.centre = Geodetic::fromDegrees(-30.89 + 40.0 * (u01(rng) - 0.5),
                                     136.56 + 120.0 * (u01(rng) - 0.5), 800.0);
  cfg.airspeed = 20.0 + 15.0 * u01(rng);
  cfg.frame_rate = 10.0;
  cfg.radius_m = 120.0 + 500.0 * u01(rng);
  cfg.revolutions = 1.0 + 3.0 * u01(rng);
  cfg.start_track = 2 * M_PI * u01(rng);
  cfg.gps_guided_track = false;
  const auto truth = generateOrbit(cfg);
  if (truth.size() < 2) { std::fprintf(stderr, "no trajectory\n"); return 1; }

  const auto cat = catalogBrighterThan(6.0);
  StarField field(cat, truth[truth.size() / 2].epoch);
  std::vector<double> vmags;
  for (size_t i = 0; i < field.size(); ++i) vmags.push_back(field.vmag(i));

  Camera cam;
  const Eigen::Matrix3d C_b_c =
      nominalCameraMount() *
      smallRotation(Eigen::Vector3d(1, 0, 0), 0.4 * kDeg2Rad);

  std::printf("seed %u: r=%.0fm %.1f rev, %zu frames available, writing %d\n",
              seed, cfg.radius_m, cfg.revolutions, truth.size(), frames);

  const size_t stride = std::max<size_t>(1, (truth.size() - 1) / frames);
  int written = 0;
  for (size_t i = 0; i + 1 < truth.size() && written < frames; i += stride) {
    SensorModel sensor;
    sensor.seed = seed * 100000u + unsigned(i);
    sensor.exposure_s = 0.02 + 0.18 * u01(rng);
    // Weather. Half the frames clear, so the network sees both.
    if (u01(rng) < 0.5) {
      sensor.cloud_amount = 0.7 * u01(rng);
      sensor.cloud_scale_px = 150.0 + 400.0 * u01(rng);
    }
    if (u01(rng) < 0.3) sensor.flare_amount = 0.6 * u01(rng);

    const RenderedFrame rf =
        renderFrame(truth[i], truth[i + 1], 1.0 / cfg.frame_rate, field, vmags,
                    cam, C_b_c, sensor, unsigned(i + 1));
    const auto mask = maskFromLabels(rf, rf.image.width, rf.image.height,
                                     sensor.psf_sigma_px);

    char base[256];
    std::snprintf(base, sizeof base, "%s/s%03u_f%05d", out.c_str(), seed,
                  written);
    if (!writePgm16(std::string(base) + "_img.pgm", rf.image)) {
      std::fprintf(stderr, "cannot write to %s -- does the directory exist?\n",
                   out.c_str());
      return 1;
    }
    writePgm8(std::string(base) + "_mask.pgm", mask, rf.image.width,
              rf.image.height);
    ++written;
    if (written % 25 == 0) {
      std::printf("\r  %d/%d", written, frames);
      std::fflush(stdout);
    }
  }
  std::printf("\r  wrote %d frame pairs to %s/\n", written, out.c_str());
  return 0;
}
