// Timing benchmark for the star pipeline.
//
// Accuracy without a frame rate is half an answer. The matched filter roughly
// triples identified stars, and it also costs ~50x the plain detector -- which
// on measurement puts it BELOW the camera frame rate. That is a deployment
// blocker, not a footnote, so it gets measured here rather than discovered on
// the aircraft.
//
//   ./build/bench            # all stages
//   ./build/bench --repeat 20
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>
#include <functional>
#include <vector>

#include "celestial/imaging.hpp"
#include "celestial/orbit.hpp"
#include "celestial/star_catalog.hpp"

using namespace celestial;

namespace {

double timeIt(int repeat, const std::function<void()>& f) {
  f();  // warm
  const auto t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < repeat; ++i) f();
  const auto t1 = std::chrono::steady_clock::now();
  return std::chrono::duration<double, std::milli>(t1 - t0).count() / repeat;
}

void row(const char* name, double ms, double budget_ms, const char* note) {
  std::printf("  %-26s %8.2f ms %8.1f Hz   %-3s  %s\n", name, ms,
              ms > 0 ? 1000.0 / ms : 0.0, ms <= budget_ms ? "ok" : "SLOW", note);
}

}  // namespace

int main(int argc, char** argv) {
  int repeat = 5;
  for (int i = 1; i < argc; ++i)
    if (!std::strcmp(argv[i], "--repeat") && i + 1 < argc)
      repeat = std::atoi(argv[++i]);

  OrbitConfig cfg;
  cfg.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  cfg.airspeed = 25.0;
  cfg.frame_rate = 10.0;
  cfg.radius_m = 150.0;
  cfg.revolutions = 0.05;
  cfg.gps_guided_track = false;
  const auto truth = generateOrbit(cfg);

  const auto cat = catalogBrighterThan(5.5);
  StarField field(cat, truth[0].epoch);
  std::vector<double> vmags;
  for (size_t i = 0; i < field.size(); ++i) vmags.push_back(field.vmag(i));

  Camera cam;
  const Eigen::Matrix3d C_b_c = nominalCameraMount();
  SensorModel sensor;
  sensor.exposure_s = 0.10;

  const RenderedFrame rf = renderFrame(truth[0], truth[1], 0.1, field, vmags,
                                       cam, C_b_c, sensor, 1u);

  const Eigen::Matrix3d Ca =
      eulerToDcm(truth[0].roll, truth[0].pitch, truth[0].yaw);
  const Eigen::Matrix3d Cb =
      eulerToDcm(truth[1].roll, truth[1].pitch, truth[1].yaw);
  const Eigen::AngleAxisd aa(Ca.transpose() * Cb);
  const Eigen::Vector3d w_cam =
      C_b_c.transpose() * (aa.axis() * aa.angle() * cfg.frame_rate);

  DetectorConfig plain, mf;
  mf.omega_cam = w_cam;
  mf.exposure_s = sensor.exposure_s;
  mf.focal_px = cam.focalPx();

  // The budget: one frame period at the camera rate. Anything slower than this
  // cannot keep up, and frames must be dropped.
  const double budget = 1000.0 / cfg.frame_rate;

  std::printf("\nStar pipeline timing -- %d x %d (%.2f Mpx), %.0f Hz camera, "
              "%.0f ms budget, %d reps\n\n",
              rf.image.width, rf.image.height,
              rf.image.width * rf.image.height / 1e6, cfg.frame_rate, budget,
              repeat);
  std::printf("  %-26s %11s %11s %5s  %s\n", "stage", "per frame", "max rate",
              "", "note");

  row("render (simulator only)",
      timeIt(repeat,
             [&] {
               renderFrame(truth[0], truth[1], 0.1, field, vmags, cam, C_b_c,
                           sensor, 2u);
             }),
      1e9, "not flown; a real camera replaces this");

  const double t_plain =
      timeIt(repeat, [&] { detectStars(rf.image, plain); });
  row("detect, plain", t_plain, budget, "median+MAD, connected components");

  const double t_mf = timeIt(repeat, [&] { detectStars(rf.image, mf); });
  char note[160];
  std::snprintf(note, sizeof note, "%.0fx plain -- SEE README.md", t_mf / t_plain);
  row("detect, matched filter", t_mf, budget, note);

  const auto dets = detectStars(rf.image, mf);
  MatcherConfig mc;
  FrameData fd;
  row("match to catalogue",
      timeIt(repeat,
             [&] {
               FrameData f;
               matchDetections(dets, truth[0].epoch, Ca, C_b_c, cam, field,
                               truth[0].pos, mc, f);
             }),
      budget, "nearest neighbour + ambiguity guard");

  matchDetections(dets, truth[0].epoch, Ca, C_b_c, cam, field, truth[0].pos, mc,
                  fd);
  fd.t = 0;
  fd.epoch = truth[0].epoch;
  fd.C_l_b_est = Ca;
  row("per-frame fix (RANSAC)",
      timeIt(repeat,
             [&] {
               Geodetic p;
               singleFrameFix(fd, C_b_c, Atmosphere::none(), p);
             }),
      budget, "3-star minimal set, 100 iterations");

  std::printf("\n  TOTAL flight path (detect + match + fix):\n");
  std::printf("    plain          %7.2f ms  -> %5.1f Hz\n", t_plain,
              1000.0 / t_plain);
  std::printf("    matched filter %7.2f ms  -> %5.1f Hz\n", t_mf,
              1000.0 / t_mf);
  std::printf("\n  Measured on this machine. A Raspberry Pi 5 is several times\n"
              "  slower, so treat anything near the budget here as over it\n"
              "  there. See the matched-filter section in README.md for the\n"
              "  polar-warp plan that makes the filter separable.\n\n");
  return 0;
}
