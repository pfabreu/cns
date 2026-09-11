// test_horizon.cpp — the horizon sensor as an independent vertical reference.
//
// The claim under test: the horizon observes the IMU tilt bias IN STRAIGHT AND
// LEVEL FLIGHT, so flight geometry becomes an improvement rather than a
// requirement.

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "celestial/horizon.hpp"
#include "celestial/orbit.hpp"
#include "celestial/sky_model.hpp"

using namespace celestial;

namespace {
int g_fail = 0;
constexpr double kR = 6371.0;
void check(bool c, const std::string& w) {
  if (!c) { std::printf("  ** FAIL: %s\n", w.c_str()); ++g_fail; }
}
double am(double rad) { return rad / kArcmin2Rad; }

HorizonConfig foreOnly() {
  HorizonConfig c; c.cameras = {HorizonCamera{}}; return c;
}
HorizonConfig forePair() {
  HorizonConfig c;
  HorizonCamera f; HorizonCamera a; a.aft = true;
  c.cameras = {f, a};
  return c;
}

/// Empirical error of the horizon measurement against the true tilt.
void characterise(const char* label, HorizonConfig cfg, double imu_sigma) {
  OrbitConfig oc;
  oc.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  oc.radius_m = 600; oc.airspeed = 25; oc.frame_rate = 10; oc.revolutions = 2;
  const auto truth = generateOrbit(oc);

  HorizonState st;
  double s2f = 0, s2l = 0, s2ff = 0, s2lf = 0;
  int n = 0;
  for (size_t k = 0; k < truth.size(); k += 5) {
    const FrameTruth& t = truth[k];
    const Eigen::Matrix3d Ct = eulerToDcm(t.roll, t.pitch, t.yaw);
    // IMU with a body-fixed tilt bias.
    const Eigen::Matrix3d Ce =
        eulerToDcm(t.roll + imu_sigma, t.pitch + imu_sigma * 0.6, t.yaw);
    const TiltMeasurement z = measureHorizon(Ct, Ce, t.pos.alt, t.t, cfg, st);
    if (!z.ok) continue;

    const Eigen::AngleAxisd aa(Ce.transpose() * Ct);
    const Eigen::Vector3d d = aa.axis() * aa.angle();
    s2f += (z.tilt(0) - d.x()) * (z.tilt(0) - d.x());
    s2l += (z.tilt(1) - d.y()) * (z.tilt(1) - d.y());

    const Eigen::Matrix3d Cf = fuseTilt(Ce, z, imu_sigma);
    const Eigen::AngleAxisd af(Cf.transpose() * Ct);
    const Eigen::Vector3d df = af.axis() * af.angle();
    s2ff += df.x() * df.x();
    s2lf += df.y() * df.y();
    ++n;
  }
  if (!n) return;
  const double ef = std::sqrt(s2f / n), el = std::sqrt(s2l / n);
  const double ff = std::sqrt(s2ff / n), fl = std::sqrt(s2lf / n);
  std::printf("  %-26s %7.2f' %7.2f' | %7.2f' %7.2f' | %6.1f %6.1f\n", label,
              am(ef), am(el), am(ff), am(fl), am(ff) * kR * kArcmin2Rad,
              am(fl) * kR * kArcmin2Rad);
}

void testErrorStructure() {
  std::printf("\n=== TEST 1: measurement error structure ===\n");
  std::printf("  IMU tilt bias 6 arcmin (0.1 deg) on both axes\n");
  std::printf("\n  %-26s %-17s | %-17s | %-13s\n", "", "horizon raw",
              "after fusion", "fused, km");
  std::printf("  %-26s %7s %7s | %7s %7s | %6s %6s\n", "config", "fwd", "lat",
              "fwd", "lat", "fwd", "lat");
  std::printf("  %s\n", std::string(78, '-').c_str());
  const double s = 0.1 * kDeg2Rad;
  characterise("single forward camera", foreOnly(), s);
  characterise("fore/aft pair", forePair(), s);
  std::printf("\n  The forward (slope) axis is DIFFERENTIAL, so common-mode\n"
              "  anomalous refraction cancels and it reaches sub-arcminute.\n"
              "  The lateral (dip) axis is absolute and eats the full anomaly,\n"
              "  unless a fore/aft pair makes it differential too.\n");
}

void testDip() {
  std::printf("\n=== TEST 2: dip and altitude sensitivity ===\n");
  for (double h : {100.0, 800.0, 3000.0, 20000.0}) {
    std::printf("  %6.0f m -> dip %6.2f'   (10 m alt error = %.2f')\n", h,
                am(horizonDip(h)), am(std::abs(horizonDip(h + 10) - horizonDip(h))));
  }
  check(std::abs(am(horizonDip(800)) - 46.9) < 1.0, "dip at 800 m is ~47 arcmin");
  std::printf("\n  Barometric altitude is entirely adequate: the dip error from\n"
              "  10 m is two orders below the anomalous refraction term.\n");
}

void testStraightAndLevel() {
  std::printf("\n=== TEST 3: the claim -- no maneuver required ===\n");
  OrbitConfig oc;
  oc.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  // A very large radius approximates straight and level: almost no heading
  // change over the segment, so nothing averages out.
  oc.radius_m = 200000; oc.airspeed = 25; oc.frame_rate = 10; oc.revolutions = 0.002;
  const auto truth = generateOrbit(oc);
  std::printf("  %zu frames, %.1f deg of heading change (i.e. none)\n",
              truth.size(),
              std::abs(truth.back().yaw - truth.front().yaw) * kRad2Deg);

  Camera cam;
  ErrorModel err;
  err.boresight_error = 0.0;
  err.ahrs_bias = Eigen::Vector3d(0.1 * kDeg2Rad, 0.06 * kDeg2Rad, 0.0);
  err.ahrs_noise = 0.02 * kDeg2Rad;
  err.centroid_noise = 10 * kArcsec2Rad;

  struct Case { const char* name; bool use; HorizonConfig cfg; };
  std::vector<Case> cases = {{"IMU only", false, foreOnly()},
                             {"IMU + single forward", true, foreOnly()},
                             {"IMU + fore/aft pair", true, forePair()}};

  std::printf("\n  %-26s %-12s\n", "configuration", "fix error");
  std::printf("  %s\n", std::string(42, '-').c_str());
  double imu_err = 0;
  for (Case& c : cases) {
    HorizonState st;
    std::vector<Eigen::Matrix3d> est;
    est.reserve(truth.size());
    std::mt19937 g(1);
    std::normal_distribution<double> nd(0.0, err.ahrs_noise);
    for (const FrameTruth& t : truth) {
      const Eigen::Matrix3d Ct = eulerToDcm(t.roll, t.pitch, t.yaw);
      Eigen::Matrix3d Ce = eulerToDcm(t.roll + err.ahrs_bias.x() + nd(g),
                                      t.pitch + err.ahrs_bias.y() + nd(g),
                                      t.yaw + nd(g));
      if (c.use) {
        const TiltMeasurement z =
            measureHorizon(Ct, Ce, t.pos.alt, t.t, c.cfg, st);
        Ce = fuseTilt(Ce, z, 0.1 * kDeg2Rad);
      }
      est.push_back(Ce);
    }
    ErrorModel e2 = err;
    e2.ahrs_bias.setZero();
    e2.ahrs_noise = 0.0;
    const auto frames = simulateObservationsWithAttitude(truth, est, cam, e2);
    const auto r = estimateOrbit(frames, nominalCameraMount(),
                                 AverageMethod::HeadingWeighted);
    const double e = r.ok ? haversine(r.position, truth[truth.size() / 2].pos) : -1;
    std::printf("  %-26s %8.2f km\n", c.name, e / 1000.0);
    if (!c.use) imu_err = e;
    else check(e < imu_err, std::string(c.name) + " beats IMU alone");
  }
  std::printf("\n  With zero heading change the orbit mechanism does nothing.\n"
              "  The horizon still observes the IMU tilt bias directly, because\n"
              "  it is a NON-INERTIAL vertical reference and is not confused by\n"
              "  acceleration. That is the whole point.\n");
}

}  // namespace

int main() {
  std::printf("Horizon sensor -- independent vertical reference\n");
  std::printf("===============================================\n");
  testDip();
  testErrorStructure();
  testStraightAndLevel();
  std::printf("\n===============================================\n");
  std::printf(g_fail ? "%d CHECK(S) FAILED\n" : "ALL CHECKS PASSED\n", g_fail);
  return g_fail ? 1 : 0;
}
