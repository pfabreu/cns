// test_orbit.cpp — MILESTONE 1.
//
// Reproduce Figure 7 of the paper: an uncalibrated strapdown camera produces a
// large circular error in latitude/longitude as the aircraft flies through a
// full revolution of compass heading, and averaging collapses it.
//
// No images. Trajectory and star geometry only.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "celestial/orbit.hpp"
#include "celestial/attitude.hpp"
#include "celestial/sky_model.hpp"

using namespace celestial;

namespace {

int g_failures = 0;
constexpr double kEarthR = 6371008.8;

void check(bool cond, const std::string& what) {
  if (!cond) {
    std::printf("  ** FAIL: %s\n", what.c_str());
    ++g_failures;
  }
}

OrbitConfig baseConfig() {
  OrbitConfig c;
  c.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  c.radius_m = 600.0;
  c.airspeed = 25.0;
  c.frame_rate = 10.0;
  c.revolutions = 1.0;
  return c;
}

const char* methodName(AverageMethod m) {
  switch (m) {
    case AverageMethod::Naive: return "naive mean";
    case AverageMethod::HeadingWeighted: return "heading-weighted";
    case AverageMethod::CircleFit: return "circle fit";
  }
  return "?";
}

// ---------------------------------------------------------------------------

void testGeometrySanity() {
  std::printf("\n=== TEST 1: trajectory and camera sanity ===\n");
  const OrbitConfig cfg = baseConfig();
  const auto truth = generateOrbit(cfg);
  check(!truth.empty(), "orbit generated");
  if (truth.empty()) return;

  const double period = 2.0 * M_PI * cfg.radius_m / cfg.airspeed;
  std::printf("  radius %.0f m, airspeed %.0f m/s\n", cfg.radius_m, cfg.airspeed);
  std::printf("  orbit period %.1f s, %zu frames at %.0f Hz\n", period,
              truth.size(), cfg.frame_rate);
  std::printf("  bank angle %.1f deg\n", truth.front().roll * kRad2Deg);

  Camera cam;
  std::printf("  camera: %dx%d, %.1f deg hfov, %.1f arcsec/px\n", cam.width,
              cam.height, cam.hfov * kRad2Deg,
              cam.pixelIfov() / kArcsec2Rad);

  ErrorModel err;
  err.boresight_error = 0.0;
  const auto frames = simulateObservations(truth, cam, err);
  size_t total = 0, min_n = 999, max_n = 0;
  for (const auto& f : frames) {
    total += f.id.size();
    min_n = std::min(min_n, f.id.size());
    max_n = std::max(max_n, f.id.size());
  }
  std::printf("  stars in frame: min %zu, mean %.1f, max %zu\n", min_n,
              static_cast<double>(total) / frames.size(), max_n);

  // With no misalignment at all, the per-frame fix should be near-perfect.
  const OrbitResult r =
      estimateOrbit(frames, nominalCameraMount(), AverageMethod::Naive);
  check(r.ok, "estimate succeeded with a perfect camera");
  if (r.ok) {
    std::printf("\n  perfectly calibrated camera -> %.1f m\n",
                haversine(r.position, cfg.centre));
    check(haversine(r.position, cfg.centre) < 200.0,
          "perfect camera gives a sub-200 m fix");
  }
}

void testCircleTraced() {
  std::printf("\n=== TEST 2: a body-fixed error traces a circle (Figure 7a) ===\n");
  const OrbitConfig cfg = baseConfig();
  const auto truth = generateOrbit(cfg);
  Camera cam;

  std::printf("  %-14s %-16s %-16s %-14s\n", "boresight", "single frame (km)",
              "circle radius (km)", "averaged (km)");
  for (double deg : {0.1, 0.2, 0.4, 0.8}) {
    ErrorModel err;
    err.boresight_error = deg * kDeg2Rad;
    const auto frames = simulateObservations(truth, cam, err);

    // Single frame, for reference.
    std::vector<FrameData> one{frames.front()};
    double single_km = 0.0;
    {
      const OrbitResult r1 =
          estimateOrbit(std::vector<FrameData>(frames.begin(), frames.begin() + 8),
                        nominalCameraMount(), AverageMethod::Naive);
      single_km = r1.ok ? haversine(r1.position, cfg.centre) / 1000.0 : -1.0;
    }

    const OrbitResult rc =
        estimateOrbit(frames, nominalCameraMount(), AverageMethod::CircleFit);
    const OrbitResult rn =
        estimateOrbit(frames, nominalCameraMount(), AverageMethod::Naive);

    std::printf("  %-14.2f %-16.1f %-16.1f %-14.2f\n", deg, single_km,
                rc.circle_radius * kEarthR / 1000.0,
                haversine(rn.position, cfg.centre) / 1000.0);

    // The circle radius must equal the misalignment: eps * 6371 km.
    const double expected_km = deg * kDeg2Rad * kEarthR / 1000.0;
    check(std::abs(rc.circle_radius * kEarthR / 1000.0 - expected_km) <
              0.15 * expected_km,
          "fitted circle radius matches the injected misalignment");
  }
  std::printf("\n  The circle radius recovers the misalignment directly:\n"
              "  0.4 deg -> 44.5 km, exactly as in Figure 7a.\n");
}

void testConvergence() {
  std::printf("\n=== TEST 3: the iteration converges (Figure 7a-d) ===\n");
  const OrbitConfig cfg = baseConfig();
  const auto truth = generateOrbit(cfg);
  Camera cam;

  ErrorModel err;
  err.boresight_error = 0.4 * kDeg2Rad;
  err.boresight_axis = Eigen::Vector3d(1, 0.3, 0);
  err.ahrs_bias = Eigen::Vector3d(0.05 * kDeg2Rad, 0.05 * kDeg2Rad, 0.0);
  err.ahrs_noise = 0.02 * kDeg2Rad;
  err.centroid_noise = 10.0 * kArcsec2Rad;
  const auto frames = simulateObservations(truth, cam, err);
  std::printf("  0.4 deg boresight + 0.05 deg AHRS bias + 0.02 deg AHRS noise\n"
              "  + 10 arcsec centroid noise\n\n");

  const Eigen::Matrix3d true_mount =
      nominalCameraMount() *
      smallRotation(err.boresight_axis, err.boresight_error);

  IterationConfig icfg;
  icfg.method = AverageMethod::Naive;
  icfg.max_iterations = 6;
  const auto steps =
      iterateOrbit(frames, nominalCameraMount(), cfg.centre, icfg, &true_mount);

  std::printf("  %-6s %-16s %-20s %-16s\n", "iter", "position err (km)",
              "mount err (deg)", "scatter (km)");
  for (size_t i = 0; i < steps.size(); ++i) {
    std::printf("  %-6zu %-16.2f %-20.4f %-16.2f\n", i + 1,
                steps[i].error_m / 1000.0, steps[i].mount_error_deg,
                steps[i].circle_radius_km);
  }
  check(!steps.empty(), "iteration produced steps");
  if (steps.size() >= 2) {
    check(steps.back().error_m < steps.front().error_m,
          "iteration improves the fix");
    check(steps.back().error_m < 4000.0, "converges to within 4 km");
    std::printf("\n  %.1f km -> %.2f km in %zu iterations.\n",
                steps.front().error_m / 1000.0, steps.back().error_m / 1000.0,
                steps.size());
  }
}

void testInitialConditions() {
  std::printf("\n=== TEST 4: insensitivity to initial mount guess (Table 3) ===\n");
  const OrbitConfig cfg = baseConfig();
  const auto truth = generateOrbit(cfg);
  Camera cam;
  ErrorModel err;
  err.boresight_error = 0.4 * kDeg2Rad;
  const auto frames = simulateObservations(truth, cam, err);
  const Eigen::Matrix3d true_mount =
      nominalCameraMount() *
      smallRotation(err.boresight_axis, err.boresight_error);

  std::printf("  %-24s %-18s %-12s\n", "initial mount error", "final err (km)",
              "iterations");
  for (double deg : {0.0, 5.0, 45.0, 60.0, 85.0, 120.0}) {
    const Eigen::Matrix3d guess =
        nominalCameraMount() *
        smallRotation(Eigen::Vector3d(0.6, 0.8, 0.0), deg * kDeg2Rad);
    IterationConfig icfg;
    icfg.max_iterations = 8;
    const auto steps = iterateOrbit(frames, guess, cfg.centre, icfg, &true_mount);
    if (steps.empty()) {
      std::printf("  %-24.0f  (no solution)\n", deg);
      continue;
    }
    std::printf("  %-24.0f %-18.2f %-12zu\n", deg,
                steps.back().error_m / 1000.0, steps.size());
    if (deg <= 85.0) {
      check(steps.back().error_m < 5000.0,
            "converges from within a hemisphere of tolerance");
    }
  }
  std::printf("\n  Within 90 deg the algorithm converges regardless of the\n"
              "  initial guess. Beyond it, the computed zenith flips hemisphere\n"
              "  and the solution lands near the antipode.\n");
}

void testWindAndSampling() {
  std::printf("\n=== TEST 5: wind, sampling uniformity, and the averaging rule ===\n");
  std::printf("  Section 4.4 of the paper: 54 km/h southerly, one orbit.\n");
  const double wind_n = 15.0;  // 54 km/h from the south -> blowing north
  Camera cam;
  ErrorModel err;
  err.boresight_error = 0.4 * kDeg2Rad;
  err.ahrs_noise = 0.02 * kDeg2Rad;
  err.centroid_noise = 10.0 * kArcsec2Rad;

  const Eigen::Matrix3d true_mount =
      nominalCameraMount() *
      smallRotation(err.boresight_axis, err.boresight_error);

  struct Case { const char* label; bool gps_track; };
  const std::vector<Case> cases = {{"GPS-guided ground track", true},
                                   {"fixed attitude (GPS denied)", false}};

  std::printf("\n  %-30s %-20s %-12s %-12s\n", "control mode", "averaging",
              "iter1 (km)", "final (km)");
  for (const Case& c : cases) {
    OrbitConfig cfg = baseConfig();
    cfg.wind_n = wind_n;
    cfg.gps_guided_track = c.gps_track;
    const auto truth = generateOrbit(cfg);
    const auto frames = simulateObservations(truth, cam, err);

    // Where the aircraft actually was, on average — the fixed-attitude case
    // drifts, so scoring against the nominal centre would be unfair.
    double lat = 0, lon = 0;
    for (const auto& f : truth) { lat += f.pos.lat; lon += f.pos.lon; }
    const Geodetic mean_pos{lat / truth.size(), lon / truth.size(),
                            cfg.centre.alt};

    if (!c.gps_track) {
      const double drift = haversine(truth.front().pos, truth.back().pos);
      std::printf("  (fixed-attitude drift over one orbit: %.2f km)\n",
                  drift / 1000.0);
    }

    for (AverageMethod m : {AverageMethod::Naive, AverageMethod::HeadingWeighted,
                            AverageMethod::CircleFit}) {
      IterationConfig icfg;
      icfg.method = m;
      icfg.max_iterations = 6;
      const auto steps =
          iterateOrbit(frames, nominalCameraMount(), mean_pos, icfg, &true_mount);
      const double e0 = steps.empty() ? -1.0 : steps.front().error_m / 1000.0;
      const double e = steps.empty() ? -1.0 : steps.back().error_m / 1000.0;
      std::printf("  %-30s %-20s %-12.2f %-12.2f\n", c.label, methodName(m), e0, e);
    }
    std::printf("\n");
  }
  std::printf("  The naive mean is biased under a GPS-guided track because\n"
              "  upwind arcs are prolonged in time, piling samples onto one\n"
              "  side of the circle. Weighting by heading increment, or fitting\n"
              "  the circle, removes that bias without changing the flight plan.\n");
}

void testWithRefraction() {
  std::printf("\n=== TEST 6: orbit + refraction, end to end ===\n");
  const OrbitConfig cfg = baseConfig();
  const auto truth = generateOrbit(cfg);
  Camera cam;

  ErrorModel err;
  err.boresight_error = 0.4 * kDeg2Rad;
  err.atmos = Atmosphere::isa(cfg.centre.alt, 10.0);
  err.ahrs_bias = Eigen::Vector3d(0.05 * kDeg2Rad, 0.05 * kDeg2Rad, 0.0);
  err.ahrs_noise = 0.02 * kDeg2Rad;
  err.centroid_noise = 10.0 * kArcsec2Rad;
  const auto frames = simulateObservations(truth, cam, err);

  std::printf("  boresight 0.4 deg, AHRS bias 0.05 deg, AHRS noise 0.02 deg,\n"
              "  centroid noise 10 arcsec, ISA refraction at 800 m\n\n");

  for (bool correct : {false, true}) {
    IterationConfig icfg;
    icfg.method = AverageMethod::CircleFit;
    icfg.max_iterations = 6;
    icfg.refraction_model =
        correct ? Atmosphere::isa(cfg.centre.alt, 10.0) : Atmosphere::none();
    const auto steps = iterateOrbit(frames, nominalCameraMount(), cfg.centre, icfg);
    std::printf("  refraction %-12s -> %8.2f km  (%zu iterations)\n",
                correct ? "CORRECTED" : "ignored",
                steps.empty() ? -1.0 : steps.back().error_m / 1000.0,
                steps.size());
  }
  std::printf("\n  The orbit removes the body-fixed misalignment. Refraction is\n"
              "  earth-fixed and survives it, so it still has to be modelled.\n");
}

/// The celestial compass survives mount calibration.
///
/// recalibrateMount solves C_b_c = C_l_b_est^T * R, so an AHRS yaw bias comes
/// back as mounting CLOCKING. Feed that mounting to celestialHeading, which
/// computes C_l_b = R * C_b_c^T, and it reconstructs the biased AHRS attitude
/// exactly -- the compass reports no error however large the bias. That is a
/// silent failure: the number looks like a converged zero.
///
/// mountTiltOnly drops the clocking. This pins both halves.
void testCompassSurvivesCalibration() {
  std::printf("\nCELESTIAL COMPASS vs MOUNT CALIBRATION\n");

  OrbitConfig cfg = baseConfig();
  cfg.revolutions = 1.0;
  const auto truth = generateOrbit(cfg);

  const double bias = 3.0 * kDeg2Rad;
  ErrorModel err;
  err.boresight_error = 0.4 * kDeg2Rad;   // a TILT, not a clocking
  err.ahrs_bias = Eigen::Vector3d(0, 0, bias);
  err.mag_limit = 5.0;
  const auto frames = simulateObservations(truth, Camera{}, err);

  const OrbitResult r =
      estimateOrbit(frames, nominalCameraMount(), AverageMethod::CircleFit);
  check(r.ok, "orbit fix converged");

  const Eigen::Matrix3d C_full = recalibrateMount(frames, r.position);
  // The alias axis is the local vertical in body axes, NOT body z. Passing
  // body z leaves bias*sin(bank) of the yaw behind as boresight tilt; that
  // case is asserted below.
  const Eigen::Matrix3d C_tilt = mountTiltOnly(C_full, meanVerticalBody(frames));

  const HeadingResult h_full = celestialHeading(frames, r.position, C_full);
  const HeadingResult h_tilt = celestialHeading(frames, r.position, C_tilt);
  check(h_full.ok && h_tilt.ok, "compass produced an estimate");

  const double res_full = (h_full.error + bias) * kRad2Deg;
  const double res_tilt = (h_tilt.error + bias) * kRad2Deg;
  std::printf("  bias %.2f deg   full mount -> residual %6.3f deg\n",
              bias * kRad2Deg, res_full);
  std::printf("                  tilt only  -> residual %6.3f deg\n", res_tilt);

  // The full mount absorbs the bias, so the compass recovers nearly none of it.
  // Asserted so the failure MODE is documented, not just the fix.
  check(std::abs(res_full) > 1.0,
        "full mount absorbs the yaw bias (guard: if this fails the aliasing is "
        "gone and mountTiltOnly may be unnecessary)");
  check(std::abs(res_tilt) < 0.25, "tilt-only mount recovers the yaw bias");

  // Clocking is discarded; the boresight DIRECTION must be untouched, since
  // that is what the position solver depends on.
  const Eigen::Matrix3d C_true = nominalCameraMount() *
      smallRotation(err.boresight_axis, err.boresight_error);
  auto bore = [&](const Eigen::Matrix3d& C) {
    return std::acos(std::clamp(C.col(2).dot(C_true.col(2)), -1.0, 1.0)) *
           kRad2Deg;
  };
  std::printf("  boresight error   full %.4f deg, tilt only %.4f deg\n",
              bore(C_full), bore(C_tilt));
  check(bore(C_tilt) < 0.1, "tilt-only mount recovers the boresight direction");

  // THE WRONG AXIS. Body z is correct only when a vertical reference re-levels
  // the attitude every frame; without one it leaves bias*sin(bank) of yaw
  // behind AS TILT, which lands in the boresight. Asserted so the two regimes
  // stay distinguished -- getting this backwards costs a factor of ~20 in
  // position error and shows up nowhere else.
  const Eigen::Matrix3d C_wrong = mountTiltOnly(C_full, Eigen::Vector3d::UnitZ());
  const double leak = bias * std::sin(h_tilt.mean_bank) * kRad2Deg;
  std::printf("  wrong axis (body z) -> boresight %.4f deg (predicted leak "
              "%.4f deg)\n", bore(C_wrong), leak);
  check(bore(C_wrong) > 4.0 * bore(C_tilt),
        "the wrong alias axis measurably corrupts the boresight");
}

/// A CONSTANT AHRS tilt bias is absorbable; a drifting one is not.
///
/// `recalibrateMount` books a body-fixed, time-invariant tilt bias into the
/// mounting, after which `C_l_b_est * C_b_c` is correct even though both
/// factors are wrong. Fixes then come out good with no vertical reference at
/// all -- which reads as a result and is an artefact of the error model. The
/// absorption is one-shot, so it survives only while the absorbed quantity
/// does not change.
///
/// This pins the distinction: absorption works with a constant bias and fails
/// with a drifting one. If the first check starts failing, the aliasing has
/// changed; if the second starts passing, the drift has stopped being injected
/// and the simulator is flattering absorption again.
void testTiltBiasDriftIsNotAbsorbable() {
  std::printf("\nCONSTANT vs DRIFTING AHRS TILT BIAS\n");

  OrbitConfig cfg = baseConfig();
  cfg.revolutions = 2.0;
  const auto truth = generateOrbit(cfg);

  auto fixErrorKm = [&](double drift_sigma_deg) {
    ErrorModel err;
    err.boresight_error = 0.0;  // isolate the AHRS bias
    err.ahrs_bias = Eigen::Vector3d(0.25 * kDeg2Rad, 0.15 * kDeg2Rad, 0.0);
    err.ahrs_drift_sigma = drift_sigma_deg * kDeg2Rad;
    err.ahrs_drift_tau = 60.0;
    err.mag_limit = 5.0;
    err.seed = 4;
    const auto frames = simulateObservations(truth, Camera{}, err);
    // Calibrate on the first half, then fix on the second -- absorption is only
    // meaningful if it still holds after the calibration window.
    const std::vector<FrameData> a(frames.begin(),
                                   frames.begin() + frames.size() / 2);
    const std::vector<FrameData> b(frames.begin() + frames.size() / 2,
                                   frames.end());
    const OrbitResult r0 =
        estimateOrbit(a, nominalCameraMount(), AverageMethod::CircleFit);
    if (!r0.ok) return -1.0;
    const Eigen::Matrix3d C = recalibrateMount(a, r0.position);
    const OrbitResult r1 = estimateOrbit(b, C, AverageMethod::CircleFit);
    if (!r1.ok) return -1.0;
    return haversine(r1.position, b.front().truth) / 1000.0;
  };

  const double constant = fixErrorKm(0.0);
  const double drifting = fixErrorKm(0.15);
  std::printf("  constant bias  -> %6.2f km  (absorbed by the mounting)\n",
              constant);
  std::printf("  drifting bias  -> %6.2f km  (not absorbable)\n", drifting);

  check(constant > 0.0 && drifting > 0.0, "both configurations produced a fix");
  check(constant < 2.0, "a constant tilt bias is absorbed by the mounting");
  check(drifting > 3.0 * constant,
        "a drifting tilt bias is NOT absorbed (guard: if this fails, the "
        "simulator has stopped injecting drift and is flattering absorption)");
}

/// The maneuver-coupled tilt error is proportional to bank and does not average
/// out over a loiter.
///
/// An accelerometer senses specific force, so in a coordinated turn its
/// vertical is pulled toward the aircraft's own down axis and the AHRS
/// under-reads the bank. Unlike a constant bias this is CORRELATED with
/// maneuvering, so a loiter -- which averages away body-fixed errors, and which
/// is where the mounting is calibrated -- does not remove it.
void testManeuverCoupledTilt() {
  std::printf("\nMANEUVER-COUPLED TILT ERROR\n");

  OrbitConfig cfg = baseConfig();
  cfg.revolutions = 1.0;
  const auto truth = generateOrbit(cfg);

  auto meanRollError = [&](double coupling) {
    ErrorModel err;
    err.ahrs_turn_coupling = coupling;
    const auto est = simulateAhrs(truth, err);
    double acc = 0;
    for (size_t i = 0; i < truth.size(); ++i)
      acc += dcmToEuler(est[i]).x() - truth[i].roll;
    return acc / truth.size() * kRad2Deg;
  };

  double mean_bank = 0;
  for (const FrameTruth& f : truth) mean_bank += f.roll;
  mean_bank = mean_bank / truth.size() * kRad2Deg;

  const double e0 = meanRollError(0.0);
  const double e2 = meanRollError(0.02);
  std::printf("  mean bank %.2f deg\n", mean_bank);
  std::printf("  coupling 0.00 -> mean roll error %+.4f deg\n", e0);
  std::printf("  coupling 0.02 -> mean roll error %+.4f deg  (expected %+.4f)\n",
              e2, -0.02 * mean_bank);

  check(std::abs(e0) < 1e-9, "no coupling gives no roll error");
  check(std::abs(e2 - (-0.02 * mean_bank)) < 1e-6,
        "roll error is proportional to bank");
  // The point of the term: it SURVIVES the loiter average, so a calibration
  // taken during one inherits it.
  check(std::abs(e2) > 0.1,
        "maneuver-coupled tilt does not average out over a full revolution");
}

/// PAPER REPLICATION. Teague & Chahl report 4 km from one orbit through 360 deg
/// of compass heading, on real flight data with a Cube Orange AHRS.
///
/// The orbit works because a body-fixed error traces a circle in the navigation
/// frame and averages away. The AHRS tilt DRIFT is not body-fixed, so it only
/// averages as sqrt(2 tau / T) -- which makes the orbit PERIOD the design
/// variable. At 25 m/s a 400 m orbit takes 100 s against a 60 s drift
/// correlation time and the drift survives; a 150 m orbit takes 38 s and it
/// does not. That is why radius matters, and it is not in the paper.
void testPaperReplication() {
  std::printf("\nPAPER REPLICATION: one orbit, no vertical reference\n");

  auto orbitFix = [](double radius, double revs, unsigned seed) {
    OrbitConfig cfg;
    cfg.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
    cfg.airspeed = 25.0;
    cfg.frame_rate = 10.0;              // the paper's rate
    cfg.radius_m = radius;
    cfg.revolutions = revs;
    cfg.wind_n = -4.0;
    cfg.wind_e = 3.0;
    cfg.gps_guided_track = false;
    const auto truth = generateOrbit(cfg);

    ErrorModel err;
    err.seed = seed;
    err.boresight_error = 0.4 * kDeg2Rad;   // mounting unknown at departure
    err.ahrs_bias = Eigen::Vector3d(0.25 * kDeg2Rad, 0.15 * kDeg2Rad, 0.0);
    err.ahrs_drift_sigma = 0.15 * kDeg2Rad;
    err.ahrs_drift_tau = 60.0;
    err.ahrs_turn_coupling = 0.02;
    err.ahrs_noise = 0.02 * kDeg2Rad;
    err.centroid_noise = 10.0 * kArcsec2Rad;
    err.mag_limit = 5.0;
    const auto frames = simulateObservations(truth, Camera{}, err);
    const OrbitResult r =
        estimateOrbit(frames, nominalCameraMount(), AverageMethod::Naive);
    return r.ok ? haversine(r.position, cfg.centre) / 1000.0 : -1.0;
  };

  std::printf("  %8s %8s %10s\n", "radius", "revs", "mean km");
  auto mean = [&](double radius, double revs) {
    double acc = 0;
    for (unsigned s = 1; s <= 6; ++s) acc += orbitFix(radius, revs, s);
    const double m = acc / 6;
    std::printf("  %7.0fm %8.0f %9.2f\n", radius, revs, m);
    return m;
  };
  const double tight = mean(150, 1);
  const double wide = mean(400, 1);
  const double many = mean(150, 4);

  check(tight > 0 && wide > 0 && many > 0, "all orbits produced a fix");
  check(tight < 6.0, "a tight single orbit reaches the paper's 4 km");
  check(wide > tight,
        "a wider orbit is worse: its period exceeds the drift correlation "
        "time, so the drift no longer averages out");
  check(many < tight, "more revolutions average more drift realisations");
}

/// The heading-weighted mean fixes a SAMPLING bias, not an attitude error.
///
/// Eq. (23) of the paper takes an arithmetic mean of the per-frame fixes, which
/// implicitly assumes uniform sampling in heading. A GPS-guided ground track in
/// wind does not provide that, and the bias scales with the residual circle
/// radius. Weighting each sample by its heading increment removes it.
void testHeadingWeightedMean() {
  std::printf("\nHEADING-WEIGHTED MEAN vs Eq. (23)\n");

  OrbitConfig cfg;
  cfg.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  cfg.airspeed = 25.0;
  cfg.frame_rate = 10.0;
  cfg.radius_m = 400.0;
  cfg.revolutions = 1.0;
  cfg.wind_n = -4.0;
  cfg.wind_e = 3.0;
  cfg.gps_guided_track = true;     // the case that breaks the arithmetic mean
  const auto truth = generateOrbit(cfg);

  auto run = [&](double drift, AverageMethod m) {
    double acc = 0;
    int n = 0;
    for (unsigned seed = 1; seed <= 6; ++seed) {
      ErrorModel err;
      err.seed = seed;
      err.boresight_error = 0.4 * kDeg2Rad;
      err.ahrs_bias = Eigen::Vector3d(0.25 * kDeg2Rad, 0.15 * kDeg2Rad, 0.0);
      err.ahrs_drift_sigma = drift * kDeg2Rad;
      err.ahrs_drift_tau = 60.0;
      err.mag_limit = 5.0;
      const auto frames = simulateObservations(truth, Camera{}, err);
      const OrbitResult r = estimateOrbit(frames, nominalCameraMount(), m);
      if (!r.ok) continue;
      acc += haversine(r.position, cfg.centre) / 1000.0;
      ++n;
    }
    return n ? acc / n : -1.0;
  };

  const double naive_clean = run(0.0, AverageMethod::Naive);
  const double wtd_clean = run(0.0, AverageMethod::HeadingWeighted);
  const double naive_drift = run(0.15, AverageMethod::Naive);
  const double wtd_drift = run(0.15, AverageMethod::HeadingWeighted);
  std::printf("  no AHRS drift : naive %6.2f km, heading-weighted %6.2f km\n",
              naive_clean, wtd_clean);
  std::printf("  with drift    : naive %6.2f km, heading-weighted %6.2f km\n",
              naive_drift, wtd_drift);

  check(wtd_clean < naive_clean / 5.0,
        "heading weighting removes the sampling bias (>5x) when it dominates");
  // Honest bound: the win does NOT survive realistic AHRS drift, because drift
  // is not body-fixed and heading weighting does nothing for it. Asserted so
  // the earlier 36x headline is not quoted outside the conditions that produced
  // it.
  check(wtd_drift > naive_drift * 0.8,
        "with drift modelled the weighting advantage is masked (guard on the "
        "retired 36x claim)");
}

/// STAR-AIDED ATTITUDE: the mechanism, in simulation only.
///
/// Three attempts to remove the AHRS tilt error from the FIX failed, all
/// against the same wall: tilt and position are one observable to a celestial
/// fix. This is a different quantity. Star directions in ECEF do not depend on
/// position, so an absolute star attitude carries none, and the SEQUENCE of
/// them observes GYRO BIAS -- a rate, not an offset.
///
/// This asserts only that the filter works: unaided gyro attitude drifts
/// without bound, aided attitude does not. It does NOT show the transit gets
/// better, because that depends on a real gyro and on what EKF3 already does.
/// See attitude.hpp. Do not quote a number from here.
void testStarAidedAttitude() {
  std::printf("\nSTAR-AIDED ATTITUDE (mechanism only, simulated gyro)\n");

  OrbitConfig cfg = baseConfig();
  cfg.radius_m = 150.0;
  cfg.revolutions = 8.0;      // ~5 min, long enough for drift to show
  cfg.frame_rate = 10.0;
  const auto truth = generateOrbit(cfg);

  ErrorModel err;
  err.seed = 3;
  err.mag_limit = 5.0;
  const auto gyro = simulateGyro(truth, err);
  check(gyro.size() == truth.size(), "gyro samples generated");

  auto trueAttitude = [&](size_t i) {
    return nedBasisEcef(truth[i].pos) *
           eulerToDcm(truth[i].roll, truth[i].pitch, truth[i].yaw);
  };

  // The star measurement: body -> ECEF, built the way the real pipeline would,
  // with noise standing in for centroid and identification error.
  std::mt19937 rng(11);
  const double star_sigma = 60.0 * kArcsec2Rad;
  std::normal_distribution<double> sn(0.0, star_sigma);

  StarAidedAttitude aided;
  StarAidedAttitude unaided;   // same filter, never updated after the seed
  aided.setAttitude(trueAttitude(0));
  unaided.setAttitude(trueAttitude(0));

  double err_aided = 0, err_unaided = 0;
  double peak_aided = 0, peak_unaided = 0;
  int n = 0;
  for (size_t i = 1; i < truth.size(); ++i) {
    const double dt = truth[i].t - truth[i - 1].t;
    aided.predict(gyro[i - 1], dt);
    unaided.predict(gyro[i - 1], dt);

    // A star fix every 2 s, as a real 10 Hz camera with matching would give.
    if (i % 20 == 0) {
      const Eigen::Vector3d nz(sn(rng), sn(rng), sn(rng));
      const Eigen::Matrix3d meas =
          trueAttitude(i) * smallRotation(nz.normalized(), nz.norm());
      aided.update(meas, star_sigma);
    }

    const double ea =
        rotationAngleBetween(aided.attitude(), trueAttitude(i)) * kRad2Deg;
    const double eu =
        rotationAngleBetween(unaided.attitude(), trueAttitude(i)) * kRad2Deg;
    err_aided += ea; err_unaided += eu;
    peak_aided = std::max(peak_aided, ea);
    peak_unaided = std::max(peak_unaided, eu);
    ++n;
  }

  std::printf("  over %.0f s of flight:\n", truth.back().t);
  std::printf("    gyro alone   mean %7.4f deg, peak %7.4f deg\n",
              err_unaided / n, peak_unaided);
  std::printf("    star-aided   mean %7.4f deg, peak %7.4f deg\n",
              err_aided / n, peak_aided);
  std::printf("    estimated gyro bias %.4f deg/hr\n",
              aided.gyroBias().norm() * kRad2Deg * 3600.0);

  check(peak_unaided > 0.05, "unaided gyro attitude drifts measurably");
  check(peak_aided < peak_unaided / 3.0,
        "star aiding bounds the attitude error");
  check(aided.attitudeSigma().norm() * kRad2Deg < 0.2,
        "aided attitude covariance converges");
}

/// The circle fit must not be iterated.
///
/// `iterateOrbit` recalibrates the mounting between passes, which is right for
/// the averaging methods and WRONG for the circle fit: the fitted radius IS the
/// misalignment, so recalibrating from the fit's own output applies the same
/// correction twice and the fixed-point iteration can walk to a wrong answer.
///
/// Measured across eight real ArduPilot SITL orbits, four GPS-aided and four
/// GNSS-denied: iteration 1 gives 6.64 km, "converged" gives 9.37 km, and
/// iterating helps in only 2 of 8. Worst observed 5.13 -> 15.19 km, and aided
/// 31.68 -> 67.53 km.
void testCircleFitIsNotIterated() {
  std::printf("\nCIRCLE FIT IS NOT ITERATED\n");

  OrbitConfig cfg = baseConfig();
  cfg.radius_m = 250.0;
  cfg.revolutions = 1.0;
  const auto truth = generateOrbit(cfg);

  int worse = 0, n = 0;
  for (unsigned seed = 1; seed <= 6; ++seed) {
    ErrorModel err;
    err.seed = seed;
    err.boresight_error = 0.4 * kDeg2Rad;
    err.ahrs_bias = Eigen::Vector3d(0.25 * kDeg2Rad, 0.15 * kDeg2Rad,
                                    0.9 * kDeg2Rad);
    err.ahrs_drift_sigma = 0.30 * kDeg2Rad;   // as measured GNSS-denied
    err.ahrs_drift_tau = 60.0;
    err.ahrs_turn_coupling = 0.004;
    err.mag_limit = 5.0;
    const auto frames = simulateObservations(truth, Camera{}, err);

    IterationConfig ic;
    ic.method = AverageMethod::CircleFit;
    ic.max_iterations = 6;
    const auto steps = iterateOrbit(frames, nominalCameraMount(), cfg.centre, ic);
    if (steps.empty()) continue;
    ++n;
    if (steps.size() > 1) ++worse;

    // What it WOULD have done, with recalibration left on.
    IterationConfig ir = ic;
    ir.recalibrate_mount = true;
    // (recalibrate_mount is overridden for CircleFit inside iterateOrbit, so
    // compare against the averaging path instead to show the loop still runs.)
    IterationConfig ih = ic;
    ih.method = AverageMethod::HeadingWeighted;
    const auto hsteps =
        iterateOrbit(frames, nominalCameraMount(), cfg.centre, ih);
    std::printf("  seed %u: circle %zu step(s), %6.2f km   |   "
                "heading-wtd %zu steps, %6.2f km\n",
                seed, steps.size(), steps.back().error_m / 1000.0,
                hsteps.size(), hsteps.empty() ? -1.0
                                              : hsteps.back().error_m / 1000.0);
  }

  check(n > 0, "produced fixes");
  check(worse == 0,
        "circle fit returns a single pass and does not iterate");
  // The averaging methods still iterate -- the change is specific to the fit
  // that already solves for misalignment.
}

}  // namespace

int main() {
  std::printf("Milestone 1 — analytic orbit, no images\n");
  std::printf("======================================\n");

  testGeometrySanity();
  testCircleTraced();
  testConvergence();
  testInitialConditions();
  testWindAndSampling();
  testWithRefraction();
  testCircleFitIsNotIterated();
  testCompassSurvivesCalibration();
  testTiltBiasDriftIsNotAbsorbable();
  testManeuverCoupledTilt();
  testPaperReplication();
  testHeadingWeightedMean();
  testStarAidedAttitude();

  std::printf("\n======================================\n");
  if (g_failures == 0) {
    std::printf("ALL CHECKS PASSED\n");
    return 0;
  }
  std::printf("%d CHECK(S) FAILED\n", g_failures);
  return 1;
}
