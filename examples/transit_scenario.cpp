// transit_scenario.cpp — the deterministic end-to-end demonstration.
//
//   ./build/transit_scenario
//
// A single process: no sockets, no MAVLink, no real time. It calls the same
// library functions the live node calls, so the same build prints the same
// numbers every time. That makes it reproducible, CI-able, and reviewable in
// thirty seconds without installing ArduPilot.
//
// ArduPilot SITL remains the VALIDATION -- it proves the system works against
// a real EKF3, with real timing and real sensor errors. This proves the
// algorithm works, repeatably. Different jobs.
//
// The scenario: a long over-water transit, GNSS-denied. One calibrating loiter
// after departure, then straight legs. Dead reckoning drifts; each celestial
// fix pulls it back.

#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "celestial/deadreckon.hpp"
#include "celestial/imaging.hpp"
#include "celestial/horizon.hpp"
#include "celestial/orbit.hpp"
#include "celestial/sky_model.hpp"
#include "celestial/star_catalog.hpp"

using namespace celestial;

namespace {

struct Args {
  double legs_km = 220.0;      // total over-water distance
  double fix_interval_s = 500.0;
  /// --imaging runs the real star pipeline: render, detect, centroid, match.
  /// Default off because it is ~50x slower and this scenario is also a
  /// regression test. Off, the star vectors are ideal, so the transit assumes
  /// perfect detection and matching -- run --imaging to find what that costs.
  bool imaging = false;
  std::string horizon;
  std::string csv_path;
  std::string dr_quality = "realistic";
  double boresight_deg = 0.4;  // uncalibrated camera at departure
  /// --exposure OVERRIDES the SensorModel default, in seconds. <=0 keeps it.
  ///
  /// WHY THIS EXISTS. celestial_node holds smear at a SETPOINT
  /// (pipeline.hpp target_smear_px) by shortening exposure, which in a turn
  /// starves the detector: measured in SITL, det 43 -> 0 as the controller
  /// wound 200 ms down to 35 ms. This scenario has no such controller -- its
  /// exposure is fixed -- so the photons/smear trade cannot be swept without
  /// a knob. Exposure drives BOTH the rendered photon count and the streak
  /// length, so this is the real trade and not a proxy for it.
  double exposure_s = 0.0;


  double drift_sigma_deg = 0.15;  ///< AHRS tilt bias drift, 1-sigma deg
  double drift_tau = 60.0;        ///< s
  int legs = 3;                   ///< transit legs, each ending in a fix orbit
  double fix_radius_m = 150.0;    ///< fix orbit radius; short period beats drift
  // 4, not 2. The fix orbit must be long enough to average several drift
  // realisations: measured on this transit, 2 revs gives 13.3 km per fix and
  // 4 gives 8.5 km. 150 s of orbiting per fix is the price of navigating
  // without a vertical reference.
  double fix_revs = 4.0;
  // 2 Hz, not the paper's 10. MEASURED: rate barely matters. Over four seeds
  // the mean per-fix error is 7.77 / 7.10 / 7.12 km at 2 / 5 / 10 Hz, and the
  // transit result does not improve at all. The fix error is a BIAS over the
  // orbit -- AHRS drift -- and more frames do not average a bias, the same
  // reason star-camera quality does not move the number either. 10 Hz costs 3x
  // the runtime for nothing, and this scenario is also a regression test.
  //
  // The paper's 10 Hz IS used where the claim is made: testPaperReplication.
  // `--rate 10` here if you want it.
  double frame_rate = 2.0;
  unsigned seed = 1;              ///< ErrorModel seed; drift makes this matter
  /// Fraction of bank leaking into the AHRS tilt estimate.
  ///
  /// MEASURED in ArduPilot SITL, not guessed: weighted fit over +5..+45 deg of
  /// bank across 100/200/400 m loiters, ~30k samples, gives 0.0009. The
  /// previous value of 0.02 was invented and was 20x too pessimistic.
  ///
  /// TWO CAVEATS. It is not linear -- implied coupling is 0.0002 at 15-20 deg
  /// rising to 0.0018 at 45 deg, so compensation DEGRADES with bank rather
  /// than leaking a fixed fraction. And it was measured with GPS ENABLED, so
  /// EKF3 had good velocity for the centripetal correction. GNSS-denied it has
  /// only airspeed and a wind estimate, and this number will be worse. Treat
  /// 0.001 as an optimistic bound until it is remeasured with GPS denied.
  double turn_coupling = 0.004;

  // Take heading from the STARS rather than the magnetometer. The compass is
  // the last non-autonomous sensor in the loop and the one an adversary can
  // spoof; the stars already contain the information.
  bool celestial_heading = true;
  double mag_bias_deg = 3.0;    // magnetometer heading bias to inject
};

/// Cascais to Cabo de Sao Vicente, offshore: one calibrating loiter after
/// departure, then straight legs south. Any long over-water route works; this
/// one is ~220 km, within the endurance of a large-battery electric platform.
std::vector<FrameTruth> buildTransit(const Args& a, const Epoch& start) {
  std::vector<FrameTruth> out;
  OrbitConfig base;
  base.centre = Geodetic::fromDegrees(38.70, -9.48, 800.0);  // Cabo Raso
  base.airspeed = 25.0;
  base.frame_rate = a.frame_rate;
  base.wind_n = -4.0;
  base.wind_e = 3.0;

  auto append = [&](OrbitConfig c) {
    c.gps_guided_track = false;   // fixed attitude: GNSS-denied
    auto seg = generateOrbit(c);
    if (seg.empty()) return;
    const double t0 = out.empty() ? 0.0 : out.back().t;
    const Geodetic anchor = out.empty() ? seg.front().pos : out.back().pos;
    const double dlat = anchor.lat - seg.front().pos.lat;
    const double dlon = anchor.lon - seg.front().pos.lon;
    for (FrameTruth f : seg) {
      f.t += t0;
      f.pos.lat += dlat;
      f.pos.lon += dlon;
      f.epoch = start.advanced(f.t);
      out.push_back(f);
    }
  };

  // 1. One calibrating loiter. The camera mounting is unknown at departure and
  //    nothing but a heading sweep can estimate it.
  OrbitConfig loiter = base;
  // 250 m, and this moved once turn coupling was MEASURED rather than guessed.
  //
  // The trade-off is real: a FIX wants a short period, because AHRS drift only
  // averages out if the orbit is quick against its correlation time; a
  // CALIBRATION wants low bank, because maneuver coupling puts a
  // bank-proportional tilt into the AHRS that recalibrateMount books into the
  // mounting. With the invented coupling of 0.02 the second term dominated and
  // 400 m won. At the SITL-measured 0.0009 it is ~20x weaker and the optimum
  // moves in: 150 m gives 9.35 km, 250 m gives 7.30 km, 400 m gives 8.84 km.
  //
  // Caveat: 0.0009 was measured with GPS aiding EKF3, so it is an optimistic
  // bound and this optimum may move again under real GNSS denial.
  //
  // A FIX wants a SHORT PERIOD: the AHRS drift only averages out if the orbit
  // is quick compared to its 60 s correlation time. 150 m is 38 s and gives
  // 4.14 km; 400 m is 100 s and gives 8.80 km (testPaperReplication).
  //
  // A CALIBRATION wants LOW BANK: maneuver coupling puts a bank-proportional
  // tilt into the AHRS, and recalibrateMount books it into the mounting, where
  // it then poisons straight flight. 150 m banks at 23 deg and leaves the
  // boresight at 0.26 deg; 400 m banks at 9 deg and reaches 0.02 deg.
  //
  // This scenario's loiter exists to CALIBRATE, so it takes the wide orbit.
  // A system that only ever fixes in orbits should take the tight one and skip
  // calibration entirely -- the averaging removes the mounting error anyway.
  loiter.radius_m = 250.0;
  loiter.revolutions = 2.0;
  append(loiter);

  // 2. Transit: straight legs, each followed by a FIX ORBIT.
  //
  // Without a vertical reference a straight-leg fix carries the full AHRS tilt
  // error and is worth ~40 km. The orbit is the only thing that removes it, so
  // the aircraft must maneuver to navigate. That is the paper's method and it
  // is the operating mode this scenario now flies: leg, orbit, leg, orbit.
  //
  // The fix orbit is TIGHT (150 m, 38 s) so the AHRS drift stays correlated
  // across it. That is the opposite of what the calibration loiter wanted, and
  // it is fine: calibration happens once, at departure, at low bank.
  const double leg_km = a.legs_km / a.legs;
  for (int k = 0; k < a.legs; ++k) {
    OrbitConfig leg = base;
    leg.radius_m = 200000.0;
    leg.revolutions = (leg_km * 1000.0) / (2 * M_PI * leg.radius_m);
    leg.start_track = 3.3 + 0.15 * k;   // roughly south, slight fan
    append(leg);

    OrbitConfig fix = base;
    fix.radius_m = a.fix_radius_m;
    fix.revolutions = a.fix_revs;
    append(fix);
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    const std::string s = argv[i];
    auto nx = [&] { return (i + 1 < argc) ? argv[++i] : ""; };
    if (s == "--csv") a.csv_path = nx();
    else if (s == "--horizon") a.horizon = nx();
    else if (s == "--imaging") a.imaging = true;
    else if (s == "--drift-sigma") a.drift_sigma_deg = std::atof(nx());
    else if (s == "--drift-tau") a.drift_tau = std::atof(nx());
    else if (s == "--seed") a.seed = (unsigned)std::atoi(nx());
    else if (s == "--legs") a.legs = std::atoi(nx());
    else if (s == "--fix-radius") a.fix_radius_m = std::atof(nx());
    else if (s == "--fix-revs") a.fix_revs = std::atof(nx());
    else if (s == "--rate") a.frame_rate = std::atof(nx());
    else if (s == "--turn-coupling") a.turn_coupling = std::atof(nx());
    else if (s == "--magnetic-heading") a.celestial_heading = false;
    else if (s == "--mag-bias") a.mag_bias_deg = std::atof(nx());
    else if (s == "--km") a.legs_km = std::atof(nx());
    else if (s == "--interval") a.fix_interval_s = std::atof(nx());
    else if (s == "--dr-quality") a.dr_quality = nx();
    else if (s == "--boresight") a.boresight_deg = std::atof(nx());
    else if (s == "--exposure") a.exposure_s = std::atof(nx());
  }

  auto epoch0 = Epoch::fromUtc(2024, 11, 15, 23, 0, 0.0);
  if (!epoch0) return 1;
  const auto truth = buildTransit(a, *epoch0);
  if (truth.size() < 100) { std::printf("no trajectory\n"); return 1; }

  double ground_km = 0;
  for (size_t k = 1; k < truth.size(); ++k)
    ground_km += haversine(truth[k - 1].pos, truth[k].pos) / 1000.0;

  std::printf("Celestial navigation, GNSS-denied transit\n");
  std::printf("=========================================\n");
  std::printf("  %.0f km over water, %.1f h, %zu frames at %.0f Hz\n",
              ground_km, truth.back().t / 3600.0, truth.size(), a.frame_rate);
  std::printf("  camera boresight %.2f deg (uncalibrated at departure)\n",
              a.boresight_deg);
  std::printf("  heading from %s\n",
              a.celestial_heading ? "the STARS (celestial compass)"
                                  : "the magnetometer");
  std::printf("  star pipeline %s\n",
              a.imaging ? "render -> detect -> centroid -> match"
                        : "ideal vectors (perfect detection and matching)");

  // --- sensors -----------------------------------------------------------
  Camera cam;
  ErrorModel err;
  err.boresight_error = a.boresight_deg * kDeg2Rad;
  // The yaw component is the MAGNETOMETER bias. GNSS-denied, EKF3 takes yaw
  // from the compass and looks up declination at the last known position, so
  // in reality this term grows as position drifts.
  err.ahrs_bias = Eigen::Vector3d(0.25 * kDeg2Rad, 0.15 * kDeg2Rad,
                                  a.mag_bias_deg * kDeg2Rad);
  err.ahrs_noise = 0.02 * kDeg2Rad;
  err.ahrs_drift_sigma = a.drift_sigma_deg * kDeg2Rad;
  err.ahrs_drift_tau = a.drift_tau;
  err.ahrs_turn_coupling = a.turn_coupling;
  err.seed = a.seed;
  err.centroid_noise = 10.0 * kArcsec2Rad;
  err.mag_limit = 5.0;



  // AHRS attitude, straight from the autopilot. Everything downstream sees
  // only this -- never the truth. There is NO vertical reference: constant
  // bias + Gauss-Markov drift + maneuver coupling + per-frame noise, exactly
  // as the paper's Cube Orange delivered it. The orbit is what removes it.
  FILE* csv = a.csv_path.empty() ? nullptr : std::fopen(a.csv_path.c_str(), "w");
  if (csv)
    std::fprintf(csv, "t,lat_true,lon_true,dr_lat,dr_lon,dr_open_lat,"
                      "dr_open_lon,dr_err_m,orbit_lat,orbit_lon,orbit_err_m\n");

  std::vector<Eigen::Matrix3d> est = simulateAhrs(truth, err);

  // Optional horizon sensor: a NON-INERTIAL vertical reference, which unlike an
  // accelerometer is not confused by the acceleration of a turn. This is the
  // deterministic place to measure it -- `fake_sitl` samples a UDP stream and
  // its runs are not reproducible, so an A/B there is meaningless.
  if (!a.horizon.empty()) {
    HorizonConfig hc;
    hc.seed = a.seed * 7919u;
    hc.cameras.clear();
    if (a.horizon == "lepton") hc.cameras.push_back(HorizonCamera::lepton25());
    else if (a.horizon == "pair") {
      hc.cameras.push_back(HorizonCamera::lepton25(false));
      hc.cameras.push_back(HorizonCamera::lepton25(true));
    } else {
      std::fprintf(stderr, "--horizon: lepton | pair\n");
      return 2;
    }
    HorizonState hs;
    for (size_t i = 0; i < truth.size(); ++i) {
      const Eigen::Matrix3d Ct =
          eulerToDcm(truth[i].roll, truth[i].pitch, truth[i].yaw);
      const TiltMeasurement z = measureHorizon(Ct, est[i], truth[i].pos.alt,
                                               truth[i].t, hc, hs);
      if (z.ok) est[i] = fuseTilt(est[i], z, 0.3 * kDeg2Rad);
    }
  }
  double tilt_raw = 0;
  for (size_t fi = 0; fi < truth.size(); ++fi) {
    const FrameTruth& f = truth[fi];
    const Eigen::Matrix3d Ct = eulerToDcm(f.roll, f.pitch, f.yaw);
    const Eigen::AngleAxisd aa(est[fi].transpose() * Ct);
    const Eigen::Vector3d d = aa.axis() * aa.angle();
    tilt_raw += std::hypot(d.x(), d.y());
  }
  tilt_raw /= truth.size();

  ErrorModel e2 = err;
  e2.ahrs_bias.setZero();
  e2.ahrs_noise = 0.0;
  std::vector<FrameData> frames;
  int n_det = 0, n_matched = 0, n_dropped = 0;
  if (!a.imaging) {
    frames = simulateObservationsWithAttitude(truth, est, cam, e2);
  } else {
    // FULL STAR PIPELINE: render -> detect -> centroid -> match, per frame.
    // Nothing here sees truth except renderFrame, which stands in for the sky.
    const Eigen::Matrix3d C_b_c_true =
        nominalCameraMount() * smallRotation(err.boresight_axis,
                                             err.boresight_error);
    const auto cat = catalogBrighterThan(err.mag_limit);
    StarField field(cat, truth[truth.size() / 2].epoch);
    std::vector<double> vmags;
    for (size_t i = 0; i < field.size(); ++i) vmags.push_back(field.vmag(i));

    SensorModel sensor;
    if (a.exposure_s > 0) sensor.exposure_s = a.exposure_s;
    double smear_sum = 0;
    DetectorConfig dcfg;
    MatcherConfig mcfg;
    const double period = 1.0 / a.frame_rate;

    // The matcher needs an assumed position. It gets the DEPARTURE point and
    // never a truth value -- a stale prior is the realistic case, and at
    // 107 arcsec/px even a few degrees of error is well inside the match
    // radius.
    const Geodetic assumed = truth.front().pos;

    for (size_t fi = 0; fi + 1 < truth.size(); ++fi) {
      const RenderedFrame rf =
          renderFrame(truth[fi], truth[fi + 1], period, field, vmags, cam,
                      C_b_c_true, sensor, unsigned(fi + 1));
      const auto dets = detectStars(rf.image, dcfg);
      FrameData fd;
      const int m =
          matchDetections(dets, truth[fi].epoch, est[fi], nominalCameraMount(),
                          cam, field, assumed, mcfg, fd);
      fd.t = truth[fi].t;
      fd.epoch = truth[fi].epoch;
      fd.yaw_est = dcmToEuler(est[fi]).z();
      fd.C_l_b_est = est[fi];
      fd.truth = truth[fi].pos;
      {
        const Eigen::Matrix3d Ca = eulerToDcm(truth[fi].roll, truth[fi].pitch,
                                              truth[fi].yaw);
        const Eigen::Matrix3d Cb = eulerToDcm(truth[fi + 1].roll,
                                              truth[fi + 1].pitch,
                                              truth[fi + 1].yaw);
        const Eigen::AngleAxisd aa(Ca.transpose() * Cb);
        const Eigen::Vector3d w = C_b_c_true.transpose() *
                                  (aa.axis() * aa.angle() * a.frame_rate);
        smear_sum += w.norm() * sensor.exposure_s * cam.focalPx();
      }
      n_det += (int)dets.size();
      n_matched += m;
      if (m < 3) ++n_dropped;
      frames.push_back(fd);
    }
    std::printf("  star pipeline: %.1f detected, %.1f matched per frame; "
                "%.1f%% of frames unusable (<3 stars)\n",
                (double)n_det / frames.size(), (double)n_matched / frames.size(),
                100.0 * n_dropped / frames.size());
    std::printf("  exposure %.0f ms, mean smear %.1f px\n",
                sensor.exposure_s * 1000.0, smear_sum / frames.size());
  }

  // --- run the estimator over the flight ---------------------------------
  const DeadReckonConfig drcfg =
      (a.dr_quality == "optimistic") ? DeadReckonConfig::optimistic()
      : (a.dr_quality == "poor")     ? DeadReckonConfig::poor()
                                     : DeadReckonConfig::realistic();
  DeadReckoner dr{drcfg}, dr_open{drcfg};

  Eigen::Matrix3d C_b_c = nominalCameraMount();
  const Eigen::Matrix3d C_true =
      nominalCameraMount() *
      smallRotation(err.boresight_axis, err.boresight_error);
  Geodetic assumed = truth.front().pos;
  bool calibrated = false;

  std::vector<FrameData> window;
  double rev_start_t = frames.front().t;
  double last_fix_t = 0.0, heading_accum = 0.0, prev_yaw = frames[0].yaw_est;
  int n_fix = 0;
  double sum_dr = 0, sum_open = 0, max_open = 0;
  int n_dr = 0;
  // Heading correction estimated from the stars, applied to dead reckoning.
  double hdg_corr = 0.0, hdg_resid = 0.0;
  int n_hdg = 0;

  std::printf("\n  %5s %8s %9s %9s %10s %9s\n", "fix", "t (s)", "fix err",
              "DR err", "DR alone", "boresight");
  std::printf("  %s\n", std::string(58, '-').c_str());

  for (size_t k = 1; k < frames.size(); ++k) {
    const FrameData& fd = frames[k];
    const double dt = fd.t - frames[k - 1].t;
    if (dr.started() && dt > 0) {
      dr.predict(dt, 25.0, fd.yaw_est + hdg_corr);
      dr_open.predict(dt, 25.0, fd.yaw_est + hdg_corr);
      // CSV for the live viewer. Written every frame the DR has, so the
      // trajectory is continuous rather than only at fixes.
      if (csv) {
        std::fprintf(csv, "%.2f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.1f,0,0,-1\n",
                     fd.t, fd.truth.lat * kRad2Deg, fd.truth.lon * kRad2Deg,
                     dr.position().lat * kRad2Deg, dr.position().lon * kRad2Deg,
                     dr_open.position().lat * kRad2Deg,
                     dr_open.position().lon * kRad2Deg,
                     haversine(dr.position(), fd.truth));
      }
      const double e1 = haversine(dr.position(), fd.truth);
      const double e2b = haversine(dr_open.position(), fd.truth);
      sum_dr += e1; sum_open += e2b; ++n_dr;
      max_open = std::max(max_open, e2b);
    }
    // Accumulate heading ONLY while actually turning. A revolution counted
    // across a straight leg produces a window that spans leg and orbit, and
    // the leg frames all share one heading, so they drag the average toward
    // the uncorrected fix. Separation is clean: a 150 m orbit at 25 m/s turns
    // at 9.5 deg/s, a 200 km transit leg at 0.007 deg/s.
    const double dyaw = std::atan2(std::sin(fd.yaw_est - prev_yaw),
                                   std::cos(fd.yaw_est - prev_yaw));
    prev_yaw = fd.yaw_est;
    const double rate = dt > 0 ? std::abs(dyaw) / dt : 0.0;
    if (rate < 1.0 * kDeg2Rad) {          // straight: not in an orbit
      heading_accum = 0.0;
      rev_start_t = fd.t;
    } else {
      heading_accum += dyaw;
    }
    if (fd.id.size() >= 5) {
      FrameData f2 = fd;
      // Transport reference: per-frame fixes are averaged as if from one
      // place, so hand estimateOrbit the dead-reckoned displacement and let it
      // remove the aircraft's own motion.
      if (dr.started()) f2.dr_ne = dr.positionNE();
      window.push_back(f2);
    }
    // BOUND THE WINDOW BY TIME. estimateOrbit averages the per-frame fixes as
    // if they came from one place, which holds in a loiter (400 m) and fails
    // badly on a straight leg: at 25 m/s a 500 s window spans 12.5 km, and an
    // unbounded one spans the whole flight. Left uncapped the "fix" lagged the
    // aircraft by 130 km.
    //
    // The proper fix is to dead-reckon each frame's fix forward to a common
    // epoch before averaging. Capping the window is the cheap approximation:
    // 60 s is 1.5 km of travel, so the averaging lags by well under a km.
    while (!window.empty() && fd.t - window.front().t > 300.0)
      window.erase(window.begin());

    // A FIX IS AN ORBIT. With no vertical reference a straight-leg fix carries
    // the whole AHRS tilt error -- 40 km -- and averaging over a window that
    // spans both a leg and an orbit is worse than either, because the leg
    // frames contribute one heading and drag the mean toward the uncorrected
    // value. So: fix only on a completed revolution, and average only the
    // frames from that revolution.
    const bool full_rev = std::abs(heading_accum) >= 2 * M_PI;
    if (!full_rev) continue;

    std::vector<FrameData> orbit;
    for (const FrameData& f : window)
      if (f.t >= rev_start_t) orbit.push_back(f);
    if (orbit.size() < 20) { heading_accum = 0; continue; }

    const AverageMethod m = calibrated ? AverageMethod::HeadingWeighted
                                       : AverageMethod::CircleFit;
    const OrbitResult r = estimateOrbit(orbit, C_b_c, m);
    if (r.ok) {
      ++n_fix;
      last_fix_t = fd.t;
      const double ferr = haversine(r.position, fd.truth);
      const double bore =
          std::acos(std::clamp(C_b_c.col(2).dot(C_true.col(2)), -1.0, 1.0)) *
          kRad2Deg;

      if (!dr.started()) { dr.start(r.position); dr_open.start(r.position); }
      else dr.update(r.position, std::clamp(r.spread * 6371000.0, 1500.0, 30000.0));

      if (csv)
        std::fprintf(csv, "%.2f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.1f,%.7f,%.7f,%.1f\n",
                     fd.t, fd.truth.lat * kRad2Deg, fd.truth.lon * kRad2Deg,
                     dr.position().lat * kRad2Deg, dr.position().lon * kRad2Deg,
                     dr_open.position().lat * kRad2Deg,
                     dr_open.position().lon * kRad2Deg,
                     haversine(dr.position(), fd.truth),
                     r.position.lat * kRad2Deg, r.position.lon * kRad2Deg, ferr);
      std::printf("  %5d %8.0f %8.2fk %8.2fk %9.2fk %9.4f\n", n_fix, fd.t,
                  ferr / 1000, haversine(dr.position(), fd.truth) / 1000,
                  haversine(dr_open.position(), fd.truth) / 1000, bore);

      // CELESTIAL COMPASS. The same Kabsch that calibrates the mounting also
      // yields the true heading, so take it rather than the magnetometer's.
      // This removes the last non-autonomous sensor from the loop.
      if (a.celestial_heading) {
        const HeadingResult h = celestialHeading(orbit, r.position, C_b_c);
        if (h.ok && std::abs(h.error) < 20 * kDeg2Rad) {
          hdg_corr = h.error;
          hdg_resid = h.error + a.mag_bias_deg * kDeg2Rad;
          ++n_hdg;
        }
      }

      assumed = r.position;
      // CALIBRATE ONCE, at the departure loiter. Every later orbit is tight
      // and banks hard, and maneuver coupling puts a bank-proportional tilt
      // into the AHRS that recalibrateMount books straight into the mounting.
      // Recalibrating at each fix orbit made the boresight wander between 0.11
      // and 0.46 deg and the fixes with it. The mounting is a bracket: it does
      // not change, so estimate it once where the bank is lowest.
      if (!calibrated && orbit.size() >= 30) {
        // Calibrate the mounting, then keep only its TILT. recalibrateMount
        // divides by the AHRS attitude, so a magnetometer yaw bias comes back
        // as mounting CLOCKING -- and a mounting that already contains the bias
        // reproduces the AHRS heading exactly, leaving the compass nothing to
        // find. Dropping clocking keeps the tilt the position solver needs and
        // leaves heading to the stars. Clocking is then assumed nominal; see
        // README.md for what that assumption buys and costs.
        const Eigen::Matrix3d C_cal = recalibrateMount(orbit, r.position);

        // mountTiltOnly is a REPAIR, not a refinement, and it must not run
        // when there is nothing to repair.
        //
        // It exists because a yaw bias gets absorbed into the mounting as
        // clocking, hiding it from the compass. That absorption happens when
        // the AHRS attitude carries a tilt error. With a horizon sensor it does
        // not: measured, 20 seeds, 3 deg yaw bias, boresight error after
        // calibration --
        //
        //                     recalibrateMount   + mountTiltOnly
        //     no horizon           0.905 deg        0.286 deg   <- repaired
        //     1x Lepton            0.056 deg        0.678 deg   <- damaged
        //
        // With a good vertical reference the calibration comes out clean, and
        // stripping the vertical component then removes real mounting tilt.
        // Identical with dip and sensor noise switched off, so this is
        // structural rather than a noise artefact.
        // KEPT ON even with a horizon, despite the damage measured above,
        // because dropping it kills the compass outright: the heading residual
        // against a 3 deg magnetometer bias goes from 0.006 to 2.875 deg -- the
        // mounting absorbs the bias as clocking and the compass has nothing
        // left to find. Better calibration is not worth losing heading.
        //
        // So this is an unresolved TRADE-OFF, not a bug with a known fix:
        //
        //     mountTiltOnly ON   compass 0.006 deg, boresight 0.678 deg
        //     mountTiltOnly OFF  compass 2.875 deg, boresight 0.056 deg
        //
        // Two hypotheses tested and both WRONG: it is not anomalous dip (the
        // damage is identical with the atmosphere off) and it is not the choice
        // of removal axis (stripping about the camera boresight instead of the
        // mean local vertical is a no-op, 0.056 vs 0.056). What is actually
        // needed is a calibration that solves for tilt and clocking JOINTLY
        // rather than solving for both and subtracting one afterwards.
        C_b_c = a.celestial_heading
                    ? mountTiltOnly(C_cal, meanVerticalBody(orbit))
                    : C_cal;
        calibrated = true;
      }
    }
    heading_accum = 0;
    rev_start_t = fd.t;
  }

  if (csv) std::fclose(csv);

  std::printf("\n  RESULT\n");
  if (a.celestial_heading && n_hdg) {
    std::printf("  heading            magnetometer bias %.2f deg -> residual"
                " %.3f deg  (%d fixes)\n",
                a.mag_bias_deg, hdg_resid * kRad2Deg, n_hdg);
  } else {
    std::printf("  heading            magnetometer, %.2f deg bias UNCORRECTED\n",
                a.mag_bias_deg);
  }
  std::printf("  AHRS tilt error      %6.2f'  (%.1f km if not averaged out)\n",
              tilt_raw / kArcmin2Rad, tilt_raw * 6371.0);
  std::printf("  celestial fixes      %d over %.0f km\n", n_fix, ground_km);
  if (n_dr) {
    std::printf("  dead reckoning ALONE mean %6.2f km, peak %6.2f km"
                "  (%.1f%% of distance)\n",
                sum_open / n_dr / 1000, max_open / 1000,
                100.0 * max_open / (ground_km * 1000));
    std::printf("  dead reckoning + FIXES mean %6.2f km"
                "                    BOUNDED\n", sum_dr / n_dr / 1000);
  }
  std::printf("\n  Unaided dead reckoning grows without limit; the celestial\n"
              "  fixes bound it. That is the whole claim.\n");
  return 0;
}
