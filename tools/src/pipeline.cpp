#include "pipeline.hpp"

#include <sys/stat.h>

#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cstdio>

#include "celestial/sky_model.hpp"
#include "celestial/star_catalog.hpp"

namespace celestial {
namespace {
double wrapPi(double a) { return std::atan2(std::sin(a), std::cos(a)); }
}  // namespace

Pipeline::Pipeline(const PipelineConfig& c) : cfg_(c) {
  int y, mo, d, h, mi;
  double sec;
  std::sscanf(cfg_.utc.c_str(), "%d-%d-%dT%d:%d:%lf", &y, &mo, &d, &h, &mi, &sec);
  if (auto e = Epoch::fromUtc(y, mo, d, h, mi, sec)) epoch0_ = *e;

  catalog_ = catalogBrighterThan(cfg_.mag_limit);
  for (const CatalogStar& s : catalog_) vmags_.push_back(s.vmag);
  // Render deep, match shallow: mutual nearest-neighbour breaks when
  // predictions vastly outnumber detections.
  const double mm = (cfg_.match_mag > 0) ? cfg_.match_mag
                                         : std::max(4.5, cfg_.mag_limit - 1.0);
  match_catalog_ = catalogBrighterThan(mm);

  sensor_.exposure_s = (cfg_.exposure_s > 0) ? cfg_.exposure_s : 0.05;
  sensor_.blur_substeps = 12;
  dcfg_.centroid = CentroidMethod::GaussianFit;
  mcfg_.radius_px = cfg_.match_radius_px;



  const DeadReckonConfig dc =
      (cfg_.dr_quality == "optimistic") ? DeadReckonConfig::optimistic()
      : (cfg_.dr_quality == "poor")     ? DeadReckonConfig::poor()
                                        : DeadReckonConfig::realistic();
  dr_ = DeadReckoner(dc);
  dr_open_ = DeadReckoner(dc);

  C_b_c_true_ = nominalCameraMount() *
                smallRotation(Eigen::Vector3d(1, 0, 0), cfg_.boresight_deg * kDeg2Rad);
  C_b_c_assumed_ = nominalCameraMount();
}

FrameReport Pipeline::step(const MavSample& s, FixReport& fix) {
  FrameReport r;
  fix = FixReport{};
  const double t = s.t;

  // Dead reckoning runs on EVERY sample, not only rendered frames: it carries
  // the estimate between fixes and wants the highest rate available.
  if (cfg_.dead_reckon && dr_.started() && s.airspeed > 1.0 && prev_t_ > 0) {
    const double dt = t - prev_t_;
    dr_.predict(dt, s.airspeed, s.yaw_est);
    dr_open_.predict(dt, s.airspeed, s.yaw_est);
  }
  if (s.airspeed > 1.0) prev_t_ = t;

  FrameTruth now;
  now.t = t;
  now.pos = s.true_pos;
  now.pos.alt = s.alt_m;
  now.roll = s.true_roll;
  now.pitch = s.true_pitch;
  now.yaw = s.true_yaw;
  now.epoch = epoch0_.advanced(t);

  static double last_render = -1e9;
  if (t - last_render < 1.0 / cfg_.render_hz) return r;

  if (!have_prev_) {
    prev_ = now; prev_C_ = s.C_l_b_est; prev_yaw_ = s.yaw_est;
    have_prev_ = true; return r;
  }
  // A stale previous pose would make the blur span the whole gap while the
  // exposure is still scaled by the nominal period: stars smear across the
  // frame and nothing is detected.
  const double gap = t - prev_.t;
  if (gap > 1.6 / cfg_.render_hz || gap <= 0) {
    prev_ = now; prev_C_ = s.C_l_b_est; prev_yaw_ = s.yaw_est;
    last_render = t; return r;
  }
  last_render = t;
  if (!have_assumed_) {
    assumed_ = now.pos;
    have_assumed_ = true;
    last_fix_t_ = t;
    // START DEAD RECKONING HERE, at the KNOWN departure position, rather than
    // at the first celestial fix.
    //
    // You always know where you took off -- GNSS denial happens in flight, not
    // on the ramp -- so beginning the dead reckoning at the first fix throws
    // that away and inherits the fix's error instead. Measured: the first fix
    // was 27 km off, so the DR track began 27 km from truth and spent the rest
    // of the flight recovering from a mistake it never needed to make.
    if (cfg_.dead_reckon) {
      dr_.start(now.pos);
      dr_open_.start(now.pos);
    }
  }

  // Auto exposure: perpendicular rates translate the field, the rate about the
  // boresight rotates it. Heavily low-passed, because gyro noise and turbulence
  // put high-frequency content in the rates that does not smear.
  if (cfg_.exposure_s <= 0) {
    const double half_diag =
        0.5 * std::hypot(cam_.hfov, cam_.hfov * cam_.height / cam_.width);
    const double inst = std::hypot(s.omega.x(), s.omega.y()) +
                        std::abs(s.omega.z()) * half_diag;
    const double a = std::exp(-(1.0 / cfg_.render_hz) / 3.0);
    smear_rate_ = a * smear_rate_ + (1 - a) * inst;
    double e = (smear_rate_ > 1e-6)
                   ? cfg_.target_smear_px * cam_.pixelIfov() / smear_rate_
                   : cfg_.max_exposure_s;
    sensor_.exposure_s = std::clamp(std::min(e, 0.9 / cfg_.render_hz),
                                    cfg_.min_exposure_s, cfg_.max_exposure_s);
  }

  const StarField field(catalog_, prev_.epoch);
  const StarField match_field(match_catalog_, prev_.epoch);

  // ---- SIMULATED SKY: the only use of the true pose --------------------
  const RenderedFrame rf = renderFrame(prev_, now, gap, field, vmags_, cam_,
                                       C_b_c_true_, sensor_, unsigned(t * 1000));
  const FrameTruth exp = prev_;
  const Eigen::Matrix3d exp_C = prev_C_;
  const double exp_yaw = prev_yaw_;
  prev_ = now; prev_C_ = s.C_l_b_est; prev_yaw_ = s.yaw_est;

  // Mid-exposure state: renderFrame integrates FORWARD, so each streak's
  // flux-weighted centroid sits at the middle of the exposure.
  const double half = 0.5 * sensor_.exposure_s;
  const Eigen::Vector3d rv = s.omega * half;
  const double ang = rv.norm();
  const Eigen::Matrix3d dC =
      (ang < 1e-12) ? Eigen::Matrix3d::Identity()
                    : Eigen::AngleAxisd(ang, rv / ang).toRotationMatrix();
  const Eigen::Matrix3d C_mid = exp_C * dC;
  const Epoch epoch_mid = exp.epoch.advanced(half);
  const Eigen::Matrix3d C_true_ = eulerToDcm(exp.roll, exp.pitch, exp.yaw);

  auto rvec = [&](const Eigen::Matrix3d& C) {
    const Eigen::AngleAxisd aa(C.transpose() * C_true_);
    return Eigen::Vector3d(aa.axis() * aa.angle());
  };

  // The AHRS attitude, optionally corrected by a horizon sensor.
  //
  // With no horizon this is a pass-through and the tilt error is whatever the
  // orbit maneuver can average away. With one, the horizon observes that tilt
  // error DIRECTLY -- it is a non-inertial vertical reference, so unlike an
  // accelerometer it is not confused by the acceleration of a turn, which is
  // exactly when the orbit needs it most.
  Eigen::Matrix3d C_nav = C_mid;
  if (!cfg_.horizon.cameras.empty()) {
    const TiltMeasurement z = measureHorizon(C_true_, C_mid, exp.pos.alt,
                                             exp.t, cfg_.horizon, hstate_);
    if (z.ok) {
      C_nav = fuseTilt(C_mid, z, cfg_.imu_tilt_sigma);
      r.horizon_used = true;
    }
  }

  // ---- NAVIGATION: pixels + AHRS attitude only --------------------------
  const auto dets = detectStars(rf.image, dcfg_);
  FrameData fd;
  const int matched = matchDetections(dets, epoch_mid, C_nav, C_b_c_assumed_,
                                      cam_, match_field, assumed_, mcfg_, fd);
  fd.t = exp.t;
  fd.yaw_est = exp_yaw;
  fd.truth = exp.pos;
  if (dr_.started()) fd.dr_ne = dr_.positionNE();

  r.rendered = true;
  if (dr_.started()) {
    r.dr_pos = dr_.position();
    r.dr_valid = true;
  }
  if (dr_open_.started()) {
    r.dr_open_pos = dr_open_.position();
    r.dr_open_valid = true;
  }
  r.detected = int(dets.size());
  r.matched = matched;
  r.exposure_s = sensor_.exposure_s;
  r.smear_px = smear_rate_ * sensor_.exposure_s / cam_.pixelIfov();
  r.tilt_imu = rvec(C_mid).head<2>() / kArcmin2Rad;
  if (matched >= 3 &&
      singleFrameFix(fd, C_b_c_assumed_, Atmosphere::none(), r.frame_fix)) {
    r.frame_err_m = haversine(r.frame_fix, exp.pos);
  }

  // ---- optional frame dump for the live viewer ---------------------------
  //
  // Writes the image as PGM and a CSV of detections alongside the catalogue
  // position each one was matched to. That pairing is the interesting part:
  // it shows what the matcher accepted and, by omission, what it rejected.
  if (!cfg_.frames_dir.empty() && (dump_count_++ % cfg_.frames_every) == 0) {
    if (!dump_ready_) {
      // Create it rather than failing silently. A missing directory used to
      // look exactly like the feature not working, because fopen just returned
      // null and the dump was skipped.
      ::mkdir(cfg_.frames_dir.c_str(), 0755);
      const std::string probe = cfg_.frames_dir + "/.w";
      if (FILE* t = std::fopen(probe.c_str(), "w")) {
        std::fclose(t);
        std::remove(probe.c_str());
        std::fprintf(stderr, "frames -> %s/ (every %d frames)\n",
                     cfg_.frames_dir.c_str(), cfg_.frames_every);
      } else {
        std::fprintf(stderr,
                     "WARNING: cannot write to %s/ -- frame dump disabled\n",
                     cfg_.frames_dir.c_str());
        cfg_.frames_dir.clear();
      }
      dump_ready_ = true;
    }
    if (cfg_.frames_dir.empty()) return r;
    char base[512];
    std::snprintf(base, sizeof base, "%s/f%06d", cfg_.frames_dir.c_str(),
                  dump_index_);
    std::string pgm = std::string(base) + ".pgm";
    std::string tmp = pgm + ".part";      // write-then-rename, so the viewer
    if (FILE* f = std::fopen(tmp.c_str(), "wb")) {   // never reads a half file
      std::fprintf(f, "P5\n%d %d\n65535\n", rf.image.width, rf.image.height);
      std::vector<uint8_t> row(size_t(rf.image.width) * 2);
      for (int y = 0; y < rf.image.height; ++y) {
        for (int x = 0; x < rf.image.width; ++x) {
          const uint16_t v = rf.image.at(x, y);
          row[2 * x] = uint8_t(v >> 8);
          row[2 * x + 1] = uint8_t(v & 0xff);
        }
        std::fwrite(row.data(), 1, row.size(), f);
      }
      std::fclose(f);
      std::rename(tmp.c_str(), pgm.c_str());
    }
    const std::string csv = std::string(base) + ".csv";
    const std::string ctmp = csv + ".part";
    if (FILE* f = std::fopen(ctmp.c_str(), "w")) {
      // Every detection, flagged with whether the matcher accepted it. Matched
      // stars are projected back to pixels from the camera-frame directions
      // FrameData carries, which is the same quantity the fix consumes -- so
      // what the viewer draws is what the estimator actually used, not a
      // parallel reconstruction of it.
      std::fprintf(f, "u,v,flux,matched,id\n");
      for (const auto& d : dets) {
        int id = -1;
        for (size_t k = 0; k < fd.id.size(); ++k) {
          double mu = 0, mv = 0;
          if (!rayToPixel(fd.v_cam[k], cam_, mu, mv)) continue;
          if (std::hypot(mu - d.u, mv - d.v) < 2.0) { id = fd.id[k]; break; }
        }
        std::fprintf(f, "%.2f,%.2f,%.1f,%d,%d\n", d.u, d.v, d.flux,
                     id >= 0 ? 1 : 0, id);
      }
      std::fclose(f);
      std::rename(ctmp.c_str(), csv.c_str());
    }
    r.dumped = base;
    ++dump_index_;
  }

  if (matched >= cfg_.min_stars) {
    window_.push_back(fd);
    window_hdg_.push_back(heading_accum_);
    window_truth_.push_back(exp.pos);
  }
  if (!yaw_init_) { prev_yaw_acc_ = exp_yaw; yaw_init_ = true; }
  heading_accum_ += wrapPi(exp_yaw - prev_yaw_acc_);
  prev_yaw_acc_ = exp_yaw;
  r.heading_deg = heading_accum_ * kRad2Deg;
  r.window = window_.size();

  // ---- fix: on a completed revolution, or on the timer -------------------
  const bool full_rev = std::abs(heading_accum_) >= 2 * M_PI;
  const bool timed = (t - last_fix_t_) >= cfg_.interval_s;
  const bool too_soon = (t - last_fix_t_) < 5.0;
  if (!(full_rev || timed) || too_soon || window_.size() < 20) return r;

  // Heading coverage from the frames that SURVIVED. heading_accum counts every
  // frame, but a hard turn drops the star count and leaves the survivors
  // clustered while the accumulator still claims a revolution.
  int bins[12] = {0};
  for (const FrameData& f : window_) {
    double hh = f.yaw_est;
    while (hh < 0) hh += 2 * M_PI;
    while (hh >= 2 * M_PI) hh -= 2 * M_PI;
    ++bins[std::min(11, int(hh / (2 * M_PI) * 12))];
  }
  int occupied = 0;
  for (int b : bins) if (b) ++occupied;
  const double hdg = std::abs(heading_accum_) * kRad2Deg;
  const bool dense = window_.size() >= 40 && hdg < 400.0 && occupied >= 9;

  fix.t = t; fix.heading_deg = hdg; fix.frames = window_.size();
  fix.bins = occupied; fix.full_rev = full_rev; fix.dense = dense;

  // ---- HOW GOOD IS THIS FIX? ---------------------------------------------
  //
  // Without a vertical reference the mounting error and the AHRS tilt bias are
  // both BODY-FIXED, and only averaging over a heading sweep removes them. So
  // fix quality is a function of HEADING COVERAGE, and it is steep. Measured,
  // mount calibrated, 8 seeds:
  //
  //     sweep      no horizon    1x Lepton 2.5
  //      36 deg      27.1 km        16.2 km
  //      90 deg      23.1 km        14.4 km
  //     180 deg       9.6 km         9.3 km
  //     270 deg       4.8 km         2.6 km
  //     360 deg       6.2 km         2.7 km
  //
  // Rather than gate hard on coverage, report the uncertainty HONESTLY and let
  // the dead-reckoning filter weight it. A 16 km fix is useless against a DR
  // estimate that is 5 km uncertain and valuable against one that is 40 km
  // uncertain, and the Kalman update already knows which is which -- it just
  // has to be told the truth about sigma.
  //
  // sigma ~ sigma_full * (12/bins)^0.8 fits the table above for both
  // configurations to within about 20%.
  const double frac = std::max(1, occupied) / 12.0;
  const double sigma_full =
      cfg_.fix_sigma_full > 0.0
          ? cfg_.fix_sigma_full
          : (cfg_.horizon.cameras.empty() ? 6000.0 : 2700.0);
  const double fix_sigma_m =
      std::clamp(sigma_full * std::pow(1.0 / frac, 0.8), 1500.0, 60000.0);
  fix.sigma_m = fix_sigma_m;

  // ---- TRIM TO WHOLE REVOLUTIONS -----------------------------------------
  //
  // A body-fixed error only cancels if every heading is sampled EQUALLY, so a
  // fractional sweep leaves the over-sampled sector weighted and the error
  // partly uncancelled. Measured, 40 seeds, no horizon:
  //
  //     1.00 rev   6.18 km      1.25 rev  10.16 km
  //     2.00 rev   3.94 km      1.50 rev   6.79 km
  //     3.00 rev   3.32 km      2.25 rev   6.39 km
  //
  // 1.25 revolutions is 64% WORSE than 1.00 despite covering more sky. Whole
  // revolutions are local minima at every period tested, so when the window
  // holds at least one, keep only the most recent whole number of them.
  //
  // This is free: it discards frames that were actively hurting.
  std::vector<FrameData> solve_window = window_;
  std::vector<Geodetic> solve_truth = window_truth_;
  if (!window_hdg_.empty()) {
    const double total = std::abs(heading_accum_);
    const int whole = int(total / (2 * M_PI));
    if (whole >= 1) {
      const double keep_from =
          heading_accum_ - (heading_accum_ < 0 ? -1 : 1) * whole * 2 * M_PI;
      size_t cut = 0;
      while (cut + 1 < window_hdg_.size() &&
             ((heading_accum_ < 0) ? window_hdg_[cut] > keep_from
                                   : window_hdg_[cut] < keep_from))
        ++cut;
      if (cut > 0 && cut < solve_window.size()) {
        solve_window.erase(solve_window.begin(), solve_window.begin() + cut);
        if (cut < solve_truth.size())
          solve_truth.erase(solve_truth.begin(), solve_truth.begin() + cut);
      }
    }
  }

  // Circle fit BOOTSTRAPS from an unknown mounting; after calibration there is
  // nothing left for it to fit and the heading-weighted mean is better.
  const AverageMethod m = calibrated_ ? AverageMethod::HeadingWeighted
                                      : AverageMethod::CircleFit;
  const OrbitResult res = estimateOrbit(solve_window, C_b_c_assumed_, m);
  const OrbitResult naive =
      estimateOrbit(solve_window, C_b_c_assumed_, AverageMethod::Naive);
  if (!res.ok) { last_fix_t_ = t; return r; }

  double la = 0, lo = 0;
  for (const Geodetic& p : solve_truth) { la += p.lat; lo += p.lon; }
  const Geodetic mean_truth{la / std::max<size_t>(1, solve_truth.size()),
                            lo / std::max<size_t>(1, solve_truth.size()),
                            s.alt_m};

  // ---- consistency gate ---------------------------------------------------
  //
  // A bad fix is not merely useless, it is CORROSIVE: it overwrites `assumed_`,
  // which the matcher and the next solve both depend on, and it can be
  // recalibrated upon. Dead reckoning is the independent check, and the test is
  // MAHALANOBIS rather than a circle because DR uncertainty is anisotropic --
  // along-track error comes from airspeed scale, cross-track from heading.
  double ref_dist = -1.0, allowed = 0.0;
  if (dr_.started()) {
    ref_dist = haversine(res.position, dr_.position());
    const double sigmas = dr_.mahalanobis(res.position);
    if (sigmas > 4.0) {
      fix.rejected = true;
      char buf[128];
      std::snprintf(buf, sizeof buf,
                    "rejected: %.1f sigma from dead reckoning (%.1f km)",
                    sigmas, ref_dist / 1000);
      fix.note = buf;
      last_fix_t_ = t;
      return r;
    }
    allowed = 60000.0;
  } else if (have_assumed_) {
    ref_dist = haversine(res.position, assumed_);
    allowed = 60000.0;   // first fix: only reject the absurd
  }
  if (ref_dist >= 0 && ref_dist > allowed) {
    fix.rejected = true;
    fix.note = "rejected: implausibly far from the prior";
    last_fix_t_ = t;
    return r;
  }

  // A hard floor remains: below ~90 deg of sweep the fix is not merely poor,
  // it is poorly CONDITIONED, and the error is a body-fixed bias rather than
  // noise -- so the filter's Gaussian assumption stops holding and reporting a
  // large sigma no longer describes it.
  if (occupied < cfg_.min_heading_bins) {
    fix.rejected = true;
    char buf[128];
    std::snprintf(buf, sizeof buf,
                  "skipped: only %d of 12 heading bins -- a body-fixed error "
                  "cannot average out without a sweep", occupied);
    fix.note = buf;
    last_fix_t_ = t;
    return r;
  }

  ++n_fix_;
  last_fix_t_ = t;
  fix.emitted = true;
  fix.index = n_fix_;
  fix.position = res.position;
  fix.err_m = haversine(res.position, mean_truth);
  fix.err_naive_m = naive.ok ? haversine(naive.position, mean_truth) : -1;
  assumed_ = res.position;

  // Recalibration must be EARNED: only dense, well-covered windows. A partial
  // arc is fine for a position and poor for a calibration.
  double nstars = 0;
  for (const FrameData& f : window_) nstars += f.id.size();
  const double spf = nstars / window_.size();
  if (window_.size() >= 30 && spf >= 4.0 && occupied >= 9) {
    const Eigen::Matrix3d C_new = recalibrateMount(window_, res.position);
    const OrbitResult before = estimateOrbit(window_, C_b_c_assumed_, m);
    const OrbitResult after = estimateOrbit(window_, C_new, m);
    if (!(after.ok && before.ok && after.spread > before.spread * 1.2)) {
      // RUNNING AVERAGE, not a fixed gain. The true mounting is static, so
      // successive estimates are repeated measurements of ONE fixed quantity
      // and should be averaged with growing confidence. A fixed gain never
      // settles: measured on a 100-minute SITL flight it took the boresight
      // 0.40 -> 0.10 deg and then let it wander back to 0.17, and since fix
      // error tracks mounting error at ~111 km/deg, that wander is worth
      // several km on every subsequent fix.
      //
      // 1/(n+1) is the rotation-average weight; floored so a long flight can
      // still follow a mounting that genuinely shifts.
      ++n_calib_;
      const double w = std::max(1.0 / (n_calib_ + 1), 0.05);
      Eigen::Quaterniond qo(C_b_c_assumed_), qn(C_new);
      if (qo.dot(qn) < 0) qn.coeffs() *= -1.0;
      C_b_c_assumed_ = qo.slerp(w, qn).toRotationMatrix();
      if (!calibrated_) { calibrated_ = true; mcfg_.radius_px = std::min(mcfg_.radius_px, 60.0); }
    } else {
      fix.note = "recalibration reverted: scatter would grow";
    }
  } else {
    fix.note = "recalibration skipped: window too thin";
  }
  fix.boresight_deg =
      std::acos(std::clamp(C_b_c_assumed_.col(2).dot(C_b_c_true_.col(2)), -1.0, 1.0)) *
      kRad2Deg;

  if (cfg_.dead_reckon) {
    if (!dr_.started()) { dr_.start(res.position); dr_open_.start(res.position); }
    // (kept for the case where dead reckoning was enabled after takeoff)
    else {
      fix.dr_drift_m = dr_.sinceFix();
      // Measurement sigma from an OBSERVABLE quantity, never the true error.
      // Honest per-fix sigma from heading coverage, not a proxy from the
      // fit spread -- the spread measures scatter about the fitted circle,
      // which says nothing about the body-fixed bias that dominates.
      dr_.update(res.position, fix_sigma_m);
    }
    fix.dr_err_m = haversine(dr_.position(), mean_truth);
    fix.dr_open_err_m = haversine(dr_open_.position(), mean_truth);
  }

  if (full_rev) {
    window_.clear(); window_truth_.clear(); window_hdg_.clear();
    heading_accum_ = 0;
  }
  const double keep = 300.0;
  size_t drop = 0;
  while (drop < window_.size() && t - window_[drop].t > keep) ++drop;
  if (drop) {
    window_.erase(window_.begin(), window_.begin() + drop);
    if (window_hdg_.size() >= size_t(drop))
      window_hdg_.erase(window_hdg_.begin(), window_hdg_.begin() + drop);
    window_truth_.erase(window_truth_.begin(), window_truth_.begin() + drop);
  }
  return r;
}

}  // namespace celestial
