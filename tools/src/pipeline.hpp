// pipeline.hpp — the NAVIGATION half of the live node.
//
// Follows the system schematic exactly, and knows nothing about MAVLink:
//
//     render (true pose)          <- simulator standing in for the sky
//       -> detect and centroid
//       -> match to catalogue
//       -> per-frame fix          <- absolute, no prior, RANSAC
//       -> orbit fix              <- departure loiter only, calibrates mounting
//       -> dead reckoning         <- consumes fixes, produces the trajectory
//
// THE TRUTH BOUNDARY. `MavSample::true_pos` and the true attitude are used for
// exactly one thing: rendering the star image. In an aircraft that comes from a
// camera pointed at the actual sky. The navigation path sees only pixels, the
// AHRS attitude, and its own previous estimate. Move this to hardware and you
// delete the one render call.

#pragma once

#include <string>
#include <vector>

#include "celestial/deadreckon.hpp"
#include "celestial/imaging.hpp"
#include "celestial/horizon.hpp"
#include "celestial/orbit.hpp"
#include "mavlink_source.hpp"

namespace celestial {

struct PipelineConfig {
  double render_hz = 2.0;
  double mag_limit = 5.0;
  double match_mag = 0.0;        ///< 0 = one magnitude brighter than mag_limit
  double exposure_s = 0.0;       ///< 0 = auto from body rates
  double target_smear_px = 8.0;
  /// FLOOR, and a load-bearing one. The auto-exposure above holds smear at
  /// target_smear_px by SHORTENING exposure, so in a turn it trades photons
  /// for sharpness without limit. Measured in SITL at a 250 m loiter
  /// (5.7 deg/s): the controller wound 200 ms down to 35 ms and detections
  /// went 43 -> 0, which cost every fix for the rest of the flight.
  ///
  /// Measured, transit_scenario --imaging, 144 runs, 10-20 seeds per cell,
  /// pooled over 150/250/400 m fix radius:
  ///
  ///     exposure   matched/frame   frames unusable   fix error
  ///        20 ms        4.2            6.03 %        13.95 km
  ///        50 ms       12.0            0.00 %         6.73 km
  ///       100 ms       24.8            0.07 %         6.61 km
  ///       200 ms       44.8            0.00 %         6.84 km
  ///
  /// There is a CLIFF below ~50 ms and a plateau above it: 50/100/200 ms are
  /// indistinguishable (spread 0.23 km against sd 2.3, SE 0.42), because past
  /// that point attitude is the constraint, not the star pipeline. So this is
  /// not tuned for the best number, it is set to stay off the cliff -- 100 ms
  /// keeps 3x margin and still sits mid-plateau.
  ///
  /// A floor rather than a bigger target_smear_px on purpose: the smear
  /// setpoint has to be re-derived for every turn rate (8 px means 35 ms at
  /// 5.7 deg/s), the floor holds regardless. Smear is NOT the thing to
  /// minimise -- 200 ms runs at 12 px of smear and has the most matched stars
  /// of any cell, because the matched filter takes the streak and photons are
  /// what is scarce.
  double min_exposure_s = 0.100;
  double max_exposure_s = 0.20;
  /// Override for the full-sweep fix sigma, metres. 0 keeps the built-in
  /// 6000 (no horizon) / 2700 (horizon).
  ///
  /// WHY IT IS A KNOB. Those constants encode a MEASURED accuracy ratio, and
  /// the dead-reckoning Kalman update weights fixes by it. If the ratio is
  /// wrong the filter mis-prices every fix. Measured on the Sagres->Porto
  /// Santo crossing: horizon fixes were 9.65 km median against 11.05 km
  /// without, a 1.15x edge, but they are priced at 6000/2700 = 2.22x, i.e.
  /// 4.9x in variance against an earned 1.3x. The horizon arm finished WORSE
  /// on the filtered track (9.05 vs 7.96 km median) despite better fixes --
  /// paying for precision it did not receive.
  double fix_sigma_full = 0.0;

  double match_radius_px = 120.0;
  int min_stars = 5;
  double boresight_deg = 0.4;    ///< the TRUE mounting error, unknown to nav
  double interval_s = 60.0;
  bool dead_reckon = true;
  std::string dr_quality = "realistic";
  std::string utc = "2024-09-15T10:30:00";

  /// Optional: dump the rendered frame, its detections and their matched
  /// catalogue predictions to this directory, for the live viewer. Empty
  /// disables, which is the default -- it costs a PGM write per dumped frame.
  ///
  /// Every Nth frame only; at 10 Hz writing all of them is 47 MB/s.
  /// Minimum heading bins (of 12) covered by the fix window before a fix is
  /// emitted at all. Without a vertical reference the dominant errors are
  /// BODY-FIXED, and only a heading sweep removes them -- a straight-leg fix
  /// carries the lot. 8 of 12 is 240 degrees, enough for the average to bite
  /// while still tolerating a partial arc.
  /// Below this many of 12 heading bins the fix is not merely poor but poorly
  /// CONDITIONED -- the error becomes a body-fixed bias rather than noise, so
  /// reporting a large sigma no longer describes it honestly. Above it, quality
  /// is reported through `FixReport::sigma_m` and the filter decides.
  int min_heading_bins = 4;

  // --- horizon sensor (optional) -------------------------------------------
  //
  // An LWIR camera looking forward gives a NON-INERTIAL vertical reference: it
  // observes the AHRS tilt error directly and is not confused by acceleration,
  // which is what an accelerometer cannot do in a turn.
  //
  // Measured, one revolution at 250 m, mount calibrated, 8 seeds:
  //
  //                          orbit fix    straight leg
  //   no horizon               6.16 km        30.87 km
  //   1x Lepton 2.5            2.46 km        16.19 km
  //   1x Boson 320             0.59 km         4.46 km
  //
  // A single cheap Lepton is worth 2.5x on the ORBIT fix. It does NOT make
  // straight-leg fixes viable -- 16 km is worse than simply orbiting -- so the
  // heading-coverage gate stays. Only a Boson-class sensor changes the flight
  // plan.
  //
  // Empty `cameras` disables it, which is the default.
  HorizonConfig horizon;
  /// Assumed 1-sigma AHRS tilt uncertainty for the Bayesian fusion.
  double imu_tilt_sigma = 0.3 * kDeg2Rad;

  std::string frames_dir;
  int frames_every = 20;
};

/// What one frame produced, for logging and display.
struct FrameReport {
  bool rendered = false;
  int detected = 0, matched = 0;
  double exposure_s = 0, smear_px = 0, heading_deg = 0;
  size_t window = 0;
  Eigen::Vector2d tilt_imu = Eigen::Vector2d::Zero();   ///< arcmin
  Geodetic frame_fix;
  double frame_err_m = -1;
  /// Dead-reckoned position at this frame, once the filter has started. This
  /// is what makes a trajectory plot possible: the fixes are sparse, the DR
  /// track is continuous, and the interesting behaviour is the second being
  /// pulled back by the first.
  Geodetic dr_pos;
  bool dr_valid = false;
  /// Dead reckoning with NO fixes applied -- the open-loop control. Plotting it
  /// alongside the corrected track is what shows the bounding actually working,
  /// because the two curves separate on every leg.
  Geodetic dr_open_pos;
  bool dr_open_valid = false;
  /// True when a horizon measurement was fused into this frame's attitude.
  bool horizon_used = false;
  /// Set when this frame was dumped to `frames_dir`; the viewer polls for it.
  std::string dumped;
};

/// What a completed fix produced. Emitted on a heading revolution or the timer.
struct FixReport {
  bool emitted = false, rejected = false, full_rev = false, dense = false;
  int index = 0, bins = 0;
  double t = 0, heading_deg = 0, err_m = 0, err_naive_m = 0;
  /// Reported 1-sigma uncertainty of this fix, from heading coverage.
  double sigma_m = -1;
  double boresight_deg = 0, dr_err_m = -1, dr_open_err_m = -1, dr_drift_m = 0;
  size_t frames = 0;
  Geodetic position;
  std::string note;
};

class Pipeline {
 public:
  explicit Pipeline(const PipelineConfig& c);
  /// Feed one MAVLink sample. Returns what this frame produced.
  FrameReport step(const MavSample& s, FixReport& fix);
  bool calibrated() const { return calibrated_; }

 private:
  PipelineConfig cfg_;
  Camera cam_;
  SensorModel sensor_;
  DetectorConfig dcfg_;
  MatcherConfig mcfg_;
  DeadReckoner dr_, dr_open_;

  std::vector<CatalogStar> catalog_, match_catalog_;
  std::vector<double> vmags_;
  Epoch epoch0_;
  Eigen::Matrix3d C_b_c_true_, C_b_c_assumed_;

  std::vector<FrameData> window_;
  std::vector<Geodetic> window_truth_;
  /// Cumulative heading at each windowed frame, for whole-revolution trimming.
  std::vector<double> window_hdg_;
  FrameTruth prev_;
  Eigen::Matrix3d prev_C_ = Eigen::Matrix3d::Identity();
  double prev_yaw_ = 0, prev_t_ = 0, last_fix_t_ = -1;
  int dump_count_ = 0, dump_index_ = 0, n_calib_ = 0;
  HorizonState hstate_;
  bool dump_ready_ = false;
  double smear_rate_ = 0, heading_accum_ = 0, prev_yaw_acc_ = 0;
  Geodetic assumed_;
  bool have_prev_ = false, have_assumed_ = false, yaw_init_ = false;
  bool calibrated_ = false;
  int n_fix_ = 0;
};

}  // namespace celestial
