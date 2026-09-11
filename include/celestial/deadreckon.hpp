// deadreckon.hpp — airspeed + heading dead reckoning, corrected by celestial.
//
// -------------------------------------------------------------------------
// WHY THIS EXISTS
// -------------------------------------------------------------------------
// A celestial fix is bias-dominated at a few km and cannot resolve a 3 km
// loiter -- averaging does not help, because the residual is a common-mode
// tilt error shared by every frame. Dead reckoning has the opposite character:
// the SHAPE is right over minutes, but it drifts without bound.
//
// Fused, you get a trajectory whose shape is right and whose absolute position
// is bounded. That is the actual product: not "reconstruct my flight path",
// but "bound absolute position error indefinitely, without emitting anything".
//
// State (4): position north/east in metres from an origin, and the wind
// vector. Wind is a state rather than a parameter because it is the dominant
// unknown in the velocity, and each celestial fix observes it indirectly --
// a position error that grows linearly in one direction IS a wind error.
//
// Deliberately NOT modelled: attitude, accelerometer bias, scale factors.
// Those belong in the autopilot's EKF. This layer answers one question -- how
// far has the aircraft moved since the last fix -- with the two signals that
// survive in a GNSS-denied aircraft: airspeed and heading.

#pragma once

#include "celestial/types.hpp"

namespace celestial {

/// Sensor error model, chosen to give the drift a cheap installation shows.
///
///   2% airspeed scale error at 25 m/s  -> 0.5 m/s  -> 1.8 km/hour
///   2 deg heading bias at 25 m/s       -> 0.9 m/s  -> 3.1 km/hour
///   20% wind error at 8 m/s            -> 1.6 m/s  -> 5.8 km/hour
///
/// i.e. of order 5-10 km/hour, which is what makes a 2 km celestial fix worth
/// having on a multi-hour flight.
struct DeadReckonConfig {
  double airspeed_scale_err = 0.02;   ///< fractional, systematic
  double airspeed_noise = 0.5;        ///< m/s, 1-sigma
  double heading_bias = 2.0 * kDeg2Rad;
  double heading_noise = 1.0 * kDeg2Rad;
  /// Wind random walk. The wind genuinely changes; this also absorbs the
  /// systematic airspeed and heading errors, which is why the filter can
  /// correct for them at all.
  double wind_walk = 0.05;            ///< m/s per sqrt(s)
  double wind_sigma0 = 5.0;           ///< m/s, initial uncertainty
  unsigned seed = 99;

  /// Named presets. The literature spread is wide and it is driven almost
  /// entirely by how well the WIND is known:
  ///
  ///   OPTIMISTIC  ~2-4% of distance. A calibrated pitot (ARSPD_RATIO
  ///               autocal), a calibrated magnetometer, and a wind estimate
  ///               kept fresh. Comparable to the ~2.2%/distance reported for
  ///               a high-end MEMS AHRS in GNSS-denied dead reckoning.
  ///
  ///   REALISTIC   ~6-8%. Wind known to about 20%, which alone contributes
  ///               6.4% at 8 m/s wind and 25 m/s airspeed -- larger than the
  ///               airspeed and heading errors combined.
  ///
  ///   POOR        ~15-30%. A stale or absent wind estimate. Note EKF3's wind
  ///               estimate normally comes from differencing GPS velocity
  ///               against airspeed, so without GPS it goes stale on its own.
  ///               A navigation patent quotes 125 nm/hr for this case.
  ///
  /// For reference, 1930s naval navigators achieved ~0.8% of distance -- but
  /// with a drift meter, actively re-estimating wind against the sea surface.
  /// That is not open-loop dead reckoning.
  static DeadReckonConfig optimistic();
  static DeadReckonConfig realistic();
  static DeadReckonConfig poor();
};

class DeadReckoner {
 public:
  DeadReckoner() = default;
  explicit DeadReckoner(const DeadReckonConfig& cfg) : cfg_(cfg) {}

  /// Anchor the filter. Call once, with the first celestial fix.
  void start(const Geodetic& origin);
  bool started() const { return started_; }

  /// Propagate. `airspeed` and `heading` are the raw sensor values; the
  /// configured scale error, bias and noise are applied here, so the caller
  /// passes truth and the model corrupts it.
  void predict(double dt, double airspeed, double heading);

  /// Fuse an absolute position fix with the given 1-sigma accuracy.
  void update(const Geodetic& fix, double sigma_m);

  Geodetic position() const;
  /// Position in metres north/east of the anchor. This is what
  /// FrameData::dr_ne wants: a common-origin displacement, so that per-frame
  /// fixes can be transported to a shared epoch before averaging.
  Eigen::Vector2d positionNE() const { return x_.head<2>(); }
  Eigen::Vector2d wind() const { return x_.tail<2>(); }
  /// 1-sigma position uncertainty, metres. Isotropic summary -- it averages
  /// the two diagonal terms, so it discards the SHAPE.
  double sigmaPos() const;

  /// Position covariance, m^2, north/east. The uncertainty is genuinely
  /// anisotropic: along-track error is driven by airspeed scale, cross-track by
  /// heading error, and they are different sizes. A gate that tests a circle
  /// therefore either admits bad fixes across the narrow axis or rejects good
  /// ones along the wide one.
  Eigen::Matrix2d covPos() const;

  /// Mahalanobis distance of a position from the current estimate, in sigmas,
  /// using the full covariance. This is the "cone" test: how many sigmas away
  /// is this fix, measured in the shape the filter actually believes.
  double mahalanobis(const Geodetic& p) const;
  /// Metres of DR travel since the last fix -- how far the estimate has been
  /// running open-loop.
  double sinceFix() const { return since_fix_m_; }

 private:
  DeadReckonConfig cfg_;
  Geodetic origin_;
  Eigen::Vector4d x_ = Eigen::Vector4d::Zero();   // n, e, wind_n, wind_e
  Eigen::Matrix4d P_ = Eigen::Matrix4d::Identity();
  bool started_ = false;
  double since_fix_m_ = 0.0;
  unsigned long step_ = 0;
};

}  // namespace celestial
