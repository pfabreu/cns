// attitude.hpp — rotation utilities in NED / FRD.
//
// Conventions (see types.hpp):
//   Local level = NED (North, East, Down)
//   Body        = FRD (Forward, Right, Down)
//   Camera      = OpenCV style: +x right in image, +y down in image,
//                 +z along the optical axis (boresight)
//
// C_l_b maps a vector from BODY to LOCAL LEVEL:  v_ned = C_l_b * v_body
// C_b_c maps a vector from CAMERA to BODY:       v_body = C_b_c * v_cam
// so the full chain of Eq. (1) is  v_ned = C_l_b * C_b_c * v_cam.

#pragma once

#include <vector>

#include "celestial/types.hpp"

namespace celestial {

/// Yaw-pitch-roll (3-2-1) Euler angles to C_l_b.
/// C_l_b = Rz(yaw) * Ry(pitch) * Rx(roll)
Eigen::Matrix3d eulerToDcm(double roll, double pitch, double yaw);

/// Inverse of eulerToDcm. Returns (roll, pitch, yaw).
Eigen::Vector3d dcmToEuler(const Eigen::Matrix3d& C);

/// Small-angle rotation about an arbitrary axis, used to inject a boresight
/// misalignment of a known magnitude and direction.
Eigen::Matrix3d smallRotation(const Eigen::Vector3d& axis, double angle);

/// Angle of the rotation taking A to B, radians. Handy for reporting how far
/// an estimated C_b_c is from the truth.
double rotationAngleBetween(const Eigen::Matrix3d& A, const Eigen::Matrix3d& B);

/// Convert azimuth/elevation to a NED unit vector.
/// NED down is positive, so an elevation above the horizon gives z < 0.
Eigen::Vector3d azElToNed(double azimuth, double elevation);

/// Elevation of a NED unit vector above the horizon.
///
/// NOTE: the paper's Eq. (22) reads  el = -atan(y / sqrt(x^2 + y^2)),
/// which mixes East into both numerator and denominator. It is a typo:
/// the Down component is z. Implemented correctly here.
double nedToElevation(const Eigen::Vector3d& v);

/// Azimuth of a NED unit vector, radians from North through East.
double nedToAzimuth(const Eigen::Vector3d& v);

/// Kabsch: the rotation R minimising sum |b_i - R a_i|^2 over paired unit
/// vectors. Algorithm 2 of the paper. Both inputs must be the same length.
///
/// Note that the paper's Algorithm 2 says to translate to the centroid; that
/// step belongs to point-cloud registration and is WRONG for direction
/// vectors, which have no translational degree of freedom. Removing the
/// centroid of a set of unit vectors destroys information. We do not do it.
Eigen::Matrix3d kabsch(const std::vector<Eigen::Vector3d>& a,
                       const std::vector<Eigen::Vector3d>& b);

/// Average of a set of rotations, via the dominant eigenvector of the summed
/// quaternion outer products (Markley's method). Used to fuse the per-frame
/// Kabsch solutions into one camera mounting estimate.
Eigen::Matrix3d averageRotation(const std::vector<Eigen::Matrix3d>& rotations);


// ===========================================================================
// STAR-AIDED ATTITUDE -- folded in from star_aided.hpp. Same domain: this
// header is attitude representation, this is attitude ESTIMATION, and the
// filter is built from the rotation helpers above.
// ===========================================================================
// ---------------------------------------------------------------------------
// STAR-AIDED ATTITUDE — killing inertial drift with an absolute reference.
//
// WHY THIS IS NOT THE THING THAT FAILED. Three attempts to remove the AHRS tilt
// error from the FIX failed, all against the same wall: a celestial fix cannot
// separate a tilt error from a position error, because they are one observable
// (1 arcmin = 1 nautical mile). Attempts to estimate the tilt after the fact
// therefore have nothing to work with. See README.md.
//
// This is a different quantity. It does not estimate absolute tilt; it
// estimates GYRO BIAS, and bias is observable from the SEQUENCE of attitudes
// rather than from any single one:
//
//   * star directions in ECEF do not depend on where the observer is, so
//     `kabsch(camera dirs, ECEF star dirs)` is an absolute attitude with NO
//     position in it anywhere;
//   * two such attitudes give the true rotation between their epochs, exactly;
//   * the gyro's integral of the same interval differs from it by the bias.
//
// So drift-as-a-RATE is observable where tilt-as-an-OFFSET is not. The two
// results do not conflict: one is about an offset that is degenerate with
// position, the other about a rate that is not.
//
// WHAT IT BUYS. An unaided gyro's attitude error grows without bound. Aided, it
// is held to the accuracy of the star measurement, so the attitude error stops
// growing and becomes very nearly CONSTANT over an orbit -- and a constant
// body-fixed error is exactly what the orbit average removes. That is the whole
// argument, and it is why this attacks the dominant error term where post-hoc
// estimation could not.
//
// STATUS: THE MECHANISM IS IMPLEMENTED AND TESTED IN SIMULATION ONLY.
// `testStarAidedAttitude` shows unaided drift growing and aided drift bounded,
// against a gyro model this project invented. That proves the filter is correct,
// NOT that it helps on your airframe: the answer depends on the real gyro's
// bias stability and on what ArduPilot's EKF3 already does internally. It is
// here so the ArduPilot experiment has something to run. Do not quote a
// simulated number for it.
//
// PRIOR ART. Sodern's Astradia does this in hardware -- a strapdown daytime
// star tracker coupled to an aviation-grade INS, reporting 1.2 Nm at 1 h and
// 2.4 Nm at 10 h, with the residual budget dominated by accelerometer bias and
// star-tracker-to-IMU alignment once attitude drift is removed. Those are the
// same two terms this project arrived at independently. What is unproven is
// doing it with a hobby-grade IMU.
//
// INTEGRATING WITH ARDUPILOT. EKF3 has no general-purpose attitude measurement
// input: it fuses external-navigation YAW (EK3_SRC*_YAW = 6) and position or
// velocity, but not roll and pitch. So there are two routes, and the cheap one
// is not the interesting one:
//
//   1. Emit yaw as external-nav yaw. Easy, and nearly pointless here --
//      `celestialHeading` already recovers heading to ~0.1 deg.
//   2. Run this filter OUTSIDE EKF3, on `RAW_IMU` / `SCALED_IMU` gyro plus the
//      star attitude, and use its output in place of EKF3's attitude in the
//      fix path. No autopilot changes. This is what the class is shaped for.
//
// A third route -- adding a tilt fusion path inside EKF3 -- is the "correct"
// one and is a much larger job.
// ---------------------------------------------------------------------------


/// Earth rotation rate, rad/s. The gyro senses rotation relative to INERTIAL
/// space; this filter tracks attitude relative to ECEF, so Earth rate has to be
/// removed. Over a 150 s orbit it is 0.6 deg -- far too large to ignore.
inline constexpr double kEarthRate = 7.292115e-5;

class StarAidedAttitude {
 public:
  struct Config {
    /// Gyro bias, 1-sigma rad/s, and its correlation time. A cheap MEMS gyro
    /// in-run bias stability is a few deg/hr; 0.5 deg/hr = 2.4e-6 rad/s.
    double bias_sigma = 0.5 * kDeg2Rad / 3600.0;
    double bias_tau = 3600.0;
    /// Angle random walk, rad/sqrt(s).
    double gyro_arw = 0.005 * kDeg2Rad;
    /// Initial attitude uncertainty, rad.
    double att_sigma0 = 2.0 * kDeg2Rad;
  };

  StarAidedAttitude() : cfg_(Config()) { reset(); }
  explicit StarAidedAttitude(const Config& cfg) : cfg_(cfg) { reset(); }

  void reset();

  /// Seed the attitude from any estimate (the AHRS at startup, say). Body to
  /// ECEF.
  void setAttitude(const Eigen::Matrix3d& C_b_e);

  /// Propagate with a gyro sample: body angular rate relative to INERTIAL, in
  /// body axes, as an IMU reports it.
  void predict(const Eigen::Vector3d& omega_ib_body, double dt);

  /// Fuse an absolute star attitude, body to ECEF, with `sigma` radians of
  /// measurement uncertainty per axis.
  ///
  /// Build it as `R * C_b_c.transpose()`, where `R = kabsch(v_cam, ecef_dirs)`
  /// -- see `recalibrateMount` for the same construction. NO POSITION is
  /// involved, which is the entire point.
  void update(const Eigen::Matrix3d& C_b_e_meas, double sigma);

  /// Current attitude, body to ECEF.
  const Eigen::Matrix3d& attitude() const { return C_b_e_; }
  /// Estimated gyro bias, rad/s, body axes.
  const Eigen::Vector3d& gyroBias() const { return bias_; }
  /// Attitude uncertainty, rad, 1-sigma per axis.
  Eigen::Vector3d attitudeSigma() const;
  bool initialised() const { return initialised_; }

 private:
  Config cfg_;
  Eigen::Matrix3d C_b_e_ = Eigen::Matrix3d::Identity();
  Eigen::Vector3d bias_ = Eigen::Vector3d::Zero();
  /// Error state [dtheta(3), dbias(3)]; the state itself is carried in
  /// `C_b_e_` and `bias_` and the error is folded back in after each update,
  /// which is the standard multiplicative formulation.
  Eigen::Matrix<double, 6, 6> P_;
  bool initialised_ = false;
};

}  // namespace celestial
