#include "celestial/attitude.hpp"

#include <Eigen/Geometry>
#include <cmath>

namespace celestial {

Eigen::Matrix3d eulerToDcm(double roll, double pitch, double yaw) {
  return (Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ()) *
          Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY()) *
          Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX()))
      .toRotationMatrix();
}

Eigen::Vector3d dcmToEuler(const Eigen::Matrix3d& C) {
  const double pitch = -std::asin(std::clamp(C(2, 0), -1.0, 1.0));
  const double roll = std::atan2(C(2, 1), C(2, 2));
  const double yaw = std::atan2(C(1, 0), C(0, 0));
  return Eigen::Vector3d(roll, pitch, yaw);
}

Eigen::Matrix3d smallRotation(const Eigen::Vector3d& axis, double angle) {
  return Eigen::AngleAxisd(angle, axis.normalized()).toRotationMatrix();
}

double rotationAngleBetween(const Eigen::Matrix3d& A, const Eigen::Matrix3d& B) {
  const Eigen::Matrix3d D = A.transpose() * B;
  const double c = std::clamp((D.trace() - 1.0) / 2.0, -1.0, 1.0);
  return std::acos(c);
}

Eigen::Vector3d azElToNed(double azimuth, double elevation) {
  const double ce = std::cos(elevation);
  return Eigen::Vector3d(ce * std::cos(azimuth), ce * std::sin(azimuth),
                         -std::sin(elevation));
}

double nedToElevation(const Eigen::Vector3d& v) {
  return -std::atan2(v.z(), std::hypot(v.x(), v.y()));
}

double nedToAzimuth(const Eigen::Vector3d& v) {
  return std::atan2(v.y(), v.x());
}

Eigen::Matrix3d kabsch(const std::vector<Eigen::Vector3d>& a,
                       const std::vector<Eigen::Vector3d>& b) {
  if (a.size() != b.size() || a.size() < 2) return Eigen::Matrix3d::Identity();

  // Cross-covariance. No centroid removal: these are directions, not points.
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  for (size_t i = 0; i < a.size(); ++i) H += b[i] * a[i].transpose();

  Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU |
                                               Eigen::ComputeFullV);
  Eigen::Matrix3d S = Eigen::Matrix3d::Identity();
  // Guard against a reflection: force det(R) = +1.
  if ((svd.matrixU() * svd.matrixV().transpose()).determinant() < 0) {
    S(2, 2) = -1.0;
  }
  return svd.matrixU() * S * svd.matrixV().transpose();
}

Eigen::Matrix3d averageRotation(const std::vector<Eigen::Matrix3d>& rotations) {
  if (rotations.empty()) return Eigen::Matrix3d::Identity();
  if (rotations.size() == 1) return rotations.front();

  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  Eigen::Quaterniond ref(rotations.front());
  for (const Eigen::Matrix3d& R : rotations) {
    Eigen::Quaterniond q(R);
    // Quaternions double-cover SO(3); align hemispheres before accumulating.
    if (q.dot(ref) < 0) q.coeffs() *= -1.0;
    const Eigen::Vector4d v = q.coeffs();
    M += v * v.transpose();
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> es(M);
  const Eigen::Vector4d best = es.eigenvectors().col(3);  // largest eigenvalue
  Eigen::Quaterniond q;
  q.coeffs() = best;
  q.normalize();
  return q.toRotationMatrix();
}


// ===========================================================================
// STAR-AIDED ATTITUDE -- folded in from star_aided.hpp. Same domain: this
// header is attitude representation, this is attitude ESTIMATION, and the
// filter is built from the rotation helpers above.
// ===========================================================================

namespace {

/// Rotation from a rotation vector, small-angle safe.
Eigen::Matrix3d expmap(const Eigen::Vector3d& r) {
  const double a = r.norm();
  if (a < 1e-12) return Eigen::Matrix3d::Identity();
  return Eigen::AngleAxisd(a, r / a).toRotationMatrix();
}

/// Rotation vector from a rotation.
Eigen::Vector3d logmap(const Eigen::Matrix3d& C) {
  const Eigen::AngleAxisd aa(C);
  return aa.axis() * aa.angle();
}

}  // namespace

void StarAidedAttitude::reset() {
  C_b_e_.setIdentity();
  bias_.setZero();
  P_.setZero();
  const double a2 = cfg_.att_sigma0 * cfg_.att_sigma0;
  const double b2 = cfg_.bias_sigma * cfg_.bias_sigma;
  for (int i = 0; i < 3; ++i) {
    P_(i, i) = a2;
    P_(3 + i, 3 + i) = b2;
  }
  initialised_ = false;
}

void StarAidedAttitude::setAttitude(const Eigen::Matrix3d& C_b_e) {
  C_b_e_ = C_b_e;
  initialised_ = true;
}

void StarAidedAttitude::predict(const Eigen::Vector3d& omega_ib_body,
                                double dt) {
  if (dt <= 0 || dt > 10.0) return;

  // The gyro senses rotation relative to INERTIAL. Attitude here is relative to
  // ECEF, so subtract Earth rate expressed in body axes. Over a 150 s orbit
  // Earth rate is 0.6 deg -- dropping it would swamp the bias being estimated.
  const Eigen::Vector3d w_ie_ecef(0.0, 0.0, kEarthRate);
  const Eigen::Vector3d w_eb_body =
      omega_ib_body - bias_ - C_b_e_.transpose() * w_ie_ecef;

  C_b_e_ = C_b_e_ * expmap(w_eb_body * dt);

  // Error-state transition. dtheta grows with the bias error; the bias itself
  // is a slow Gauss-Markov process.
  Eigen::Matrix<double, 6, 6> F = Eigen::Matrix<double, 6, 6>::Identity();
  F.block<3, 3>(0, 0) = expmap(w_eb_body * dt).transpose();
  F.block<3, 3>(0, 3) = -Eigen::Matrix3d::Identity() * dt;
  const double a_b = std::exp(-dt / std::max(1e-6, cfg_.bias_tau));
  F.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity() * a_b;

  Eigen::Matrix<double, 6, 6> Q = Eigen::Matrix<double, 6, 6>::Zero();
  const double q_att = cfg_.gyro_arw * cfg_.gyro_arw * dt;
  const double q_bias =
      cfg_.bias_sigma * cfg_.bias_sigma * (1.0 - a_b * a_b);
  for (int i = 0; i < 3; ++i) {
    Q(i, i) = q_att;
    Q(3 + i, 3 + i) = q_bias;
  }

  P_ = F * P_ * F.transpose() + Q;
}

void StarAidedAttitude::update(const Eigen::Matrix3d& C_b_e_meas,
                               double sigma) {
  if (!initialised_) {
    setAttitude(C_b_e_meas);
    for (int i = 0; i < 3; ++i) P_(i, i) = sigma * sigma;
    return;
  }

  // Residual: the small rotation taking the estimate onto the measurement,
  // expressed as a rotation vector in body axes.
  const Eigen::Vector3d y = logmap(C_b_e_.transpose() * C_b_e_meas);

  Eigen::Matrix<double, 3, 6> H = Eigen::Matrix<double, 3, 6>::Zero();
  H.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();

  const Eigen::Matrix3d R = Eigen::Matrix3d::Identity() * (sigma * sigma);
  const Eigen::Matrix3d S = H * P_ * H.transpose() + R;
  const Eigen::Matrix<double, 6, 3> K = P_ * H.transpose() * S.inverse();
  const Eigen::Matrix<double, 6, 1> dx = K * y;

  // Multiplicative correction: fold the attitude error into the quaternion /
  // DCM rather than carrying it in the state.
  C_b_e_ = C_b_e_ * expmap(dx.head<3>());
  bias_ += dx.tail<3>();

  const Eigen::Matrix<double, 6, 6> I =
      Eigen::Matrix<double, 6, 6>::Identity();
  P_ = (I - K * H) * P_;
  P_ = 0.5 * (P_ + P_.transpose());  // keep it symmetric
}

Eigen::Vector3d StarAidedAttitude::attitudeSigma() const {
  return Eigen::Vector3d(std::sqrt(std::max(0.0, P_(0, 0))),
                         std::sqrt(std::max(0.0, P_(1, 1))),
                         std::sqrt(std::max(0.0, P_(2, 2))));
}

}  // namespace celestial
