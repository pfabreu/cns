#include "celestial/deadreckon.hpp"

#include <cmath>
#include <random>

namespace celestial {
namespace {
constexpr double kEarthR = 6371008.8;

double gauss(unsigned long k, unsigned seed) {
  // Deterministic per-step noise: reproducible runs, no shared generator.
  std::mt19937 g(seed ^ static_cast<unsigned>(k * 2654435761u));
  std::normal_distribution<double> n(0.0, 1.0);
  return n(g);
}
}  // namespace

DeadReckonConfig DeadReckonConfig::optimistic() {
  DeadReckonConfig c;
  c.airspeed_scale_err = 0.01;
  c.heading_bias = 1.0 * kDeg2Rad;
  c.wind_walk = 0.03;
  c.wind_sigma0 = 3.0;
  return c;
}

DeadReckonConfig DeadReckonConfig::realistic() {
  DeadReckonConfig c;             // the defaults
  c.airspeed_scale_err = 0.03;
  c.heading_bias = 3.0 * kDeg2Rad;
  c.wind_walk = 0.08;
  c.wind_sigma0 = 6.0;
  return c;
}

DeadReckonConfig DeadReckonConfig::poor() {
  DeadReckonConfig c;
  c.airspeed_scale_err = 0.05;
  c.airspeed_noise = 1.0;
  c.heading_bias = 6.0 * kDeg2Rad;
  c.heading_noise = 2.0 * kDeg2Rad;
  c.wind_walk = 0.20;             // wind estimate going stale
  c.wind_sigma0 = 10.0;
  return c;
}

void DeadReckoner::start(const Geodetic& origin) {
  origin_ = origin;
  x_.setZero();
  P_.setZero();
  P_(0, 0) = P_(1, 1) = 4000.0 * 4000.0;                 // 4 km, one fix
  P_(2, 2) = P_(3, 3) = cfg_.wind_sigma0 * cfg_.wind_sigma0;
  started_ = true;
  since_fix_m_ = 0.0;
}

void DeadReckoner::predict(double dt, double airspeed, double heading) {
  if (!started_ || dt <= 0 || dt > 5.0) return;
  ++step_;

  // Corrupt the sensors: a systematic scale error and bias, plus noise. This
  // is what makes the drift realistic -- a perfect airspeed and heading would
  // dead-reckon indefinitely and there would be nothing for the celestial fix
  // to correct.
  const double va = airspeed * (1.0 + cfg_.airspeed_scale_err) +
                    cfg_.airspeed_noise * gauss(step_, cfg_.seed);
  const double psi = heading + cfg_.heading_bias +
                     cfg_.heading_noise * gauss(step_, cfg_.seed ^ 0xABCDu);

  // Ground velocity = air velocity + wind.
  const double vn = va * std::cos(psi) + x_(2);
  const double ve = va * std::sin(psi) + x_(3);
  x_(0) += vn * dt;
  x_(1) += ve * dt;
  since_fix_m_ += std::hypot(vn, ve) * dt;

  // F is identity plus dt coupling wind into position.
  Eigen::Matrix4d F = Eigen::Matrix4d::Identity();
  F(0, 2) = dt;
  F(1, 3) = dt;

  // Process noise. The velocity uncertainty is dominated by the airspeed and
  // heading errors, which are SYSTEMATIC -- so position uncertainty grows
  // roughly linearly with time, not as sqrt(t). Modelling it as white noise
  // would make the filter overconfident between fixes.
  const double sv_a = cfg_.airspeed_noise;
  const double sv_h = airspeed * cfg_.heading_noise;
  const double sv2 = sv_a * sv_a + sv_h * sv_h;
  Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();
  Q(0, 0) = Q(1, 1) = sv2 * dt * dt;
  Q(2, 2) = Q(3, 3) = cfg_.wind_walk * cfg_.wind_walk * dt;

  P_ = F * P_ * F.transpose() + Q;
}

void DeadReckoner::update(const Geodetic& fix, double sigma_m) {
  if (!started_) { start(fix); return; }

  const double dn = (fix.lat - origin_.lat) * kEarthR;
  const double de = (fix.lon - origin_.lon) * kEarthR * std::cos(origin_.lat);
  const Eigen::Vector2d z(dn, de);

  Eigen::Matrix<double, 2, 4> H = Eigen::Matrix<double, 2, 4>::Zero();
  H(0, 0) = 1.0;
  H(1, 1) = 1.0;
  const Eigen::Matrix2d R =
      Eigen::Matrix2d::Identity() * (sigma_m * sigma_m);

  const Eigen::Vector2d y = z - H * x_;
  const Eigen::Matrix2d S = H * P_ * H.transpose() + R;
  const Eigen::Matrix<double, 4, 2> K = P_ * H.transpose() * S.inverse();
  x_ += K * y;
  P_ = (Eigen::Matrix4d::Identity() - K * H) * P_;
  since_fix_m_ = 0.0;
}

Geodetic DeadReckoner::position() const {
  Geodetic p = origin_;
  p.lat += x_(0) / kEarthR;
  p.lon += x_(1) / (kEarthR * std::cos(origin_.lat));
  return p;
}

Eigen::Matrix2d DeadReckoner::covPos() const {
  return P_.topLeftCorner<2, 2>();
}

double DeadReckoner::mahalanobis(const Geodetic& p) const {
  if (!started_) return 0.0;
  // Offset in the same north/east frame the filter's state uses.
  const double dn = (p.lat - origin_.lat) * kEarthR - x_(0);
  const double de = (p.lon - origin_.lon) * kEarthR *
                        std::cos(origin_.lat) - x_(1);
  const Eigen::Vector2d d(dn, de);
  Eigen::Matrix2d C = covPos();
  // Floor the eigenvalues: a filter that has just been updated can be
  // arbitrarily confident, and an over-confident gate rejects everything.
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(C);
  Eigen::Vector2d ev = es.eigenvalues();
  const double floor_m = 1500.0;
  ev(0) = std::max(ev(0), floor_m * floor_m);
  ev(1) = std::max(ev(1), floor_m * floor_m);
  C = es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();
  const double m2 = (d.transpose() * C.inverse() * d)(0, 0);
  return std::sqrt(std::max(0.0, m2));
}

double DeadReckoner::sigmaPos() const {
  return std::sqrt(std::max(0.0, 0.5 * (P_(0, 0) + P_(1, 1))));
}

}  // namespace celestial
