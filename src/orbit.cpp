#include "celestial/orbit.hpp"

#include "celestial/attitude.hpp"
#include "celestial/sky_model.hpp"

#include <algorithm>
#include <cmath>
#include <random>

#include "celestial/star_catalog.hpp"

namespace celestial {
namespace {

constexpr double kG = 9.80665;
constexpr double kEarthR = 6371008.8;  // WGS-84 mean radius

/// Offset a geodetic position by a local NED displacement, small-angle.
Geodetic offsetNed(const Geodetic& p, double dn, double de) {
  Geodetic q = p;
  q.lat += dn / (kEarthR + p.alt);
  q.lon += de / ((kEarthR + p.alt) * std::cos(p.lat));
  return q;
}

/// Ground speed achievable on a given track with a given wind.
/// Solving |V_g * u - w| = V_a for V_g gives
///     V_g = (u.w) + sqrt(V_a^2 - |w|^2 + (u.w)^2)
double groundSpeedOnTrack(double track, double airspeed, double wn, double we) {
  const double ux = std::cos(track), uy = std::sin(track);
  const double udotw = ux * wn + uy * we;
  const double w2 = wn * wn + we * we;
  const double disc = airspeed * airspeed - w2 + udotw * udotw;
  if (disc <= 0.0) return 0.0;  // wind too strong to hold this track
  return udotw + std::sqrt(disc);
}

}  // namespace

/// Interpret a zenith unit vector as a geodetic position. See
/// position_solver.hpp: the latitude is already geodetic.
Geodetic fromZenith(const Eigen::Vector3d& z, double alt) {
  const Eigen::Vector3d n = z.normalized();
  return Geodetic{std::asin(std::clamp(n.z(), -1.0, 1.0)),
                  std::atan2(n.y(), n.x()), alt};
}

// ---------------------------------------------------------------------------

bool Camera::inFov(const Eigen::Vector3d& v_cam) const {
  if (v_cam.z() <= 1e-6) return false;  // behind the sensor
  const double f = focalPx();
  const double u = f * v_cam.x() / v_cam.z() + width / 2.0;
  const double v = f * v_cam.y() / v_cam.z() + height / 2.0;
  return u >= 0.0 && u < width && v >= 0.0 && v < height;
}

Eigen::Matrix3d nominalCameraMount() {
  Eigen::Matrix3d C;
  // Columns are the camera axes expressed in body FRD.
  C.col(0) = Eigen::Vector3d(0, 1, 0);   // cam +x -> starboard
  C.col(1) = Eigen::Vector3d(1, 0, 0);   // cam +y -> forward
  C.col(2) = Eigen::Vector3d(0, 0, -1);  // cam +z (boresight) -> up
  return C;
}

// ---------------------------------------------------------------------------

std::vector<FrameTruth> generateOrbit(const OrbitConfig& cfg) {
  std::vector<FrameTruth> out;
  auto epoch0 = Epoch::fromUtc(cfg.year, cfg.month, cfg.day, cfg.hour,
                               cfg.minute, cfg.second);
  if (!epoch0) return out;

  const double dt = 1.0 / cfg.frame_rate;
  const double yaw_rate = cfg.airspeed / cfg.radius_m;
  const double bank = std::atan2(cfg.airspeed * yaw_rate, kG);

  if (cfg.gps_guided_track) {
    // Circular GROUND track. Track angle theta advances at V_g(track)/r, so
    // upwind arcs take longer and produce more frames per unit heading.
    double theta = cfg.start_track;
    const double theta_end = cfg.start_track + 2.0 * M_PI * cfg.revolutions;
    double t = 0.0;
    while (theta < theta_end) {
      const double heading = theta + M_PI_2;
      // Holding a fixed GROUND radius in wind means the turn rate in the NED
      // frame varies with ground speed, so the autopilot must vary bank to
      // suit. This is the mechanism behind Section 4.4: it is not just that
      // heading is sampled non-uniformly in time, it is that the aircraft
      // attitude itself is modulated at exactly the orbit frequency.
      const double vg_here =
          groundSpeedOnTrack(heading, cfg.airspeed, cfg.wind_n, cfg.wind_e);
      const double psi_dot = vg_here / cfg.radius_m;
      FrameTruth f;
      f.t = t;
      f.pos = offsetNed(cfg.centre, cfg.radius_m * std::cos(theta),
                        cfg.radius_m * std::sin(theta));
      f.pos.alt = cfg.centre.alt;
      f.roll = std::atan2(cfg.airspeed * psi_dot, kG);
      f.pitch = 0.0;
      f.yaw = heading;
      f.epoch = epoch0->advanced(t);
      out.push_back(f);

      const double vg =
          groundSpeedOnTrack(heading, cfg.airspeed, cfg.wind_n, cfg.wind_e);
      if (vg <= 0.1) break;
      theta += (vg / cfg.radius_m) * dt;
      t += dt;
    }
  } else {
    // Fixed attitude: constant bank, constant yaw rate. Heading is uniform in
    // time; position drifts downwind.
    double heading = cfg.start_track + M_PI_2;
    Geodetic pos = offsetNed(cfg.centre, cfg.radius_m * std::cos(cfg.start_track),
                             cfg.radius_m * std::sin(cfg.start_track));
    pos.alt = cfg.centre.alt;
    const int n = static_cast<int>(2.0 * M_PI * cfg.revolutions / (yaw_rate * dt));
    for (int k = 0; k < n; ++k) {
      const double t = k * dt;
      FrameTruth f;
      f.t = t;
      f.pos = pos;
      f.roll = bank;
      f.pitch = 0.0;
      f.yaw = heading;
      f.epoch = epoch0->advanced(t);
      out.push_back(f);

      // Ground velocity = air velocity + wind.
      const double vn = cfg.airspeed * std::cos(heading) + cfg.wind_n;
      const double ve = cfg.airspeed * std::sin(heading) + cfg.wind_e;
      pos = offsetNed(pos, vn * dt, ve * dt);
      pos.alt = cfg.centre.alt;
      heading += yaw_rate * dt;
    }
  }
  return out;
}

// ---------------------------------------------------------------------------

std::vector<Eigen::Matrix3d> simulateAhrs(const std::vector<FrameTruth>& truth,
                                          const ErrorModel& err) {
  std::mt19937 rng(err.seed ^ 0x9e3779b9u);
  std::normal_distribution<double> n(0.0, std::max(err.ahrs_noise, 1e-12));
  std::normal_distribution<double> u(0.0, 1.0);
  const bool jitter = err.ahrs_noise > 0.0;
  const bool drift = err.ahrs_drift_sigma > 0.0 && err.ahrs_drift_tau > 0.0;

  std::vector<Eigen::Matrix3d> est;
  est.reserve(truth.size());
  // Gauss-Markov tilt drift, roll and pitch. Started at a steady-state draw
  // rather than zero, so the run does not begin with an artificially good
  // vertical that then decays.
  double d_roll = drift ? err.ahrs_drift_sigma * u(rng) : 0.0;
  double d_pitch = drift ? err.ahrs_drift_sigma * u(rng) : 0.0;
  double t_prev = truth.empty() ? 0.0 : truth.front().t;

  for (const FrameTruth& ft : truth) {
    if (drift) {
      const double dt = std::max(0.0, ft.t - t_prev);
      const double a = std::exp(-dt / err.ahrs_drift_tau);
      const double q = err.ahrs_drift_sigma * std::sqrt(std::max(0.0, 1 - a * a));
      d_roll = a * d_roll + q * u(rng);
      d_pitch = a * d_pitch + q * u(rng);
    }
    t_prev = ft.t;
    // Coordinated-turn leakage: the accelerometer's vertical is pulled toward
    // the aircraft's own down axis, so the estimate under-reads the bank.
    const double turn = -err.ahrs_turn_coupling * ft.roll;
    est.push_back(eulerToDcm(
        ft.roll + err.ahrs_bias.x() + d_roll + turn + (jitter ? n(rng) : 0.0),
        ft.pitch + err.ahrs_bias.y() + d_pitch + (jitter ? n(rng) : 0.0),
        ft.yaw + err.ahrs_bias.z() + (jitter ? n(rng) : 0.0)));
  }
  return est;
}

std::vector<Eigen::Vector3d> simulateGyro(const std::vector<FrameTruth>& truth,
                                          const ErrorModel& err) {
  std::vector<Eigen::Vector3d> out;
  if (truth.size() < 2) return out;
  out.reserve(truth.size());

  std::mt19937 rng(err.seed ^ 0x5bf03635u);
  std::normal_distribution<double> u(0.0, 1.0);

  Eigen::Vector3d bias(err.gyro_bias_sigma * u(rng),
                       err.gyro_bias_sigma * u(rng),
                       err.gyro_bias_sigma * u(rng));

  const Eigen::Vector3d w_ie_ecef(0.0, 0.0, kEarthRate);

  for (size_t i = 0; i < truth.size(); ++i) {
    const size_t j = (i + 1 < truth.size()) ? i + 1 : i;
    const size_t k = (i + 1 < truth.size()) ? i : i - 1;
    const double dt = std::max(1e-6, truth[j].t - truth[k].t);

    const Eigen::Matrix3d Ca = eulerToDcm(truth[k].roll, truth[k].pitch,
                                          truth[k].yaw);
    const Eigen::Matrix3d Cb = eulerToDcm(truth[j].roll, truth[j].pitch,
                                          truth[j].yaw);
    // Body rate relative to NED, from the truth attitude sequence.
    const Eigen::AngleAxisd aa(Ca.transpose() * Cb);
    Eigen::Vector3d w_nb_body = aa.axis() * aa.angle() / dt;

    // A gyro senses rotation relative to INERTIAL, so add Earth rate. The
    // transport rate (NED rotating as the aircraft moves) is ~1e-6 rad/s at
    // 25 m/s and is neglected.
    const Eigen::Matrix3d C_ecef_ned = nedBasisEcef(truth[i].pos);
    const Eigen::Matrix3d C_b_e = C_ecef_ned * eulerToDcm(truth[i].roll,
                                                          truth[i].pitch,
                                                          truth[i].yaw);
    const Eigen::Vector3d w_ib_body =
        w_nb_body + C_b_e.transpose() * w_ie_ecef;

    // Gauss-Markov bias walk plus white noise.
    const double a = std::exp(-dt / std::max(1e-6, err.gyro_bias_tau));
    const double q =
        err.gyro_bias_sigma * std::sqrt(std::max(0.0, 1 - a * a));
    for (int c = 0; c < 3; ++c) bias(c) = a * bias(c) + q * u(rng);

    const double wn = err.gyro_arw / std::sqrt(dt);
    out.push_back(w_ib_body + bias +
                  Eigen::Vector3d(wn * u(rng), wn * u(rng), wn * u(rng)));
  }
  return out;
}

std::vector<FrameData> simulateObservations(const std::vector<FrameTruth>& truth,
                                            const Camera& cam,
                                            const ErrorModel& err) {
  return simulateObservationsWithAttitude(truth, simulateAhrs(truth, err), cam,
                                          err);
}

std::vector<FrameData> simulateObservationsWithAttitude(
    const std::vector<FrameTruth>& truth,
    const std::vector<Eigen::Matrix3d>& C_l_b_est, const Camera& cam,
    const ErrorModel& err) {
  std::vector<FrameData> out;
  if (truth.empty() || C_l_b_est.size() != truth.size()) return out;
  out.reserve(truth.size());

  const auto catalog = catalogBrighterThan(err.mag_limit);
  const StarField field(catalog, truth[truth.size() / 2].epoch);

  const Eigen::Matrix3d C_b_c_nominal = nominalCameraMount();
  // The TRUE mounting differs from nominal by the boresight error.
  const Eigen::Matrix3d C_b_c_true =
      C_b_c_nominal * smallRotation(err.boresight_axis, err.boresight_error);

  std::mt19937 rng(err.seed);
  std::normal_distribution<double> cent_n(0.0, std::max(err.centroid_noise, 1e-12));

  std::vector<Eigen::Vector3d> ecef;
  for (size_t fi = 0; fi < truth.size(); ++fi) {
    const FrameTruth& ft = truth[fi];
    const Eigen::Matrix3d C_l_b_true = eulerToDcm(ft.roll, ft.pitch, ft.yaw);
    const Eigen::Matrix3d C_ecef_ned = nedBasisEcef(ft.pos);

    FrameData fd;
    fd.t = ft.t;
    fd.epoch = ft.epoch;
    fd.C_l_b_est = C_l_b_est[fi];
    fd.yaw_est = dcmToEuler(fd.C_l_b_est).z();
    fd.truth = ft.pos;

    field.ecefDirections(ft.epoch, ecef);
    const Eigen::Matrix3d cam_from_ned = (C_l_b_true * C_b_c_true).transpose();

    for (size_t i = 0; i < ecef.size(); ++i) {
      // Star direction in the local NED frame at the TRUE position.
      Eigen::Vector3d v_ned = C_ecef_ned.transpose() * ecef[i];
      const double el_true = nedToElevation(v_ned);
      if (el_true < err.min_elevation) continue;

      // Forward refraction: the star APPEARS higher than it is.
      if (err.atmos.pressure_hpa > 0.0) {
        const double az = nedToAzimuth(v_ned);
        const double r = refraction::saemundsson(el_true, err.atmos.pressure_hpa,
                                                 err.atmos.temperature_c);
        v_ned = azElToNed(az, el_true + r);
      }

      Eigen::Vector3d v_cam = cam_from_ned * v_ned;
      if (!cam.inFov(v_cam)) continue;

      if (err.centroid_noise > 0) {
        v_cam += Eigen::Vector3d(cent_n(rng), cent_n(rng), 0.0);
        v_cam.normalize();
      }
      fd.id.push_back(field.id(i));
      fd.v_cam.push_back(v_cam);
    }
    out.push_back(std::move(fd));
  }
  return out;
}

// ---------------------------------------------------------------------------

namespace {

/// Per-frame fix using an assumed camera mounting.
bool frameFixZenith(const FrameData& fd, const StarField& field,
                    const Eigen::Matrix3d& C_b_c,
                    const Atmosphere& refraction_model,
                    Eigen::Vector3d& zenith_out) {
  if (fd.id.size() < 3) return false;
  const Eigen::Matrix3d R = fd.C_l_b_est * C_b_c;

  std::vector<StarSight> sights;
  sights.reserve(fd.id.size());
  for (size_t i = 0; i < fd.id.size(); ++i) {
    // Recover the NED direction using the ESTIMATED attitude and mounting.
    // This is Eq. (1) and it is where the misalignment enters.
    const Eigen::Vector3d v_ned = R * fd.v_cam[i];
    const double el = nedToElevation(v_ned);

    const int k = field.indexOf(fd.id[i]);
    if (k < 0) continue;

    StarSight ss;
    ss.id = fd.id[i];
    // The sub-stellar point IS the star's ECEF direction.
    ss.gp = field.ecefDirection(static_cast<size_t>(k), fd.epoch);
    ss.zenith_angle = M_PI_2 - el;
    sights.push_back(ss);
  }
  if (sights.size() < 3) return false;

  if (refraction_model.pressure_hpa > 0.0) {
    correctRefraction(sights, refraction_model);
    weightForRefraction(sights, refraction_model);
  }

  // RANSAC. A graduated non-convexity solver was tried here and lost by an
  // order of magnitude under misidentified stars; see position_solver.hpp.
  // RANSAC also preserves the refraction weights applied above.
  RansacConfig rc;
  rc.tolerance = 0.05 * kDeg2Rad;
  const Fix f = solveFixRansac(sights, rc);
  if (!f.ok) return false;
  zenith_out = f.zenith_ecef;
  return true;
}

/// Fit a small circle to unit vectors: find axis c and radius rho with
/// c . p_i = cos(rho). This is a plane fit, so the axis is the eigenvector of
/// the smallest eigenvalue of the scatter matrix about the centroid.
void fitSmallCircle(const std::vector<Eigen::Vector3d>& pts,
                    Eigen::Vector3d& axis, double& radius) {
  Eigen::Vector3d mean = Eigen::Vector3d::Zero();
  for (const auto& p : pts) mean += p;
  mean /= static_cast<double>(pts.size());

  Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
  for (const auto& p : pts) {
    const Eigen::Vector3d d = p - mean;
    C += d * d.transpose();
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
  Eigen::Vector3d n = es.eigenvectors().col(0);  // smallest eigenvalue
  if (n.dot(mean) < 0) n = -n;
  n.normalize();

  axis = n;
  radius = std::acos(std::clamp(n.dot(mean.normalized()) * 1.0, -1.0, 1.0));
  // More directly: cos(rho) = n . p_i, averaged.
  double cs = 0.0;
  for (const auto& p : pts) cs += n.dot(p);
  cs /= static_cast<double>(pts.size());
  radius = std::acos(std::clamp(cs, -1.0, 1.0));
}

}  // namespace

bool singleFrameFix(const FrameData& frame, const Eigen::Matrix3d& C_b_c,
                    const Atmosphere& refraction_model, Geodetic& out) {
  const StarField field(brightStarCatalog(), frame.epoch);
  Eigen::Vector3d z;
  if (!frameFixZenith(frame, field, C_b_c, refraction_model, z)) return false;
  out = fromZenith(z, frame.truth.alt);
  return true;
}

OrbitResult estimateOrbit(const std::vector<FrameData>& frames,
                          const Eigen::Matrix3d& C_b_c_assumed,
                          AverageMethod method,
                          const Atmosphere& refraction_model) {
  std::vector<Eigen::Vector3d> pts;
  std::vector<double> headings;
  if (frames.empty()) return OrbitResult{};
  pts.reserve(frames.size());

  const StarField field(brightStarCatalog(),
                        frames[frames.size() / 2].epoch);

  // Transport every fix to the LAST frame's epoch, so that averaging removes
  // the misalignment circle rather than smearing in the aircraft's own motion.
  // Zero when dr_ne is unset, which reduces to the original behaviour.
  const Eigen::Vector2d ref_ne = frames.back().dr_ne;

  for (const FrameData& fd : frames) {
    Eigen::Vector3d z;
    if (!frameFixZenith(fd, field, C_b_c_assumed, refraction_model, z)) continue;

    const Eigen::Vector2d d = ref_ne - fd.dr_ne;
    if (d.squaredNorm() > 1.0) {
      // Small-angle displacement of the zenith direction: move the point by
      // d north/east on the sphere.
      const Eigen::Matrix3d C = nedBasisEcef(fromZenith(z, 0.0));
      z += (C.col(0) * d.x() + C.col(1) * d.y()) / kEarthR;
      z.normalize();
    }
    pts.push_back(z);
    headings.push_back(fd.yaw_est);
  }

  OrbitResult res;
  if (pts.size() < 8) return res;

  // ROBUST OUTLIER REJECTION.
  //
  // A per-frame fix from the minimum three stars is exactly determined: one
  // bad identification or centroid goes straight into the answer with no
  // redundancy and no residual to detect it. Such fixes land thousands of km
  // away and, being averaged in, drag the whole estimate.
  //
  // Drop points far from the median direction. The threshold is generous --
  // this is for catastrophes, not for trimming the distribution, and it must
  // not eat the legitimate circle traced by an uncalibrated mounting.
  {
    Eigen::Vector3d c = Eigen::Vector3d::Zero();
    for (const auto& p : pts) c += p;
    if (c.norm() > 1e-9) {
      c.normalize();
      std::vector<double> ang;
      ang.reserve(pts.size());
      for (const auto& p : pts) {
        ang.push_back(std::acos(std::clamp(p.dot(c), -1.0, 1.0)));
      }
      std::vector<double> srt = ang;
      std::nth_element(srt.begin(), srt.begin() + srt.size() / 2, srt.end());
      const double med = srt[srt.size() / 2];
      // 5x the median deviation, floored at 5 deg (~550 km) so a genuine
      // 44 km circle is never trimmed.
      const double lim = std::max(5.0 * med, 5.0 * kDeg2Rad);
      std::vector<Eigen::Vector3d> kp;
      std::vector<double> kh;
      for (size_t i = 0; i < pts.size(); ++i) {
        if (ang[i] <= lim) { kp.push_back(pts[i]); kh.push_back(headings[i]); }
      }
      if (kp.size() >= 8) { pts.swap(kp); headings.swap(kh); }
    }
  }

  res.n_frames = static_cast<int>(pts.size());

  Eigen::Vector3d mean = Eigen::Vector3d::Zero();

  if (method == AverageMethod::CircleFit) {
    double rho = 0.0;
    fitSmallCircle(pts, mean, rho);
    res.circle_radius = rho;
  } else if (method == AverageMethod::HeadingWeighted) {
    // Weight each sample by the heading interval it represents. Equivalent to
    // resampling uniformly in heading, which is what the naive mean assumes
    // but does not get when the aircraft holds a ground track in wind.
    double wsum = 0.0;
    const size_t n = pts.size();
    for (size_t i = 0; i < n; ++i) {
      const double prev = headings[(i + n - 1) % n];
      const double next = headings[(i + 1) % n];
      double dpsi = std::atan2(std::sin(next - prev), std::cos(next - prev)) / 2.0;
      dpsi = std::abs(dpsi);
      if (dpsi < 1e-9) dpsi = 1e-9;
      mean += dpsi * pts[i];
      wsum += dpsi;
    }
    if (wsum > 0) mean /= wsum;
  } else {
    for (const auto& p : pts) mean += p;
    mean /= static_cast<double>(pts.size());
  }

  if (mean.norm() < 1e-9) return res;
  const Eigen::Vector3d axis = mean.normalized();

  double s2 = 0.0;
  for (const auto& p : pts) {
    const double a = std::acos(std::clamp(p.dot(axis), -1.0, 1.0));
    s2 += a * a;
  }
  res.spread = std::sqrt(s2 / pts.size());

  res.ok = true;
  res.position = fromZenith(axis, frames.front().truth.alt);
  return res;
}

HeadingResult celestialHeading(const std::vector<FrameData>& frames,
                               const Geodetic& pos,
                               const Eigen::Matrix3d& C_b_c) {
  HeadingResult out;
  if (frames.empty()) return out;
  const StarField field(brightStarCatalog(), frames[frames.size() / 2].epoch);
  const Eigen::Matrix3d C_ecef_ned = nedBasisEcef(pos);

  // Accumulate as unit vectors so the circular mean is well defined.
  double sx = 0, sy = 0, bank = 0;
  std::vector<double> errs;
  errs.reserve(frames.size());

  for (const FrameData& fd : frames) {
    if (fd.id.size() < 3) continue;

    std::vector<Eigen::Vector3d> a, b;
    a.reserve(fd.id.size());
    b.reserve(fd.id.size());
    for (size_t i = 0; i < fd.id.size(); ++i) {
      const int k = field.indexOf(fd.id[i]);
      if (k < 0) continue;
      const Eigen::Vector3d e =
          field.ecefDirection(static_cast<size_t>(k), fd.epoch);
      a.push_back(fd.v_cam[i]);
      b.push_back(C_ecef_ned.transpose() * e);
    }
    if (a.size() < 3) continue;

    // Kabsch gives camera -> NED. With the mounting known, that is the true
    // body attitude, and its yaw is the heading the stars actually see.
    const Eigen::Matrix3d R = kabsch(a, b);
    const Eigen::Matrix3d C_l_b = R * C_b_c.transpose();
    const double yaw = std::atan2(C_l_b(1, 0), C_l_b(0, 0));
    const double err = std::atan2(std::sin(yaw - fd.yaw_est),
                                  std::cos(yaw - fd.yaw_est));
    sx += std::cos(err);
    sy += std::sin(err);
    errs.push_back(err);
    bank += std::abs(std::atan2(C_l_b(2, 1), C_l_b(2, 2)));
    ++out.n;
  }
  if (out.n < 3) return out;

  out.error = std::atan2(sy / out.n, sx / out.n);
  out.mean_bank = bank / out.n;
  double v = 0;
  for (double e : errs) {
    const double d = std::atan2(std::sin(e - out.error), std::cos(e - out.error));
    v += d * d;
  }
  out.spread = std::sqrt(v / out.n);
  out.ok = true;
  return out;
}

Eigen::Matrix3d recalibrateMount(const std::vector<FrameData>& frames,
                                 const Geodetic& pos) {
  if (frames.empty()) return nominalCameraMount();
  const StarField field(brightStarCatalog(), frames[frames.size() / 2].epoch);
  const Eigen::Matrix3d C_ecef_ned = nedBasisEcef(pos);

  std::vector<Eigen::Matrix3d> mounts;
  mounts.reserve(frames.size());

  for (const FrameData& fd : frames) {
    if (fd.id.size() < 3) continue;

    std::vector<Eigen::Vector3d> a;  // observed, camera frame
    std::vector<Eigen::Vector3d> b;  // theoretical, NED at the assumed position
    a.reserve(fd.id.size());
    b.reserve(fd.id.size());
    for (size_t i = 0; i < fd.id.size(); ++i) {
      const int k = field.indexOf(fd.id[i]);
      if (k < 0) continue;
      const Eigen::Vector3d e =
          field.ecefDirection(static_cast<size_t>(k), fd.epoch);
      a.push_back(fd.v_cam[i]);
      b.push_back(C_ecef_ned.transpose() * e);
    }
    if (a.size() < 3) continue;

    // R = C_l_b * C_b_c maps camera -> NED, so C_b_c = C_l_b^T * R.
    const Eigen::Matrix3d R = kabsch(a, b);
    mounts.push_back(fd.C_l_b_est.transpose() * R);
  }
  if (mounts.empty()) return nominalCameraMount();
  return averageRotation(mounts);
}

Eigen::Vector3d meanVerticalBody(const std::vector<FrameData>& frames) {
  Eigen::Vector3d u = Eigen::Vector3d::Zero();
  for (const FrameData& fd : frames) {
    // C_l_b maps body -> NED, so C_l_b^T * z_ned is the third ROW of C_l_b:
    // the local vertical in body axes.
    u += fd.C_l_b_est.row(2).transpose();
  }
  if (u.norm() < 1e-9) return Eigen::Vector3d::UnitZ();
  return u.normalized();
}

Eigen::Matrix3d mountTiltOnly(const Eigen::Matrix3d& C_b_c,
                              const Eigen::Vector3d& vertical_body) {
  const Eigen::Matrix3d nominal = nominalCameraMount();
  Eigen::Vector3d u = vertical_body;
  if (u.norm() < 1e-9) u = Eigen::Vector3d::UnitZ();
  u.normalize();

  // The yaw bias enters as a BODY-frame rotation on the left, so take the
  // deviation from nominal in the body frame and project out the vertical.
  const Eigen::AngleAxisd aa(C_b_c * nominal.transpose());
  Eigen::Vector3d d = aa.axis() * aa.angle();
  d -= u * d.dot(u);
  const double n = d.norm();
  if (n < 1e-12) return nominal;
  return Eigen::AngleAxisd(n, d / n).toRotationMatrix() * nominal;
}

std::vector<IterationStep> iterateOrbit(const std::vector<FrameData>& frames,
                                        Eigen::Matrix3d C_b_c,
                                        const Geodetic& truth,
                                        const IterationConfig& cfg,
                                        const Eigen::Matrix3d* true_mount) {
  std::vector<IterationStep> steps;
  Geodetic prev;
  bool have_prev = false;
  double last_move = -1.0;

  for (int it = 0; it < cfg.max_iterations; ++it) {
    const OrbitResult r =
        estimateOrbit(frames, C_b_c, cfg.method, cfg.refraction_model);
    if (!r.ok) break;

    IterationStep step;
    step.position = r.position;
    step.error_m = haversine(r.position, truth);
    step.circle_radius_km =
        (cfg.method == AverageMethod::CircleFit)
            ? r.circle_radius * kEarthR / 1000.0
            : r.spread * kEarthR / 1000.0;
    step.mount_error_deg =
        true_mount ? rotationAngleBetween(C_b_c, *true_mount) * kRad2Deg : 0.0;
    steps.push_back(step);

    if (have_prev && haversine(r.position, prev) < cfg.convergence_m) break;

    // DIVERGENCE GUARD -- a CONTRACTION test.
    //
    // The loop is a fixed-point iteration: estimate position, recalibrate the
    // mounting on it, re-estimate. That converges only if each pass moves the
    // estimate LESS than the last. When it does not, the mounting is absorbing
    // the position error and feeding it back, and the iteration walks away to
    // a wrong fixed point.
    //
    // Observed in ArduPilot SITL: a circle fit went 5.13 km at iteration 1 to
    // 15.19 km after six, and aided, 31.68 km to 67.53 km. Iterating was
    // actively destroying the best estimate available.
    //
    // A first attempt bounded the step against the fitted circle radius. That
    // was useless: the radius is the misalignment expressed as a distance,
    // tens of km for a fraction of a degree of mounting error, so the bound sat
    // far above any real step and never fired.
    //
    // The contraction test needs no scale at all. A convergent iteration has
    // shrinking steps; the moment one grows, stop and keep what came before.
    const double moved = haversine(r.position, prev);
    if (last_move > 0 && moved > last_move) {
      steps.pop_back();   // the step that grew is not trusted
      break;
    }
    last_move = moved;

    prev = r.position;
    have_prev = true;

    // The circle fit already estimates the misalignment; recalibrating from
    // its own output double-counts it. See IterationConfig::recalibrate_mount.
    if (!cfg.recalibrate_mount || cfg.method == AverageMethod::CircleFit) break;
    C_b_c = recalibrateMount(frames, r.position);
  }
  return steps;
}

}  // namespace celestial
