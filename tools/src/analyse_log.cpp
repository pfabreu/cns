// analyse_log.cpp — MILESTONE 3.
//
//   ./build/analyse_log flight.csv [--imaging] [--utc YYYY-MM-DDTHH:MM:SS]
//
// Three things, in order of importance:
//
//   1. Decompose the EKF3 attitude error into body-frame and NED components
//      around each orbit. This tests the assumption Milestone 1 baked in and
//      the paper never verified: that the AHRS error is BODY-FIXED and
//      therefore averages away over a heading revolution.
//
//   2. Run the position estimate on the real EKF3 attitude with IDEAL star
//      vectors (the Milestone 1 path).
//
//   3. Optionally re-run it through the full imaging chain (Milestone 2 path)
//      with --imaging, so the two error sources stay separable.

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>


#include "celestial/imaging.hpp"
#include "celestial/orbit.hpp"
#include "celestial/sky_model.hpp"
#include "celestial/star_catalog.hpp"

#include <Eigen/Geometry>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <unordered_map>

// ---------------------------------------------------------------------------
// Flight-log loading and attitude-error decomposition. Inlined here because
// analyse_log is its only caller.
// ---------------------------------------------------------------------------
namespace celestial {

/// One resampled row from an ArduPilot log. Produced by tools/sitl/dump_log.py.
struct LogRow {
  double t = 0.0;  ///< seconds since the first row
  Geodetic pos_true;
  double roll_true = 0, pitch_true = 0, yaw_true = 0;  ///< SIM message
  double roll_ekf = 0, pitch_ekf = 0, yaw_ekf = 0;     ///< ATT message
  double gps_week = -1, gps_ms = -1;
  double airspeed = 0.0;
  /// Raw gyro, body rates relative to INERTIAL, rad/s. Zero if the log carried
  /// no IMU messages, which disables --star-aided.
  Eigen::Vector3d gyro = Eigen::Vector3d::Zero();

  Eigen::Matrix3d C_true() const {
    return eulerToDcm(roll_true, pitch_true, yaw_true);
  }
  Eigen::Matrix3d C_ekf() const {
    return eulerToDcm(roll_ekf, pitch_ekf, yaw_ekf);
  }
};

/// Load the CSV emitted by tools/sitl/dump_log.py. Angles in the file are degrees;
/// they are converted to radians here.
std::vector<LogRow> loadLogCsv(const std::string& path);

/// UTC epoch from GPS week + milliseconds of week.
/// leap_seconds is GPS-UTC (18 as of 2017; SITL logs are usually consistent).
Epoch epochFromGps(double gps_week, double gps_ms, int leap_seconds = 18);

// ---------------------------------------------------------------------------

/// Attitude error at one frame.
struct AttitudeError {
  double t = 0.0;
  double yaw = 0.0;         ///< true heading, for binning
  Eigen::Vector3d rv_body;  ///< rotation vector of C_true^T C_ekf, BODY frame
  Eigen::Vector3d rv_ned;   ///< the same rotation vector in NED
  double tilt = 0.0;        ///< magnitude of the horizontal (NED) part
};

std::vector<AttitudeError> attitudeErrors(const std::vector<LogRow>& rows);

/// Summary of how the attitude error behaves around a heading revolution.
///
/// `body_mean` is the average body-frame error; if the error is genuinely
/// body-fixed this is the whole story and `body_std_vs_heading` is small.
/// `ned_mean` is the average in NED; a body-fixed error averages to ~0 there,
/// so a LARGE ned_mean is the signature of an earth-fixed component that the
/// orbit will not remove.
struct ErrorDecomposition {
  int n = 0;
  Eigen::Vector3d body_mean = Eigen::Vector3d::Zero();
  Eigen::Vector3d body_std = Eigen::Vector3d::Zero();
  Eigen::Vector3d ned_mean = Eigen::Vector3d::Zero();
  Eigen::Vector3d ned_std = Eigen::Vector3d::Zero();
  /// Predicted residual position error from the NON-averaging (NED-mean) part,
  /// metres. This is the number that matters.
  double predicted_residual_m = 0.0;
};

ErrorDecomposition decompose(const std::vector<AttitudeError>& errs);

/// Per-heading-bin means of the body-frame error, for the flat-vs-varying test.
struct HeadingBin {
  double heading_deg = 0.0;
  int n = 0;
  Eigen::Vector3d body_mean = Eigen::Vector3d::Zero();
};
std::vector<HeadingBin> binByHeading(const std::vector<AttitudeError>& errs,
                                     int n_bins = 12);

// ---------------------------------------------------------------------------

/// A contiguous stretch of log covering one full revolution of heading.
struct OrbitSegment {
  size_t begin = 0;
  size_t end = 0;  ///< exclusive
  double total_heading_change = 0.0;  ///< signed, radians
  bool clockwise = false;
  Geodetic centre;  ///< mean position, used as the truth for scoring
  double mean_radius_m = 0.0;
};

/// Split a log into full heading revolutions. Straight legs are discarded.
std::vector<OrbitSegment> findOrbits(const std::vector<LogRow>& rows);

/// Build FrameTruth records (true pose + epoch) for a segment.
std::vector<FrameTruth> toFrameTruth(const std::vector<LogRow>& rows,
                                     const OrbitSegment& seg,
                                     const Epoch& start_epoch);

/// The EKF3 attitude matrices for the same segment, to pass to
/// simulateObservationsWithAttitude().
std::vector<Eigen::Matrix3d> ekfAttitudes(const std::vector<LogRow>& rows,
                                          const OrbitSegment& seg);


namespace {

constexpr double kEarthR = 6371008.8;

double wrapPi(double a) { return std::atan2(std::sin(a), std::cos(a)); }

/// Rotation vector (axis * angle) of a rotation matrix.
Eigen::Vector3d rotationVector(const Eigen::Matrix3d& R) {
  const Eigen::AngleAxisd aa(R);
  return aa.axis() * aa.angle();
}

}  // namespace

std::vector<LogRow> loadLogCsv(const std::string& path) {
  std::vector<LogRow> rows;
  std::ifstream f(path);
  if (!f) return rows;

  std::string line;
  if (!std::getline(f, line)) return rows;

  // Header-driven so column order does not matter.
  std::unordered_map<std::string, int> col;
  {
    std::stringstream ss(line);
    std::string name;
    int i = 0;
    while (std::getline(ss, name, ',')) {
      name.erase(std::remove_if(name.begin(), name.end(),
                                [](char c) { return c == '\r' || c == ' '; }),
                 name.end());
      col[name] = i++;
    }
  }
  auto has = [&](const char* k) { return col.count(k) > 0; };

  while (std::getline(f, line)) {
    std::vector<double> v;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
      try { v.push_back(std::stod(cell)); } catch (...) { v.push_back(0.0); }
    }
    auto get = [&](const char* k, double dflt = 0.0) {
      if (!has(k)) return dflt;
      const int i = col[k];
      return (i < static_cast<int>(v.size())) ? v[i] : dflt;
    };

    LogRow r;
    r.t = get("t");
    r.pos_true.lat = get("lat_true") * kDeg2Rad;
    r.pos_true.lon = get("lon_true") * kDeg2Rad;
    r.pos_true.alt = get("alt_true");
    r.roll_true = get("roll_true") * kDeg2Rad;
    r.pitch_true = get("pitch_true") * kDeg2Rad;
    r.yaw_true = get("yaw_true") * kDeg2Rad;
    r.roll_ekf = get("roll_ekf") * kDeg2Rad;
    r.pitch_ekf = get("pitch_ekf") * kDeg2Rad;
    r.yaw_ekf = get("yaw_ekf") * kDeg2Rad;
    r.gps_week = get("gps_week", -1);
    r.gps_ms = get("gps_ms", -1);
    r.airspeed = get("airspeed");
    r.gyro = Eigen::Vector3d(get("gyro_x"), get("gyro_y"), get("gyro_z"));
    rows.push_back(r);
  }
  return rows;
}

Epoch epochFromGps(double gps_week, double gps_ms, int leap_seconds) {
  // GPS epoch: 1980-01-06 00:00:00 UTC.
  auto base = Epoch::fromUtc(1980, 1, 6, 0, 0, 0.0);
  const double secs =
      gps_week * 7.0 * 86400.0 + gps_ms / 1000.0 - double(leap_seconds);
  return base ? base->advanced(secs) : Epoch();
}

// ---------------------------------------------------------------------------

std::vector<AttitudeError> attitudeErrors(const std::vector<LogRow>& rows) {
  std::vector<AttitudeError> out;
  out.reserve(rows.size());
  for (const LogRow& r : rows) {
    const Eigen::Matrix3d Ct = r.C_true();
    const Eigen::Matrix3d Ce = r.C_ekf();
    // dC maps body-true into body-estimated: the error as seen in body axes.
    const Eigen::Matrix3d dC = Ct.transpose() * Ce;
    AttitudeError e;
    e.t = r.t;
    e.yaw = r.yaw_true;
    e.rv_body = rotationVector(dC);
    e.rv_ned = Ct * e.rv_body;
    // Only the horizontal part of the NED rotation vector tilts the zenith;
    // the vertical part is a heading error, which the position fix ignores.
    e.tilt = std::hypot(e.rv_ned.x(), e.rv_ned.y());
    out.push_back(e);
  }
  return out;
}

ErrorDecomposition decompose(const std::vector<AttitudeError>& errs) {
  ErrorDecomposition d;
  d.n = static_cast<int>(errs.size());
  if (errs.empty()) return d;

  for (const AttitudeError& e : errs) {
    d.body_mean += e.rv_body;
    d.ned_mean += e.rv_ned;
  }
  d.body_mean /= d.n;
  d.ned_mean /= d.n;

  for (const AttitudeError& e : errs) {
    const Eigen::Vector3d db = e.rv_body - d.body_mean;
    const Eigen::Vector3d dn = e.rv_ned - d.ned_mean;
    d.body_std += db.cwiseProduct(db);
    d.ned_std += dn.cwiseProduct(dn);
  }
  d.body_std = (d.body_std / d.n).cwiseSqrt();
  d.ned_std = (d.ned_std / d.n).cwiseSqrt();

  // A body-fixed error averages to ~0 in NED over a full revolution. Whatever
  // survives in the NED mean is the part the orbit CANNOT remove, and its
  // horizontal component maps to position at ~111 km/deg.
  d.predicted_residual_m =
      std::hypot(d.ned_mean.x(), d.ned_mean.y()) * kEarthR;
  return d;
}

std::vector<HeadingBin> binByHeading(const std::vector<AttitudeError>& errs,
                                     int n_bins) {
  std::vector<HeadingBin> bins(n_bins);
  for (int i = 0; i < n_bins; ++i) {
    bins[i].heading_deg = 360.0 * i / n_bins;
  }
  for (const AttitudeError& e : errs) {
    double h = e.yaw;
    while (h < 0) h += 2 * M_PI;
    while (h >= 2 * M_PI) h -= 2 * M_PI;
    const int b = std::min(n_bins - 1, int(h / (2 * M_PI) * n_bins));
    bins[b].body_mean += e.rv_body;
    bins[b].n++;
  }
  for (HeadingBin& b : bins) {
    if (b.n) b.body_mean /= b.n;
  }
  return bins;
}

// ---------------------------------------------------------------------------

std::vector<OrbitSegment> findOrbits(const std::vector<LogRow>& rows) {
  std::vector<OrbitSegment> out;
  if (rows.size() < 20) return out;

  size_t start = 0;
  double accum = 0.0;
  for (size_t i = 1; i < rows.size(); ++i) {
    const double d = wrapPi(rows[i].yaw_true - rows[i - 1].yaw_true);
    // A big jump means a discontinuity, not a turn: restart.
    if (std::abs(d) > 0.5) { start = i; accum = 0.0; continue; }
    // Sign change of more than a few degrees means the turn reversed.
    if (accum != 0.0 && d * accum < 0 && std::abs(d) > 0.02 * std::abs(accum)) {
      // keep accumulating; a little wobble is normal in wind
    }
    accum += d;

    if (std::abs(accum) >= 2.0 * M_PI) {
      OrbitSegment seg;
      seg.begin = start;
      seg.end = i + 1;
      seg.total_heading_change = accum;
      seg.clockwise = accum > 0;

      double lat = 0, lon = 0, alt = 0;
      for (size_t k = seg.begin; k < seg.end; ++k) {
        lat += rows[k].pos_true.lat;
        lon += rows[k].pos_true.lon;
        alt += rows[k].pos_true.alt;
      }
      const double n = double(seg.end - seg.begin);
      seg.centre = Geodetic{lat / n, lon / n, alt / n};

      double rsum = 0;
      for (size_t k = seg.begin; k < seg.end; ++k) {
        rsum += haversine(rows[k].pos_true, seg.centre);
      }
      seg.mean_radius_m = rsum / n;

      out.push_back(seg);
      start = i + 1;
      accum = 0.0;
    }
  }
  return out;
}

std::vector<FrameTruth> toFrameTruth(const std::vector<LogRow>& rows,
                                     const OrbitSegment& seg,
                                     const Epoch& start_epoch) {
  std::vector<FrameTruth> out;
  out.reserve(seg.end - seg.begin);
  const double t0 = rows[seg.begin].t;
  for (size_t k = seg.begin; k < seg.end; ++k) {
    FrameTruth f;
    f.t = rows[k].t - t0;
    f.pos = rows[k].pos_true;
    f.roll = rows[k].roll_true;
    f.pitch = rows[k].pitch_true;
    f.yaw = rows[k].yaw_true;
    f.epoch = start_epoch.advanced(f.t);
    out.push_back(f);
  }
  return out;
}

std::vector<Eigen::Matrix3d> ekfAttitudes(const std::vector<LogRow>& rows,
                                          const OrbitSegment& seg) {
  std::vector<Eigen::Matrix3d> out;
  out.reserve(seg.end - seg.begin);
  for (size_t k = seg.begin; k < seg.end; ++k) out.push_back(rows[k].C_ekf());
  return out;
}

}  // namespace celestial

using namespace celestial;

namespace {

constexpr double kEarthR = 6371008.8;

void printDecomposition(const std::vector<AttitudeError>& errs) {
  const ErrorDecomposition d = decompose(errs);
  std::printf("    attitude error, %d samples (degrees)\n", d.n);
  std::printf("      body-frame mean  roll %+7.4f  pitch %+7.4f  yaw %+7.4f\n",
              d.body_mean.x() * kRad2Deg, d.body_mean.y() * kRad2Deg,
              d.body_mean.z() * kRad2Deg);
  std::printf("      body-frame std   roll %7.4f  pitch %7.4f  yaw %7.4f\n",
              d.body_std.x() * kRad2Deg, d.body_std.y() * kRad2Deg,
              d.body_std.z() * kRad2Deg);
  std::printf("      NED mean         N    %+7.4f  E     %+7.4f  D   %+7.4f\n",
              d.ned_mean.x() * kRad2Deg, d.ned_mean.y() * kRad2Deg,
              d.ned_mean.z() * kRad2Deg);

  std::printf("\n      body-frame mean magnitude : %.4f deg  (%.1f km if it did"
              " NOT average)\n",
              d.body_mean.head<2>().norm() * kRad2Deg,
              d.body_mean.head<2>().norm() * kEarthR / 1000.0);
  std::printf("      NED-mean horizontal       : %.4f deg  -> PREDICTED"
              " RESIDUAL %.2f km\n",
              std::hypot(d.ned_mean.x(), d.ned_mean.y()) * kRad2Deg,
              d.predicted_residual_m / 1000.0);
  std::printf("\n      A body-fixed error averages to ~0 in NED over a full\n"
              "      revolution. Whatever survives in the NED mean is the part\n"
              "      the orbit CANNOT remove.\n");
}

void printHeadingBins(const std::vector<AttitudeError>& errs) {
  const auto bins = binByHeading(errs, 12);
  std::printf("\n    body-frame error vs heading (deg) -- FLAT means body-fixed\n");
  std::printf("      %-10s %6s %9s %9s %9s\n", "heading", "n", "roll", "pitch",
              "yaw");
  double mn[3] = {1e9, 1e9, 1e9}, mx[3] = {-1e9, -1e9, -1e9};
  for (const HeadingBin& b : bins) {
    if (!b.n) continue;
    std::printf("      %-10.0f %6d %+9.4f %+9.4f %+9.4f\n", b.heading_deg, b.n,
                b.body_mean.x() * kRad2Deg, b.body_mean.y() * kRad2Deg,
                b.body_mean.z() * kRad2Deg);
    for (int i = 0; i < 3; ++i) {
      mn[i] = std::min(mn[i], b.body_mean[i] * kRad2Deg);
      mx[i] = std::max(mx[i], b.body_mean[i] * kRad2Deg);
    }
  }
  std::printf("      %-10s %6s %+9.4f %+9.4f %+9.4f   <- peak-to-peak\n",
              "range", "", mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]);
  const double worst = std::max(mx[0] - mn[0], mx[1] - mn[1]);
  std::printf("\n      Worst tilt swing %.4f deg = %.1f km of position, which\n"
              "      is the part of the error that does NOT behave as the\n"
              "      paper's method assumes.\n",
              worst, worst * kDeg2Rad * kEarthR / 1000.0);
}

}  // namespace

static int runSynth(double body_deg, double ned_deg) {

  OrbitConfig cfg;
  cfg.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  cfg.radius_m = 600; cfg.airspeed = 25; cfg.frame_rate = 10; cfg.revolutions = 2;
  const auto truth = generateOrbit(cfg);

  std::printf("t,lat_true,lon_true,alt_true,roll_true,pitch_true,yaw_true,"
              "roll_ekf,pitch_ekf,yaw_ekf,gps_week,gps_ms,airspeed\n");
  for (const FrameTruth& f : truth) {
    const Eigen::Matrix3d Ct = eulerToDcm(f.roll, f.pitch, f.yaw);
    // Body-fixed pitch error, plus an earth-fixed tilt about NED north.
    const Eigen::Matrix3d dB =
        Eigen::AngleAxisd(body_deg * kDeg2Rad, Eigen::Vector3d::UnitY()).toRotationMatrix();
    const Eigen::Matrix3d dN =
        Eigen::AngleAxisd(ned_deg * kDeg2Rad, Eigen::Vector3d::UnitX()).toRotationMatrix();
    const Eigen::Matrix3d Ce = dN * Ct * dB;
    const Eigen::Vector3d e = dcmToEuler(Ce);
    std::printf("%.6f,%.9f,%.9f,%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,-1,-1,25\n",
                f.t, f.pos.lat * kRad2Deg, f.pos.lon * kRad2Deg, f.pos.alt,
                f.roll * kRad2Deg, f.pitch * kRad2Deg, f.yaw * kRad2Deg,
                e.x() * kRad2Deg, e.y() * kRad2Deg, e.z() * kRad2Deg);
  }
  return 0;
}

int main(int argc, char** argv) {
  // Synthetic mode: generate a log with a KNOWN error and analyse it, to
  // validate the decomposition before trusting it on real data.
  if (argc > 1 && !std::strcmp(argv[1], "--synth")) {
    const double b = (argc > 2) ? std::atof(argv[2]) : 0.10;
    const double n = (argc > 3) ? std::atof(argv[3]) : 0.03;
    return runSynth(b, n);
  }

  if (argc < 2) {
    std::printf("usage: analyse_log flight.csv [--imaging] [--utc ISO8601]\n");
    return 1;
  }
  bool imaging = false, star_aided = false;
  std::string utc;
  for (int i = 2; i < argc; ++i) {
    if (!std::strcmp(argv[i], "--imaging")) imaging = true;
    else if (!std::strcmp(argv[i], "--star-aided")) star_aided = true;
    else if (!std::strcmp(argv[i], "--utc") && i + 1 < argc) utc = argv[++i];
  }

  const auto rows = loadLogCsv(argv[1]);
  if (rows.size() < 50) {
    std::printf("only %zu rows loaded -- check the CSV\n", rows.size());
    return 1;
  }
  std::printf("loaded %zu rows, %.1f s\n", rows.size(), rows.back().t);

  // Absolute time. Prefer the log's GPS time; --utc overrides.
  Epoch start;
  if (!utc.empty()) {
    int y, mo, d, h, mi;
    double sec;
    if (std::sscanf(utc.c_str(), "%d-%d-%dT%d:%d:%lf", &y, &mo, &d, &h, &mi,
                    &sec) == 6) {
      auto e = Epoch::fromUtc(y, mo, d, h, mi, sec);
      if (e) start = *e;
    }
  } else if (rows.front().gps_week > 0) {
    start = epochFromGps(rows.front().gps_week, rows.front().gps_ms);
    std::printf("epoch from GPS week %.0f\n", rows.front().gps_week);
  } else {
    auto e = Epoch::fromUtc(2024, 9, 15, 10, 30, 0.0);
    if (e) start = *e;
    std::printf("no GPS time in log; assuming 2024-09-15T10:30:00Z\n");
  }

  const auto orbits = findOrbits(rows);
  std::printf("found %zu full heading revolutions\n", orbits.size());
  if (orbits.empty()) {
    std::printf("\nNo complete 360 deg heading sweep in this log. Fly a\n"
                "LOITER (or a fixed-bank FBWA turn) for at least one full\n"
                "revolution.\n");
    return 1;
  }

  Camera cam;
  const auto catalog = catalogBrighterThan(4.0);
  std::vector<double> vmags;
  for (const CatalogStar& s : catalog) vmags.push_back(s.vmag);

  ErrorModel err;
  err.boresight_error = 0.4 * kDeg2Rad;  // uncalibrated camera, as in Fig. 7a
  err.ahrs_bias.setZero();               // the REAL EKF3 error comes from the log
  err.ahrs_noise = 0.0;
  err.mag_limit = 4.0;
  const Eigen::Matrix3d true_mount =
      nominalCameraMount() *
      smallRotation(err.boresight_axis, err.boresight_error);

  for (size_t oi = 0; oi < orbits.size(); ++oi) {
    const OrbitSegment& seg = orbits[oi];
    std::printf("\n================ orbit %zu ================\n", oi + 1);
    std::printf("  %s, radius %.0f m, %zu samples, %.1f s\n",
                seg.clockwise ? "CW" : "CCW", seg.mean_radius_m,
                seg.end - seg.begin, rows[seg.end - 1].t - rows[seg.begin].t);

    // --- 1. attitude error decomposition -----------------------------------
    std::vector<LogRow> sub(rows.begin() + seg.begin, rows.begin() + seg.end);
    const auto errs = attitudeErrors(sub);
    std::printf("\n  [1] EKF3 ATTITUDE ERROR\n");
    printDecomposition(errs);
    printHeadingBins(errs);

    // --- 2. ideal-vector position estimate (Milestone 1 path) --------------
    const auto truth = toFrameTruth(rows, seg, start);
    const auto est = ekfAttitudes(rows, seg);
    const auto frames = simulateObservationsWithAttitude(truth, est, cam, err);

    size_t nstars = 0;
    for (const auto& f : frames) nstars += f.id.size();
    std::printf("\n  [2] POSITION, IDEAL STAR VECTORS  (%.1f stars/frame)\n",
                frames.empty() ? 0.0 : double(nstars) / frames.size());

    IterationConfig icfg;
    icfg.max_iterations = 6;
    for (AverageMethod m : {AverageMethod::Naive, AverageMethod::HeadingWeighted,
                            AverageMethod::CircleFit}) {
      icfg.method = m;
      const auto steps = iterateOrbit(frames, nominalCameraMount(), seg.centre,
                                      icfg, &true_mount);
      const char* nm = (m == AverageMethod::Naive) ? "naive"
                       : (m == AverageMethod::HeadingWeighted) ? "heading-wtd"
                                                               : "circle fit";
      std::printf("      %-14s iter1 %8.2f km   final %8.2f km   (%zu it)\n",
                  nm, steps.empty() ? -1 : steps.front().error_m / 1000.0,
                  steps.empty() ? -1 : steps.back().error_m / 1000.0,
                  steps.size());
    }

    // --- 2b. STAR-AIDED ATTITUDE -------------------------------------------
    //
    // Does correcting EKF3's attitude with an absolute star reference improve
    // the fix? This is the experiment the filter was written for, and until now
    // it had only ever been run against a synthetic gyro in a unit test.
    //
    // The star attitude carries NO POSITION. `kabsch(v_cam, ecef_dirs)` is
    // camera-to-ECEF, and star directions in ECEF do not depend on where the
    // observer is. So the SEQUENCE of them observes gyro bias -- a rate --
    // even though a single fix cannot separate tilt from position.
    //
    // Converting back to NED for the fix does need a position, but only weakly
    // and only constantly: over one orbit the aircraft moves a few hundred
    // metres, and even a 10 km error in the assumed position rotates the NED
    // frame by 0.09 deg CONSTANTLY, which is a body-fixed error and exactly
    // what the orbit average already removes.
    if (star_aided) {
      double gyro_mag = 0;
      for (const LogRow& r : sub) gyro_mag += r.gyro.norm();
      if (gyro_mag < 1e-9) {
        std::printf("\n  [2b] STAR-AIDED ATTITUDE: no gyro in the log "
                    "(re-run dump_log.py; needs IMU messages)\n");
      } else {
        StarAidedAttitude saa;
        StarField field(catalog, truth[truth.size() / 2].epoch);
        std::vector<Eigen::Vector3d> ecef;

        // NOT seg.centre. The star attitude is in ECEF and the fix needs it in
        // NED, and that conversion takes a position -- so using the TRUE centre
        // here would feed truth into the correction and the comparison would be
        // meaningless. It gave 0.01 km, which is what a circular test looks
        // like.
        //
        // Use the fix this orbit produces from EKF3's own attitude: an estimate
        // the system actually has. The conversion is weakly and CONSTANTLY
        // dependent on it -- a 10 km position error rotates the NED frame by
        // 0.09 deg, the same amount all orbit, which is a body-fixed error and
        // exactly what the orbit average removes.
        IterationConfig ref_cfg;
        ref_cfg.max_iterations = 6;
        ref_cfg.method = AverageMethod::HeadingWeighted;
        const auto ref_steps =
            iterateOrbit(frames, nominalCameraMount(), seg.centre, ref_cfg,
                         &true_mount);
        Geodetic ref_pos = frames.front().truth;   // fallback: departure prior
        if (!ref_steps.empty()) ref_pos = ref_steps.back().position;
        const Eigen::Matrix3d C_ecef_ned = nedBasisEcef(ref_pos);

        std::vector<Eigen::Matrix3d> corrected;
        corrected.reserve(frames.size());
        int n_upd = 0;
        double err_ekf = 0, err_saa = 0;
        int n_cmp = 0;

        for (size_t i = 0; i < frames.size(); ++i) {
          if (i > 0) {
            const double dt = frames[i].t - frames[i - 1].t;
            const size_t ri = std::min(seg.begin + i, rows.size() - 1);
            saa.predict(rows[ri].gyro, dt);
          }
          // Absolute star attitude, body -> ECEF, position-free.
          if (frames[i].id.size() >= 3) {
            field.ecefDirections(frames[i].epoch, ecef);
            std::vector<Eigen::Vector3d> a, b;
            for (size_t k = 0; k < frames[i].id.size(); ++k) {
              const int j = field.indexOf(frames[i].id[k]);
              if (j < 0) continue;
              a.push_back(frames[i].v_cam[k]);
              b.push_back(ecef[size_t(j)]);
            }
            if (a.size() >= 3) {
              const Eigen::Matrix3d R = kabsch(a, b);
              saa.update(R * nominalCameraMount().transpose(),
                         60.0 * kArcsec2Rad);
              ++n_upd;
            }
          }
          // Back to NED for the fix.
          corrected.push_back(saa.initialised()
                                  ? C_ecef_ned.transpose() * saa.attitude()
                                  : frames[i].C_l_b_est);
          const size_t ri = std::min(seg.begin + i, rows.size() - 1);
          const Eigen::Matrix3d Ct = rows[ri].C_true();
          err_ekf += rotationAngleBetween(frames[i].C_l_b_est, Ct);
          err_saa += rotationAngleBetween(corrected.back(), Ct);
          ++n_cmp;
        }

        std::vector<FrameData> f2 = frames;
        for (size_t i = 0; i < f2.size(); ++i) {
          f2[i].C_l_b_est = corrected[i];
          f2[i].yaw_est = dcmToEuler(corrected[i]).z();
        }

        // --- DE-DRIFT: correct only the TIME-VARYING part -------------------
        //
        // Replacing EKF3's attitude wholesale (above) fails because converting
        // an absolute star attitude into NED needs a position, and that error
        // is a constant NED rotation the orbit cannot average away.
        //
        // RELATIVE rotations escape both problems. Write the star attitude as
        // S_i = T_i M (true body->ECEF times an unknown CONSTANT mounting) and
        // EKF3's as C_be,i = T_i B_i (times its body-frame error). Then
        //
        //     dS_i = S_i S_0^T   = T_i T_0^T            <- M cancels EXACTLY
        //     dE_i = C_i C_0^T   = T_i B_i B_0^T T_0^T
        //     D_i  = dE_i^T dS_i = T_0 (B_0 B_i^T) T_0^T
        //
        // which is pure DRIFT. The assumed position survives only as a
        // similarity transform of a small rotation, i.e. second order.
        //
        // Conjugating D_i back into the body frame gives the correction that
        // takes B_i to B_0. Subtracting its orbit MEAN freezes the error at the
        // mean instead of at frame zero -- lower variance, and either way what
        // is left is CONSTANT and body-fixed, which is exactly what the orbit
        // average removes.
        std::vector<Eigen::Matrix3d> corr(frames.size(),
                                          Eigen::Matrix3d::Identity());
        {
          std::vector<Eigen::Matrix3d> Sv(frames.size());
          std::vector<bool> ok(frames.size(), false);
          for (size_t i = 0; i < frames.size(); ++i) {
            if (frames[i].id.size() < 3) continue;
            field.ecefDirections(frames[i].epoch, ecef);
            std::vector<Eigen::Vector3d> a, b;
            for (size_t k = 0; k < frames[i].id.size(); ++k) {
              const int j = field.indexOf(frames[i].id[k]);
              if (j < 0) continue;
              a.push_back(frames[i].v_cam[k]);
              b.push_back(ecef[size_t(j)]);
            }
            if (a.size() < 3) continue;
            Sv[i] = kabsch(a, b) * nominalCameraMount().transpose();
            ok[i] = true;
          }
          size_t i0 = 0;
          while (i0 < ok.size() && !ok[i0]) ++i0;
          if (i0 < ok.size()) {
            // Want X_i with C_ekf,i X_i = C_true,i B_0: freeze the attitude
            // error at its frame-0 value, removing everything time-varying.
            // Since C_ekf,i = C_true,i B_i,
            //
            //     X_i = B_i^T B_0 = C_ekf,i^T (C_true,i C_true,0^T) C_ekf,0
            //
            // and the middle term -- the TRUE rotation between the two frames
            // -- is what the stars give, mount-free, as P^T (S_i S_0^T) P.
            //
            // Applying the correction in the body frame at frame i is the part
            // an earlier attempt got wrong: it conjugated into the body frame
            // at frame ZERO and applied it at frame i, which are different
            // frames once the aircraft has turned.
            const Eigen::Matrix3d S0 = Sv[i0];
            // CNS_ORACLE_P substitutes the TRUE position here. It is a
            // diagnostic and not a method -- it proves the correction is
            // correct and isolates what limits it:
            //
            //   reference position   residual scatter   fix
            //   estimated (honest)       0.0521 deg     5.79 km
            //   TRUE (oracle)            0.0011 deg     0.00 km
            //
            // The mechanism works. It simply inherits the position error it is
            // given, because P^T (S_i S_0^T) P conjugates a LARGE rotation --
            // up to 360 deg round the orbit -- by the position error, and
            // conjugating a large rotation by a small angle perturbs it to
            // FIRST order. (Second order only for small rotations; that was
            // the mistaken step.) 0.0521 deg x 111 km = 5.8 km, which is the
            // reference error that went in.
            const Eigen::Matrix3d P = getenv("CNS_ORACLE_P")
                                          ? nedBasisEcef(seg.centre)
                                          : C_ecef_ned;
            const Eigen::Matrix3d Ce0 = frames[i0].C_l_b_est;
            for (size_t i = 0; i < frames.size(); ++i) {
              if (!ok[i]) continue;
              const Eigen::Matrix3d dS = Sv[i] * S0.transpose();      // ECEF
              const Eigen::Matrix3d dTrue = P.transpose() * dS * P;   // NED
              corr[i] = frames[i].C_l_b_est.transpose() * dTrue * Ce0;
            }
            // Carry the last good correction across frames with too few stars,
            // so a thin frame does not snap back to the uncorrected attitude.
            Eigen::Matrix3d last = Eigen::Matrix3d::Identity();
            for (size_t i = 0; i < frames.size(); ++i) {
              if (ok[i]) last = corr[i];
              else corr[i] = last;
            }
          }
        }
        std::vector<FrameData> f3 = frames;
        double err_dd = 0;
        for (size_t i = 0; i < f3.size(); ++i) {
          f3[i].C_l_b_est = frames[i].C_l_b_est * corr[i];
          f3[i].yaw_est = dcmToEuler(f3[i].C_l_b_est).z();
          const size_t ri = std::min(seg.begin + i, rows.size() - 1);
          err_dd += rotationAngleBetween(f3[i].C_l_b_est, rows[ri].C_true());
        }

        std::printf("\n  [2b] STAR-AIDED ATTITUDE  (%d star updates, "
                    "gyro bias %.3f deg/hr)\n", n_upd,
                    saa.gyroBias().norm() * kRad2Deg * 3600.0);
        std::printf("      attitude error  EKF3 %.4f deg   star-aided %.4f deg\n",
                    n_cmp ? err_ekf / n_cmp * kRad2Deg : -1.0,
                    n_cmp ? err_saa / n_cmp * kRad2Deg : -1.0);
        // TILT only. Yaw rotates about the local vertical and does not move
        // the zenith, so it does not move the fix -- reporting total rotation
        // angle hides whether the thing that matters improved.
        auto tiltErr = [&](const std::vector<FrameData>& fs) {
          double acc = 0;
          int n = 0;
          for (size_t i = 0; i < fs.size(); ++i) {
            const size_t ri = std::min(seg.begin + i, rows.size() - 1);
            const Eigen::AngleAxisd aa(fs[i].C_l_b_est.transpose() *
                                       rows[ri].C_true());
            const Eigen::Vector3d d = aa.axis() * aa.angle();
            acc += std::hypot(d.x(), d.y());
            ++n;
          }
          return n ? acc / n * kRad2Deg : -1.0;
        };
        std::printf("      attitude error  de-drifted %.4f deg\n",
                    n_cmp ? err_dd / n_cmp * kRad2Deg : -1.0);
        std::printf("      TILT error      EKF3 %.4f   star-aided %.4f   "
                    "de-drifted %.4f deg\n",
                    tiltErr(frames), tiltErr(f2), tiltErr(f3));
        {   // DIAGNOSTIC: is the correction non-trivial, and is the residual
            // actually CONSTANT in the body frame (which is what the orbit
            // average needs)?
          double cmag = 0;
          std::vector<Eigen::Matrix3d> resid;
          for (size_t i = 0; i < f3.size(); ++i) {
            const Eigen::AngleAxisd aa(corr[i]);
            cmag += std::abs(aa.angle());
            const size_t ri = std::min(seg.begin + i, rows.size() - 1);
            resid.push_back(rows[ri].C_true().transpose() * f3[i].C_l_b_est);
          }
          const Eigen::Matrix3d rbar = averageRotation(resid);
          double scatter = 0;
          for (const auto& R : resid)
            scatter += rotationAngleBetween(R, rbar);
          std::printf("      [diag] mean |correction| %.4f deg, residual "
                      "scatter about its mean %.4f deg\n",
                      cmag / f3.size() * kRad2Deg,
                      scatter / resid.size() * kRad2Deg);
        }
        std::printf("      %-14s %10s %12s %12s\n", "", "EKF3", "star-aided",
                    "DE-DRIFTED");
        for (AverageMethod m : {AverageMethod::HeadingWeighted,
                                AverageMethod::CircleFit}) {
          icfg.method = m;
          const auto s0 = iterateOrbit(frames, nominalCameraMount(), seg.centre,
                                       icfg, &true_mount);
          const auto s1 = iterateOrbit(f2, nominalCameraMount(), seg.centre,
                                       icfg, &true_mount);
          const auto s2 = iterateOrbit(f3, nominalCameraMount(), seg.centre,
                                       icfg, &true_mount);
          const char* nm = (m == AverageMethod::HeadingWeighted)
                               ? "heading-wtd" : "circle fit";
          std::printf("      %-14s %7.2f km %9.2f km %9.2f km\n", nm,
                      s0.empty() ? -1 : s0.back().error_m / 1000.0,
                      s1.empty() ? -1 : s1.back().error_m / 1000.0,
                      s2.empty() ? -1 : s2.back().error_m / 1000.0);
        }
      }
    }

    // --- 3. full imaging chain (Milestone 2 path) --------------------------
    if (imaging) {
      std::printf("\n  [3] POSITION, FULL IMAGING CHAIN\n");
      const StarField field(catalog, truth[truth.size() / 2].epoch);
      SensorModel sensor;
      sensor.exposure_s = 0.05;
      sensor.blur_substeps = 16;
      DetectorConfig dcfg;
      dcfg.centroid = CentroidMethod::GaussianFit;
      MatcherConfig mcfg;

      // Sub-sample: rendering every log row is far too slow, and a full
      // revolution is what matters, not the frame rate.
      const int stride = std::max(1, int((seg.end - seg.begin) / 120));
      const double period = (truth.size() > 1) ? (truth[1].t - truth[0].t) : 0.1;

      std::vector<std::vector<Detection>> dets;
      std::vector<size_t> idx;
      for (size_t k = 0; k + 1 < truth.size(); k += stride) {
        const RenderedFrame rf =
            renderFrame(truth[k], truth[k + 1], period * stride, field, vmags,
                        cam, true_mount, sensor, static_cast<unsigned>(k));
        dets.push_back(detectStars(rf.image, dcfg));
        idx.push_back(k);
      }

      Eigen::Matrix3d C_b_c = nominalCameraMount();
      Geodetic assumed = seg.centre;
      double last = -1.0;
      int nit = 0;
      for (int it = 0; it < 4; ++it) {
        std::vector<FrameData> ff;
        for (size_t j = 0; j < dets.size(); ++j) {
          FrameData fd;
          const int m = matchDetections(dets[j], truth[idx[j]].epoch,
                                        est[idx[j]], C_b_c, cam, field, assumed,
                                        mcfg, fd);
          fd.t = truth[idx[j]].t;
          fd.yaw_est = dcmToEuler(est[idx[j]]).z();
          fd.truth = truth[idx[j]].pos;
          if (m >= 3) ff.push_back(fd);
        }
        if (ff.size() < 15) break;
        const OrbitResult r =
            estimateOrbit(ff, C_b_c, AverageMethod::HeadingWeighted);
        if (!r.ok) break;
        assumed = r.position;
        last = haversine(r.position, seg.centre);
        nit = it + 1;
        C_b_c = recalibrateMount(ff, r.position);
      }
      std::printf("      %-14s final %8.2f km   (%d it, %zu frames)\n",
                  "heading-wtd", last / 1000.0, nit, dets.size());
      std::printf("\n      Difference against [2] is the cost of the imaging\n"
                  "      chain, with the attitude held identical.\n");
    }
  }
  return 0;
}
