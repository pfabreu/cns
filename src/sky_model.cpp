#include "celestial/sky_model.hpp"

#include <erfa.h>
#include <algorithm>
#include <erfam.h>  // ERFA_DJ00, ERFA_DJY and friends live here

#include <cmath>
#include <unordered_map>

#include <numeric>

#include <random>

namespace celestial {

std::optional<Epoch> Epoch::fromUtc(int year, int month, int day, int hour,
                                    int minute, double second, double dut1) {
  Epoch e;
  if (eraDtf2d("UTC", year, month, day, hour, minute, second, &e.utc1_,
               &e.utc2_) != 0) {
    return std::nullopt;
  }
  double tai1, tai2;
  if (eraUtctai(e.utc1_, e.utc2_, &tai1, &tai2) != 0) return std::nullopt;
  if (eraTaitt(tai1, tai2, &e.tt1_, &e.tt2_) != 0) return std::nullopt;
  if (eraUtcut1(e.utc1_, e.utc2_, dut1, &e.ut11_, &e.ut12_) != 0) {
    return std::nullopt;
  }
  e.dut1_ = dut1;
  return e;
}

Epoch Epoch::advanced(double seconds) const {
  Epoch e;
  e.utc1_ = utc1_;
  e.utc2_ = utc2_ + seconds / 86400.0;
  e.dut1_ = dut1_;
  double tai1, tai2;
  if (eraUtctai(e.utc1_, e.utc2_, &tai1, &tai2) == 0) {
    eraTaitt(tai1, tai2, &e.tt1_, &e.tt2_);
  }
  eraUtcut1(e.utc1_, e.utc2_, e.dut1_, &e.ut11_, &e.ut12_);
  return e;
}

double Epoch::earthRotationAngle() const { return eraEra00(ut11_, ut12_); }

double Epoch::julianYearsSinceJ2000() const {
  return ((tt1_ - ERFA_DJ00) + tt2_) / ERFA_DJY;
}

Atmosphere Atmosphere::isa(double altitude_m, double temperature_c) {
  // Simple ISA troposphere. Only used for refraction, which is a
  // sub-arcminute effect above ~30 deg elevation, so this is plenty.
  const double p = 1013.25 * std::pow(1.0 - 2.25577e-5 * altitude_m, 5.25588);
  return Atmosphere{p, temperature_c, 0.5, 0.55};
}

Observation observe(const CatalogStar& star, const Geodetic& site,
                    const Epoch& epoch, const Atmosphere& atmos) {
  double aob, zob, hob, dob, rob, eo;
  eraAtco13(star.ra, star.dec, star.pm_ra_star, star.pm_dec,
            /*px arcsec=*/0.0, /*rv km/s=*/0.0,
            epoch.utc1(), epoch.utc2(), epoch.dut1(),
            site.lon, site.lat, site.alt,
            /*xp=*/0.0, /*yp=*/0.0,
            atmos.pressure_hpa, atmos.temperature_c, atmos.relative_humidity,
            atmos.wavelength_um,
            &aob, &zob, &hob, &dob, &rob, &eo);

  Observation o;
  o.id = star.id;
  o.azimuth = aob;
  o.zenith_angle = zob;
  return o;
}

std::vector<Observation> observeVisible(const std::vector<CatalogStar>& catalog,
                                        const Geodetic& site,
                                        const Epoch& epoch,
                                        double min_elevation,
                                        const Atmosphere& atmos) {
  std::vector<Observation> out;
  out.reserve(catalog.size());
  for (const CatalogStar& s : catalog) {
    Observation o = observe(s, site, epoch, atmos);
    if (o.elevation() >= min_elevation) out.push_back(o);
  }
  return out;
}

Eigen::Vector3d subStellarPoint(const CatalogStar& star, const Epoch& epoch) {
  // ICRS -> CIRS apparent place: proper motion, parallax, light deflection,
  // annual aberration, precession-nutation. NOT diurnal aberration or
  // refraction, both of which are properties of the observer, not the star.
  double ri, di, eo;
  eraAtci13(star.ra, star.dec, star.pm_ra_star, star.pm_dec,
            /*px=*/0.0, /*rv=*/0.0, epoch.tt1(), epoch.tt2(), &ri, &di, &eo);

  // CIO-based hour angle: h = ERA - ri, so the sub-stellar longitude
  // (east positive) is ri - ERA. This is the CIO form of the paper's
  // lambda = alpha - GHA_Aries.
  const double lon = eraAnpm(ri - epoch.earthRotationAngle());
  const double lat = di;

  return Eigen::Vector3d(std::cos(lon) * std::cos(lat),
                         std::sin(lon) * std::cos(lat), std::sin(lat));
}

StarField::StarField(const std::vector<CatalogStar>& catalog,
                     const Epoch& reference) {
  ids_.reserve(catalog.size());
  vmag_.reserve(catalog.size());
  cirs_.reserve(catalog.size());

  // One astrometry context for the whole field: precession-nutation, annual
  // aberration and light deflection are all folded in here, once.
  eraASTROM astrom;
  double eo;
  eraApci13(reference.tt1(), reference.tt2(), &astrom, &eo);

  for (const CatalogStar& s : catalog) {
    double ri, di;
    eraAtciq(s.ra, s.dec, s.pm_ra_star, s.pm_dec, 0.0, 0.0, &astrom, &ri, &di);
    ids_.push_back(s.id);
    vmag_.push_back(s.vmag);
    cirs_.emplace_back(std::cos(di) * std::cos(ri), std::cos(di) * std::sin(ri),
                       std::sin(di));
  }
}

Eigen::Vector3d StarField::ecefDirection(size_t i, const Epoch& epoch) const {
  const double era = epoch.earthRotationAngle();
  const double c = std::cos(era), s = std::sin(era);
  const Eigen::Vector3d& v = cirs_[i];
  // Rz(-ERA)
  return Eigen::Vector3d(c * v.x() + s * v.y(), -s * v.x() + c * v.y(), v.z());
}

void StarField::ecefDirections(const Epoch& epoch,
                               std::vector<Eigen::Vector3d>& out) const {
  const double era = epoch.earthRotationAngle();
  const double c = std::cos(era), s = std::sin(era);
  out.resize(cirs_.size());
  for (size_t i = 0; i < cirs_.size(); ++i) {
    const Eigen::Vector3d& v = cirs_[i];
    out[i] = Eigen::Vector3d(c * v.x() + s * v.y(), -s * v.x() + c * v.y(),
                             v.z());
  }
}

int StarField::indexOf(int id) const {
  for (size_t i = 0; i < ids_.size(); ++i) {
    if (ids_[i] == id) return static_cast<int>(i);
  }
  return -1;
}

Eigen::Matrix3d nedBasisEcef(const Geodetic& p) {
  const double sp = std::sin(p.lat), cp = std::cos(p.lat);
  const double sl = std::sin(p.lon), cl = std::cos(p.lon);
  Eigen::Matrix3d C;
  C.col(0) = Eigen::Vector3d(-sp * cl, -sp * sl, cp);  // north
  C.col(1) = Eigen::Vector3d(-sl, cl, 0.0);            // east
  C.col(2) = Eigen::Vector3d(-cp * cl, -cp * sl, -sp); // down
  return C;
}

std::vector<StarSight> makeSights(const std::vector<Observation>& obs,
                                  const std::vector<CatalogStar>& catalog,
                                  const Epoch& epoch) {
  std::unordered_map<int, const CatalogStar*> by_hip;
  by_hip.reserve(catalog.size());
  for (const CatalogStar& s : catalog) by_hip[s.id] = &s;

  std::vector<StarSight> sights;
  sights.reserve(obs.size());
  for (const Observation& o : obs) {
    auto it = by_hip.find(o.id);
    if (it == by_hip.end()) continue;  // unidentified detection
    StarSight s;
    s.id = o.id;
    s.gp = subStellarPoint(*it->second, epoch);
    s.zenith_angle = o.zenith_angle;
    s.weight = 1.0;
    sights.push_back(s);
  }
  return sights;
}

namespace refraction {
namespace {

/// Pressure and temperature scaling, common to both models.
/// Reference conditions are 1010 hPa and 10 C.
double ptScale(double pressure_hpa, double temperature_c) {
  return (pressure_hpa / 1010.0) * (283.0 / (273.0 + temperature_c));
}

constexpr double kArcmin2Rad = M_PI / (180.0 * 60.0);

}  // namespace

double bennett(double apparent_el, double pressure_hpa, double temperature_c) {
  if (pressure_hpa <= 0.0) return 0.0;
  const double el_deg = apparent_el * kRad2Deg;
  // Below the horizon the formula is meaningless; clamp rather than produce
  // nonsense, and let the caller's elevation cut do the real filtering.
  if (el_deg < -1.0) return 0.0;
  const double arg = (el_deg + 7.31 / (el_deg + 4.4)) * kDeg2Rad;
  const double r_arcmin = 1.0 / std::tan(arg);
  return r_arcmin * kArcmin2Rad * ptScale(pressure_hpa, temperature_c);
}

double saemundsson(double true_el, double pressure_hpa, double temperature_c) {
  if (pressure_hpa <= 0.0) return 0.0;
  const double el_deg = true_el * kRad2Deg;
  if (el_deg < -1.0) return 0.0;
  const double arg = (el_deg + 10.3 / (el_deg + 5.11)) * kDeg2Rad;
  const double r_arcmin = 1.02 / std::tan(arg);
  return r_arcmin * kArcmin2Rad * ptScale(pressure_hpa, temperature_c);
}

}  // namespace refraction

void correctRefraction(std::vector<StarSight>& sights,
                       const Atmosphere& atmos) {
  if (atmos.pressure_hpa <= 0.0) return;
  for (StarSight& s : sights) {
    const double apparent_el = M_PI_2 - s.zenith_angle;
    const double r = refraction::bennett(apparent_el, atmos.pressure_hpa,
                                         atmos.temperature_c);
    // The star appeared too high, so the true zenith angle is LARGER.
    s.zenith_angle += r;
  }
}

void weightForRefraction(std::vector<StarSight>& sights,
                         const Atmosphere& atmos, double centroid_sigma_rad,
                         double model_error_frac) {
  for (StarSight& s : sights) {
    // Use the (already corrected or not) zenith angle to estimate elevation.
    // The difference between apparent and true elevation is at most a few
    // arcminutes, which is irrelevant for choosing a weight.
    const double el = M_PI_2 - s.zenith_angle;
    const double r = (atmos.pressure_hpa > 0.0)
                         ? refraction::bennett(el, atmos.pressure_hpa,
                                               atmos.temperature_c)
                         : 0.0;
    const double model_sigma = model_error_frac * r;
    const double var = centroid_sigma_rad * centroid_sigma_rad +
                       model_sigma * model_sigma;
    s.weight = (var > 0.0) ? 1.0 / var : 1.0;
  }

  // Normalise so the mean weight is 1. Purely cosmetic — it keeps the
  // reported residual_rms comparable across weighting schemes.
  double sum = 0.0;
  for (const StarSight& s : sights) sum += s.weight;
  if (sum > 0.0 && !sights.empty()) {
    const double scale = static_cast<double>(sights.size()) / sum;
    for (StarSight& s : sights) s.weight *= scale;
  }
}



// ===========================================================================
// POSITION SOLVER -- folded in from position_solver.hpp. It operates on the
// StarSight values this header defines, so the two were always one module.
// ===========================================================================

namespace {

/// Core weighted least-squares solve of A x = p. Returns false if the geometry
/// is degenerate (e.g. all stars in a line, or fewer than 3 sights).
bool leastSquares(const std::vector<StarSight>& sights,
                  const std::vector<int>& idx, Eigen::Vector3d& x_out,
                  double& residual_rms, double& condition) {
  const int n = static_cast<int>(idx.size());
  if (n < 3) return false;

  Eigen::MatrixXd A(n, 3);
  Eigen::VectorXd p(n);
  Eigen::VectorXd w(n);
  for (int i = 0; i < n; ++i) {
    const StarSight& s = sights[idx[i]];
    A.row(i) = s.gp.transpose();
    p(i) = std::cos(s.zenith_angle);
    w(i) = s.weight;
  }

  // Weighted normal equations. Solved by SVD rather than by forming
  // (A^T W A)^-1 explicitly: same answer, far better behaved when the stars
  // are clustered, which they always are in a 53 deg field.
  const Eigen::MatrixXd Aw = w.cwiseSqrt().asDiagonal() * A;
  const Eigen::VectorXd pw = w.cwiseSqrt().asDiagonal() * p;

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Aw, Eigen::ComputeThinU |
                                                Eigen::ComputeThinV);
  const auto& sv = svd.singularValues();
  if (sv(2) < 1e-12) return false;
  condition = (sv(0) / sv(2)) * (sv(0) / sv(2));  // cond(A^T A) = cond(A)^2

  const Eigen::Vector3d x = svd.solve(pw);
  const Eigen::VectorXd r = A * x - p;
  residual_rms = std::sqrt(r.squaredNorm() / n);
  x_out = x;
  return true;
}

Fix finalise(const Eigen::Vector3d& x, double residual_rms, double condition,
             int n_used) {
  Fix f;
  const double norm = x.norm();
  if (norm < 1e-9) return f;  // ok stays false

  // NOTE: x is the observer's ZENITH direction, so asin of its Z component is
  // GEODETIC latitude directly. See the header comment. Do not apply a
  // geocentric-to-geodetic conversion here.
  const Eigen::Vector3d z = x / norm;
  f.ok = true;
  f.zenith_ecef = z;
  f.lat = std::asin(std::clamp(z.z(), -1.0, 1.0));
  f.lon = std::atan2(z.y(), z.x());
  f.residual_rms = residual_rms;
  f.condition = condition;
  f.n_used = n_used;
  return f;
}

}  // namespace

Fix solveFix(const std::vector<StarSight>& sights) {
  std::vector<int> idx(sights.size());
  std::iota(idx.begin(), idx.end(), 0);

  Eigen::Vector3d x;
  double rms = 0.0, cond = 0.0;
  if (!leastSquares(sights, idx, x, rms, cond)) return Fix{};
  return finalise(x, rms, cond, static_cast<int>(sights.size()));
}

Fix solveFixRansac(const std::vector<StarSight>& sights,
                   const RansacConfig& cfg) {
  const int n = static_cast<int>(sights.size());
  if (n < 6) return solveFix(sights);

  // The inlier test compares |p_k - a_k . x| against a tolerance. Both sides
  // are cosines of angles near the measurement, so converting an angular
  // tolerance to a cosine residual needs a scale factor: near zero,
  // d(cos zeta) = -sin(zeta) d(zeta). sin(zeta) varies across the field, so a
  // flat cosine tolerance is tight near the zenith and loose near the horizon
  // -- and a camera pointed up sees mostly near-zenith stars, i.e. the tight
  // end. Normalising makes `tolerance` mean the angle it claims to.
  //
  // The floor of 0.2 (zenith angle ~11.5 deg) bounds the amplification; without
  // it a star directly overhead divides by ~0. Worth 0.368 -> 0.305 km at a 40%
  // outlier rate and neutral below ~25%, so this is about the tolerance being
  // honest rather than about the numbers.
  const double tol = cfg.tolerance;

  std::mt19937 rng(cfg.seed);
  std::uniform_int_distribution<int> pick(0, n - 1);

  double best_error = std::numeric_limits<double>::infinity();
  std::vector<int> best_inliers;

  for (int iter = 0; iter < cfg.iterations; ++iter) {
    // Three distinct random sights.
    int a = pick(rng), b = pick(rng), c = pick(rng);
    if (a == b || b == c || a == c) continue;
    const std::vector<int> sample{a, b, c};

    Eigen::Vector3d x;
    double rms = 0.0, cond = 0.0;
    if (!leastSquares(sights, sample, x, rms, cond)) continue;
    if (x.norm() < 1e-9) continue;
    const Eigen::Vector3d xn = x / x.norm();

    std::vector<int> inliers;
    double inlier_error = 0.0;
    for (int k = 0; k < n; ++k) {
      const double e = std::abs(std::cos(sights[k].zenith_angle) -
                                sights[k].gp.dot(xn)) /
                       std::max(0.2, std::sin(sights[k].zenith_angle));
      if (e < tol) {
        inliers.push_back(k);
        inlier_error += e;
      }
    }

    if (static_cast<double>(inliers.size()) / n < cfg.min_inlier_ratio) continue;
    const double mean_error = inlier_error / inliers.size();
    if (mean_error < best_error) {
      best_error = mean_error;
      best_inliers = std::move(inliers);
    }
  }

  if (best_inliers.size() < 3) return Fix{};

  // Refit on the full inlier set — the 3-star model that won is only a seed.
  Eigen::Vector3d x;
  double rms = 0.0, cond = 0.0;
  if (!leastSquares(sights, best_inliers, x, rms, cond)) return Fix{};
  return finalise(x, rms, cond, static_cast<int>(best_inliers.size()));
}

double haversine(const Geodetic& a, const Geodetic& b) {
  // WGS-84 mean radius (IUGG): (2a + b) / 3.
  constexpr double kR = (2.0 * kWgs84A + kWgs84A * (1.0 - kWgs84F)) / 3.0;
  const double dlat = b.lat - a.lat;
  const double dlon = b.lon - a.lon;
  const double h = std::sin(dlat / 2) * std::sin(dlat / 2) +
                   std::cos(a.lat) * std::cos(b.lat) * std::sin(dlon / 2) *
                       std::sin(dlon / 2);
  return 2.0 * kR * std::asin(std::sqrt(std::clamp(h, 0.0, 1.0)));
}

}  // namespace celestial
