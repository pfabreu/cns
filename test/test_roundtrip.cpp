// test_roundtrip.cpp — MILESTONE 0.
//
// Generate the star sights that SHOULD be seen from a known position and time,
// feed them to the position solver, and check we get the position back.
//
// If this test does not pass to well under 100 m, nothing built on top of it
// can possibly work. No images, no ROS, no simulator — just the geometry.

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "celestial/sky_model.hpp"
#include "celestial/star_catalog.hpp"
#include "celestial/orbit.hpp"

using namespace celestial;

namespace {

int g_failures = 0;

void check(bool cond, const std::string& what) {
  if (!cond) {
    std::printf("  ** FAIL: %s\n", what.c_str());
    ++g_failures;
  }
}

struct Scenario {
  const char* name;
  Geodetic site;
  int y, mo, d, h, mi;
  double s;
};

// The first scenario is roughly the paper's flight: outback South Australia,
// 800 m AGL, shortly after astronomical twilight. The others check that
// nothing depends on hemisphere, longitude sign, or being far from a pole.
// Woomera-ish, 800 m, as in the paper's flight -- used by the refraction tests.
const Geodetic kSite = Geodetic::fromDegrees(-30.89, 136.56, 800.0);

const std::vector<Scenario> kScenarios = {
    {"Woomera SA (paper-like)", Geodetic::fromDegrees(-30.89, 136.56, 800.0),
     2024, 9, 15, 10, 30, 0.0},
    {"North Sea",               Geodetic::fromDegrees(56.20, 3.10, 500.0),
     2024, 1, 20, 22, 15, 30.0},
    {"Mid-Pacific",             Geodetic::fromDegrees(-8.50, -150.75, 1200.0),
     2025, 6, 3, 8, 5, 12.5},
    {"High latitude (Tromso)",  Geodetic::fromDegrees(69.65, 18.95, 300.0),
     2024, 11, 11, 1, 45, 0.0},
    {"Equator / prime meridian",Geodetic::fromDegrees(0.10, 0.20, 50.0),
     2026, 3, 1, 3, 0, 0.0},
};

/// Apply a common-mode tilt of the estimated vertical to a set of sights.
/// tau_n / tau_e are the NED-frame tilt components in radians.
///
/// Derivation: with zenith z, star s at zenith angle zeta and azimuth az,
///   cos(zeta') = (z + tau_n*n + tau_e*e) . s = cos(zeta) + sin(zeta)*(...)
/// so to first order  d(zeta) = -(tau_n*cos(az) + tau_e*sin(az)).
///
/// This is the dominant error term in the real system: it is what an AHRS
/// pitch/roll bias or a camera boresight misalignment looks like.
void applyTilt(std::vector<StarSight>& sights,
               const std::vector<Observation>& obs, double tau_n,
               double tau_e) {
  for (size_t i = 0; i < sights.size(); ++i) {
    const double az = obs[i].azimuth;
    sights[i].zenith_angle -= (tau_n * std::cos(az) + tau_e * std::sin(az));
  }
}

double fixError(const Fix& f, const Geodetic& truth);

void reportFix(const char* label, const Fix& fix, const Geodetic& truth) {
  if (!fix.ok) {
    std::printf("  %-28s SOLVER FAILED\n", label);
    ++g_failures;
    return;
  }
  const Geodetic est{fix.lat, fix.lon, truth.alt};
  const double err = haversine(est, truth);
  std::printf("  %-28s n=%2d  err=%10.2f m  rms=%.2e  cond=%7.1f\n", label,
              fix.n_used, err, fix.residual_rms, fix.condition);
}

double fixError(const Fix& fix, const Geodetic& truth) {
  if (!fix.ok) return std::numeric_limits<double>::infinity();
  return haversine(Geodetic{fix.lat, fix.lon, truth.alt}, truth);
}

/// Initial great-circle bearing from a to b, degrees from North through East.
double bearingDeg(const Geodetic& a, const Geodetic& b) {
  const double dlon = b.lon - a.lon;
  const double y = std::sin(dlon) * std::cos(b.lat);
  const double x = std::cos(a.lat) * std::sin(b.lat) -
                   std::sin(a.lat) * std::cos(b.lat) * std::cos(dlon);
  double brg = std::atan2(y, x) * kRad2Deg;
  if (brg < 0) brg += 360.0;
  return brg;
}

// ---------------------------------------------------------------------------

void testCleanRoundTrip() {
  std::printf("\n=== TEST 1: clean round trip (no refraction, no noise) ===\n");
  const auto catalog = catalogBrighterThan(3.0);

  for (const Scenario& sc : kScenarios) {
    auto epoch = Epoch::fromUtc(sc.y, sc.mo, sc.d, sc.h, sc.mi, sc.s);
    check(epoch.has_value(), std::string("epoch parse for ") + sc.name);
    if (!epoch) continue;

    // 15 deg elevation cut: below that refraction and the horizon get ugly.
    const auto obs = observeVisible(catalog, sc.site, *epoch,
                                    15.0 * kDeg2Rad, Atmosphere::none());
    const auto sights = makeSights(obs, catalog, *epoch);

    std::printf("\n %s  (%zu stars above 15 deg)\n", sc.name, sights.size());
    check(sights.size() >= 4,
          std::string("enough stars visible at ") + sc.name);
    if (sights.size() < 4) continue;

    const Fix fix = solveFix(sights);
    reportFix("least squares", fix, sc.site);
    check(fixError(fix, sc.site) < 100.0,
          std::string("round trip < 100 m at ") + sc.name);
  }
}

void testTiltSensitivity() {
  std::printf("\n=== TEST 2: tilt sensitivity (expect ~111 km per degree) ===\n");
  const Scenario& sc = kScenarios[0];
  auto epoch = Epoch::fromUtc(sc.y, sc.mo, sc.d, sc.h, sc.mi, sc.s);
  if (!epoch) return;

  const auto catalog = catalogBrighterThan(3.0);
  const auto obs =
      observeVisible(catalog, sc.site, *epoch, 15.0 * kDeg2Rad, Atmosphere::none());
  const auto clean = makeSights(obs, catalog, *epoch);

  std::printf("  %-12s %-14s %-14s\n", "tilt (deg)", "error (km)", "km per deg");
  for (double tilt_deg : {0.05, 0.1, 0.2, 0.5, 1.0}) {
    auto sights = clean;
    applyTilt(sights, obs, tilt_deg * kDeg2Rad, 0.0);
    const Fix fix = solveFix(sights);
    const double err_km = fixError(fix, sc.site) / 1000.0;
    std::printf("  %-12.2f %-14.2f %-14.1f\n", tilt_deg, err_km,
                err_km / tilt_deg);
    // The paper's rule of thumb: 1 deg of attitude error ~ 100 km of position.
    check(std::abs(err_km / tilt_deg - 111.2) < 5.0,
          "tilt-to-position scale is ~111 km/deg");
  }
  std::printf("\n  This is the entire problem the paper solves: a tilt bias is\n"
              "  indistinguishable from a position offset in a single frame.\n");
}

void testGeodeticTrap() {
  std::printf("\n=== TEST 3: the geodetic/geocentric trap ===\n");
  const Scenario& sc = kScenarios[0];
  auto epoch = Epoch::fromUtc(sc.y, sc.mo, sc.d, sc.h, sc.mi, sc.s);
  if (!epoch) return;

  const auto catalog = catalogBrighterThan(3.0);
  const auto obs =
      observeVisible(catalog, sc.site, *epoch, 15.0 * kDeg2Rad, Atmosphere::none());
  const auto sights = makeSights(obs, catalog, *epoch);
  const Fix fix = solveFix(sights);
  if (!fix.ok) return;

  // What the solver correctly returns:
  const double correct_err = fixError(fix, sc.site);

  // What you get if you assume the solve returned a GEOCENTRIC latitude and
  // "helpfully" convert it to geodetic. This is wrong: x is the zenith
  // direction, so it is already geodetic.
  const double e2 = kWgs84F * (2.0 - kWgs84F);
  const double bogus_lat = std::atan(std::tan(fix.lat) / (1.0 - e2));
  const double bogus_err =
      haversine(Geodetic{bogus_lat, fix.lon, sc.site.alt}, sc.site);

  std::printf("  correct (no conversion) : %10.2f m\n", correct_err);
  std::printf("  after bogus geo->geodetic: %10.2f m\n", bogus_err);
  std::printf("\n  A 'harmless looking' coordinate conversion costs ~%.0f km.\n",
              bogus_err / 1000.0);
  check(correct_err < 100.0, "correct interpretation is accurate");
  check(bogus_err > 10000.0, "bogus conversion is demonstrably large");
}

void testRansacWithOutliers() {
  std::printf("\n=== TEST 4: RANSAC rejects false detections ===\n");
  const Scenario& sc = kScenarios[0];
  auto epoch = Epoch::fromUtc(sc.y, sc.mo, sc.d, sc.h, sc.mi, sc.s);
  if (!epoch) return;

  const auto catalog = catalogBrighterThan(3.0);
  const auto obs =
      observeVisible(catalog, sc.site, *epoch, 15.0 * kDeg2Rad, Atmosphere::none());
  auto sights = makeSights(obs, catalog, *epoch);
  const size_t n_good = sights.size();

  // Corrupt three sights the way a satellite streak or a misidentified star
  // would: plausible geometry, badly wrong zenith angle.
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> bad(-0.3, 0.3);
  for (int k = 0; k < 3 && k < static_cast<int>(sights.size()); ++k) {
    sights[k * 3 % sights.size()].zenith_angle += bad(rng);
  }

  const Fix plain = solveFix(sights);
  RansacConfig cfg;
  cfg.iterations = 400;
  const Fix robust = solveFixRansac(sights, cfg);

  std::printf("  %zu sights, 3 corrupted\n", n_good);
  reportFix("plain least squares", plain, sc.site);
  reportFix("RANSAC", robust, sc.site);
  check(fixError(robust, sc.site) < fixError(plain, sc.site),
        "RANSAC beats plain least squares with outliers present");
  check(fixError(robust, sc.site) < 1000.0, "RANSAC recovers a good fix");
}

void testNoiseAndRefraction() {
  std::printf("\n=== TEST 5: measurement noise and unmodelled refraction ===\n");
  const Scenario& sc = kScenarios[0];
  auto epoch = Epoch::fromUtc(sc.y, sc.mo, sc.d, sc.h, sc.mi, sc.s);
  if (!epoch) return;
  const auto catalog = catalogBrighterThan(3.0);

  // Centroiding noise. The paper's optics are ~100 arcsec/pixel; good
  // subpixel centroiding gives maybe 0.1 px, so ~10 arcsec is realistic.
  std::mt19937 rng(7);
  for (double sigma_arcsec : {1.0, 10.0, 60.0}) {
    const auto obs = observeVisible(catalog, sc.site, *epoch, 15.0 * kDeg2Rad,
                                    Atmosphere::none());
    auto sights = makeSights(obs, catalog, *epoch);
    std::normal_distribution<double> n(0.0, sigma_arcsec * kArcsec2Rad);
    for (auto& s : sights) s.zenith_angle += n(rng);
    const Fix fix = solveFix(sights);
    std::printf("  noise %5.1f arcsec  ->  %8.1f m\n", sigma_arcsec,
                fixError(fix, sc.site));
  }

  // Now: observations WITH refraction, solver WITHOUT a refraction model.
  // This is what happens if you forget it, and it is a bias, not noise.
  const auto obs_ref = observeVisible(catalog, sc.site, *epoch, 15.0 * kDeg2Rad,
                                      Atmosphere::isa(sc.site.alt, 10.0));
  const auto sights_ref = makeSights(obs_ref, catalog, *epoch);
  const Fix fix_ref = solveFix(sights_ref);
  std::printf("\n  unmodelled refraction (15 deg cut) -> %8.1f m\n",
              fixError(fix_ref, sc.site));

  const auto obs_ref40 = observeVisible(catalog, sc.site, *epoch,
                                        40.0 * kDeg2Rad,
                                        Atmosphere::isa(sc.site.alt, 10.0));
  const auto sights_ref40 = makeSights(obs_ref40, catalog, *epoch);
  if (sights_ref40.size() >= 4) {
    const Fix f40 = solveFix(sights_ref40);
    std::printf("  unmodelled refraction (40 deg cut) -> %8.1f m  (%zu stars)\n",
                fixError(f40, sc.site), sights_ref40.size());
  }
  std::printf("\n  Refraction is a systematic error: raising the elevation cut\n"
              "  reduces it but costs you stars and weakens the geometry.\n");
}

void testModelAgainstErfa() {
  std::printf("\n=== TEST 1: Bennett vs ERFA's refraction model ===\n");
  auto epoch = Epoch::fromUtc(2024, 9, 15, 10, 30, 0.0);
  if (!epoch) return;

  const Atmosphere atmos = Atmosphere::isa(kSite.alt, 10.0);
  std::printf("  atmosphere: %.1f hPa, %.1f C  (ISA at %.0f m)\n",
              atmos.pressure_hpa, atmos.temperature_c, kSite.alt);
  std::printf("\n  %-10s %-12s %-12s %-12s %-10s\n", "HR", "el (deg)",
              "ERFA (\")", "Bennett (\")", "diff (\")");

  const auto catalog = catalogBrighterThan(3.0);
  int shown = 0;
  double worst = 0.0;
  for (const CatalogStar& s : catalog) {
    const Observation dry = observe(s, kSite, *epoch, Atmosphere::none());
    if (dry.elevation() < 15.0 * kDeg2Rad) continue;
    const Observation wet = observe(s, kSite, *epoch, atmos);

    // ERFA's refraction, from the change in zenith angle.
    const double r_erfa = dry.zenith_angle - wet.zenith_angle;
    // Bennett, evaluated at the APPARENT elevation, as it would be in flight.
    const double r_bennett = refraction::bennett(
        wet.elevation(), atmos.pressure_hpa, atmos.temperature_c);
    const double diff = (r_bennett - r_erfa) / kArcsec2Rad;
    worst = std::max(worst, std::abs(diff));

    if (shown++ < 12) {
      std::printf("  %-10s %-12.2f %-12.2f %-12.2f %-+10.3f\n", std::to_string(s.id).c_str(),
                  dry.elevation() * kRad2Deg, r_erfa / kArcsec2Rad,
                  r_bennett / kArcsec2Rad, diff);
    }
  }
  std::printf("\n  worst |Bennett - ERFA| above 15 deg: %.3f arcsec"
              "  (~%.0f m of position)\n", worst, worst * kArcsec2Rad * 6371000.0);
  check(worst < 5.0, "Bennett agrees with ERFA to better than 5 arcsec");
}

struct AblationRow {
  double cut_deg;
  size_t n;
  double none;         // no refraction anywhere (geometry floor)
  double uncorrected;  // refracted observations, no correction
  double bennett;      // Bennett correction, unit weights
  double weighted;     // Bennett correction + variance weighting
};

AblationRow runAblation(const Epoch& epoch, double cut_deg,
                        const Atmosphere& atmos) {
  const auto catalog = catalogBrighterThan(3.0);
  const double cut = cut_deg * kDeg2Rad;
  AblationRow row{cut_deg, 0, 0, 0, 0, 0};

  // Geometry floor: no refraction in the forward model at all.
  {
    const auto obs = observeVisible(catalog, kSite, epoch, cut, Atmosphere::none());
    const auto sights = makeSights(obs, catalog, epoch);
    row.n = sights.size();
    if (sights.size() < 3) return row;
    row.none = fixError(solveFix(sights), kSite);
  }

  // Refracted observations from here on.
  const auto obs = observeVisible(catalog, kSite, epoch, cut, atmos);
  const auto base = makeSights(obs, catalog, epoch);
  if (base.size() < 3) return row;

  row.uncorrected = fixError(solveFix(base), kSite);

  auto corrected = base;
  correctRefraction(corrected, atmos);
  row.bennett = fixError(solveFix(corrected), kSite);

  auto weighted = base;
  correctRefraction(weighted, atmos);
  weightForRefraction(weighted, atmos);
  row.weighted = fixError(solveFix(weighted), kSite);

  return row;
}

void testAblation() {
  std::printf("\n=== TEST 2: ablation vs elevation cut ===\n");
  auto epoch = Epoch::fromUtc(2024, 9, 15, 10, 30, 0.0);
  if (!epoch) return;
  const Atmosphere atmos = Atmosphere::isa(kSite.alt, 10.0);

  std::printf("\n  %-6s %-5s %-12s %-14s %-14s %-14s\n", "cut", "n",
              "no refr (m)", "uncorr (m)", "Bennett (m)", "+weights (m)");
  for (double cut : {10.0, 15.0, 20.0, 30.0, 40.0}) {
    const AblationRow r = runAblation(*epoch, cut, atmos);
    if (r.n < 3) {
      std::printf("  %-6.0f %-5zu  (too few stars)\n", cut, r.n);
      continue;
    }
    std::printf("  %-6.0f %-5zu %-12.1f %-14.1f %-14.1f %-14.1f\n", cut, r.n,
                r.none, r.uncorrected, r.bennett, r.weighted);
  }

  const AblationRow r15 = runAblation(*epoch, 15.0, atmos);
  check(r15.bennett < r15.uncorrected / 5.0,
        "Bennett correction improves the fix by at least 5x");
  std::printf("\n  At a 15 deg cut the correction takes %.0f m -> %.0f m,\n"
              "  a factor of %.1f, for about forty lines of code.\n",
              r15.uncorrected, r15.bennett, r15.uncorrected / r15.bennett);
  std::printf("\n  Note the uncorrected column barely improves as the cut rises:\n"
              "  you remove the worst-refracted stars but also lose geometry.\n"
              "  Cutting is the wrong trade. Correct and weight instead.\n");
}

void testPressureSensitivity() {
  std::printf("\n=== TEST 3: does the 800 m pressure correction matter? ===\n");
  auto epoch = Epoch::fromUtc(2024, 9, 15, 10, 30, 0.0);
  if (!epoch) return;

  const Atmosphere truth = Atmosphere::isa(kSite.alt, 10.0);
  const auto catalog = catalogBrighterThan(3.0);
  const auto obs =
      observeVisible(catalog, kSite, *epoch, 15.0 * kDeg2Rad, truth);
  const auto base = makeSights(obs, catalog, *epoch);

  struct Case { const char* label; Atmosphere a; };
  const std::vector<Case> cases = {
      {"correct (ISA @ 800 m)", truth},
      {"sea level assumed", Atmosphere{1013.25, 10.0, 0.5, 0.55}},
      {"temp off by +15 K", Atmosphere{truth.pressure_hpa, 25.0, 0.5, 0.55}},
      {"temp off by -15 K", Atmosphere{truth.pressure_hpa, -5.0, 0.5, 0.55}},
  };

  for (const Case& c : cases) {
    auto s = base;
    correctRefraction(s, c.a);
    std::printf("  %-24s %8.1f hPa %6.1f C  ->  %8.1f m\n", c.label,
                c.a.pressure_hpa, c.a.temperature_c, fixError(solveFix(s), kSite));
  }
  std::printf("\n  This is a second, independent reason to have a real static\n"
              "  source and an OAT probe rather than the autopilot's internal\n"
              "  barometer temperature.\n");
}

void testSkyDependence() {
  std::printf("\n=== TEST 4: refraction bias depends on WHICH stars are up ===\n");
  std::printf("  Sweeping through a night. The aircraft never moves.\n");
  const Atmosphere atmos = Atmosphere::isa(kSite.alt, 10.0);
  const auto catalog = catalogBrighterThan(3.0);

  std::printf("\n  %-8s %-5s %-14s %-12s %-14s\n", "UTC", "n",
              "uncorr err (m)", "bearing", "Bennett (m)");

  double min_err = 1e30, max_err = 0.0;
  for (int hour = 8; hour <= 20; hour += 2) {
    auto epoch = Epoch::fromUtc(2024, 9, 15, hour, 0, 0.0);
    if (!epoch) continue;
    const auto obs =
        observeVisible(catalog, kSite, *epoch, 15.0 * kDeg2Rad, atmos);
    auto sights = makeSights(obs, catalog, *epoch);
    if (sights.size() < 4) continue;

    const Fix raw = solveFix(sights);
    auto corrected = sights;
    correctRefraction(corrected, atmos);
    const Fix fixed = solveFix(corrected);

    const double e = fixError(raw, kSite);
    min_err = std::min(min_err, e);
    max_err = std::max(max_err, e);
    const double brg =
        raw.ok ? bearingDeg(kSite, Geodetic{raw.lat, raw.lon, kSite.alt}) : 0.0;

    std::printf("  %02d:00    %-5zu %-14.1f %-12.0f %-14.1f\n", hour,
                sights.size(), e, brg, fixError(fixed, kSite));
  }

  std::printf("\n  Uncorrected bias ranges %.0f m to %.0f m and the BEARING\n"
              "  swings right around the compass, purely because a different\n"
              "  part of the sky is up.\n", min_err, max_err);
  std::printf("\n  This matters for reading the paper. Table 2 reports 1.73 to\n"
              "  9.89 km across ten orbits, attributed to pitch/roll variance\n"
              "  and sample count (their Figure 8). But refraction bias is\n"
              "  earth-fixed, so orbit averaging does NOT remove it, and it\n"
              "  changes with sky region. Some of that orbit-to-orbit spread\n"
              "  may be refraction rather than attitude variance.\n");
  check(max_err / std::max(min_err, 1.0) > 1.5,
        "refraction bias varies materially with sky region");
}

void testDoesNotAverageOut() {
  std::printf("\n=== TEST 5: proof that the orbit does NOT remove refraction ===\n");
  auto epoch = Epoch::fromUtc(2024, 9, 15, 10, 30, 0.0);
  if (!epoch) return;
  const Atmosphere atmos = Atmosphere::isa(kSite.alt, 10.0);
  const auto catalog = catalogBrighterThan(3.0);
  const auto obs =
      observeVisible(catalog, kSite, *epoch, 15.0 * kDeg2Rad, atmos);
  const auto base = makeSights(obs, catalog, *epoch);

  // Simulate one orbit: sweep heading through 360 deg with a body-fixed
  // boresight misalignment of 0.4 deg, and average the ECEF zenith vectors
  // (Eq. 23). Refraction stays in the measurements throughout.
  const double eps = 0.4 * kDeg2Rad;
  const int n_frames = 360;

  auto runOrbit = [&](bool correct) {
    Eigen::Vector3d acc = Eigen::Vector3d::Zero();
    int used = 0;
    for (int k = 0; k < n_frames; ++k) {
      const double psi = 2.0 * M_PI * k / n_frames;
      // Body-fixed tilt of magnitude eps, rotated into NED by heading.
      const double tau_n = eps * std::cos(psi);
      const double tau_e = eps * std::sin(psi);

      auto s = base;
      for (size_t i = 0; i < s.size(); ++i) {
        const double az = obs[i].azimuth;
        s[i].zenith_angle -= (tau_n * std::cos(az) + tau_e * std::sin(az));
      }
      if (correct) correctRefraction(s, atmos);
      const Fix f = solveFix(s);
      if (!f.ok) continue;
      acc += f.zenith_ecef;
      ++used;
    }
    if (used == 0) return std::numeric_limits<double>::infinity();
    acc /= used;
    acc.normalize();
    const Geodetic est{std::asin(std::clamp(acc.z(), -1.0, 1.0)),
                       std::atan2(acc.y(), acc.x()), kSite.alt};
    return haversine(est, kSite);
  };

  // Single frame, heading zero, for reference.
  auto single = base;
  for (size_t i = 0; i < single.size(); ++i) {
    single[i].zenith_angle -= eps * std::cos(obs[i].azimuth);
  }
  const double single_err = fixError(solveFix(single), kSite);

  std::printf("  boresight misalignment: %.2f deg  (%.1f km if uncancelled)\n",
              eps * kRad2Deg, eps * 6371.0);
  std::printf("\n  single frame, refraction uncorrected : %9.1f m\n", single_err);
  std::printf("  orbit-averaged, refraction uncorrected: %9.1f m\n",
              runOrbit(false));
  std::printf("  orbit-averaged, refraction corrected  : %9.1f m\n",
              runOrbit(true));

  std::printf("\n  The orbit removes the 44 km misalignment completely. It\n"
              "  leaves the refraction bias essentially untouched, because\n"
              "  refraction is fixed in NED and does not rotate with heading.\n");

  check(runOrbit(false) < single_err,
        "orbit averaging removes the body-fixed misalignment");
  check(runOrbit(true) < runOrbit(false) / 3.0,
        "correction still needed after orbit averaging");
}

/// How much does robustness buy, and where does it stop working?
///
/// A graduated non-convexity solver used to sit on the fix path and was removed
/// after losing to RANSAC by an order of magnitude; the numbers are in
/// position_solver.hpp. What remains worth pinning is the RANSAC-vs-plain gap
/// and the star count below which robustness is impossible -- with a 3-star
/// minimal set you need 6+ sights before an outlier can be outvoted.
///
/// Note the simulator produces NO misidentifications, so the outlier-free row
/// is the only one the transit scenario exercises. That is why a bad solver
/// choice survived here unnoticed.
void testRobustSolverChoice() {
  std::printf("\n=== TEST: robust solver under misidentified stars ===\n");

  auto epoch = Epoch::fromUtc(2024, 11, 15, 14, 0, 0.0);
  const auto cat = catalogBrighterThan(5.5);
  Atmosphere atmos;
  atmos.pressure_hpa = 1013.25;
  atmos.temperature_c = 10.0;

  // A 53 deg field pointed up sees only near-zenith stars. Sampling the whole
  // sky flatters GNC, because the residual metric divides by sin(zenith).
  std::vector<Observation> obs;
  for (const auto& o : observeVisible(cat, kSite, *epoch, 20.0 * kDeg2Rad))
    if (o.elevation() > 64.0 * kDeg2Rad) obs.push_back(o);
  check(obs.size() >= 20, "enough near-zenith stars to sample");

  const double sigma = 20.0 * kArcsec2Rad;
  auto meanErr = [&](double outlier_frac, int which) {
    double acc = 0;
    int n = 0;
    for (unsigned seed = 1; seed <= 60; ++seed) {
      std::mt19937 rng(seed);
      std::normal_distribution<double> noise(0.0, sigma);
      std::uniform_real_distribution<double> u(0.0, 1.0);
      std::uniform_real_distribution<double> gross(0.3 * kDeg2Rad,
                                                   3.0 * kDeg2Rad);
      std::vector<Observation> pick = obs;
      std::shuffle(pick.begin(), pick.end(), rng);
      pick.resize(20);
      auto sights = makeSights(pick, cat, *epoch);
      for (auto& sg : sights) {
        sg.zenith_angle += noise(rng);
        if (u(rng) < outlier_frac) sg.zenith_angle += gross(rng);
      }
      correctRefraction(sights, atmos);
      Fix f;
      if (which == 0) f = solveFix(sights);
      else { RansacConfig rc; rc.seed = seed; f = solveFixRansac(sights, rc); }
      if (!f.ok) continue;
      acc += fixError(f, kSite) / 1000.0;
      ++n;
    }
    return n ? acc / n : -1.0;
  };

  std::printf("  %10s %10s %10s\n", "outliers", "plain LS", "RANSAC");
  double ls_o = 0, rs_o = 0;
  for (double of : {0.0, 0.10, 0.25}) {
    const double a = meanErr(of, 0), c = meanErr(of, 2);
    std::printf("  %9.0f%% %9.2fk %9.2fk\n", of * 100, a, c);
    if (of > 0.05) { ls_o = a; rs_o = c; }
  }

  check(rs_o > 0 && rs_o < 5.0, "RANSAC survives 25% misidentification");
  check(rs_o < ls_o / 20.0, "RANSAC beats plain least squares by >20x");
}

}  // namespace

int main() {
  std::printf("Milestone 0 — geometry round trip and refraction\n");
  std::printf("==========================================\n");

  testCleanRoundTrip();
  testTiltSensitivity();
  testGeodeticTrap();
  testRansacWithOutliers();
  testRobustSolverChoice();
  testNoiseAndRefraction();
  testModelAgainstErfa();
  testAblation();
  testPressureSensitivity();
  testSkyDependence();
  testDoesNotAverageOut();

  std::printf("\n==========================================\n");
  if (g_failures == 0) {
    std::printf("ALL CHECKS PASSED\n");
    return 0;
  }
  std::printf("%d CHECK(S) FAILED\n", g_failures);
  return 1;
}
