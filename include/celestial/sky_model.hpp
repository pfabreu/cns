// sky_model.hpp — astrometry, wrapping ERFA (the BSD-licensed fork of IAU SOFA).
//
// Two directions are provided and they are deliberately INDEPENDENT code paths:
//
//   observe()          FORWARD  : known site + star  -> observed az / zenith angle
//                                 (via eraAtco13: proper motion, parallax,
//                                  annual + diurnal aberration, light
//                                  deflection, precession-nutation, Earth
//                                  rotation, optional refraction)
//
//   subStellarPoint()  INVERSE  : star -> sub-stellar point on Earth
//                                 (via eraAtci13 + eraEra00)
//
// The round-trip test drives observe() to synthesise sights, then feeds them to
// subStellarPoint() + solveFix(). Because the two paths share only ERFA's
// low-level models and not the geometry, a clean round trip is a genuine test
// of the frames, the time handling and the solver — not a tautology.

#pragma once

#include <optional>
#include <vector>

#include "celestial/types.hpp"

namespace celestial {

/// A UTC epoch held as ERFA's two-part quasi-Julian date, plus the derived
/// TT and UT1 scales.
class Epoch {
 public:
  /// Construct from a UTC calendar date/time. Returns nullopt if ERFA rejects
  /// the date (e.g. an invalid month or a bad leap-second boundary).
  ///
  /// dut1 = UT1 - UTC in seconds. |dut1| < 0.9 s always. Neglecting it costs
  /// up to 0.9 s of Earth rotation = 13.5 arcsec = ~415 m of longitude, so for
  /// flight use fetch it from IERS Bulletin A. For this test 0 is fine.
  static std::optional<Epoch> fromUtc(int year, int month, int day, int hour,
                                      int minute, double second,
                                      double dut1 = 0.0);

  double utc1() const { return utc1_; }
  double utc2() const { return utc2_; }
  double tt1() const { return tt1_; }
  double tt2() const { return tt2_; }
  double ut11() const { return ut11_; }
  double ut12() const { return ut12_; }
  double dut1() const { return dut1_; }

  /// A new epoch `seconds` later. Use this rather than rebuilding from a
  /// calendar time with an out-of-range seconds field: eraDtf2d rejects
  /// seconds >= 60, which silently truncates any run longer than a minute.
  Epoch advanced(double seconds) const;

  /// Earth Rotation Angle at this epoch, radians. IAU 2000.
  double earthRotationAngle() const;

  /// Julian years since J2000 (TT), used for proper-motion bookkeeping.
  double julianYearsSinceJ2000() const;

 private:
  double utc1_ = 0, utc2_ = 0;
  double tt1_ = 0, tt2_ = 0;
  double ut11_ = 0, ut12_ = 0;
  double dut1_ = 0;
};

/// Atmospheric state for refraction. Set pressure_hpa = 0 to disable
/// refraction entirely (the correct choice when validating pure geometry).
struct Atmosphere {
  double pressure_hpa = 0.0;
  double temperature_c = 15.0;
  double relative_humidity = 0.5;  ///< 0..1
  double wavelength_um = 0.55;

  static Atmosphere none() { return Atmosphere{0.0, 15.0, 0.5, 0.55}; }
  /// Rough ISA pressure at a given geometric height, for refraction modelling.
  static Atmosphere isa(double altitude_m, double temperature_c = 15.0);
};

/// Result of the forward observation model.
struct Observation {
  int id = 0;
  double azimuth = 0.0;      ///< radians, from North through East
  double zenith_angle = 0.0; ///< radians
  double elevation() const { return M_PI_2 - zenith_angle; }
};

/// FORWARD model. Where does `star` appear in the sky, as seen from `site` at
/// `epoch`? Uses eraAtco13.
Observation observe(const CatalogStar& star, const Geodetic& site,
                    const Epoch& epoch,
                    const Atmosphere& atmos = Atmosphere::none());

/// All stars in `catalog` that are above `min_elevation` as seen from `site`.
std::vector<Observation> observeVisible(const std::vector<CatalogStar>& catalog,
                                        const Geodetic& site,
                                        const Epoch& epoch,
                                        double min_elevation,
                                        const Atmosphere& atmos = Atmosphere::none());

/// INVERSE model. The unit vector, in ECEF, of the point on Earth where `star`
/// is at the zenith at `epoch`.
///
/// Uses eraAtci13 to get CIRS apparent place, then
///     longitude_east = ra_cirs - ERA
///     latitude       = dec_cirs
/// which is the CIO-based form of the paper's lambda = alpha - GHA_Aries.
Eigen::Vector3d subStellarPoint(const CatalogStar& star, const Epoch& epoch);

// ---------------------------------------------------------------------------
// StarField — the fast path
// ---------------------------------------------------------------------------
//
// Calling eraAtci13/eraAtco13 per star per frame is far too slow for a 10 Hz
// pipeline: each call redoes the whole precession-nutation-aberration chain.
//
// Over a single orbit (~150 s) that chain is constant to well under 0.01
// arcsec. So reduce ICRS -> CIRS ONCE, cache the unit vectors, and per frame
// apply only Earth rotation:
//
//     v_ecef = Rz(-ERA) * v_cirs
//
// Two useful consequences:
//   * the sub-stellar point of a star IS its ECEF direction, so the GP needed
//     by the position solver falls straight out with no extra work;
//   * the inner loop becomes pure linear algebra, which is what you want when
//     this eventually runs on a Raspberry Pi.
//
// This is the representation the real-time implementation should use.
class StarField {
 public:
  /// Reduce a catalogue to CIRS at a reference epoch (use the orbit mid-time).
  StarField(const std::vector<CatalogStar>& catalog, const Epoch& reference);

  size_t size() const { return ids_.size(); }
  int id(size_t i) const { return ids_[i]; }
  double vmag(size_t i) const { return vmag_[i]; }
  const Eigen::Vector3d& cirs(size_t i) const { return cirs_[i]; }

  /// ECEF unit directions of every star at `epoch`. Also the sub-stellar
  /// point unit vectors, which is what the position solver wants.
  void ecefDirections(const Epoch& epoch,
                      std::vector<Eigen::Vector3d>& out) const;

  /// ECEF direction of one star at `epoch`.
  Eigen::Vector3d ecefDirection(size_t i, const Epoch& epoch) const;

  /// Index of a star by catalogue id, or -1.
  int indexOf(int id) const;

 private:
  std::vector<int> ids_;
  std::vector<double> vmag_;
  std::vector<Eigen::Vector3d> cirs_;
};

/// ECEF basis of the local NED frame at a geodetic position.
/// Columns are (north, east, down) expressed in ECEF, i.e. C_ecef_ned.
Eigen::Matrix3d nedBasisEcef(const Geodetic& p);

// ---------------------------------------------------------------------------
// Refraction
// ---------------------------------------------------------------------------
//
// Lives here rather than in its own module because it is only ever applied
// alongside the sky model: refraction is a property of the OBSERVER's
// atmosphere, and the only thing that needs it is the reduction of an observed
// direction to a true one.

namespace refraction {

/// Bennett (1982). Takes APPARENT (refracted) elevation — i.e. what you
/// actually measure — and returns the refraction in radians. This is the one
/// you want on the measurement side: no iteration required.
///
/// Accurate to a few arcseconds above ~15 deg elevation.
double bennett(double apparent_el, double pressure_hpa, double temperature_c);

/// Saemundsson (1986). Takes TRUE (unrefracted) elevation and returns the
/// refraction in radians. This is the inverse of Bennett and is what you want
/// in the FORWARD direction, e.g. inside the image simulator when you know the
/// true geometry and need to place the star where it will appear.
double saemundsson(double true_el, double pressure_hpa, double temperature_c);

}  // namespace refraction

/// OPTION A: correct the measurements in place.
///
/// Each sight's zenith angle is INCREASED by the refraction, undoing the
/// upward displacement. This keeps subStellarPoint() observer-independent,
/// which is the right separation of concerns: refraction is a property of the
/// observer's atmosphere, not of the star.
///
/// Uses the sight's own zenith angle to derive apparent elevation, so it is
/// self-contained and needs no extra bookkeeping.
///
/// Pass an Atmosphere built with Atmosphere::isa(altitude) — at 800 m the
/// pressure is ~921 hPa and refraction is ~9% weaker than at sea level.
/// Ignoring that over-corrects by ~130 m.
void correctRefraction(std::vector<StarSight>& sights, const Atmosphere& atmos);

/// Set per-sight weights from the expected measurement variance.
///
/// After correction, the residual error in a sight is dominated by the
/// refraction MODEL error, which scales with the size of the correction:
///
///     sigma_i^2 = sigma_centroid^2 + (model_error_frac * R(el_i))^2
///     w_i       = 1 / sigma_i^2
///
/// This is strictly better than a hard elevation cut. A cut throws away stars
/// and weakens the geometry (see the 15 deg vs 40 deg comparison in the
/// refraction ablation); weighting keeps low stars contributing what they are
/// worth without letting them dominate.
///
/// model_error_frac ~ 0.05 is a reasonable default: 5% of the correction,
/// covering Bennett model error plus pressure and temperature uncertainty.
void weightForRefraction(std::vector<StarSight>& sights,
                         const Atmosphere& atmos,
                         double centroid_sigma_rad = 10.0 * kArcsec2Rad,
                         double model_error_frac = 0.05);


/// Convenience: turn forward observations into solver input.
std::vector<StarSight> makeSights(const std::vector<Observation>& obs,
                                  const std::vector<CatalogStar>& catalog,
                                  const Epoch& epoch);


// ===========================================================================
// POSITION SOLVER -- folded in from position_solver.hpp. It operates on the
// StarSight values this header defines, so the two were always one module.
// ===========================================================================


/// Least-squares fix from three or more sights. Eq. (7)-(10) of the paper.
///
/// Weights, if set on the sights, are applied as a diagonal weighting of the
/// normal equations: x = (A^T W A)^-1 A^T W p.
Fix solveFix(const std::vector<StarSight>& sights);

struct RansacConfig {
  int iterations = 100;
  /// Inlier threshold on the plane residual |p_k - a_k . x|, expressed as an
  /// angle. 0.05 deg at the surface is ~5.5 km, a sane starting point.
  double tolerance = 0.05 * kDeg2Rad;
  double min_inlier_ratio = 0.5;
  unsigned seed = 1;
};

/// RANSAC position estimate, Algorithm 1 of the paper. Rejects misidentified
/// stars, satellites, aircraft and hot pixels. Requires >= 6 sights to be
/// meaningful; falls through to solveFix() below that.
///
/// This is the fix path. A Graduated Non-Convexity / truncated-least-squares
/// solver was tried here and REMOVED: it is deterministic, uses every
/// measurement each iteration and needs no initial guess, all of which is true
/// and none of which made it competitive. Measured on a near-zenith field
/// (a 53 deg camera pointed up sees nothing below ~64 deg elevation), 20
/// stars, 60 trials, mean fix error:
///
///   misidentified   plain LS      GNC     RANSAC
///        0%          0.34 km    0.34 km   0.34 km
///       10%         25.48 km    6.22 km   0.34 km
///       25%         45.06 km   31.69 km   0.42 km
///
/// A misidentified star here is GROSS -- 0.3 to 3 deg -- which is the regime a
/// minimal-seed hypothesis test handles best, and the corrupted least-squares
/// seed a graduated method starts from is exactly what it handles worst.
/// See README.md before reaching for a robust solver again.
Fix solveFixRansac(const std::vector<StarSight>& sights,
                   const RansacConfig& cfg = {});

}  // namespace celestial
