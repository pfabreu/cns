// types.hpp — core value types for the celestial navigation library.
//
// Conventions used EVERYWHERE in this library:
//   * Angles are radians unless a name ends in _deg.
//   * Longitude is EAST-POSITIVE.
//   * Local level frame is NED (North, East, Down).
//   * Body frame is FRD (Forward, Right, Down).
//   * "ECEF" here means the Terrestrial Intermediate Reference System with
//     polar motion neglected (xp = yp = 0). The error from neglecting polar
//     motion is < 0.5 arcsec (~15 m), far below our target accuracy.
//
// If you are wiring this to MAVROS later: convert ENU/FLU to NED/FRD at the
// MAVROS boundary and never let ENU leak into this library.

#pragma once

#include <Eigen/Dense>
#include <string>

namespace celestial {

inline constexpr double kDeg2Rad = M_PI / 180.0;
inline constexpr double kRad2Deg = 180.0 / M_PI;
inline constexpr double kArcsec2Rad = M_PI / (180.0 * 3600.0);
inline constexpr double kArcmin2Rad = M_PI / (180.0 * 60.0);
inline constexpr double kMas2Rad = kArcsec2Rad / 1000.0;

// WGS-84
inline constexpr double kWgs84A = 6378137.0;
inline constexpr double kWgs84F = 1.0 / 298.257223563;

/// Geodetic position on the WGS-84 ellipsoid.
struct Geodetic {
  double lat = 0.0;  ///< geodetic latitude, radians
  double lon = 0.0;  ///< longitude, radians, east positive
  double alt = 0.0;  ///< height above ellipsoid, metres

  static Geodetic fromDegrees(double lat_deg, double lon_deg, double alt_m) {
    return Geodetic{lat_deg * kDeg2Rad, lon_deg * kDeg2Rad, alt_m};
  }
};

/// A star as stored in the catalogue: ICRS / J2000, with proper motion.
struct CatalogStar {
  int id = 0;             ///< catalogue identifier (HR number for the built-in BSC5)
  std::string name;        ///< optional common name; empty for BSC5 rows
  double ra = 0.0;         ///< ICRS right ascension at J2000, radians
  double dec = 0.0;        ///< ICRS declination at J2000, radians
  double pm_ra_star = 0.0; ///< mu_alpha* = (dRA/dt)*cos(dec), rad/year
  double pm_dec = 0.0;     ///< dDec/dt, rad/year
  double vmag = 0.0;       ///< visual magnitude
};

/// One reduced star observation, ready for the position solver.
///
/// `gp` is the sub-stellar point ("geographical position" in the classical
/// nautical sense): the unit vector, in ECEF, of the place on Earth where this
/// star is exactly overhead at the observation epoch.
///
/// `zenith_angle` is the measured angle between the star and the observer's
/// vertical. In the real system this comes from the star tracker centroid
/// rotated through C_l/b * C_b/c; here it comes from the forward model.
struct StarSight {
  int id = 0;
  Eigen::Vector3d gp = Eigen::Vector3d::UnitZ();
  double zenith_angle = 0.0;
  double weight = 1.0;
};

/// Output of the position solve.
struct Fix {
  bool ok = false;
  double lat = 0.0;  ///< radians
  double lon = 0.0;  ///< radians, east positive
  /// The unit vector solved for. See position_solver.hpp for why this is the
  /// observer's ZENITH direction and not a position vector.
  Eigen::Vector3d zenith_ecef = Eigen::Vector3d::UnitZ();
  double residual_rms = 0.0;  ///< rms of (p - A x), dimensionless (cos units)
  int n_used = 0;
  /// Condition number of A^T A. Large values mean the stars were poorly
  /// distributed and the fix is geometrically weak (the celestial equivalent
  /// of GDOP).
  double condition = 0.0;
};

/// Great-circle distance between two geodetic positions, metres.
/// Haversine on a sphere of the WGS-84 mean radius — matches the error metric
/// used in the Teague & Chahl paper.
double haversine(const Geodetic& a, const Geodetic& b);

}  // namespace celestial
