// horizon.hpp — optical horizon sensor as an independent vertical reference.
//
// -------------------------------------------------------------------------
// WHY THIS EXISTS
// -------------------------------------------------------------------------
// The stars are exact; the vertical is the error. Everything else in this
// library is about supplying or cleaning up the vertical:
//
//   IMU          supplies it. Required, always. Every zenith angle needs
//                C_l_b on every frame.
//   orbit/circle REMOVES the body-fixed part of whatever error the IMU has,
//                by exploiting the fact that it rotates with heading.
//   horizon      supplies a SECOND, INDEPENDENT measurement of it.
//
// The horizon is non-inertial, so unlike the accelerometers it is not confused
// by acceleration. Differencing it against the IMU therefore observes the IMU
// tilt bias IN STRAIGHT AND LEVEL FLIGHT, with no maneuver. That is what makes
// flight geometry an improvement rather than a requirement.
//
// -------------------------------------------------------------------------
// ERROR STRUCTURE (this is the interesting part)
// -------------------------------------------------------------------------
// A forward-looking camera measures two things from one horizon line:
//
//   line OFFSET -> dip   -> tilt about the body LATERAL axis
//   line SLOPE           -> tilt about the body FORWARD axis
//
// They have completely different error behaviour, because anomalous refraction
// near the tangent point is largely COMMON-MODE across the field:
//
//   * a common-mode shift moves the line up or down  -> corrupts the dip
//   * a common-mode shift does NOT rotate the line   -> slope is immune
//
// So the slope axis is differential and reaches sub-arcminute, while the dip
// axis eats the full 2-5 arcmin of anomalous refraction. The measurement
// covariance is strongly anisotropic and must be reported as such.
//
// A FORE/AFT PAIR makes the dip differential too:
//     d_fore = d + tau_lat + delta_fore
//     d_aft  = d - tau_lat + delta_aft
// so (d_fore - d_aft)/2 recovers the tilt with the common-mode delta removed.
// Same cancellation as the orbit maneuver, done geometrically.
//
// -------------------------------------------------------------------------
// SCOPE
// -------------------------------------------------------------------------
// Over water only. Over terrain the "horizon" is a ridgeline: 100 m of relief
// at 100 km range is 3.4 arcmin of false dip that varies with azimuth. This
// model does not attempt terrain.
//
// The optics are modelled; the computer vision is not. We assume the horizon
// line is extracted at the stated per-pixel accuracy, which for an LWIR camera
// over water at night (a 20-40 K sea/sky step) is a reasonable idealisation.

#pragma once

#include <deque>
#include <vector>

#include "celestial/attitude.hpp"
#include "celestial/imaging.hpp"
#include "celestial/types.hpp"

namespace celestial {

/// One camera. Boresight is body +x (fore) or body -x (aft).
struct HorizonCamera {
  bool aft = false;
  double hfov = 32.0 * kDeg2Rad;  ///< FLIR Boson 640 with a 32 deg lens
  int width = 640;
  double edge_sigma_px = 0.5;  ///< subpixel line-fit accuracy per column

  // --- optical imperfections (used by the renderer; the CV does not know
  // about them, which is the point) ---------------------------------------

  /// Radial distortion of the TRUE optics. The extractor assumes a pinhole,
  /// so a distorted horizon is a CURVE fitted by a straight line, giving a
  /// systematic error that depends on where the line sits in the frame --
  /// i.e. on pitch. LWIR lenses run 1-3% at the field edge; calibration gets
  /// you to ~0.1% residual.
  ///
  /// Note this error is COMMON-MODE between fore and aft (the line sits at
  /// the same image row in both), so the pair cancels it. Another reason the
  /// pair earns its keep beyond anomalous refraction.
  double k1 = 0.0;

  /// Fractional focal length error. Scales the dip measurement directly:
  /// 1% on a 47 arcmin dip is 0.5 arcmin, which is not negligible.
  double focal_err = 0.0;

  /// Column fixed-pattern noise, in counts. The dominant FPN mode in
  /// microbolometers. It biases the per-column edge estimate systematically,
  /// and unlike a mounting offset it DRIFTS with sensor temperature -- so it
  /// is the horizon sensor's own bias term, the honest counterweight to
  /// "the horizon's error is calibratable and the IMU's is not".
  /// Modelled simply: a fixed per-column offset, scaled by fpn_drift.
  double fpn_counts = 0.0;

  /// Static mounting misalignment, radians. Unlike an accelerometer bias this
  /// is a mechanical offset on a rigid bracket: it does not drift with
  /// temperature or vibration, so it is calibratable once on the ground.
  /// It is also body-fixed, so an orbit removes whatever is left of it.
  double mount_lat = 0.0;
  double mount_fwd = 0.0;

  double ifov() const { return hfov / width; }

  /// FLIR Lepton 2.5 -- 80x60, ~51 deg HFOV, uncooled microbolometer, a
  /// couple of hundred dollars against thousands for a research-grade core.
/// The question is
  /// not whether it resolves the horizon but whether the LINE FIT across 80
  /// columns is good enough, and the answer is set by the atmosphere rather
  /// than the sensor: anomalous dip is ~3 arcmin and irreducible, so a sensor
  /// better than that is wasted.
  ///
  /// edge_sigma_px is 1.2 rather than a research core's ~0.5: fewer, larger
  /// pixels, a cheaper lens and more FPN per column. That number is an
  /// ESTIMATE, not a measurement, and it is the one to check first against
  /// real hardware.
  static HorizonCamera lepton25(bool aft = false) {
    HorizonCamera c;
    c.aft = aft;
    c.width = 80;
    c.hfov = 51.0 * kDeg2Rad;
    c.edge_sigma_px = 1.2;
    c.k1 = 0.02;          // cheap silicon optics, ~2% at the field edge
    c.focal_err = 0.005;  // 0.5% focal error after a rough calibration
    c.fpn_counts = 2.0;   // column FPN, worse than a research-grade core
    return c;
  }
};

struct HorizonConfig {
  std::vector<HorizonCamera> cameras{HorizonCamera{}};  ///< 1 = fore, 2 = fore/aft

  /// Anomalous dip: refraction anomaly at the tangent point. Modelled as a
  /// Gauss-Markov process SHARED between cameras -- that sharing is what makes
  /// the fore/aft cancellation work, so do not randomise it per camera.
  ///
  /// CORRELATION TIME MATTERS: at ~600 s this is effectively constant over an
  /// orbit, so it behaves as a bias, not noise. Inflating R instead of
  /// modelling it will make the covariance dishonest.
  double anom_dip_sigma = 3.0 * kArcmin2Rad;
  double anom_dip_tau = 600.0;
  /// Across-field gradient of the anomalous dip -> apparent line rotation.
  /// This, not edge noise, usually dominates the slope axis.
  double anom_grad_sigma = 0.5 * kArcmin2Rad;

  /// Fraction of the anomalous dip that is NOT common-mode between cameras.
  ///
  /// The fore/aft cancellation is only as good as the correlation between two
  /// tangent points ~200 km apart. Setting this to 0 would make the pair
  /// perfect, which is the easiest way to fool yourself with this model.
  double anom_decorrelated_frac = 0.33;

  /// Slowly varying scale on the column FPN, standing in for NUC drift with
  /// sensor temperature. 1.0 = as calibrated.
  double fpn_drift = 1.0;

  /// Sea-state blur of the limb, pixels. The sea/sky transition is not a step:
  /// whitecaps, spray and atmospheric blur spread it over several pixels, and
  /// this is where the real edge uncertainty comes from.
  double limb_blur_px = 1.5;

  /// Moving-average window for the tilt measurement, seconds. 0 disables.
  ///
  /// The original code argued no filter was needed, because "the orbit
  /// machinery downstream already performs the time-averaging". That holds for
  /// a research-grade core at 0.5 px of edge noise. A Lepton is 1.2 px and is
  /// NOISE-limited -- halving its edge noise is worth 2.2x on the fix -- so the
  /// argument no longer applies and a short average pays.
  ///
  /// Measured, 1x Lepton, 30 seeds, one revolution at 250 m:
  ///
  ///     averaging   tilt error   orbit fix
  ///       none        0.200 deg    2.31 km
  ///       1.0 s       0.104 deg    1.48 km
  ///       2.5 s       0.084 deg    1.38 km
  ///       4.0 s       0.075 deg    1.93 km   <- tilt still improving, fix worse
  ///
  /// Note the LAST row: past ~3 s the tilt error keeps falling while the fix
  /// gets worse. The average lags, and during an orbit a lag is a
  /// heading-correlated error -- which is exactly the kind the orbit average
  /// cannot remove. Optimising the tilt number past this point makes the
  /// navigation worse.
  double tilt_average_s = 2.5;

  unsigned seed = 4242;
};

/// Evolving state of the atmosphere model. Keep one per flight.
struct HorizonState {
  /// Recent tilt measurements, for the moving average. See
  /// HorizonConfig::tilt_average_s.
  std::deque<Eigen::Vector2d> tilt_hist;
  double last_meas_t = -1.0;
  double anom_dip = 0.0;
  double anom_grad = 0.0;
  double last_t = -1.0;
};

/// A measurement of the IMU's TILT ERROR, in the body frame.
///
/// tilt(0) = error about the body FORWARD axis (roll-like), from line slope
/// tilt(1) = error about the body LATERAL axis (pitch-like), from dip
struct TiltMeasurement {
  bool ok = false;
  Eigen::Vector2d tilt = Eigen::Vector2d::Zero();
  Eigen::Matrix2d R = Eigen::Matrix2d::Identity();
};

/// Geometric dip of the horizon, reduced by terrestrial refraction.
/// ~46.9 arcmin at 800 m. A 10 m altitude error costs only 0.29 arcmin, so
/// barometric altitude is entirely adequate here.
double horizonDip(double altitude_m, double refraction_k = 0.13);

/// Simulate the sensor. Compares the TRUE attitude against the IMU's estimate
/// and returns a noisy measurement of the difference.
///
/// In hardware this function is replaced by horizon extraction from a real
/// camera; nothing downstream changes.
TiltMeasurement measureHorizon(const Eigen::Matrix3d& C_true,
                               const Eigen::Matrix3d& C_est, double altitude_m,
                               double t, const HorizonConfig& cfg,
                               HorizonState& state);

/// A horizon line extracted from an image, in pixel coordinates.
struct HorizonLine {
  bool ok = false;
  double v_centre = 0.0;   ///< row of the line at the image centre column
  double slope = 0.0;      ///< d(row)/d(column)
  double sigma_v = 0.0;    ///< from the fit residuals, not assumed
  double sigma_slope = 0.0;
  int n_inliers = 0;
  double residual_mad = 0.0;  ///< robust scatter, px -- a cloud/limb quality gate
};

/// Extract the horizon by per-column gradient peak with parabolic subpixel
/// refinement, then a THEIL-SEN robust line fit.
///
/// Theil-Sen (median of pairwise slopes) rather than least squares because
/// cloud tops, ships and whitecaps produce outlier columns, and a single bad
/// column drags an ordinary regression. Its 29% breakdown point is ample here.
///
/// The covariance comes from the fit RESIDUALS, so the measurement degrades
/// gracefully and self-reports: a ragged limb widens sigma instead of quietly
/// biasing the answer.
HorizonLine detectHorizon(const Image& img, const HorizonCamera& cam);

/// Where the horizon line SHOULD be for a given attitude, assuming ideal
/// pinhole optics. Differencing this against detectHorizon() gives the tilt.
HorizonLine predictHorizonLine(const Eigen::Matrix3d& C, double altitude_m,
                               const HorizonCamera& cam);

/// The CV path: render -> detect -> tilt measurement.
///
/// Structurally identical to the star path (render -> detect -> match), and it
/// MEASURES the edge accuracy rather than assuming it.
TiltMeasurement measureHorizonCV(const Eigen::Matrix3d& C_true,
                                 const Eigen::Matrix3d& C_est,
                                 double altitude_m, double t,
                                 const HorizonConfig& cfg, HorizonState& state);



/// Static Bayesian fusion of the horizon measurement with the IMU prior.
///
/// NO FILTER IS NEEDED FOR THIS. The horizon measures the IMU tilt error at
/// the same instant, so a single update is already correct. A filter would add
/// two things later -- averaging the white per-frame noise, and separating the
/// anomalous-dip bias from genuine tilt -- but the orbit machinery downstream
/// already performs the time-averaging for the body-fixed part.
///
/// Returns a corrected attitude. Everything downstream is unchanged: pass the
/// result wherever you previously passed the raw AHRS attitude.
Eigen::Matrix3d fuseTilt(const Eigen::Matrix3d& C_est, const TiltMeasurement& z,
                         double imu_tilt_sigma,
                         Eigen::Matrix2d* posterior = nullptr);

/// Synthetic LWIR view of the horizon. VISUALISATION ONLY -- the measurement
/// path models the optics analytically and never touches pixels.
///
/// Over water at night the sea/sky boundary is a 20-40 K step, which is why
/// thermal makes this easy where visible light does not. Each pixel is
/// classified by whether its line of sight falls below the depressed horizon.
Image renderHorizonView(const Eigen::Matrix3d& C_true, double altitude_m,
                        const HorizonCamera& cam, const HorizonState& state,
                        unsigned seed = 0);

}  // namespace celestial
