// orbit.hpp — MILESTONE 1: analytic orbit, no images.
//
// Two halves:
//
//   1. A trajectory + camera simulator that produces, per frame, the set of
//      star direction vectors AS SEEN IN THE CAMERA FRAME, together with the
//      AHRS attitude the autopilot would report. No pixels, no noise model
//      beyond angular jitter — the image simulator is Milestone 2.
//
//   2. The orbit-level estimator: per-frame fixes, averaging over a heading
//      revolution (Eq. 23), and the Kabsch recalibration loop (Algorithm 2).
//
// The point of the whole exercise: a body-fixed error traces a circle in the
// navigation frame as heading sweeps 360 deg, and averaging removes it.

#pragma once

#include <vector>

#include "celestial/attitude.hpp"
#include "celestial/sky_model.hpp"
#include "celestial/types.hpp"

namespace celestial {

// ---------------------------------------------------------------------------
// Camera
// ---------------------------------------------------------------------------

/// Pinhole camera. Defaults match the paper's payload: Alvium 1800 U-240,
/// 1936 x 1216, 6 mm f/1.4 lens giving 53.5 deg horizontal field of view.
struct Camera {
  int width = 1936;
  int height = 1216;
  double hfov = 53.5 * kDeg2Rad;

  double focalPx() const {
    return (width / 2.0) / std::tan(hfov / 2.0);
  }
  /// Angular size of one pixel at the centre of the field, radians.
  double pixelIfov() const { return 1.0 / focalPx(); }
  /// True if a camera-frame direction lands on the sensor.
  bool inFov(const Eigen::Vector3d& v_cam) const;
};

/// The nominal "camera looking straight up" mounting, in FRD body axes:
/// camera +z (boresight) -> body -z (up), camera +x -> body +y (starboard),
/// camera +y -> body +x (forward).
Eigen::Matrix3d nominalCameraMount();

// ---------------------------------------------------------------------------
// Trajectory
// ---------------------------------------------------------------------------

struct OrbitConfig {
  Geodetic centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  double radius_m = 600.0;
  double airspeed = 25.0;      ///< true airspeed, m/s
  double frame_rate = 10.0;    ///< Hz, as in the paper
  double start_track = 0.0;    ///< initial bearing from centre, radians
  double revolutions = 1.0;

  double wind_n = 0.0;  ///< m/s
  double wind_e = 0.0;

  /// true  : hold a circular GROUND track (what the paper's flight did, using
  ///         GPS). In wind this makes upwind arcs take longer, so heading is
  ///         sampled NON-UNIFORMLY in time and the naive mean is biased.
  /// false : hold a fixed attitude, i.e. constant bank and constant yaw rate.
  ///         Heading is sampled uniformly; the aircraft drifts downwind.
  bool gps_guided_track = true;

  /// Epoch of the first frame.
  int year = 2024, month = 9, day = 15, hour = 10, minute = 30;
  double second = 0.0;
};

/// True state at one frame.
struct FrameTruth {
  double t = 0.0;
  Geodetic pos;
  double roll = 0.0, pitch = 0.0, yaw = 0.0;  ///< TRUE attitude
  Epoch epoch;
};

std::vector<FrameTruth> generateOrbit(const OrbitConfig& cfg);

// ---------------------------------------------------------------------------
// Observation generation
// ---------------------------------------------------------------------------

/// Errors injected between truth and what the navigation code gets to see.
struct ErrorModel {
  /// Boresight misalignment: the TRUE camera mounting differs from the assumed
  /// one by this angle about this body axis. This is the paper's dominant term.
  double boresight_error = 0.4 * kDeg2Rad;
  Eigen::Vector3d boresight_axis = Eigen::Vector3d(1, 0, 0);

  /// Constant AHRS attitude bias, body frame (roll, pitch, yaw), radians.
  /// Aliased with the boresight error — that is the point.
  Eigen::Vector3d ahrs_bias = Eigen::Vector3d::Zero();

  /// Zero-mean AHRS noise, radians, 1-sigma per axis per frame.
  double ahrs_noise = 0.0;

  /// SLOWLY DRIFTING part of the AHRS tilt bias: a first-order Gauss-Markov
  /// process added to `ahrs_bias` on the roll and pitch axes, 1-sigma steady
  /// state in radians, with correlation time `ahrs_drift_tau` seconds.
  ///
  /// This exists because a purely CONSTANT tilt bias is absorbable. It is
  /// body-fixed and time-invariant, so `recalibrateMount` books it into the
  /// mounting at the calibration loiter and the product `C_l_b_est * C_b_c`
  /// comes out right ever after -- a wrong attitude times a compensating wrong
  /// mounting gives correct fixes. With a constant bias and NO horizon sensor
  /// the transit scores 0.66 km, better than the 2.68 km baseline WITH one,
  /// which is not a result about celestial navigation but an artefact of the
  /// error model. A one-shot absorption survives only if the absorbed quantity
  /// never changes.
  ///
  /// OFF by default, following every other term here: an ErrorModel should
  /// contain exactly the errors a caller asked for, and the analytic tests rely
  /// on that. The transit scenario opts in at 0.15 deg / 60 s, matching
  /// `TiltFilter::Config` so the simulator generates the process the filter
  /// claims to model. `--drift-sigma 0` restores the constant-bias behaviour.
  double ahrs_drift_sigma = 0.0;
  double ahrs_drift_tau = 60.0;

  /// MANEUVER-DEPENDENT tilt error: the fraction of the centripetal tilt that
  /// survives the AHRS's compensation, applied about the roll axis in
  /// proportion to bank.
  ///
  /// An accelerometer senses specific force, not gravity. In a coordinated turn
  /// the net specific force lies along the aircraft's own down axis, so an
  /// accelerometer-only vertical reads WINGS LEVEL however hard the aircraft is
  /// banked. An EKF removes most of this using velocity and turn rate, but
  /// GNSS-denied it has only airspeed and an imperfect wind estimate, so a
  /// residual proportional to bank survives. 0.02 means 2 percent of the bank
  /// leaks into the tilt estimate: 0.18 deg at 9 deg of bank, comparable to the
  /// constant bias.
  ///
  /// This is the term that makes the horizon camera worth carrying. It is
  /// correlated with maneuvering, so it does NOT average out over a loiter --
  /// and the loiter is where the mounting is calibrated, so it contaminates the
  /// calibration rather than being averaged away by it. A horizon camera is
  /// immune to it, being geometric rather than inertial.
  ///
  /// OFF by default, like every other term here.
  double ahrs_turn_coupling = 0.0;

  /// GYRO model, used only by `simulateGyro` for the star-aided attitude
  /// experiment. Nothing in the fix path reads these: the AHRS attitude is
  /// synthesised directly by `simulateAhrs`, which is why this project could
  /// not test attitude aiding until the gyro was added.
  ///
  /// Defaults are a cheap MEMS part: 0.5 deg/hr in-run bias stability with a
  /// 1 hr correlation time, and 0.005 deg/sqrt(s) angle random walk.
  double gyro_bias_sigma = 0.5 * kDeg2Rad / 3600.0;   ///< rad/s
  double gyro_bias_tau = 3600.0;                     ///< s
  double gyro_arw = 0.005 * kDeg2Rad;                ///< rad/sqrt(s)

  /// Centroiding noise applied to the camera-frame direction, radians.
  double centroid_noise = 0.0;

  /// Atmosphere used to generate the observations. Set pressure to 0 for
  /// no refraction.
  Atmosphere atmos = Atmosphere::none();

  /// Elevation below which stars are discarded, in addition to the FOV test.
  double min_elevation = 10.0 * kDeg2Rad;

  /// Detector sensitivity: stars fainter than this are not seen. The paper's
  /// f/1.4 6 mm optics at 10 Hz tracked ~27 stars in frame, which corresponds
  /// to roughly V = 4.0 for this field of view.
  double mag_limit = 4.0;

  unsigned seed = 1;
};

/// One frame's worth of data, as the navigation code would receive it.
struct FrameData {
  double t = 0.0;
  Epoch epoch;
  double yaw_est = 0.0;                ///< AHRS heading (for binning)
  Eigen::Matrix3d C_l_b_est;           ///< AHRS attitude estimate
  std::vector<int> id;                ///< identified star IDs
  std::vector<Eigen::Vector3d> v_cam;  ///< unit directions in the CAMERA frame
  Geodetic truth;                      ///< carried through for scoring only

  /// Cumulative dead-reckoned position at this frame, metres north/east from
  /// an arbitrary common origin.
  ///
  /// estimateOrbit averages the per-frame fixes as if they all came from ONE
  /// PLACE. That holds in a 400 m loiter and fails badly on a straight leg: at
  /// 25 m/s a 150 s window spans 3.75 km, so the "fix" lags the aircraft by
  /// half of that, and an unbounded window lags by the whole flight.
  ///
  /// Supplying this lets estimateOrbit TRANSPORT each fix forward to the
  /// window's reference epoch before averaging, removing the aircraft's own
  /// motion. Leave at zero to disable (the loiter-only behaviour).
  ///
  /// It is a DEAD-RECKONED quantity -- airspeed, heading and estimated wind.
  /// No truth is involved.
  Eigen::Vector2d dr_ne = Eigen::Vector2d::Zero();
};

/// Synthesise the raw AHRS attitude sequence from an error model: truth plus
/// the constant bias, plus the Gauss-Markov tilt drift, plus per-frame noise.
///
/// Shared by `simulateObservations` and by callers that need the raw estimate
/// before fusing a vertical reference onto it (the transit scenario does this).
/// Deterministic in `ErrorModel::seed`.
std::vector<Eigen::Matrix3d> simulateAhrs(const std::vector<FrameTruth>& truth,
                                          const ErrorModel& err);

/// Synthesise gyro samples: the TRUE body rate relative to inertial, plus bias
/// and noise. One sample per frame interval, valid over the interval starting
/// at that frame.
///
/// The true rate is differentiated from the truth attitude sequence and has
/// Earth rate added, since a real gyro senses rotation relative to inertial
/// space rather than relative to ECEF.
///
/// Exists so the star-aided attitude filter can be exercised. See
/// attitude.hpp for why aiding is a different problem from the tilt
/// estimation attempts that failed.
std::vector<Eigen::Vector3d> simulateGyro(const std::vector<FrameTruth>& truth,
                                          const ErrorModel& err);

std::vector<FrameData> simulateObservations(const std::vector<FrameTruth>& truth,
                                            const Camera& cam,
                                            const ErrorModel& err);

/// As above, but the AHRS attitude is SUPPLIED per frame rather than
/// synthesised from ErrorModel::ahrs_bias / ahrs_noise.
///
/// This is the entry point for Milestone 3: the estimated attitude comes from
/// a real EKF3 (SITL or flight log) instead of an assumed error model, so the
/// assumption that the attitude error is body-fixed is tested rather than
/// baked in.
std::vector<FrameData> simulateObservationsWithAttitude(
    const std::vector<FrameTruth>& truth,
    const std::vector<Eigen::Matrix3d>& C_l_b_est, const Camera& cam,
    const ErrorModel& err);

// ---------------------------------------------------------------------------
// Orbit-level estimation
// ---------------------------------------------------------------------------

/// How to combine the per-frame fixes into one position.
enum class AverageMethod {
  /// Eq. (23): plain mean of the ECEF zenith unit vectors. What the paper does.
  Naive,
  /// Mean weighted by heading increment. Equivalent to resampling uniformly in
  /// heading, so it is immune to the non-uniform sampling that ruins the naive
  /// mean under a GPS-guided track in wind.
  HeadingWeighted,
  /// Fit a small circle to the point cloud. The axis is the position and the
  /// angular radius is the misalignment magnitude, so this also yields a
  /// direct estimate of the boresight error.
  CircleFit,
};

struct OrbitResult {
  bool ok = false;
  Geodetic position;
  int n_frames = 0;
  /// Angular radius of the fitted circle, radians. Only set for CircleFit.
  /// Multiply by 6371 km to read it as the misalignment in position terms.
  double circle_radius = 0.0;
  /// Scatter of the per-frame fixes about the estimate, radians.
  double spread = 0.0;
};

/// Fix from a SINGLE frame. Returns false if fewer than three stars are
/// identified or the solve is degenerate. This is what the paper's Figure 6
/// plots frame by frame, and each such fix is one blue dot in Figure 7.
bool singleFrameFix(const FrameData& frame, const Eigen::Matrix3d& C_b_c,
                    const Atmosphere& refraction_model, Geodetic& out);

/// One pass: per-frame fixes with the given camera mount, then averaged.
OrbitResult estimateOrbit(const std::vector<FrameData>& frames,
                          const Eigen::Matrix3d& C_b_c_assumed,
                          AverageMethod method,
                          const Atmosphere& refraction_model = Atmosphere::none());

/// Recalibrate the camera mounting from an assumed position (Algorithm 2).
/// Per frame: build the theoretical NED star vectors at `pos`, Kabsch against
/// the observed camera vectors to get R = C_l_b * C_b_c, then C_b_c =
/// C_l_b^T * R. The per-frame results are averaged.
/// Heading estimated from the stars alone — a celestial compass.
///
/// WHY THIS EXISTS. Dead reckoning takes heading from `ATTITUDE.yaw`, which
/// GNSS-denied is magnetometer-derived. That is the last non-autonomous sensor
/// in the loop, and it is the weakest: a magnetometer can be spoofed, suffers
/// hard and soft iron error, and — worst — EKF3 looks up magnetic declination
/// at the LAST KNOWN POSITION, so heading error grows with position error,
/// which grows heading error. On a 1400 km leg declination changes by degrees.
///
/// The information was already being computed and discarded. `recalibrateMount`
/// runs Kabsch, which solves the FULL 3-DOF rotation between the observed and
/// catalogue directions; the mounting takes what it needs and the heading is
/// thrown away — worse, an EKF3 yaw error is absorbed as if it were mechanical
/// misalignment. Given a calibrated mounting, that same rotation yields the
/// true heading directly.
///
/// OBSERVABILITY. Camera mounting CLOCKING and vehicle YAW are rotations about
/// the body z and the local vertical respectively. In level flight those axes
/// coincide and the two are perfectly aliased.
///
/// A bank does NOT rescue this when the mounting is calibrated from the same
/// flight. The bias enters the per-frame mounting as a rotation about
/// `C_l_b^T * z_ned`, which depends on roll and pitch but not on yaw, so a
/// constant-bank loiter pins that axis however far heading sweeps and every
/// candidate bias fits equally well. Separating them needs bank DIVERSITY.
///
/// So this function requires the CLOCKING to be known, not solved for: pass a
/// mounting from `mountTiltOnly`, never one straight out of `recalibrateMount`.
/// Given that, heading is recovered well in straight and level flight -- see the
/// measured table in README.md -- because the residual mounting error is
/// then a TILT, which leaks into heading only when banked.
struct HeadingResult {
  bool ok = false;
  double heading = 0.0;        ///< estimated true heading, radians
  double error = 0.0;          ///< estimated minus AHRS yaw, radians
  double spread = 0.0;         ///< frame-to-frame scatter, radians
  double mean_bank = 0.0;      ///< the observability driver
  int n = 0;
};

HeadingResult celestialHeading(const std::vector<FrameData>& frames,
                               const Geodetic& pos,
                               const Eigen::Matrix3d& C_b_c);

Eigen::Matrix3d recalibrateMount(const std::vector<FrameData>& frames,
                                 const Geodetic& pos);

/// Mean local vertical, expressed in BODY axes, over a window of frames.
///
/// This is the axis an AHRS yaw bias rotates the estimated mounting about, so
/// it is the direction `mountTiltOnly` must remove. Level it is body z; banked
/// it tilts by the bank angle, and using body z instead would leave
/// `bias * sin(bank)` of the yaw behind AS TILT -- which lands in the boresight
/// and is far more damaging there. Returns body z if the window is empty.
Eigen::Vector3d meanVerticalBody(const std::vector<FrameData>& frames);



/// Keep the TILT of a calibrated mounting and discard the component aliased
/// with vehicle yaw.
///
/// `recalibrateMount` solves `C_b_c = C_l_b_est^T * R`, so any yaw error in the
/// supplied attitude is returned as part of the mounting. Substituting
/// `E = Rz(b) * T` and `R = T * M` gives
///
///     E^T R  =  [rotation by -b about C_l_b^T * z_ned] * M
///
/// so an AHRS yaw bias `b` appears as a body-frame rotation about the LOCAL
/// VERTICAL in body axes. A magnetometer bias and a rotated camera bracket are
/// then the same object and the calibration cannot prefer one; it returns the
/// mounting, because that is all it can return. Feed that mounting back to
/// `celestialHeading` and it reconstructs the biased attitude exactly, so the
/// compass reports no error however large the bias.
///
/// This projects the mounting's deviation from nominal onto the plane
/// perpendicular to `vertical_body` and drops the rest. What survives is the
/// TILT of the boresight, which is what the position solver needs and which is
/// not aliased with yaw; what goes is the clocking, leaving heading for the
/// stars to determine.
///
/// PASS THE MEASURED VERTICAL, not body z. Level the two coincide, but in a
/// bank they differ by the bank angle, and using body z leaves `b * sin(bank)`
/// of the yaw bias behind AS TILT -- which lands directly in the boresight and
/// is far more damaging there than the clocking error it was removing. At 9 deg
/// of bank and a 3 deg bias that is 0.47 deg of boresight error, against a
/// calibration that otherwise reaches 0.02 deg. Use `meanVerticalBody`.
///
/// The cost is that clocking is no longer estimated in flight: it is assumed
/// nominal. That is an assumption about the airframe, not a result, and heading
/// accuracy is conditional on it. See README.md.
Eigen::Matrix3d mountTiltOnly(const Eigen::Matrix3d& C_b_c,
                              const Eigen::Vector3d& vertical_body);

struct IterationConfig {
  /// Recalibrate the mounting between passes.
  ///
  /// TURN THIS OFF FOR CircleFit, which is what `iterateOrbit` now does by
  /// default. The circle fit already solves for the misalignment -- the fitted
  /// radius IS the mounting error -- so recalibrating from its output and
  /// re-fitting applies the same correction twice, and the fixed-point
  /// iteration can walk to a wrong answer.
  ///
  /// Measured across eight real ArduPilot SITL orbits, four GPS-aided and four
  /// denied: iteration 1 gives 6.64 km, "converged" gives 9.37 km, and
  /// iterating helps in only 2 of 8. Worst case observed 5.13 -> 15.19 km.
  bool recalibrate_mount = true;
  int max_iterations = 8;
  /// Stop when successive position estimates move less than this, metres.
  double convergence_m = 50.0;
  AverageMethod method = AverageMethod::Naive;
  Atmosphere refraction_model = Atmosphere::none();
};

struct IterationStep {
  Geodetic position;
  double error_m = 0.0;
  double mount_error_deg = 0.0;  ///< vs the true mounting, if known
  double circle_radius_km = 0.0;
};

/// The full converging loop of Figure 7: fix -> recalibrate -> fix -> ...
std::vector<IterationStep> iterateOrbit(const std::vector<FrameData>& frames,
                                        Eigen::Matrix3d C_b_c_assumed,
                                        const Geodetic& truth_for_scoring,
                                        const IterationConfig& cfg,
                                        const Eigen::Matrix3d* true_mount = nullptr);

}  // namespace celestial
