// imaging.hpp — MILESTONE 2, first version.
//
// Render -> detect -> centroid -> match -> the Milestone 1 solver.
//
// The point of this stage is NOT to produce pretty pictures. It is to put a
// number on the cost of the imaging chain: everything downstream of the
// centroids is already validated, so any degradation against the Milestone 1
// ideal-vector result is attributable to the front end alone.
//
// Deliberately dependency-free (no OpenCV yet). The connected-component
// labeller here is ~40 lines; swap in cv::connectedComponentsWithStats when
// the pipeline moves into a ROS 2 node, where OpenCV arrives free with
// cv_bridge anyway.

#pragma once

#include <cstdint>
#include <vector>

#include "celestial/orbit.hpp"
#include "celestial/types.hpp"

namespace celestial {

// ---------------------------------------------------------------------------
// Sensor
// ---------------------------------------------------------------------------

/// Photometric and noise model of the imaging chain.
///
/// Defaults approximate the paper's payload: Alvium 1800 U-240 with a 6 mm
/// f/1.4 lens on a moonless night.
struct SensorModel {
  // --- optics / photometry -------------------------------------------------
  /// Electrons per second from a V = 0 star. Derived from
  ///   1000 photons/s/cm^2/Angstrom at V=0
  ///   x aperture area (6 mm / f1.4 -> D = 4.29 mm -> 0.144 cm^2)
  ///   x bandwidth 2000 A x QE 0.5 x transmission 0.8
  /// Gives ~1.15e5 e/s at V=0, so ~2900 e/s at V=4.
  double zeropoint_e_per_s = 1.15e5;

  /// Gaussian PSF sigma in pixels.
  ///
  /// NOTE: star trackers are deliberately DEFOCUSED. At 107 arcsec/px a
  /// focused 6 mm lens puts the whole star inside one pixel, which destroys
  /// all subpixel information. Spreading it over ~2-3 px is what makes
  /// centroiding to 0.1 px possible.
  double psf_sigma_px = 1.2;

  // --- timing --------------------------------------------------------------
  /// Exposure, seconds. The paper never states this, yet it sets the motion
  /// blur completely: at ~2.4 deg/s yaw rate and 107 arcsec/px, 100 ms is
  /// ~8 px of smear.
  double exposure_s = 0.05;
  /// Sub-steps used to integrate the smear across the exposure. The streak
  /// CURVES during a turn, so a directional kernel is not equivalent.
  int blur_substeps = 12;

  /// DIAGNOSTIC ONLY. If > 0, photon collection (stars, sky, dark) uses this
  /// exposure while the SMEAR still uses exposure_s.
  ///
  /// Exposure physically sets both signal and blur at once, so sweeping it
  /// measures the combined trade rather than blur itself. Pinning the signal
  /// here isolates the blur mechanism. Leave at 0 for physical behaviour.
  double photometric_exposure_s = 0.0;

  // --- noise ---------------------------------------------------------------
  double read_noise_e = 6.0;
  double dark_e_per_s = 5.0;
  /// Sky background. A moonless dark site is ~21.5 mag/arcsec^2, which over a
  /// 107 arcsec pixel is ~11.3 mag/px, i.e. only a few e/s. Raise this to
  /// explore twilight or daylight.
  double sky_mag_per_arcsec2 = 21.5;
  bool shot_noise = true;
  double prnu_frac = 0.005;      ///< pixel response non-uniformity, 1-sigma
  double hot_pixel_rate = 2e-5;  ///< fraction of pixels with elevated dark
  double hot_pixel_e_per_s = 400.0;

  // --- digitisation --------------------------------------------------------
  double gain_e_per_adu = 0.5;
  double full_well_e = 15000.0;
  double bias_adu = 100.0;
  int bit_depth = 12;

  // --- weather (for the learned-detector experiment) ------------------------
  //
  // Motion blur is a KNOWN signal shape and can be handled deterministically.
  // Cloud cannot: it is structured, non-stationary and looks
  // like signal at some scales -- and cloud is precisely where Teague & Chahl
  // (2026) report a UNet holding F1 0.77 while every classical baseline drops
  // below 0.26. Nothing can be evaluated on that claim until the simulator has
  // weather, so here it is. Deliberately crude; see README.md.
  //
  // BOTH DEFAULT TO ZERO, so every existing number is unchanged.

  /// Cloud opacity, 0-1. Cloud is not an overlay: it ATTENUATES starlight
  /// multiplicatively and RAISES the background additively, because it both
  /// blocks what is behind it and scatters ground light back at you.
  double cloud_amount = 0.0;
  /// Spatial scale of the cloud field, pixels.
  double cloud_scale_px = 300.0;
  /// Background lift through fully opaque cloud, e/s per pixel.
  double cloud_glow_e_per_s = 60.0;

  /// DIAGNOSTIC. Cloud does two separable things: it ATTENUATES starlight
  /// (photons that never arrive, unrecoverable) and it ADDS glow (a background
  /// gradient, in principle removable). Turning each off independently is what
  /// decides whether background estimation is worth building -- see
  /// README.md.
  bool cloud_attenuates = true;

  /// Lens flare, 0-1: broad radial gradients from a bright off-axis source,
  /// fixed in the CAMERA frame so they do not move with the stars. This is
  /// what destroys single-frame adaptive thresholding, and what a temporal
  /// median removes almost for free.
  double flare_amount = 0.0;

  unsigned seed = 12345;

  /// Everything off: no noise, no blur, no fixed pattern. The starting point
  /// for the closure test.
  static SensorModel ideal();
};

/// A rendered frame, 16-bit.
struct Image {
  int width = 0;
  int height = 0;
  std::vector<uint16_t> data;
  uint16_t at(int x, int y) const { return data[y * width + x]; }
};

/// Ground truth for one star in one frame. Emit these from day one: they are
/// the labels for any learned centroider or cloud masker, and they cost
/// nothing to produce now.
struct StarLabel {
  int id = 0;
  double vmag = 0.0;
  double u = 0.0, v = 0.0;    ///< flux-weighted true centroid over the exposure
  double u0 = 0, v0 = 0;      ///< streak start (exposure open)
  double u1 = 0, v1 = 0;      ///< streak end (exposure close)
  double total_e = 0.0;
  bool truncated = false;     ///< streak partly outside the sensor
};

struct RenderedFrame {
  Image image;
  std::vector<StarLabel> labels;
};

/// Render one frame. Motion blur is integrated between poses `a` and `b`,
/// covering the fraction of that interval given by the exposure.
///
/// `field` must be built from a catalogue already limited to the detector's
/// magnitude depth.
RenderedFrame renderFrame(const FrameTruth& a, const FrameTruth& b,
                          double frame_period_s, const StarField& field,
                          const std::vector<double>& vmags, const Camera& cam,
                          const Eigen::Matrix3d& C_b_c_true,
                          const SensorModel& sensor, unsigned frame_seed);

// ---------------------------------------------------------------------------
// Detection
// ---------------------------------------------------------------------------

/// How to turn a detected blob into a subpixel position.
enum class CentroidMethod {
  /// Background-subtracted intensity-weighted mean over the thresholded
  /// component. Cheap, unbiased for a symmetric PSF that is fully above
  /// threshold, but biased when the threshold cuts into the wings.
  WholeComponentCoG,
  /// Gauss-Newton fit of A*exp(-((x-x0)^2+(y-y0)^2)/2s^2) + B over a fixed
  /// window, with A, x0, y0, s, B all free.
  ///
  /// If the PSF really is Gaussian and the noise model is right, this is the
  /// maximum-likelihood estimator and no method -- learned or otherwise -- can
  /// beat it, because it attains the Cramer-Rao bound.
  GaussianFit,
};

struct DetectorConfig {
  /// Threshold at background + k * sigma. The paper uses k = 5 with the frame
  /// MEAN and standard deviation; we default to median + MAD, which is robust
  /// to the stars themselves and to a bright corner from vignetting.
  double threshold_k = 5.0;
  bool robust_background = true;
  int min_pixels = 2;
  int max_pixels = 400;
  /// GaussianFit by default: it is 7-10x better than the paper's 3x3 CoG and
  /// reaches the Cramer-Rao bound (testCramerRao). Every real caller already
  /// set it explicitly; the default was the one that loses.
  CentroidMethod centroid = CentroidMethod::GaussianFit;
  /// Half-width of the fixed fit window, pixels. Must comfortably contain the
  /// PSF wings; too small reintroduces truncation bias.
  int fit_halfwidth = 5;

  // --- mesh background ------------------------------------------------------
  //
  // Cell size in pixels for a SExtractor-style background mesh: a sigma-clipped
  // estimate per cell, median-filtered across nodes so one bright cell cannot
  // poison its neighbourhood, bilinearly interpolated and subtracted.
  //
  // This exists because the global median+MAD threshold assumes a FLAT
  // background. A smooth gradient inflates the global MAD, the threshold rises
  // across the whole frame, and faint stars are lost even in the clean parts.
  // That loss is self-inflicted and a local background removes it.
  //
  // Measured, matched stars per frame (see testMeshBackground):
  //
  //   condition      plain   mesh    gain
  //   clear           38.3   36.8   0.96x
  //   cloud 0.6       25.8   26.2   1.01x
  //   flare 0.5       19.3   34.0   1.76x
  //   cloud + flare   13.8   24.2   1.75x
  //
  // Note WHICH contaminant it fixes. Lens flare is a steep additive gradient
  // and the mesh recovers nearly all of it. Cloud it does not help with at all,
  // because cloud's damage is not threshold bias: it ATTENUATES starlight, and
  // its glow raises the SHOT NOISE floor. Subtracting a background cannot
  // un-add Poisson noise or recover photons that never arrived.
  //
  // 0 disables. Costs ~1 ms and is slightly negative in clean conditions, so
  // it is off by default and worth enabling when flare or sky gradients are
  // expected.
  int bg_mesh_px = 0;

  // --- detection prior (for a learned detector) -----------------------------
  //
  // A segmentation network can HALLUCINATE a star and can DELETE a real one.
  // Neither is acceptable in a navigation sensor, and neither should be
  // guarded against by trusting the training. So the network does not get to
  // decide anything: it supplies a PRIOR, and the threshold becomes spatially
  // varying,
  //
  //     k(x,y) = threshold_k - (threshold_k - threshold_k_low) * prior(x,y)
  //
  // Where the network is confident the bar drops to `threshold_k_low`; where
  // it is not, the bar is unchanged. The detection still has to clear a real
  // significance floor computed from ACTUAL PHOTONS in the original image, so
  // a confident prediction over empty sky yields nothing. The network can
  // lower the bar. It cannot remove it.
  //
  // The mirror failure -- the network erasing a real star -- is handled by
  // running the ordinary detector too and taking the UNION, so the prior can
  // only ever ADD candidates.
  //
  // Both guarantees are structural. Neither depends on the model being good.
  // `testHostilePrior` asserts them against a deliberately adversarial prior.
  //
  // Empty disables. Size must match the image, row-major, 0-1.
  std::vector<float> prior;
  double threshold_k_low = 3.0;
};

struct Detection {
  double u = 0.0, v = 0.0;  ///< background-subtracted weighted centroid
  double flux = 0.0;        ///< summed background-subtracted counts
  int n_pixels = 0;
  double elongation = 1.0;  ///< major/minor axis ratio; > ~2 means a streak
};

/// Threshold, label connected components, centroid each one.
std::vector<Detection> detectStars(const Image& img, const DetectorConfig& cfg);

// ---------------------------------------------------------------------------
// Identification (tracking mode)
// ---------------------------------------------------------------------------

/// Match detections to catalogue stars by projecting the catalogue through the
/// ESTIMATED attitude and camera mounting.
///
/// This is the mode that runs 99.9% of the time. A lost-in-space solver
/// (log-polar, pyramid, tetra3) is only needed to bootstrap and to recover:
/// once you have any attitude estimate, even a 3 deg error is ~100 px at
/// 107 arcsec/px, which nearest-neighbour matching handles easily.
///
/// Uses no ground truth.
struct MatcherConfig {
  double radius_px = 40.0;
  /// Reject a detection if a second candidate lies within this factor of the
  /// best one; ambiguous matches are worse than no match.
  double ambiguity_ratio = 2.0;
};

/// Match a set of already-computed detections. Detections do not depend on the
/// assumed mounting, so cache them and call this once per iteration rather
/// than re-detecting.
int matchDetections(const std::vector<Detection>& dets, const Epoch& epoch,
                    const Eigen::Matrix3d& C_l_b_est,
                    const Eigen::Matrix3d& C_b_c_assumed, const Camera& cam,
                    const StarField& field, const Geodetic& assumed_pos,
                    const MatcherConfig& mcfg, FrameData& out);

/// Build solver input from an image, given the AHRS attitude and an assumed
/// mounting. Returns the number of stars matched.
int buildFrameData(const Image& img, const Epoch& epoch,
                   const Eigen::Matrix3d& C_l_b_est,
                   const Eigen::Matrix3d& C_b_c_assumed, const Camera& cam,
                   const StarField& field, const Geodetic& assumed_pos,
                   const DetectorConfig& dcfg, const MatcherConfig& mcfg,
                   FrameData& out);

/// Camera-frame unit vector for a pixel.
Eigen::Vector3d pixelToRay(double u, double v, const Camera& cam);
/// Pixel coordinates of a camera-frame direction. Returns false if behind.
bool rayToPixel(const Eigen::Vector3d& v_cam, const Camera& cam, double& u,
                double& v);

}  // namespace celestial
