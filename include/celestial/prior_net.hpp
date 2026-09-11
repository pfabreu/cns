// ---------------------------------------------------------------------------
// PRIOR NET — a learned detection prior, with the classical guarantee intact.
//
// WHAT THIS IS NOT. It is not a star detector. Teague & Chahl (2026) segment
// stars with a UNet and use the segmentation as the detection; that works, and
// it means a bad model can invent a star that never existed or erase one that
// did. In a navigation sensor neither is acceptable, and neither should be
// guarded against by trusting the training.
//
// WHAT IT IS. The network outputs a PRIOR, which is fed to
// `DetectorConfig::prior` and lowers the detection threshold between
// `threshold_k` and `threshold_k_low`:
//
//     k(x,y) = k_high - (k_high - k_low) * p(x,y)
//
// The detection still has to clear a real significance floor computed from
// actual photons, and the centroid is still measured on the original image. So
// the model CANNOT hallucinate -- confidence over empty sky yields nothing,
// there being no flux to exceed even k_low -- and CANNOT delete, since a zero
// prior leaves the nominal threshold untouched. `testPriorIsBounded` asserts
// the exact equivalence: prior=1 is identical to hand-setting k_low, prior=0 is
// identical to the ordinary detector. The worst an adversarial model can do is
// move you to an operating point you could have picked yourself.
//
// That bound is only worth anything because the worst case is ACCEPTABLE, so
// treat `threshold_k_low` as a safety budget rather than a tuning knob: set it
// to the most permissive threshold you would accept with no network at all.
//
// WHY THIS DECOMPOSITION IS ALSO THE STATISTICALLY CORRECT ONE. Detection is a
// hypothesis test, and the optimal test compares a likelihood ratio against a
// threshold set by the prior odds. A fixed global threshold is optimal only if
// the prior is spatially uniform, and it is not -- detectability varies with
// cloud, flare and local noise. So:
//
//     prior       <- the network. P(star | context). Wide receptive field,
//                    structured background, "is this region attenuated?"
//     likelihood  <- the matched filter. P(data | star). Provably optimal for
//                    a known signal shape in Gaussian noise; there is nothing
//                    for a network to add.
//     decision    <- the threshold, combining them.
//
// A network that outputs detections directly is using the prior AS the
// posterior and discarding the likelihood. That is simultaneously the unsafe
// design and the wrong one.
//
// TWO HONEST CAVEATS. The model is not calibrated -- Dice+BCE does not produce
// probabilities, so the linear map from p to k is a heuristic wearing Bayesian
// clothes. And the prior is computed FROM the image, including the pixel under
// test, so evidence is partly double counted and the product is not strictly a
// posterior. Neither weakens the safety bound, which is mechanical; they mean
// the optimality story is looser than the safety story.
//
// RESOLUTION. Inference runs on a DECIMATED image and the prior is upsampled.
// That is not only for speed: a coarse prior cannot target an individual pixel,
// so the model is structurally incapable of fine-grained mischief even within
// its permitted range.
//
// BUILD. Requires OpenCV (`cv::dnn`), which is optional. Without it this
// compiles to a stub that reports unavailable, and every other build target is
// unaffected. See scripts/README.md for training and export.
// ---------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

#include "celestial/imaging.hpp"

namespace celestial {

class PriorNet {
 public:
  struct Config {
    /// Inference runs at 1/decimate resolution. 4 keeps a 2.35 Mpx frame down
    /// to ~150 kpx, which is what makes this affordable inside the frame
    /// budget -- check ./build/bench.
    int decimate = 4;
    /// Clamp on the prior, so a saturated model cannot pin the whole frame to
    /// k_low. Belt and braces: the threshold floor already bounds the damage.
    float max_prior = 1.0f;
  };

  PriorNet() = default;
  explicit PriorNet(const Config& cfg) : cfg_(cfg) {}

  /// True if this build has OpenCV and can run a model at all.
  static bool available();

  /// Load an ONNX model. Returns false if unavailable or the file will not
  /// load; callers should carry on without a prior rather than fail.
  bool load(const std::string& onnx_path);
  bool loaded() const { return loaded_; }

  /// Produce a full-resolution prior for `img`, sized width*height, values
  /// 0-1, ready to assign to `DetectorConfig::prior`. Returns false if no
  /// model is loaded, in which case `out` is cleared and the caller simply
  /// gets the ordinary detector.
  bool infer(const Image& img, std::vector<float>& out);

 private:
  Config cfg_;
  bool loaded_ = false;
  void* net_ = nullptr;  ///< cv::dnn::Net, type-erased to keep OpenCV out of
                         ///< this header and out of every other translation
                         ///< unit.
};

}  // namespace celestial
