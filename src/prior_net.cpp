#include "celestial/prior_net.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#ifdef CNS_WITH_OPENCV
#include <opencv2/core.hpp>
#include <opencv2/dnn.hpp>
#include <opencv2/imgproc.hpp>
#endif

namespace celestial {

#ifndef CNS_WITH_OPENCV

// Stub build. Everything still compiles and links; callers get "no prior" and
// fall through to the ordinary detector, which is the correct degradation.
bool PriorNet::available() { return false; }
bool PriorNet::load(const std::string&) { return false; }
bool PriorNet::infer(const Image&, std::vector<float>& out) {
  out.clear();
  return false;
}

#else

bool PriorNet::available() { return true; }

bool PriorNet::load(const std::string& path) {
  delete static_cast<cv::dnn::Net*>(net_);
  net_ = nullptr;
  loaded_ = false;
  try {
    // NOTE the deliberate absence of setPreferableBackend/setPreferableTarget.
    //
    // OpenCV 5 ships three DNN engines behind one API and defaults to
    // ENGINE_AUTO: the new graph-based engine first, falling back to the
    // classic 4.x engine if a model will not load. The new engine is where the
    // 80%+ ONNX operator coverage, shape inference, constant folding and
    // operator fusion live, and it is CPU-only for now. Pinning a backend
    // pushes you onto the CLASSIC engine -- which is the right move only if
    // you want CUDA or OpenVINO, and is a straight loss on CPU.
    //
    // Leaving both unset means: best available engine on OpenCV 5, and
    // identical behaviour on 4.x where there is only one.
    auto* n = new cv::dnn::Net(cv::dnn::readNetFromONNX(path));
    if (n->empty()) {
      delete n;
      return false;
    }
    net_ = n;
    loaded_ = true;

    // WARM UP. cv::dnn initialises layers, allocates buffers and picks
    // implementations lazily on the first forward pass, so frame one costs
    // ~165 ms against ~7 ms steady state. On an aircraft that lands on
    // whichever frame happens to be first, which is not a cost anyone budgeted
    // for -- so pay it here, at load, where there is time.
    try {
      cv::Mat dummy = cv::Mat::zeros(64, 64, CV_32F);
      n->setInput(cv::dnn::blobFromImage(dummy));
      n->forward();
    } catch (const cv::Exception&) {
      // A model that cannot run a 64x64 probe may still be fine at the real
      // size; the warm-up is an optimisation, not a validation.
    }
  } catch (const cv::Exception& e) {
    std::fprintf(stderr, "PriorNet: cannot load %s: %s\n", path.c_str(),
                 e.what());
    return false;
  }
  return true;
}

bool PriorNet::infer(const Image& img, std::vector<float>& out) {
  out.clear();
  if (!loaded_ || !net_) return false;

  const int d = std::max(1, cfg_.decimate);
  const int sw = img.width / d, sh = img.height / d;
  if (sw < 32 || sh < 32) return false;

  // Wrap the 16-bit frame without copying, then decimate by area averaging --
  // the reduction that preserves total flux, which is what the network was
  // trained on.
  const cv::Mat raw(img.height, img.width, CV_16UC1,
                    const_cast<uint16_t*>(img.data.data()));
  cv::Mat small;
  cv::resize(raw, small, cv::Size(sw, sh), 0, 0, cv::INTER_AREA);
  small.convertTo(small, CV_32F);

  // SAME normalisation the trainer uses: per-image median and MAD. The
  // background level moves with cloud, flare and exposure, so a model handed
  // raw ADU would be learning offsets rather than what a star looks like. If
  // this ever diverges from scripts/train_unet.py the model silently degrades.
  // nth_element, not cv::sort: the median needs a selection, not an ordering,
  // and sorting 150k floats twice was a measurable slice of the inference time.
  std::vector<float> v(small.begin<float>(), small.end<float>());
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  const float med = v[v.size() / 2];
  for (float& x : v) x = std::fabs(x - med);
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  const float mad = std::max(1e-6f, v[v.size() / 2]);
  small = (small - med) / (1.4826f * mad);

  cv::Mat prob;
  try {
    auto* net = static_cast<cv::dnn::Net*>(net_);
    net->setInput(cv::dnn::blobFromImage(small));
    const cv::Mat logits = net->forward();
    // NCHW with N=C=1; view it as a plain 2-D image.
    const cv::Mat l(logits.size[2], logits.size[3], CV_32F,
                    const_cast<float*>(logits.ptr<float>()));
    cv::exp(-l, prob);
    prob = 1.0f / (1.0f + prob);  // sigmoid; the model exports raw logits
  } catch (const cv::Exception& e) {
    std::fprintf(stderr, "PriorNet: inference failed: %s\n", e.what());
    return false;
  }

  cv::Mat full;
  cv::resize(prob, full, cv::Size(img.width, img.height), 0, 0,
             cv::INTER_LINEAR);

  out.assign(size_t(img.width) * img.height, 0.f);
  for (int y = 0; y < img.height; ++y) {
    const float* row = full.ptr<float>(y);
    for (int x = 0; x < img.width; ++x)
      out[size_t(y) * img.width + x] =
          std::clamp(row[x], 0.0f, cfg_.max_prior);
  }
  return true;
}

#endif

}  // namespace celestial
