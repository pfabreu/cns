// test_imaging.cpp — MILESTONE 2, first version.
//
// The closure test: render -> detect -> centroid -> match -> the SAME solver
// used in Milestone 1, on the SAME trajectory. Everything downstream of the
// centroids is already validated, so any degradation is attributable to the
// imaging chain alone.
//
// Start from a perfect sensor (must reproduce Milestone 1 exactly), then turn
// on one effect at a time.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <utility>
#include <cstdlib>
#include <string>
#include <vector>

#include "celestial/imaging.hpp"
#include "celestial/prior_net.hpp"
#include "celestial/orbit.hpp"
#include "celestial/sky_model.hpp"
#include "celestial/star_catalog.hpp"

using namespace celestial;

namespace {

int g_failures = 0;
constexpr double kEarthR = 6371008.8;

void check(bool cond, const std::string& what) {
  if (!cond) {
    std::printf("  ** FAIL: %s\n", what.c_str());
    ++g_failures;
  }
}

struct Rig {
  OrbitConfig cfg;
  std::vector<FrameTruth> truth;
  Camera cam;
  Eigen::Matrix3d C_b_c_true;
  std::vector<CatalogStar> catalog;
  std::vector<double> vmags;
  double boresight_deg = 0.4;
};

Rig makeRig(double mag_limit = 4.0) {
  Rig r;
  r.cfg.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  r.cfg.radius_m = 600.0;
  r.cfg.airspeed = 25.0;
  // 2 Hz rather than 10 Hz: one full revolution is then ~300 frames instead of
  // 1508, which keeps the sweep tractable. Blur is set by EXPOSURE, not frame
  // rate, so nothing physical is lost.
  r.cfg.frame_rate = 2.0;
  r.cfg.revolutions = 1.0;
  r.truth = generateOrbit(r.cfg);
  r.C_b_c_true =
      nominalCameraMount() *
      smallRotation(Eigen::Vector3d(1, 0, 0), r.boresight_deg * kDeg2Rad);
  r.catalog = catalogBrighterThan(mag_limit);
  for (const CatalogStar& s : r.catalog) r.vmags.push_back(s.vmag);
  return r;
}

struct RunStats {
  double err_km = -1.0;
  double mean_detections = 0.0;
  double mean_matched = 0.0;
  /// Fraction of in-frame catalogue stars actually detected.
  ///
  /// REPORT THIS ALONGSIDE ANY CENTROID METRIC. A centroid RMS computed only
  /// over detected stars is optimistically biased: at low completeness the
  /// survivors are exactly those that caught a favourable subpixel phase or a
  /// positive noise fluctuation, so the metric measures the SELECTION rather
  /// than the estimator. Ignoring this is how you appear to beat the
  /// Cramer-Rao bound.
  double completeness = 0.0;
  double centroid_rms_px = 0.0;       ///< all label-associated, 10 px gate
  double centroid_rms_gated = -1.0;   ///< the same, truncated at 2 px
  double centroid_rms_bright = -1.0;  ///< V <= 3
  double centroid_rms_faint = -1.0;   ///< V >  3
  int iterations = 0;
};

/// One full pass: render every frame, detect, match, solve, iterate.
RunStats runPipeline(const Rig& rig, const SensorModel& sensor,
                     const DetectorConfig& dcfg, const MatcherConfig& mcfg) {
  const StarField field(rig.catalog, rig.truth[rig.truth.size() / 2].epoch);
  const double period = 1.0 / rig.cfg.frame_rate;

  // Render and detect ONCE. Neither depends on the assumed mounting, so the
  // iteration loop below only re-runs the matcher.
  const int kStride = 4;  // ~75 frames over a full revolution
  std::vector<std::vector<Detection>> detections;
  std::vector<size_t> frame_idx;
  RunStats st;
  double s2 = 0.0, s2b = 0.0, s2f = 0.0, s2g = 0.0;
  int nc = 0, ndet = 0, nb = 0, nf = 0, ng = 0;
  int n_labels = 0, n_found = 0;

  for (size_t k = 0; k + 1 < rig.truth.size(); k += kStride) {
    const RenderedFrame rf =
        renderFrame(rig.truth[k], rig.truth[k + 1], period, field, rig.vmags,
                    rig.cam, rig.C_b_c_true, sensor, static_cast<unsigned>(k));
    auto dets = detectStars(rf.image, dcfg);
    ndet += static_cast<int>(dets.size());
    // LABEL-driven association, so that undetected stars are counted as
    // misses rather than silently vanishing from the statistic. A generous
    // 10 px gate keeps the error tail visible; the 2 px version is reported
    // separately to show how much the tighter gate flatters the result.
    for (const StarLabel& L : rf.labels) {
      if (L.truncated) continue;
      if (L.u < 20 || L.v < 20 || L.u > rig.cam.width - 20 ||
          L.v > rig.cam.height - 20) continue;
      ++n_labels;
      double bd = 1e18;
      for (const Detection& d : dets) {
        const double r2 = (L.u - d.u) * (L.u - d.u) + (L.v - d.v) * (L.v - d.v);
        if (r2 < bd) bd = r2;
      }
      if (bd > 100.0) continue;  // not detected
      ++n_found;
      s2 += bd; ++nc;
      if (bd < 4.0) { s2g += bd; ++ng; }
      if (L.vmag <= 3.0) { s2b += bd; ++nb; } else { s2f += bd; ++nf; }
    }
    detections.push_back(std::move(dets));
    frame_idx.push_back(k);
  }
  st.mean_detections = detections.empty() ? 0 : double(ndet) / detections.size();
  st.centroid_rms_px = nc ? std::sqrt(s2 / nc) : -1.0;
  st.centroid_rms_gated = ng ? std::sqrt(s2g / ng) : -1.0;
  st.centroid_rms_bright = nb ? std::sqrt(s2b / nb) : -1.0;
  st.centroid_rms_faint = nf ? std::sqrt(s2f / nf) : -1.0;
  st.completeness = n_labels ? double(n_found) / n_labels : 0.0;

  // Iterate: assumed mounting starts at nominal, i.e. 0.4 deg wrong.
  Eigen::Matrix3d C_b_c = nominalCameraMount();
  Geodetic assumed = rig.cfg.centre;
  double last_err = -1.0;
  int total_matched = 0, frames_used = 0;

  for (int it = 0; it < 4; ++it) {
    std::vector<FrameData> frames;
    frames.reserve(detections.size());
    total_matched = 0;
    for (size_t j = 0; j < detections.size(); ++j) {
      const FrameTruth& ft = rig.truth[frame_idx[j]];
      // AHRS attitude taken as truth here: attitude error is explored in
      // Milestone 1, and this test isolates the IMAGING chain.
      const Eigen::Matrix3d C_l_b_est = eulerToDcm(ft.roll, ft.pitch, ft.yaw);
      FrameData fd;
      const int m = matchDetections(detections[j], ft.epoch, C_l_b_est, C_b_c,
                                    rig.cam, field, assumed, mcfg, fd);
      fd.t = ft.t;
      fd.yaw_est = ft.yaw;
      fd.truth = ft.pos;
      total_matched += m;
      if (m >= 3) frames.push_back(fd);
    }
    frames_used = static_cast<int>(frames.size());
    if (frames.size() < 15) break;

    const OrbitResult r =
        estimateOrbit(frames, C_b_c, AverageMethod::HeadingWeighted);
    if (!r.ok) break;
    assumed = r.position;
    last_err = haversine(r.position, rig.cfg.centre);
    st.iterations = it + 1;
    C_b_c = recalibrateMount(frames, r.position);
  }

  st.mean_matched =
      detections.empty() ? 0 : double(total_matched) / detections.size();
  st.err_km = last_err / 1000.0;
  (void)frames_used;
  return st;
}

void report(const char* label, const RunStats& s) {
  std::printf("  %-32s %5.1f %5.1f %6.0f%% %8.3f %8.3f %8.3f %3d\n", label,
              s.mean_detections, s.mean_matched, 100.0 * s.completeness,
              s.centroid_rms_px, s.centroid_rms_gated, s.err_km, s.iterations);
}

void header() {
  std::printf("\n  %-32s %5s %5s %7s %8s %8s %8s %3s\n", "configuration", "det",
              "match", "compl", "cent px", "<2px", "err km", "it");
  std::printf("  %s\n", std::string(90, '-').c_str());
}

// ---------------------------------------------------------------------------

void testIdeal() {
  std::printf("\n=== TEST 1: perfect sensor must reproduce Milestone 1 ===\n");
  const Rig rig = makeRig();
  std::printf("  %zu frames, %.1f deg boresight error, %zu catalogue stars\n",
              rig.truth.size(), rig.boresight_deg, rig.catalog.size());

  SensorModel s = SensorModel::ideal();
  DetectorConfig d;
  MatcherConfig m;
  header();
  const RunStats st = runPipeline(rig, s, d, m);
  report("ideal (no noise, no blur)", st);

  check(st.err_km >= 0 && st.err_km < 0.5,
        "ideal sensor reproduces the Milestone 1 result");
  std::printf("\n  centroid RMS by brightness:  V<=3  %.3f px   V>3  %.3f px\n",
              st.centroid_rms_bright, st.centroid_rms_faint);
  std::printf("  completeness %.0f%% -- with a perfect sensor every catalogue\n"
              "  star is found, so this metric is unbiased here.\n",
              100.0 * st.completeness);
  check(st.centroid_rms_bright >= 0 && st.centroid_rms_bright < 0.08,
        "bright-star centroids are better than 0.08 px with a perfect sensor");
  std::printf("\n  With NO noise at all the residual is threshold-truncation\n"
              "  bias: the k-sigma cut keeps a different fraction of the PSF\n"
              "  depending on the star's brightness and its subpixel phase, so\n"
              "  faint stars are systematically worse. This is a bias, not\n"
              "  noise, and it is exactly what a learned centroider trained on\n"
              "  these labels would remove.\n");
  std::printf("\n  A defocused PSF (sigma = %.1f px) is what makes subpixel\n"
              "  centroiding possible at all: a focused 6 mm lens at 107\n"
              "  arcsec/px would put the whole star inside one pixel.\n",
              s.psf_sigma_px);
}

void testNoiseLadder() {
  std::printf("\n=== TEST 2: turn on one noise source at a time ===\n");
  const Rig rig = makeRig();
  DetectorConfig d;
  MatcherConfig m;
  header();

  SensorModel s = SensorModel::ideal();
  report("ideal", runPipeline(rig, s, d, m));

  s = SensorModel::ideal(); s.read_noise_e = 6.0;
  report("+ read noise 6 e", runPipeline(rig, s, d, m));

  s.shot_noise = true;
  report("+ shot noise", runPipeline(rig, s, d, m));

  s.dark_e_per_s = 5.0; s.sky_mag_per_arcsec2 = 21.5;
  report("+ dark & sky", runPipeline(rig, s, d, m));

  s.prnu_frac = 0.005;
  report("+ PRNU 0.5%", runPipeline(rig, s, d, m));

  s.hot_pixel_rate = 2e-5; s.hot_pixel_e_per_s = 400.0;
  report("+ hot pixels (2e-5)", runPipeline(rig, s, d, m));

  std::printf("\n  Hot pixels are a FIXED pattern: the same pixels fire every\n"
              "  frame and track perfectly. They are exactly what the matcher's\n"
              "  ambiguity guard and the solver's RANSAC exist to reject.\n");
}

void testBlurLadder() {
  std::printf("\n=== TEST 3: motion blur, the effect the paper never states ===\n");
  const Rig rig = makeRig();
  DetectorConfig d;
  MatcherConfig m;

  const double yaw_rate = rig.cfg.airspeed / rig.cfg.radius_m;  // rad/s
  const double px_per_s = yaw_rate / rig.cam.pixelIfov();
  std::printf("  yaw rate %.2f deg/s -> %.1f px/s of image motion\n",
              yaw_rate * kRad2Deg, px_per_s);
  header();

  std::printf("\n  (a) BLUR ISOLATED: photon count pinned at 50 ms, smear varied.\n");
  header();
  for (double exp_ms : {0.0, 25.0, 50.0, 100.0, 200.0}) {
    SensorModel s = SensorModel::ideal();
    s.read_noise_e = 6.0;
    s.shot_noise = true;
    s.dark_e_per_s = 5.0;
    s.exposure_s = exp_ms / 1000.0;
    s.photometric_exposure_s = 0.05;  // pinned
    s.blur_substeps = (exp_ms > 0) ? 16 : 1;
    char buf[80];
    std::snprintf(buf, sizeof buf, "smear %5.0f ms  (%4.1f px)", exp_ms,
                  px_per_s * exp_ms / 1000.0);
    report(buf, runPipeline(rig, s, d, m));
  }

  std::printf("\n  (b) PHYSICAL: exposure sets photon count AND smear together.\n");
  header();
  for (double exp_ms : {10.0, 25.0, 50.0, 100.0, 200.0}) {
    SensorModel s = SensorModel::ideal();
    s.read_noise_e = 6.0;
    s.shot_noise = true;
    s.dark_e_per_s = 5.0;
    s.exposure_s = exp_ms / 1000.0;
    s.blur_substeps = 16;
    char buf[80];
    std::snprintf(buf, sizeof buf, "exposure %5.0f ms  (%4.1f px)", exp_ms,
                  px_per_s * exp_ms / 1000.0);
    report(buf, runPipeline(rig, s, d, m));
  }

  std::printf("\n  Ladder (a) shows blur alone COSTS detections: the same photons\n"
              "  spread over more pixels drop below a k-sigma threshold.\n"
              "  Ladder (b) shows the real trade, where the extra integration\n"
              "  time more than pays for the smear at these rates. That is why\n"
              "  the paper can run a 21x21 tracking box at all -- but it also\n"
              "  means exposure is a tuning parameter they never report.\n");
}

void testCentroidBiasUnderBlur() {
  std::printf("\n=== TEST 4: centroid methods compared ===\n");
  std::printf("  Eq. (12) uses a fixed 3x3 window with no background\n"
              "  subtraction. Under blur the star is a streak several pixels\n"
              "  long, so a 3x3 window does not contain it.\n");
  std::printf("\n  Label-driven association, 10 px gate, completeness reported.\n");
  const Rig rig = makeRig();
  const StarField field(rig.catalog, rig.truth[rig.truth.size() / 2].epoch);
  const double period = 1.0 / rig.cfg.frame_rate;

  std::printf("\n  %-10s %-8s %-8s %-11s %-11s %-11s\n", "exposure",
              "smear", "compl", "3x3 CoG", "whole-comp", "Gauss fit");
  std::printf("  %s\n", std::string(60, '-').c_str());

  for (double exp_ms : {0.0, 25.0, 50.0, 100.0}) {
    SensorModel s = SensorModel::ideal();
    s.read_noise_e = 6.0;
    s.shot_noise = true;
    s.exposure_s = exp_ms / 1000.0;
    s.photometric_exposure_s = 0.05;  // pin photons: isolate the CENTROIDER
    s.blur_substeps = (exp_ms > 0) ? 16 : 1;

    DetectorConfig d_cog;
    d_cog.centroid = CentroidMethod::WholeComponentCoG;
    DetectorConfig d_fit;
    d_fit.centroid = CentroidMethod::GaussianFit;

    double s2_cog = 0, s2_3x3 = 0, s2_fit = 0;
    int n = 0, n_labels = 0;

    for (size_t k = 0; k + 1 < rig.truth.size(); k += 7) {
      const RenderedFrame rf =
          renderFrame(rig.truth[k], rig.truth[k + 1], period, field, rig.vmags,
                      rig.cam, rig.C_b_c_true, s, static_cast<unsigned>(k));
      const auto dc = detectStars(rf.image, d_cog);
      const auto df = detectStars(rf.image, d_fit);

      for (const StarLabel& L : rf.labels) {
        if (L.truncated) continue;
        if (L.u < 20 || L.v < 20 || L.u > rig.cam.width - 20 ||
            L.v > rig.cam.height - 20) continue;
        ++n_labels;

        auto nearest = [&](const std::vector<Detection>& v, double& du,
                           double& dv) {
          double bd = 1e18;
          for (const Detection& d : v) {
            const double r2 =
                (L.u - d.u) * (L.u - d.u) + (L.v - d.v) * (L.v - d.v);
            if (r2 < bd) { bd = r2; du = d.u; dv = d.v; }
          }
          return bd;
        };
        double cu = 0, cv = 0, fu = 0, fv = 0;
        const double bd_c = nearest(dc, cu, cv);
        const double bd_f = nearest(df, fu, fv);
        if (bd_c > 100.0 || bd_f > 100.0) continue;  // not detected
        ++n;

        // 3x3 CoG about the local peak, as in Eq. (12): raw counts, no
        // background subtraction.
        int px = std::clamp(int(cu), 1, rf.image.width - 2);
        int py = std::clamp(int(cv), 1, rf.image.height - 2);
        uint16_t peak = 0;
        for (int yy = py - 4; yy <= py + 4; ++yy)
          for (int xx = px - 4; xx <= px + 4; ++xx) {
            if (xx < 1 || yy < 1 || xx >= rf.image.width - 1 ||
                yy >= rf.image.height - 1) continue;
            if (rf.image.at(xx, yy) > peak) {
              peak = rf.image.at(xx, yy); px = xx; py = yy;
            }
          }
        double w = 0, uw = 0, vw = 0;
        for (int yy = py - 1; yy <= py + 1; ++yy)
          for (int xx = px - 1; xx <= px + 1; ++xx) {
            const double val = rf.image.at(xx, yy);
            w += val; uw += val * (xx + 0.5); vw += val * (yy + 0.5);
          }
        if (w <= 0) { --n; continue; }

        s2_3x3 += (L.u - uw / w) * (L.u - uw / w) +
                  (L.v - vw / w) * (L.v - vw / w);
        s2_cog += bd_c;
        s2_fit += bd_f;
      }
    }
    if (!n) continue;
    const double smear = (rig.cfg.airspeed / rig.cfg.radius_m) /
                         rig.cam.pixelIfov() * exp_ms / 1000.0;
    char ebuf[16], sbuf[16];
    std::snprintf(ebuf, sizeof ebuf, "%.0f ms", exp_ms);
    std::snprintf(sbuf, sizeof sbuf, "%.1f px", smear);
    std::printf("  %-10s %-8s %6.0f%%    %-11.3f %-11.3f %-11.3f\n", ebuf, sbuf,
                100.0 * n / n_labels, std::sqrt(s2_3x3 / n),
                std::sqrt(s2_cog / n), std::sqrt(s2_fit / n));
    check(std::sqrt(s2_fit / n) < std::sqrt(s2_cog / n),
          "Gaussian fit beats whole-component CoG");
    check(std::sqrt(s2_cog / n) < std::sqrt(s2_3x3 / n),
          "whole-component CoG beats the 3x3 window");
  }
  std::printf("\n  Photon count is pinned at 50 ms across all rows, so only the\n"
              "  smear varies and this isolates the CENTROIDER rather than the\n"
              "  exposure trade.\n");
}

void testDepth() {
  std::printf("\n=== TEST 5: detection depth vs threshold ===\n");
  const Rig rig = makeRig(5.5);
  std::printf("  catalogue to V=5.5 (%zu stars); how many does the detector\n"
              "  actually recover at 50 ms with realistic noise?\n",
              rig.catalog.size());
  SensorModel s;
  s.exposure_s = 0.05;
  s.blur_substeps = 16;
  MatcherConfig m;
  header();
  for (double k : {3.0, 5.0, 8.0}) {
    DetectorConfig d;
    d.threshold_k = k;
    char buf[64];
    std::snprintf(buf, sizeof buf, "threshold mad + %.0f sigma", k);
    report(buf, runPipeline(rig, s, d, m));
  }
  std::printf("\n  A lower threshold finds more stars but admits noise peaks;\n"
              "  the matcher's ambiguity guard is what keeps those out of the\n"
              "  fix rather than the threshold itself.\n");
}

void testRealistic() {
  std::printf("\n=== TEST 6: everything on ===\n");
  const Rig rig = makeRig(5.0);
  SensorModel s;  // defaults: 50 ms, read 6e, dark, sky, PRNU, hot pixels
  DetectorConfig d;
  MatcherConfig m;
  header();
  const RunStats st = runPipeline(rig, s, d, m);
  report("full sensor model, V<=5.0", st);
  check(st.err_km >= 0 && st.err_km < 4.0,
        "full imaging chain still lands inside the paper's 4 km");
  std::printf("\n  Compare against Milestone 1's ideal-vector result on the\n"
              "  same trajectory: the difference is the cost of the imaging\n"
              "  chain, with everything else held fixed.\n");
}

void testCramerRao() {
  std::printf("\n=== TEST 7: centroid accuracy vs the Cramer-Rao bound ===\n");
  std::printf("  Per magnitude bin, with DETECTION COMPLETENESS exposed.\n");
  const Rig rig = makeRig(4.5);
  const StarField field(rig.catalog, rig.truth[rig.truth.size() / 2].epoch);
  const double period = 1.0 / rig.cfg.frame_rate;

  SensorModel s = SensorModel::ideal();
  s.read_noise_e = 6.0;
  s.shot_noise = true;
  DetectorConfig d;
  d.centroid = CentroidMethod::GaussianFit;

  const double lo[] = {1.0, 2.0, 2.5, 3.0, 3.5, 4.0};
  const double hi[] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5};
  const int NB = 6;
  std::vector<std::vector<double>> err(NB);
  std::vector<int> nlab(NB, 0), nfound(NB, 0);
  std::vector<double> flux(NB, 0.0);

  for (size_t k = 0; k + 1 < rig.truth.size(); k += 8) {
    const RenderedFrame rf =
        renderFrame(rig.truth[k], rig.truth[k + 1], period, field, rig.vmags,
                    rig.cam, rig.C_b_c_true, s, static_cast<unsigned>(k));
    const auto dets = detectStars(rf.image, d);
    for (const StarLabel& L : rf.labels) {
      if (L.truncated) continue;
      if (L.u < 20 || L.v < 20 || L.u > rig.cam.width - 20 ||
          L.v > rig.cam.height - 20) continue;
      int b = -1;
      for (int i = 0; i < NB; ++i)
        if (L.vmag >= lo[i] && L.vmag < hi[i]) b = i;
      if (b < 0) continue;
      ++nlab[b];
      flux[b] += L.total_e;
      double bd = 1e18;
      for (const Detection& dd : dets) {
        const double r2 =
            (L.u - dd.u) * (L.u - dd.u) + (L.v - dd.v) * (L.v - dd.v);
        if (r2 < bd) bd = r2;
      }
      if (bd < 100.0) { ++nfound[b]; err[b].push_back(std::sqrt(bd)); }
    }
  }

  std::printf("\n  %-10s %6s %7s %9s %10s %10s %9s\n", "V bin", "Nlab",
              "det%", "mean e-", "RMS all", "RMS<2px", "CRB(2par)");
  std::printf("  %s\n", std::string(72, '-').c_str());
  double lim_mag = -1.0;
  for (int i = 0; i < NB; ++i) {
    if (!nlab[i]) continue;
    const double N = flux[i] / nlab[i];
    const double sg = 1.2, rn = 6.0;
    const double crb =
        std::sqrt(sg * sg / N + 8 * M_PI * std::pow(sg, 4) * rn * rn / (N * N));
    double a = 0.0, g = 0.0;
    int ng = 0;
    for (double e : err[i]) {
      a += e * e;
      if (e < 2.0) { g += e * e; ++ng; }
    }
    const double rms_all = err[i].empty() ? -1 : std::sqrt(a / err[i].size());
    const double rms_g = ng ? std::sqrt(g / ng) : -1;
    const double frac = double(nfound[i]) / nlab[i];
    if (lim_mag < 0 && frac < 0.5) lim_mag = lo[i];
    char bin[16];
    std::snprintf(bin, sizeof bin, "%.1f-%.1f", lo[i], hi[i]);
    std::printf("  %-10s %6d %7.1f %9.0f %10.3f %10.3f %9.3f\n", bin, nlab[i],
                100.0 * frac, N, rms_all, rms_g, crb);
    if (frac > 0.95) {
      check(rms_all > crb,
            std::string("bin ") + bin + " does not beat the Cramer-Rao bound");
    }
  }

  std::printf("\n  50%% detection limit is near V = %.1f at %.0f ms, NOT the\n"
              "  V=4.0 catalogue depth. That gap is why the ladders above show\n"
              "  ~14 detections rather than ~25.\n", lim_mag, 50.0);
  std::printf("\n  On COMPLETE bins (det%% ~ 100) the fit sits about 1.5x above\n"
              "  the bound, which is what a 5-parameter fit (A,x,y,sigma,B) over\n"
              "  a finite window should cost against a 2-parameter idealised\n"
              "  bound with sigma and background known.\n");
  std::printf("\n  On INCOMPLETE bins the measured RMS appears to beat the\n"
              "  bound. It does not: only the luckiest realisations of those\n"
              "  stars were detected at all. Compare RMS all against RMS<2px in\n"
              "  the 3.5-4.0 row to see how much a tight association gate\n"
              "  flatters the number on top of that.\n");
}


/// MESH BACKGROUND: which contaminant it actually fixes.
///
/// The global median+MAD threshold assumes a flat background. A gradient
/// inflates the global MAD, the threshold rises across the WHOLE frame, and
/// faint stars are lost even where the sky is clean. A local background removes
/// that self-inflicted loss.
///
/// The interesting part is what it does NOT fix. Cloud attenuates starlight and
/// its glow raises the shot-noise floor; neither is recoverable by subtracting
/// anything. Measured separately: flare costs 50% of matched stars, cloud glow
/// 18%, cloud attenuation 14%. The mesh recovers the flare and none of the
/// cloud, which is the whole reason a learned background estimator is NOT
/// obviously worth building -- see scripts/README.md.
void testMeshBackground() {
  std::printf("\n=== TEST: mesh background under flare and cloud ===\n");

  OrbitConfig cfg;
  cfg.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  cfg.airspeed = 25.0;
  cfg.frame_rate = 10.0;
  cfg.radius_m = 150.0;
  cfg.revolutions = 0.2;
  cfg.gps_guided_track = false;
  const auto truth = generateOrbit(cfg);

  const auto cat = catalogBrighterThan(5.5);
  StarField field(cat, truth[truth.size() / 2].epoch);
  std::vector<double> vmags;
  for (size_t i = 0; i < field.size(); ++i) vmags.push_back(field.vmag(i));
  Camera cam;
  const Eigen::Matrix3d C_b_c = nominalCameraMount();

  struct Case { const char* name; double cloud, flare; };
  const Case cases[] = {{"clear", 0.0, 0.0},
                        {"cloud 0.6", 0.6, 0.0},
                        {"flare 0.5", 0.0, 0.5}};

  std::printf("  %-14s %9s %9s %8s\n", "condition", "plain", "mesh", "gain");
  double flare_gain = 0, clear_gain = 0;
  for (const Case& c : cases) {
    double a = 0, b = 0;
    int nf = 0;
    for (size_t i = 3; i + 1 < truth.size() && nf < 4; i += 9, ++nf) {
      SensorModel sensor;
      sensor.exposure_s = 0.10;
      sensor.seed = 100u + unsigned(i);
      sensor.cloud_amount = c.cloud;
      sensor.flare_amount = c.flare;
      const RenderedFrame rf =
          renderFrame(truth[i], truth[i + 1], 0.1, field, vmags, cam, C_b_c,
                      sensor, unsigned(i + 1));
      const Eigen::Matrix3d Ca =
          eulerToDcm(truth[i].roll, truth[i].pitch, truth[i].yaw);
      DetectorConfig d;
      DetectorConfig dm = d;
      dm.bg_mesh_px = 64;

      MatcherConfig mc;
      FrameData f1, f2;
      a += matchDetections(detectStars(rf.image, d), truth[i].epoch, Ca, C_b_c,
                           cam, field, truth[0].pos, mc, f1);
      b += matchDetections(detectStars(rf.image, dm), truth[i].epoch, Ca, C_b_c,
                           cam, field, truth[0].pos, mc, f2);
    }
    const double g = a > 0 ? b / a : 0.0;
    std::printf("  %-14s %9.1f %9.1f %7.2fx\n", c.name, a / nf, b / nf, g);
    if (!std::strcmp(c.name, "flare 0.5")) flare_gain = g;
    if (!std::strcmp(c.name, "clear")) clear_gain = g;
  }

  // 1.15, not 1.3. The 1.76x this once asserted was measured with the matched
  // filter enabled in this test; when the filter was removed from the library
  // the same measurement gave 1.20x. The mesh still recovers most of the flare
  // loss, but a plain detector has less to recover.
  check(flare_gain > 1.15, "mesh background recovers the flare loss");
  check(clear_gain > 0.85,
        "and costs little in clean conditions (guard: it is not free, so it is "
        "off by default)");
}

/// THE GUARANTEE: a prior can only move the operating point.
///
/// A segmentation network can hallucinate a star and can erase a real one, and
/// neither may be guarded against by trusting the training. So the prior does
/// not decide anything. It lowers the DETECTION THRESHOLD between
/// `threshold_k` and `threshold_k_low`, and the detection still has to clear a
/// real significance floor computed from actual photons.
///
/// The precise, checkable statement is an EQUIVALENCE:
///
///   prior = 1 everywhere  ==  running the ordinary detector at threshold_k_low
///   prior = 0 everywhere  ==  running the ordinary detector at threshold_k
///
/// So the worst a network can do -- however wrong, however adversarial -- is
/// move you to an operating point you could have selected by hand. It cannot
/// reach outside that range, cannot manufacture a detection where there is no
/// flux, and cannot delete one, because a zero prior leaves the nominal
/// threshold untouched.
///
/// That is a bound on the DAMAGE, not a claim the model is good, and it holds
/// without any assumption about training.
void testPriorIsBounded() {
  std::printf("\n=== TEST: detection prior is bounded by construction ===\n");

  OrbitConfig cfg;
  cfg.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  cfg.airspeed = 25.0;
  cfg.frame_rate = 10.0;
  cfg.radius_m = 150.0;
  cfg.revolutions = 0.15;
  cfg.gps_guided_track = false;
  const auto truth = generateOrbit(cfg);

  const auto cat = catalogBrighterThan(5.5);
  StarField field(cat, truth[truth.size() / 2].epoch);
  std::vector<double> vmags;
  for (size_t i = 0; i < field.size(); ++i) vmags.push_back(field.vmag(i));
  Camera cam;
  const Eigen::Matrix3d C_b_c = nominalCameraMount();

  SensorModel sensor;
  sensor.exposure_s = 0.05;
  sensor.seed = 4242;
  const RenderedFrame rf = renderFrame(truth[2], truth[3], 0.1, field, vmags,
                                       cam, C_b_c, sensor, 7u);
  const size_t npix = size_t(rf.image.width) * rf.image.height;
  const Eigen::Matrix3d Ca =
      eulerToDcm(truth[2].roll, truth[2].pitch, truth[2].yaw);
  MatcherConfig mc;

  auto run = [&](const DetectorConfig& d, const char* name) {
    const auto dets = detectStars(rf.image, d);
    int real = 0;
    for (const auto& det : dets) {
      double best = 1e18;
      for (const auto& L : rf.labels) {
        if (L.truncated) continue;
        best = std::min(best, std::hypot(det.u - L.u, det.v - L.v));
      }
      if (best < 5.0) ++real;
    }
    FrameData f;
    const int m = matchDetections(dets, truth[2].epoch, Ca, C_b_c, cam, field,
                                  truth[0].pos, mc, f);
    std::printf("  %-28s %4zu det, %4d real, %4d matched\n", name, dets.size(),
                real, m);
    return dets.size();
  };

  DetectorConfig nominal;
  DetectorConfig low;
  low.threshold_k = nominal.threshold_k_low;   // hand-picked low operating point
  DetectorConfig p1 = nominal, p0 = nominal;
  p1.prior.assign(npix, 1.0f);
  p0.prior.assign(npix, 0.0f);

  const size_t n_nom = run(nominal, "no prior (k_high)");
  const size_t n_low = run(low, "hand-set k_low, no prior");
  const size_t n_p1 = run(p1, "prior = 1 EVERYWHERE");
  const size_t n_p0 = run(p0, "prior = 0 EVERYWHERE");

  check(n_p0 == n_nom, "a prior of zero cannot delete anything");
  check(n_p1 == n_low,
        "a prior of one is EXACTLY the hand-set low threshold -- an adversarial "
        "network cannot reach outside the operating range");
  check(n_low >= n_nom, "the low threshold is the more permissive end");
}

/// PriorNet degrades correctly when there is no model.
///
/// The learned prior is OPTIONAL in every sense: optional dependency, optional
/// model file, optional at runtime. Every failure path must leave the caller
/// with the ordinary detector rather than an error, because a navigation
/// sensor that stops working when a file is missing is worse than one that
/// never had the file.
void testPriorNetFallback() {
  std::printf("\n=== TEST: PriorNet fallback ===\n");
  std::printf("  built with OpenCV: %s\n",
              PriorNet::available() ? "yes" : "no (stub)");

  OrbitConfig cfg;
  cfg.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  cfg.airspeed = 25.0;
  cfg.frame_rate = 10.0;
  cfg.radius_m = 150.0;
  cfg.revolutions = 0.05;
  cfg.gps_guided_track = false;
  const auto truth = generateOrbit(cfg);
  const auto cat = catalogBrighterThan(5.5);
  StarField field(cat, truth[0].epoch);
  std::vector<double> vmags;
  for (size_t i = 0; i < field.size(); ++i) vmags.push_back(field.vmag(i));
  Camera cam;
  SensorModel sensor;
  const RenderedFrame rf = renderFrame(truth[0], truth[1], 0.1, field, vmags,
                                       cam, nominalCameraMount(), sensor, 1u);

  PriorNet pn;
  check(!pn.loaded(), "starts unloaded");
  check(!pn.load("/nonexistent/model.onnx"), "missing model reports failure");

  std::vector<float> prior;
  check(!pn.infer(rf.image, prior), "inference without a model reports failure");
  check(prior.empty(), "and leaves the prior empty");

  // An empty prior must be a no-op, so a missing model costs nothing.
  DetectorConfig plain, with_empty;
  with_empty.prior = prior;
  check(detectStars(rf.image, plain).size() ==
            detectStars(rf.image, with_empty).size(),
        "an empty prior is exactly the ordinary detector");
}

}  // namespace

int main(int argc, char** argv) {
  const int only = (argc > 1) ? std::atoi(argv[1]) : 0;
  std::printf("Milestone 2 — image simulator, first version\n");
  std::printf("===========================================\n");

  if (!only || only == 1) testIdeal();
  if (!only || only == 2) testNoiseLadder();
  if (!only || only == 3) testBlurLadder();
  if (!only || only == 4) testCentroidBiasUnderBlur();
  if (!only || only == 5) testDepth();
  if (!only || only == 6) testRealistic();
  if (!only || only == 7) testCramerRao();
  if (!only || only == 9) testMeshBackground();
  if (!only || only == 10) testPriorIsBounded();
  if (!only || only == 11) testPriorNetFallback();

  std::printf("\n===========================================\n");
  if (g_failures == 0) {
    std::printf("ALL CHECKS PASSED\n");
    return 0;
  }
  std::printf("%d CHECK(S) FAILED\n", g_failures);
  return 1;
}
