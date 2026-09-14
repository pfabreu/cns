// visualise.cpp — figures, in two modes.
//
//   ./build/visualise frames  out/ [n]      star frames as PGM
//   ./build/visualise match   out/ [bore]   star identification as SVG
//
// One tool because they answer the same question -- what does the sensor
// actually see -- and neither is big enough to deserve its own binary.
//
// `match` is the more useful of the two. It renders detections, catalogue
// predictions and accepted matches as a self-contained SVG, so the systematic
// offset caused by an uncalibrated mounting is visible in a single frame: at
// 0.4 deg it is ~14 px at 107 arcsec/px, coherent across every star. That is
// precisely the quantity the orbit maneuver spends a revolution measuring.

#include "celestial/imaging.hpp"
#include "celestial/orbit.hpp"
#include "celestial/sky_model.hpp"
#include "celestial/star_catalog.hpp"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
using namespace celestial;

static void writePgm(const std::string& path, const Image& im, int shift) {
  FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) return;
  std::fprintf(f, "P5\n%d %d\n255\n", im.width, im.height);
  std::vector<unsigned char> row(im.width);
  for (int y = 0; y < im.height; ++y) {
    for (int x = 0; x < im.width; ++x)
      row[x] = (unsigned char)std::min(255, im.data[size_t(y)*im.width+x] >> shift);
    std::fwrite(row.data(), 1, row.size(), f);
  }
  std::fclose(f);
}



static int runFrames(const std::string& dir, int nframes);
static int runMatch(const std::string& dir, double bore);

int main(int argc, char** argv) {
  const std::string mode = (argc > 1) ? argv[1] : "match";
  const std::string dir = (argc > 2) ? argv[2] : ".";
  if (mode == "frames") return runFrames(dir, (argc > 3) ? std::atoi(argv[3]) : 60);
  if (mode == "match") return runMatch(dir, (argc > 3) ? std::atof(argv[3]) : 0.4);
  std::printf("usage: visualise frames|match <dir> [n|boresight_deg]\n");
  return 1;
}

static int runFrames(const std::string& dir, int nframes) {

  OrbitConfig cfg;
  cfg.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  cfg.radius_m = 600; cfg.airspeed = 25; cfg.frame_rate = 10; cfg.revolutions = 1;
  const auto truth = generateOrbit(cfg);
  const auto cat = catalogBrighterThan(4.5);
  std::vector<double> vm; for (const auto& s : cat) vm.push_back(s.vmag);
  const StarField field(cat, truth[truth.size()/2].epoch);

  Camera cam;
  SensorModel sensor; sensor.exposure_s = 0.05; sensor.blur_substeps = 12;
  const Eigen::Matrix3d mount =
      nominalCameraMount() * smallRotation(Eigen::Vector3d(1,0,0), 0.4*kDeg2Rad);

  const int stride = int(truth.size()) / nframes;
  for (int i = 0; i < nframes; ++i) {
    const size_t k = size_t(i) * stride;
    if (k + 1 >= truth.size()) break;
    char buf[256];

    const RenderedFrame rf = renderFrame(truth[k], truth[k+1], stride/10.0,
                                         field, vm, cam, mount, sensor, unsigned(k));
    std::snprintf(buf, sizeof buf, "%s/star_%03d.pgm", dir.c_str(), i);
    writePgm(buf, rf.image, 0);   // star field is faint: no shift

    std::printf("\rframe %d/%d ", i+1, nframes); std::fflush(stdout);
  }
  std::printf("\ndone: %d frames in %s\n", nframes, dir.c_str());
  return 0;
}

static int runMatch(const std::string& dir, double bore) {

  OrbitConfig cfg;
  cfg.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  cfg.radius_m = 600; cfg.airspeed = 25; cfg.frame_rate = 10; cfg.revolutions = 1;
  const auto truth = generateOrbit(cfg);
  const auto cat = catalogBrighterThan(4.5);
  std::vector<double> vm; for (const auto& s : cat) vm.push_back(s.vmag);
  const StarField field(cat, truth[truth.size()/2].epoch);

  Camera cam;
  SensorModel sensor; sensor.exposure_s = 0.05; sensor.blur_substeps = 12;
  DetectorConfig dcfg; dcfg.centroid = CentroidMethod::GaussianFit;
  MatcherConfig mcfg;
  const Eigen::Matrix3d mount_true =
      nominalCameraMount() * smallRotation(Eigen::Vector3d(1,0,0), bore*kDeg2Rad);

  // One frame to draw, and a full orbit to recalibrate from.
  const size_t K = truth.size()/3;
  const RenderedFrame rf = renderFrame(truth[K], truth[K+1], 0.1, field, vm,
                                       cam, mount_true, sensor, 7u);
  const auto dets = detectStars(rf.image, dcfg);
  const Eigen::Matrix3d C_l_b =
      eulerToDcm(truth[K].roll, truth[K].pitch, truth[K].yaw);

  ErrorModel em; em.boresight_error = bore*kDeg2Rad;
  em.centroid_noise = 10*kArcsec2Rad;
  const auto frames = simulateObservations(truth, cam, em);

  Eigen::Matrix3d Cbc = nominalCameraMount();
  for (int iter = 0; iter < 2; ++iter) {
    // Predicted pixel positions with the currently assumed mounting.
    std::vector<Eigen::Vector3d> ecef;
    field.ecefDirections(truth[K].epoch, ecef);
    const Eigen::Matrix3d C_en = nedBasisEcef(truth[K].pos);
    const Eigen::Matrix3d cam_from_ned = (C_l_b * Cbc).transpose();

    std::string svg =
      "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 1936 1216' "
      "width='100%'><rect width='1936' height='1216' fill='#080a12'/>\n";
    int nmatch = 0; double resid = 0;
    // One-to-one: without this, two predictions can claim the same detection
    // and the match count exceeds the detection count.
    std::vector<char> used(dets.size(), 0);

    for (size_t i = 0; i < ecef.size(); ++i) {
      const Eigen::Vector3d v = cam_from_ned * (C_en.transpose()*ecef[i]);
      double pu, pv;
      if (!rayToPixel(v, cam, pu, pv)) continue;
      if (pu < -40 || pv < -40 || pu > cam.width+40 || pv > cam.height+40) continue;
      char b[512];
      std::snprintf(b, sizeof b,
        "<circle cx='%.1f' cy='%.1f' r='11' fill='none' stroke='#f0a500' "
        "stroke-width='1.6' opacity='0.85'/>\n", pu, pv);
      svg += b;

      // Label. Named stars get their common name; the rest fall back to the
      // HR number, which is what you want on a plot anyway.
      const char* nm = starName(field.id(i));
      char lbl[64];
      if (nm) std::snprintf(lbl, sizeof lbl, "%s", nm);
      else std::snprintf(lbl, sizeof lbl, "HR %d", field.id(i));
      // Offset the text away from the frame edge so it stays readable.
      const bool left = pu > cam.width * 0.75;
      std::snprintf(b, sizeof b,
        "<text x='%.1f' y='%.1f' fill='#f0a500' font-family='monospace' "
        "font-size='19' opacity='0.9' text-anchor='%s'>%s</text>\n",
        pu + (left ? -16 : 16), pv - 14, left ? "end" : "start", lbl);
      svg += b;
      // nearest detection
      double bd = 1e18; int bi = -1;
      for (size_t k = 0; k < dets.size(); ++k) {
        if (used[k]) continue;
        const double r2 = (dets[k].u-pu)*(dets[k].u-pu) +
                          (dets[k].v-pv)*(dets[k].v-pv);
        if (r2 < bd) { bd = r2; bi = int(k); }
      }
      const Detection* bdet = (bi >= 0) ? &dets[bi] : nullptr;
      if (bdet && bd < mcfg.radius_px*mcfg.radius_px) {
        used[bi] = 1;
        std::snprintf(b, sizeof b,
          "<line x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f' stroke='#3ad39a' "
          "stroke-width='1.4' opacity='0.9'/>\n", pu, pv, bdet->u, bdet->v);
        svg += b;
        ++nmatch; resid += std::sqrt(bd);
      }
    }
    for (const auto& d : dets) {
      const double r = std::clamp(2.0 + std::log10(std::max(1.0,d.flux))*1.6, 2.0, 9.0);
      char b[256];
      std::snprintf(b, sizeof b,
        "<circle cx='%.1f' cy='%.1f' r='%.1f' fill='#dfe6ff'/>\n", d.u, d.v, r);
      svg += b;
    }
    char hdr[512];
    std::snprintf(hdr, sizeof hdr,
      "<text x='24' y='44' fill='#dfe6ff' font-family='monospace' "
      "font-size='30'>iteration %d   %zu detections   %d matched   "
      "mean offset %.1f px</text>\n"
      "<text x='24' y='84' fill='#8d97b8' font-family='monospace' "
      "font-size='22'>white = detected   orange = catalogue prediction   "
      "green = accepted match</text>\n",
      iter, dets.size(), nmatch, nmatch ? resid/nmatch : 0.0);
    svg = svg.substr(0, svg.find('\n')+1) + hdr + svg.substr(svg.find('\n')+1);
    svg += "</svg>\n";

    char path[512];
    std::snprintf(path, sizeof path, "%s/match_iter%d.svg", dir.c_str(), iter);
    FILE* f = std::fopen(path, "w");
    if (f) { std::fwrite(svg.data(), 1, svg.size(), f); std::fclose(f); }
    std::printf("iter %d: %zu detections, %d matched, mean offset %.2f px "
                "(%.3f deg)\n", iter, dets.size(), nmatch,
                nmatch ? resid/nmatch : 0.0,
                nmatch ? (resid/nmatch)*cam.pixelIfov()*kRad2Deg : 0.0);

    // Recalibrate from the whole orbit, as the real pipeline does.
    const auto r = estimateOrbit(frames, Cbc, AverageMethod::HeadingWeighted);
    if (r.ok) Cbc = recalibrateMount(frames, r.position);
  }
  return 0;
}
