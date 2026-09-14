// celestial_node.cpp — the live node: argument parsing, the loop, and output.
//
// The work is in two places, split along the system schematic:
//
//   mavlink_source.{hpp,cpp}   TRANSPORT -- getting a consistent (true pose,
//                              AHRS attitude) pair out of a MAVLink stream
//   pipeline.{hpp,cpp}         NAVIGATION -- render, detect, fuse, match, fix,
//                              calibrate, dead reckon
//
// This file owns neither. It wires them together and prints.

#include <cstdio>
#include <cstring>
#include <memory>
#include <string>

#include "mavlink_source.hpp"
#include "pipeline.hpp"

using namespace celestial;

int main(int argc, char** argv) {
  MavlinkSource::Config mc;
  PipelineConfig pc;
  // NO HORIZON UNLESS ASKED. HorizonConfig default-constructs with ONE camera
  // (a Boson 640), not an empty vector, so without this every run silently had
  // a horizon sensor -- and a research-grade one. Consequences when this was
  // found: a SITL "no horizon vs horizon" comparison was really Boson 640 vs
  // Lepton 2.5, --horizon was a DOWNGRADE (0.06 -> 5.13 arcmin tilt floor, so
  // the per-frame fix went 10.0 -> 16.5 km), and pipeline.cpp's
  // `cameras.empty() ? 6000 : 2700` fix sigma never took the 6000 branch.
  // pipeline.hpp says "Empty `cameras` disables it, which is the default" --
  // that comment describes the intent, this line implements it.
  pc.horizon.cameras.clear();
  std::string csv_path;
  bool inject = true;

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto nx = [&] { return (i + 1 < argc) ? argv[++i] : ""; };
    if (a == "--port") mc.port = std::atoi(nx());
    else if (a == "--max-skew") mc.max_skew_s = std::atof(nx());
    else if (a == "--rate") pc.render_hz = std::atof(nx());
    else if (a == "--mag") pc.mag_limit = std::atof(nx());
    else if (a == "--match-mag") pc.match_mag = std::atof(nx());
    else if (a == "--exposure") pc.exposure_s = std::atof(nx());
    else if (a == "--target-smear") pc.target_smear_px = std::atof(nx());
    else if (a == "--match-radius") pc.match_radius_px = std::atof(nx());
    else if (a == "--min-stars") pc.min_stars = std::atoi(nx());
    else if (a == "--boresight") pc.boresight_deg = std::atof(nx());
    else if (a == "--fix-sigma") pc.fix_sigma_full = std::atof(nx());
    else if (a == "--interval") pc.interval_s = std::atof(nx());
    else if (a == "--dr-quality") pc.dr_quality = nx();
    else if (a == "--no-dr") pc.dead_reckon = false;
    else if (a == "--utc") pc.utc = nx();
    else if (a == "--horizon") {
      // Single forward-facing Lepton 2.5: measured 2.5x on the orbit fix.
      pc.horizon.cameras.assign(1, HorizonCamera::lepton25());
    } else if (a == "--horizon-pair") {
      pc.horizon.cameras = {HorizonCamera::lepton25(false),
                            HorizonCamera::lepton25(true)};
    } else if (a == "--frames-dir") pc.frames_dir = nx();
    else if (a == "--frames-every") pc.frames_every = std::atoi(nx());
    else if (a == "--csv") csv_path = nx();
    else if (a == "--no-inject") inject = false;
  }

  MavlinkSource src(mc);
  if (!src.open()) { std::printf("failed to bind UDP %d\n", mc.port); return 1; }
  Pipeline pipe(pc);



  std::printf("listening on UDP %d, render %.1f Hz, V<=%.1f\n", mc.port,
              pc.render_hz, pc.mag_limit);
  std::printf("dead reckoning %s (%s)\n",
              pc.dead_reckon ? "on" : "off", pc.dr_quality.c_str());

  FILE* csv = csv_path.empty() ? nullptr : std::fopen(csv_path.c_str(), "w");
  if (csv) {
    std::fprintf(csv,
                 "t,lat_true,lon_true,lat_fix,lon_fix,frame_err_m,"
                 "tilt_imu_fwd_am,tilt_imu_lat_am,stars,"
                 "yaw_deg,orbit_lat,orbit_lon,orbit_err_m,det,exposure_ms,"
                 "boresight_deg,calibrated,dr_lat,dr_lon,dr_err_m,"
                 "dr_sigma_m,dr_since_fix_m,dr_open_err_m,dr_open_lat,dr_open_lon,"
                 "\n");
  }

  MavSample s;
  FixReport pend;
  bool have_pend = false;

  while (true) {
    if (!src.poll(s)) {
      const auto r = src.rates();
      if (r.warn) {
        std::printf("\n  WARNING: SIMSTATE only %.1f Hz (ATTITUDE %.1f Hz).\n"
                    "  SIMSTATE is in no stream for Plane, so no parameter"
                    " helps -- SET_MESSAGE_INTERVAL\n  is the only route, and"
                    " MAVProxy's own streamrate overrides it. In MAVProxy:\n"
                    "      set streamrate 20\n", r.simstate_hz, r.attitude_hz);
      }
      if (r.skew_rejects)
        std::printf("\n  %d frame(s) rejected on skew > %.0f ms\n",
                    r.skew_rejects, r.skew_limit_s * 1000);
      continue;
    }

    FixReport fix;
    const FrameReport fr = pipe.step(s, fix);
    if (!fr.rendered) continue;

    std::printf("\rt=%7.1f s  det %2d  matched %2d  exp %3.0f ms  smear %4.1f px"
                "  hdg %+5.0f  frames %3zu   ",
                s.t, fr.detected, fr.matched, fr.exposure_s * 1000, fr.smear_px,
                fr.heading_deg, fr.window);
    std::fflush(stdout);

    if (fix.rejected) {
      std::printf("\n[skip] %.0f deg of heading, %zu frames, %d/12 bins: %s\n",
                  fix.heading_deg, fix.frames, fix.bins, fix.note.c_str());
    }
    if (fix.emitted) {
      std::printf("\n[fix %d] %.0f deg of heading, %zu frames, %d/12 bins%s\n",
                  fix.index, fix.heading_deg, fix.frames, fix.bins,
                  fix.dense ? "" : "  <- SPARSE");
      std::printf("          %-16s %.2f km\n",
                  pipe.calibrated() ? "heading-weighted" : "circle fit",
                  fix.err_m / 1000);
      std::printf("          %-16s %.2f km   (paper Eq. 23)\n", "naive mean",
                  fix.err_naive_m / 1000);
      std::printf("          boresight %.4f deg from truth\n", fix.boresight_deg);
      if (!fix.note.empty()) std::printf("          %s\n", fix.note.c_str());
      if (fix.dr_err_m >= 0)
        std::printf("          DR %.2f km, DR alone %.2f km (%.2f km open-loop)\n",
                    fix.dr_err_m / 1000, fix.dr_open_err_m / 1000,
                    fix.dr_drift_m / 1000);
      if (inject) src.sendFix(fix.position, 4000.0, s.t);
      pend = fix;
      have_pend = true;
    }

    if (csv) {
      std::fprintf(csv,
                   "%.3f,%.7f,%.7f,%.7f,%.7f,%.1f,%.3f,%.3f,"
                   "%d,%.2f,%.7f,%.7f,%.1f,%d,%.1f,%.5f,%d,"
                   "%.7f,%.7f,%.1f,0,%.1f,%.1f,%.7f,%.7f\n",
                   s.t, s.true_pos.lat * kRad2Deg, s.true_pos.lon * kRad2Deg,
                   fr.frame_err_m >= 0 ? fr.frame_fix.lat * kRad2Deg : 0.0,
                   fr.frame_err_m >= 0 ? fr.frame_fix.lon * kRad2Deg : 0.0,
                   fr.frame_err_m, fr.tilt_imu.x(), fr.tilt_imu.y(),
                   fr.matched,
                   s.yaw_est * kRad2Deg,
                   have_pend ? pend.position.lat * kRad2Deg : 0.0,
                   have_pend ? pend.position.lon * kRad2Deg : 0.0,
                   have_pend ? pend.err_m : -1.0, fr.detected,
                   fr.exposure_s * 1000,
                   have_pend ? pend.boresight_deg : pc.boresight_deg,
                   pipe.calibrated() ? 1 : 0,
                   fr.dr_valid ? fr.dr_pos.lat * kRad2Deg : 0.0,
                   fr.dr_valid ? fr.dr_pos.lon * kRad2Deg : 0.0,
                   have_pend ? pend.dr_err_m : -1.0,
                   have_pend ? pend.dr_drift_m : 0.0,
                   have_pend ? pend.dr_open_err_m : -1.0,
                   fr.dr_open_valid ? fr.dr_open_pos.lat * kRad2Deg : 0.0,
                   fr.dr_open_valid ? fr.dr_open_pos.lon * kRad2Deg : 0.0);
      std::fflush(csv);
      have_pend = false;
    }
  }
  return 0;
}
