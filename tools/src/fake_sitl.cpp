// fake_sitl.cpp — replay a generated orbit as MAVLink, to exercise
// celestial_node without ArduPilot running.
//
//   ./build/fake_sitl --port 14556 --speed 20
//   ./build/celestial_node --port 14556 --no-inject
//
// Emits HEARTBEAT, SIMSTATE (true pose) and ATTITUDE (true pose plus a
// body-fixed bias and noise, standing in for EKF3 error). Not a substitute for
// SITL -- it is a loopback test that the node's MAVLink plumbing, timing and
// orbit detection all work before you plug in the real autopilot.
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <ardupilotmega/mavlink.h>
#include "celestial/orbit.hpp"
#include "celestial/sky_model.hpp"

using namespace celestial;

int main(int argc, char** argv) {
  int port = 14556; double speed = 20.0; double bias_deg = 0.05;
  std::string pattern = "orbit";
  double revs = 3.0;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto nx = [&]{ return (i + 1 < argc) ? argv[++i] : ""; };
    if (a == "--port") port = std::atoi(nx());
    else if (a == "--speed") speed = std::atof(nx());
    else if (a == "--bias") bias_deg = std::atof(nx());
    else if (a == "--pattern") pattern = nx();   // orbit | straight | transit
    else if (a == "--revs") revs = std::atof(nx());
  }

  OrbitConfig cfg;
  cfg.centre = Geodetic::fromDegrees(-30.89, 136.56, 800.0);
  cfg.airspeed = 25; cfg.frame_rate = 20;
  std::vector<FrameTruth> truth;
  if (pattern == "straight") {
    // A 200 km turn radius over a small arc is straight and level to within
    // a degree of heading: the orbit mechanism has nothing to work with.
    cfg.radius_m = 200000; cfg.revolutions = 0.004;
    truth = generateOrbit(cfg);
  } else if (pattern == "transit") {
    // The mission this system is actually for: long straight legs with an
    // occasional loiter. Dead reckoning drifts on the legs; each loiter
    // produces a celestial fix that pulls it back. A short flight cannot show
    // that -- the drift has to exceed the fix error before the fix is worth
    // anything, which at ~8 km/hour and a 2 km fix takes ~15 minutes.
    OrbitConfig leg = cfg, loi = cfg;
    leg.radius_m = 200000; leg.revolutions = 0.02;   // ~17 min of straight
    loi.radius_m = 600;    loi.revolutions = 1.2;    // one calibrating orbit
    for (int k = 0; k < int(revs); ++k) {
      OrbitConfig a2 = (k % 2) ? loi : leg;
      a2.start_track = 1.1 * k;
      auto seg = generateOrbit(a2);
      const double t0 = truth.empty() ? 0.0 : truth.back().t;
      Geodetic anchor = truth.empty() ? seg.front().pos : truth.back().pos;
      const double dlat = anchor.lat - seg.front().pos.lat;
      const double dlon = anchor.lon - seg.front().pos.lon;
      for (auto f : seg) {
        f.t += t0;
        f.pos.lat += dlat;
        f.pos.lon += dlon;
        truth.push_back(f);
      }
    }
  } else {
    cfg.radius_m = 600; cfg.revolutions = revs;
    truth = generateOrbit(cfg);
  }
  if (truth.empty()) { std::printf("no trajectory\n"); return 1; }
  double gt = 0;
  for (size_t k = 1; k < truth.size(); ++k) gt += haversine(truth[k-1].pos, truth[k].pos);
  std::printf("ground track %.1f km\n", gt / 1000.0);
  const double hdg = std::abs(truth.back().yaw - truth.front().yaw) * kRad2Deg;
  std::printf("pattern '%s': %zu poses, %.0f s, %.1f deg of heading change,"
              " replayed at %.0fx into UDP %d\n",
              pattern.c_str(), truth.size(), truth.back().t, hdg, speed, port);

  int fd = socket(AF_INET, SOCK_DGRAM, 0);
  sockaddr_in me{}; me.sin_family = AF_INET;
  me.sin_addr.s_addr = htonl(INADDR_ANY); me.sin_port = htons(port + 1);
  bind(fd, (sockaddr*)&me, sizeof me);
  sockaddr_in to{}; to.sin_family = AF_INET;
  to.sin_addr.s_addr = inet_addr("127.0.0.1"); to.sin_port = htons(port);

  auto send = [&](const mavlink_message_t& m) {
    uint8_t b[MAVLINK_MAX_PACKET_LEN];
    sendto(fd, b, mavlink_msg_to_send_buffer(b, &m), 0, (sockaddr*)&to, sizeof to);
  };

  std::mt19937 rng(7);
  std::normal_distribution<double> n(0.0, 0.02 * kDeg2Rad);
  const double bias = bias_deg * kDeg2Rad;

  for (size_t k = 0; k < truth.size(); ++k) {
    const FrameTruth& f = truth[k];
    mavlink_message_t m;
    if (k % 20 == 0) {
      mavlink_msg_heartbeat_pack(1, 1, &m, MAV_TYPE_FIXED_WING,
                                 MAV_AUTOPILOT_ARDUPILOTMEGA, 0, 0,
                                 MAV_STATE_ACTIVE);
      send(m);
    }
    mavlink_msg_simstate_pack(1, 1, &m, f.roll, f.pitch, f.yaw, 0, 0, -9.81,
                              0, 0, 0, int32_t(f.pos.lat * kRad2Deg * 1e7),
                              int32_t(f.pos.lon * kRad2Deg * 1e7));
    send(m);
    // "EKF3": truth plus a body-fixed pitch bias plus noise.
    //
    // Body rates MUST be populated: the node uses them to bridge the gap
    // between ATTITUDE and SIMSTATE, and to choose the exposure. Sending zeros
    // makes the loopback test silently skip both paths -- it reported 0.0 px
    // of smear and pinned the exposure at its ceiling.
    double p_rate = 0, q_rate = 0, r_rate = 0;
    if (k + 1 < truth.size()) {
      const FrameTruth& g = truth[k + 1];
      const double dt = std::max(1e-3, g.t - f.t);
      auto wrap = [](double d) { return std::atan2(std::sin(d), std::cos(d)); };
      // Euler rates -> body rates, small-angle: adequate at these attitudes.
      const double dphi = wrap(g.roll - f.roll) / dt;
      const double dth = wrap(g.pitch - f.pitch) / dt;
      const double dpsi = wrap(g.yaw - f.yaw) / dt;
      p_rate = dphi - dpsi * std::sin(f.pitch);
      q_rate = dth * std::cos(f.roll) + dpsi * std::sin(f.roll) * std::cos(f.pitch);
      r_rate = -dth * std::sin(f.roll) + dpsi * std::cos(f.roll) * std::cos(f.pitch);
    }
    mavlink_msg_attitude_pack(1, 1, &m, uint32_t(f.t * 1000),
                              f.roll + n(rng), f.pitch + bias + n(rng),
                              f.yaw + n(rng), float(p_rate), float(q_rate),
                              float(r_rate));
    send(m);
    // VFR_HUD carries airspeed, which the dead-reckoning filter needs. A stub
    // that omits a message the real source sends does not fail -- it silently
    // disables whatever consumes it. Same trap as the zero body rates.
    mavlink_msg_vfr_hud_pack(1, 1, &m, float(cfg.airspeed),
                             float(cfg.airspeed), int16_t(f.yaw * kRad2Deg),
                             50, float(f.pos.alt), 0.0f);
    send(m);
    mavlink_msg_global_position_int_pack(1, 1, &m, uint32_t(f.t * 1000),
                                         int32_t(f.pos.lat * kRad2Deg * 1e7),
                                         int32_t(f.pos.lon * kRad2Deg * 1e7),
                                         int32_t(f.pos.alt * 1000),
                                         int32_t(f.pos.alt * 1000), 0, 0, 0, 0);
    send(m);
    usleep(useconds_t(1e6 / cfg.frame_rate / speed));
  }
  std::printf("replay finished\n");
  return 0;
}
