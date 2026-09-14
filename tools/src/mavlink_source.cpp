#include "mavlink_source.hpp"

#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

#include <Eigen/Geometry>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

#include <ardupilotmega/mavlink.h>

namespace celestial {
namespace {
// Identify as a GCS, NOT as system 1. A COMMAND_LONG arriving with the
// vehicle's own system id is liable to be dropped by MAVLink routing, which
// silently kills the SET_MESSAGE_INTERVAL request.
constexpr int kSysId = 255;
constexpr int kCompId = MAV_COMP_ID_MISSIONPLANNER;
constexpr int kTargetSys = 1;
constexpr int kTargetComp = 1;

double nowSeconds() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch()).count();
}
}  // namespace

struct MavlinkSource::Impl {
  int fd = -1;
  sockaddr_in peer{};
  bool have_peer = false;
  bool requested = false;
  bool have_sim = false, have_att = false;

  MavSample cur;
  Eigen::Matrix3d C_latest = Eigen::Matrix3d::Identity();
  double yaw_latest = 0.0;
  Eigen::Vector3d omega_latest = Eigen::Vector3d::Zero();
  double t_sim = -1e9, t_att = -1e9;

  int n_sim = 0, n_att = 0, n_reject = 0;
  double t_report = nowSeconds();
  double hz_att = 20.0, skew_limit = 0.08;
  bool warn = false;

  void send(const mavlink_message_t& m) {
    if (!have_peer) return;
    uint8_t buf[MAVLINK_MAX_PACKET_LEN];
    const uint16_t n = mavlink_msg_to_send_buffer(buf, &m);
    sendto(fd, buf, n, 0, (sockaddr*)&peer, sizeof peer);
  }

  void requestAll(double hz) {
    for (auto id : {MAVLINK_MSG_ID_SIMSTATE, MAVLINK_MSG_ID_ATTITUDE,
                    MAVLINK_MSG_ID_VFR_HUD, MAVLINK_MSG_ID_GLOBAL_POSITION_INT}) {
      mavlink_message_t m;
      mavlink_msg_command_long_pack(kSysId, kCompId, &m, kTargetSys, kTargetComp,
                                    MAV_CMD_SET_MESSAGE_INTERVAL, 0, float(id),
                                    float(1e6 / hz), 0, 0, 0, 0, 0);
      send(m);
    }
  }
};

MavlinkSource::~MavlinkSource() {
  if (impl_) { if (impl_->fd >= 0) ::close(impl_->fd); delete impl_; }
}

bool MavlinkSource::open() {
  impl_ = new Impl();
  impl_->fd = socket(AF_INET, SOCK_DGRAM, 0);
  if (impl_->fd < 0) return false;
  sockaddr_in a{};
  a.sin_family = AF_INET;
  a.sin_addr.s_addr = htonl(INADDR_ANY);
  a.sin_port = htons(cfg_.port);
  if (bind(impl_->fd, (sockaddr*)&a, sizeof a) < 0) return false;
  timeval tv{0, 20000};
  setsockopt(impl_->fd, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof tv);
  // Rendering blocks the loop for ~100 ms; the socket must hold a burst.
  int rcvbuf = 1 << 20;
  setsockopt(impl_->fd, SOL_SOCKET, SO_RCVBUF, &rcvbuf, sizeof rcvbuf);
  impl_->skew_limit = (cfg_.max_skew_s > 0) ? cfg_.max_skew_s : 0.08;
  return true;
}

bool MavlinkSource::poll(MavSample& out) {
  Impl& d = *impl_;

  // Drain everything queued, not one datagram: rendering blocks long enough
  // that a backlog builds, and draining one per pass leaves arrival stamps
  // bunched and mutually skewed.
  uint8_t buf[2048];
  for (int pkt = 0; pkt < 256; ++pkt) {
    sockaddr_in from{};
    socklen_t flen = sizeof from;
    const ssize_t n = recvfrom(d.fd, buf, sizeof buf, pkt ? MSG_DONTWAIT : 0,
                               (sockaddr*)&from, &flen);
    if (n <= 0) break;
    d.peer = from;
    d.have_peer = true;
    mavlink_message_t msg;
    mavlink_status_t st;
    for (ssize_t i = 0; i < n; ++i) {
      if (!mavlink_parse_char(MAVLINK_COMM_0, buf[i], &msg, &st)) continue;

      if (msg.msgid == MAVLINK_MSG_ID_HEARTBEAT && !d.requested) {
        d.requestAll(cfg_.request_hz);
        d.requested = true;
      } else if (msg.msgid == MAVLINK_MSG_ID_ATTITUDE) {
        mavlink_attitude_t a;
        mavlink_msg_attitude_decode(&msg, &a);
        d.C_latest = eulerToDcm(a.roll, a.pitch, a.yaw);
        d.yaw_latest = a.yaw;
        d.omega_latest = Eigen::Vector3d(a.rollspeed, a.pitchspeed, a.yawspeed);
        d.cur.t = a.time_boot_ms * 1e-3;
        d.t_att = nowSeconds();
        d.have_att = true;
        ++d.n_att;
      } else if (msg.msgid == MAVLINK_MSG_ID_SIMSTATE) {
        mavlink_simstate_t s;
        mavlink_msg_simstate_decode(&msg, &s);
        d.cur.true_roll = s.roll;
        d.cur.true_pitch = s.pitch;
        d.cur.true_yaw = s.yaw;
        d.cur.true_pos.lat = s.lat * 1e-7 * kDeg2Rad;
        d.cur.true_pos.lon = s.lng * 1e-7 * kDeg2Rad;
        d.cur.true_pos.alt = d.cur.alt_m;
        d.t_sim = nowSeconds();
        ++d.n_sim;

        // PAIR AT ARRIVAL, extrapolating the AHRS over the gap with its own
        // reported body rates. Reconciling timestamps later cannot work:
        // SIMSTATE has none.
        d.cur.skew_s = d.have_att ? (d.t_sim - d.t_att) : 1e9;
        const double dt = std::clamp(d.cur.skew_s, 0.0, 0.2);
        const Eigen::Vector3d rv = d.omega_latest * dt;
        const double ang = rv.norm();
        const Eigen::Matrix3d dC =
            (ang < 1e-9) ? Eigen::Matrix3d::Identity()
                         : Eigen::AngleAxisd(ang, rv / ang).toRotationMatrix();
        d.cur.C_l_b_est = d.C_latest * dC;
        d.cur.yaw_est = d.yaw_latest + d.omega_latest.z() * dt;
        d.cur.omega = d.omega_latest;
        d.have_sim = true;
      } else if (msg.msgid == MAVLINK_MSG_ID_VFR_HUD) {
        mavlink_vfr_hud_t v;
        mavlink_msg_vfr_hud_decode(&msg, &v);
        d.cur.airspeed = v.airspeed;
      } else if (msg.msgid == MAVLINK_MSG_ID_GLOBAL_POSITION_INT) {
        mavlink_global_position_int_t g;
        mavlink_msg_global_position_int_decode(&msg, &g);
        d.cur.alt_m = g.alt * 1e-3;
      }
    }
  }

  if (!d.have_sim || !d.have_att) return false;
  if (d.cur.skew_s > d.skew_limit) { ++d.n_reject; return false; }
  d.cur.ok = true;
  out = d.cur;
  return true;
}

MavlinkSource::Rates MavlinkSource::rates() {
  Impl& d = *impl_;
  Rates r;
  const double now = nowSeconds();
  const double dt = now - d.t_report;
  if (dt < 5.0) { r.skew_limit_s = d.skew_limit; return r; }

  r.simstate_hz = d.n_sim / dt;
  r.attitude_hz = d.n_att / dt;
  r.skew_rejects = d.n_reject;
  if (r.attitude_hz > 0.5) d.hz_att = r.attitude_hz;
  if (cfg_.max_skew_s <= 0) {
    // Pairing happens at arrival, so residual skew is bounded by one ATTITUDE
    // period. This threshold only has to catch stalls.
    d.skew_limit = std::clamp(1.5 / d.hz_att, 0.05, 0.50);
  }
  r.skew_limit_s = d.skew_limit;
  r.warn = r.simstate_hz < 5.0;
  if (r.warn) d.requestAll(cfg_.request_hz);
  d.n_sim = d.n_att = d.n_reject = 0;
  d.t_report = now;
  return r;
}

void MavlinkSource::sendFix(const Geodetic& p, double h_acc_m, double t) {
  mavlink_message_t m;
  const uint16_t ignore =
      GPS_INPUT_IGNORE_FLAG_VEL_HORIZ | GPS_INPUT_IGNORE_FLAG_VEL_VERT |
      GPS_INPUT_IGNORE_FLAG_SPEED_ACCURACY | GPS_INPUT_IGNORE_FLAG_VDOP;
  mavlink_msg_gps_input_pack(
      kSysId, kCompId, &m, uint64_t(t * 1e6), 0, ignore, 0, 0, 3,
      int32_t(p.lat * kRad2Deg * 1e7), int32_t(p.lon * kRad2Deg * 1e7),
      float(p.alt), 1.0f, 1.0f, 0, 0, 0, 0.0f, float(h_acc_m), 10.0f, 10, 0);
  impl_->send(m);
}

}  // namespace celestial
