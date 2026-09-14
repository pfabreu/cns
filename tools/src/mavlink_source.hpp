// mavlink_source.hpp — the TRANSPORT half of the live node.
//
// Everything to do with getting a consistent (true pose, AHRS attitude) pair
// out of a MAVLink stream, and nothing to do with navigation. See the system
// schematic: this is the box feeding the sensors in, and it exists separately
// because almost every bug found during SITL bring-up lived here rather than
// in the estimator.
//
// The four that shaped this file:
//
//   * SIMSTATE is in NO stream for ArduPlane -- only SET_MESSAGE_INTERVAL
//     reaches it, and the request must come from a GCS system id or MAVLink
//     routing drops it.
//   * MAVProxy re-sends REQUEST_DATA_STREAM at its own rate and overrides that
//     request, so both streams pin to the same low rate.
//   * SIMSTATE carries no timestamp, so the AHRS must be PAIRED AT ARRIVAL.
//     Comparing arrival times after the fact always shows one message period
//     of skew no matter how fast the streams run.
//   * The residual skew is bridged with the body rates ATTITUDE already
//     carries, which are estimated quantities -- the truth boundary holds.

#pragma once

#include <cstdint>
#include <vector>

#include "celestial/attitude.hpp"
#include "celestial/types.hpp"

namespace celestial {

/// One consistent sample: the simulator's true pose, and the AHRS attitude as
/// of the same instant.
struct MavSample {
  bool ok = false;
  double t = 0.0;                 ///< autopilot boot time, seconds
  Geodetic true_pos;              ///< SIMSTATE -- RENDERING AND SCORING ONLY
  double true_roll = 0, true_pitch = 0, true_yaw = 0;
  Eigen::Matrix3d C_l_b_est = Eigen::Matrix3d::Identity();  ///< ATTITUDE
  double yaw_est = 0.0;
  Eigen::Vector3d omega = Eigen::Vector3d::Zero();  ///< body rates
  double airspeed = 0.0;
  double alt_m = 800.0;
  double skew_s = 0.0;            ///< how far apart the two arrived
};

class MavlinkSource {
 public:
  struct Config {
    int port = 14556;
    double request_hz = 20.0;
    /// 0 = adaptive: 1.5 x the measured ATTITUDE period.
    double max_skew_s = 0.0;
  };

  explicit MavlinkSource(const Config& c) : cfg_(c) {}
  ~MavlinkSource();

  bool open();
  /// Read whatever is queued and update the current sample. Returns true when
  /// a fresh, consistently-paired sample is available.
  bool poll(MavSample& out);

  /// Send a position fix back as GPS_INPUT.
  void sendFix(const Geodetic& p, double h_acc_m, double t);

  /// Diagnostics, printed by the caller every few seconds.
  struct Rates { double simstate_hz = 0, attitude_hz = 0; int skew_rejects = 0;
                 double skew_limit_s = 0; bool warn = false; };
  Rates rates();

 private:
  struct Impl;
  Impl* impl_ = nullptr;
  Config cfg_;
};

}  // namespace celestial
