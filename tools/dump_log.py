#!/usr/bin/env python3
"""ArduPilot .BIN -> CSV for Milestone 3.

    pip install pymavlink
    python3 tools/dump_log.py flight.BIN > flight.csv

Emits one row per ATT message, carrying forward the most recent SIM, GPS and
ARSP values.

Messages used:
    ATT   EKF3 attitude -- what the navigation algorithm actually gets
    SIM   TRUE attitude and position -- SITL ground truth
    GPS   GWk / GMS for absolute time
    ARSP  airspeed
    IMU   raw gyro, for the star-aided attitude experiment. Body rates relative
          to INERTIAL space, rad/s, which is what an IMU reports and what
          StarAidedAttitude expects. Zero if the log has no IMU messages, in
          which case analyse_log --star-aided is unavailable.

SIM is only present in SITL logs. On a real flight log there is no attitude
truth, so the attitude-error decomposition is unavailable and you can only run
the position estimate itself.
"""
import sys

try:
    from pymavlink import mavutil
except ImportError:
    sys.exit("pip install pymavlink")


def get(msg, *names, default=0.0):
    for n in names:
        if hasattr(msg, n):
            return getattr(msg, n)
    return default


def main(path):
    mlog = mavutil.mavlink_connection(path)

    cols = ["t", "lat_true", "lon_true", "alt_true",
            "roll_true", "pitch_true", "yaw_true",
            "roll_ekf", "pitch_ekf", "yaw_ekf",
            "gps_week", "gps_ms", "airspeed",
            "gyro_x", "gyro_y", "gyro_z"]
    print(",".join(cols))

    sim = None
    gps = None
    arsp = None
    imu = None
    t0 = None
    n = 0

    while True:
        m = mlog.recv_match(type=["ATT", "SIM", "GPS", "ARSP", "IMU"])
        if m is None:
            break
        t = getattr(m, "TimeUS", None)
        if t is None:
            continue
        t = t * 1e-6

        typ = m.get_type()
        if typ == "SIM":
            sim = m
            continue
        if typ == "GPS":
            gps = m
            continue
        if typ == "ARSP":
            arsp = m
            continue
        if typ == "IMU":
            # First IMU only; later instances are the same rates from redundant
            # sensors and mixing them would alias.
            if imu is None or getattr(m, "I", 0) == getattr(imu, "I", 0):
                imu = m
            continue

        # ATT
        if sim is None:
            continue  # no truth yet
        if t0 is None:
            t0 = t

        row = [
            t - t0,
            get(sim, "Lat"), get(sim, "Lng"), get(sim, "Alt"),
            get(sim, "Roll"), get(sim, "Pitch"), get(sim, "Yaw"),
            get(m, "Roll"), get(m, "Pitch"), get(m, "Yaw"),
            get(gps, "GWk", default=-1) if gps else -1,
            get(gps, "GMS", default=-1) if gps else -1,
            get(arsp, "Airspeed", "AS", default=0.0) if arsp else 0.0,
            get(imu, "GyrX", default=0.0) if imu else 0.0,
            get(imu, "GyrY", default=0.0) if imu else 0.0,
            get(imu, "GyrZ", default=0.0) if imu else 0.0,
        ]
        print(",".join("%.9g" % v for v in row))
        n += 1

    print("wrote %d rows" % n, file=sys.stderr)
    if n == 0:
        print("No ATT+SIM pairs found. Is this a SITL log?", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1])
