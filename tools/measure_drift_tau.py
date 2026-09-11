#!/usr/bin/env python3
"""Measure ahrs_drift_tau from a flight log instead of assuming it.

    python3 tools/dump_log.py logs/NNNNNNNN.BIN > flight.csv
    python3 tools/measure_drift_tau.py flight.csv

WHY. `ahrs_drift_tau` and `ahrs_drift_sigma` are the two parameters this project
never measured -- they come from the source paper's Figure 5. NOTES-eskf.md shows
the entire case for an error-state filter turns on tau: at 60 s nothing an orbit
learns about tilt survives to the next orbit, at 600 s a fifth of it does.

WHAT IT ACTUALLY MEASURES. Not the gyro. `ahrs_drift_tau` is the correlation
time of the EKF3 TILT ERROR, which is a closed-loop property of the filter --
how fast its accelerometer aiding pulls attitude back toward vertical against
gyro integration. A SITL log gives the real EKF3 algorithm running on the IMU
model in ardupilot/params/cns_sitl.parm, so the answer is as good as that model
and much better than a number read off someone else's figure. It is NOT a real
airframe; state it that way.

MEASURE IN THE DENIED PHASE. With GPS the EKF has position aiding and its tilt
error behaves differently. `--after` skips the aided portion.

STRAIGHT LEGS ONLY, BY DEFAULT. The real EKF3 tilt error carries a
one-per-revolution component during orbits (see CLAUDE.md), which shows up as
ringing in the autocorrelation and biases a naive exponential fit. `--max-rate`
keeps only low-turn-rate samples so what is fitted is drift, not the orbit.
"""
import argparse, csv, math

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("csv")
ap.add_argument("--after", type=float, default=None,
                help="ignore samples before this time, s (skip the GPS-aided phase)")
ap.add_argument("--max-rate", type=float, default=1.0,
                help="deg/s of heading change; above this counts as a turn")
ap.add_argument("--max-lag", type=float, default=1200.0, help="seconds")
a = ap.parse_args()

t, dr, dp, yaw = [], [], [], []
for r in csv.DictReader(open(a.csv)):
    try:
        tt = float(r["t"])
        t.append(tt)
        dr.append(float(r["roll_ekf"]) - float(r["roll_true"]))
        dp.append(float(r["pitch_ekf"]) - float(r["pitch_true"]))
        yaw.append(float(r["yaw_ekf"]))
    except (ValueError, KeyError):
        continue
if len(t) < 500:
    raise SystemExit(f"only {len(t)} rows with both ATT and SIM in {a.csv}")

t = np.array(t); dr = np.array(dr); dp = np.array(dp); yaw = np.unwrap(np.array(yaw))
print(f"  {len(t)} samples, {t[0]:.0f}..{t[-1]:.0f} s, mean rate "
      f"{len(t)/(t[-1]-t[0]):.1f} Hz")

# radians in the log? roll/pitch are small; if the spread looks like degrees, say so
unit = "rad" if np.percentile(np.abs(dr), 99) < 0.2 else "deg"
scale = (180.0 / math.pi * 60.0) if unit == "rad" else 60.0   # -> arcmin
dr *= scale; dp *= scale
print(f"  attitude columns look like {unit}; tilt error converted to arcmin")

keep = np.ones(len(t), bool)
if a.after is not None:
    keep &= t >= a.after
rate = np.abs(np.gradient(np.degrees(yaw) if unit == "rad" else yaw, t))
keep &= rate < a.max_rate
print(f"  keeping {keep.sum()} samples "
      f"({100*keep.sum()/len(t):.0f}%) after --after and --max-rate {a.max_rate} deg/s")

# Resample to a uniform grid; autocorrelation needs even spacing.
dt = 0.25
grid = np.arange(t[keep][0], t[keep][-1], dt)
out = {}
for name, series in (("roll (lateral tilt)", dr), ("pitch (forward tilt)", dp)):
    v = np.interp(grid, t[keep], series[keep])
    v = v - v.mean()
    n = len(v)
    ac = np.correlate(v, v, mode="full")[n - 1:]
    ac /= ac[0]
    lags = np.arange(len(ac)) * dt
    m = lags <= a.max_lag
    ac, lags = ac[m], lags[m]

    # DO NOT fit a single exponential. Measured on SITL the autocorrelation has
    # a fast component (vibration, gyro noise) on top of a slow one (bias
    # drift), and a one-exponential fit returns whichever the fitting window
    # happened to see -- 92 s and 3793 s from the same series. Report the
    # empirical curve at the lags that matter instead: the ESKF question is
    # literally "what fraction survives the gap between orbits", and that is
    # read off directly, with no model assumed.
    def at(lag):
        i = int(round(lag / dt))
        return ac[i] if i < len(ac) else float("nan")

    out[name] = (v.std(), at(950.0))
    print(f"\n  {name}")
    print(f"    1-sigma  {v.std():7.2f} arcmin  ({v.std()/60:.3f} deg)")
    print(f"    autocorrelation:", "  ".join(
        f"{L:g}s {at(L):+.3f}" for L in (1, 10, 60, 120, 300, 630, 950)))
    print(f"    -> {100*max(0.0, at(950.0)):.1f} % of the tilt error survives the "
          f"950 s gap between orbits")

print(f"\n  model defaults: ahrs_drift_sigma 0.15 deg = 9.0 arcmin, "
      f"ahrs_drift_tau 60 s")
print("  NOTE: SITL's IMU model, real EKF3. Not a real airframe.")
