#!/usr/bin/env python3
"""Plots from a celestial_node CSV. Writes PNGs, no display needed.

    python3 tools/plot_trajectory.py live.csv trajectory.png error.png

trajectory.png  true track and per-frame fixes, in a fixed box (--box, 5 km)
                Fixes outside the box are dropped rather than rescaling it:
                before the camera mounting is calibrated they sit tens of km
                out, and letting them set the scale makes the ~1 km track a
                single pixel.

error.png       position error over time. The line is a rolling mean; the
                faint dots are individual per-frame fixes; red dots are the
                orbit-averaged fixes.

WHAT THE PER-FRAME FIX IS: a single-frame celestial position, computed from
the horizon-fused attitude (IMU corrected by the horizon sensor) and the
current camera mounting estimate. It therefore carries the boresight error
until the first recalibration, which is why it starts tens of km out.
"""
import argparse
import csv
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = 6371.0
NEED = ("t", "lat_true", "lon_true", "lat_fix", "lon_fix", "frame_err_m",
        "tilt_imu_fwd_am", "tilt_imu_lat_am", "tilt_fused_fwd_am",
        "tilt_fused_lat_am", "yaw_deg")


def load(path):
    raw = list(csv.DictReader(open(path)))
    if not raw:
        sys.exit(f"{path}: empty")
    missing = [k for k in NEED if k not in raw[0]]
    if missing:
        sys.exit(f"{path}: missing column(s) {', '.join(missing)}")
    rows = [x for x in raw if all(x.get(k) not in (None, "") for k in NEED)]
    if len(rows) < len(raw):
        print(f"  dropped {len(raw)-len(rows)} incomplete row(s)", file=sys.stderr)
    if not rows:
        sys.exit(f"{path}: no complete rows")
    return rows


def rolling(t, v, win):
    """Mean of v over a trailing window of `win` seconds."""
    out, j = [], 0
    for i in range(len(t)):
        while t[i] - t[j] > win:
            j += 1
        seg = v[j:i + 1]
        out.append(sum(seg) / len(seg))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("traj", nargs="?", default="trajectory.png")
    ap.add_argument("err", nargs="?", default="error.png")
    ap.add_argument("--box", type=float, default=5.0, help="km, full width")
    ap.add_argument("--window", type=float, default=20.0,
                    help="seconds, rolling mean for the error line")
    ap.add_argument("--center", choices=("start", "mean"), default="start",
                    help="box centre: the starting position (default) or the "
                         "mean of the whole track. 'mean' drifts off the start "
                         "whenever the flight includes a cruise leg.")
    a = ap.parse_args()

    r = load(a.csv)
    t = [float(x["t"]) for x in r]
    half = a.box / 2.0

    if a.center == "start":
        lat0, lon0 = float(r[0]["lat_true"]), float(r[0]["lon_true"])
    else:
        lat0 = sum(float(x["lat_true"]) for x in r) / len(r)
        lon0 = sum(float(x["lon_true"]) for x in r) / len(r)
    kx = math.cos(math.radians(lat0)) * math.pi / 180 * R
    ky = math.pi / 180 * R
    xy = lambda la, lo: ((lo - lon0) * kx, (la - lat0) * ky)

    tx, ty, fx, fy, fa = [], [], [], [], []
    dropped = 0
    for x in r:
        px, py = xy(float(x["lat_true"]), float(x["lon_true"]))
        tx.append(px); ty.append(py)
        if float(x["lat_fix"]) != 0:
            px, py = xy(float(x["lat_fix"]), float(x["lon_fix"]))
            if abs(px) <= half and abs(py) <= half:
                fx.append(px); fy.append(py); fa.append(float(x["t"]))
            else:
                dropped += 1

    # Smoothed fix track: a rolling mean of the per-frame fixes, so it can be
    # compared against the true track point for point. Individual fixes are too
    # scattered to read; the mean is what the system actually knows.
    sm_t, sm_x, sm_y = [], [], []
    j = 0
    for i in range(len(fa)):
        while fa[i] - fa[j] > a.window:
            j += 1
        seg = slice(j, i + 1)
        sm_t.append(fa[i])
        sm_x.append(sum(fx[seg]) / (i - j + 1))
        sm_y.append(sum(fy[seg]) / (i - j + 1))

    # Turn rate, to mark where the aircraft was loitering vs flying straight.
    yaw = [float(x["yaw_deg"]) for x in r]
    turning = []
    for i in range(len(r)):
        k0, k1 = max(0, i - 2), min(len(r) - 1, i + 2)
        dt = t[k1] - t[k0]
        d = abs(((yaw[k1] - yaw[k0] + 180) % 360 - 180) / dt) if dt > 1e-3 else 0
        turning.append(d > 2.0)

    # ---- trajectory + error scatter -------------------------------------
    # Two panels, because they answer different questions at scales that
    # cannot share an axis. A typical run covers 3 km of ground while the
    # fixes sit tens of km out, so an absolute map shows one or the other.
    #
    # LEFT   what the aircraft did, auto-scaled to the track.
    # RIGHT  where each fix landed RELATIVE TO TRUTH AT THAT INSTANT. This is
    #        the scale-free view: the aircraft's own motion drops out, so the
    #        circle traced by a mounting error is visible directly and the
    #        picture stays meaningful whether the error is 40 km or 400 m.
    fig, axs = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)
    ax = axs[0]
    if fx:
        sc = ax.scatter(fx, fy, c=fa, cmap="viridis", s=11, alpha=0.45, lw=0,
                        zorder=3, label="per-frame fix")
        cb = fig.colorbar(sc, ax=ax, shrink=0.82)
        cb.set_label("time (s)", fontsize=9)
    if sm_x:
        ax.plot(sm_x, sm_y, lw=1.8, color="#4C3FBF", zorder=5,
                label=f"fix, {a.window:.0f} s mean")
    # True track, dashed where straight and solid where turning, so the loiter
    # and cruise phases are distinguishable.
    for i in range(1, len(tx)):
        ax.plot(tx[i-1:i+1], ty[i-1:i+1], lw=2.2, zorder=4,
                color="#111" if turning[i] else "#999",
                ls="-" if turning[i] else "--")
    # Time markers every 30 s on BOTH tracks, so a moment on the trajectory can
    # be matched to the fix at that moment.
    step = 30.0
    nxt = t[0] + step
    for i in range(len(r)):
        if t[i] >= nxt:
            ax.plot(tx[i], ty[i], "|", color="#111", ms=9, zorder=7)
            ax.annotate(f"{t[i]:.0f}", (tx[i], ty[i]), fontsize=7,
                        color="#111", textcoords="offset points",
                        xytext=(4, 4), zorder=8)
            nxt += step
    nxt = (sm_t[0] + step) if sm_t else 1e18
    for i in range(len(sm_t)):
        if sm_t[i] >= nxt:
            ax.plot(sm_x[i], sm_y[i], "|", color="#4C3FBF", ms=8, zorder=7)
            ax.annotate(f"{sm_t[i]:.0f}", (sm_x[i], sm_y[i]), fontsize=7,
                        color="#4C3FBF", textcoords="offset points",
                        xytext=(4, -9), zorder=8)
            nxt += step
    ax.plot([], [], lw=2.2, color="#111", label="true track (turning)")
    ax.plot([], [], lw=2.2, color="#999", ls="--", label="true track (straight)")
    # Dead-reckoned track: airspeed + heading, corrected by each celestial fix.
    # This is the trajectory the system actually produces -- the fixes bound it
    # absolutely, the DR gives it shape between them.
    if "dr_lat" in r[0]:
        dxy = [xy(float(x["dr_lat"]), float(x["dr_lon"])) for x in r
               if x.get("dr_lat") not in (None, "") and float(x["dr_lat"]) != 0]
        if dxy:
            ax.plot([p[0] for p in dxy], [p[1] for p in dxy], lw=1.6,
                    color="#E24B4A", zorder=5, label="dead reckoning + fixes")
    ax.plot(tx[0], ty[0], "o", color="#1D9E75", ms=8, zorder=6, label="start")
    ax.plot(tx[-1], ty[-1], "s", color="#E24B4A", ms=7, zorder=6, label="end")
    tspan = max(0.3, max(abs(v) for v in tx + ty) * 1.25)
    ax.set_xlim(-tspan, tspan); ax.set_ylim(-tspan, tspan)
    ax.set_aspect("equal")
    ax.set_xlabel("east (km)"); ax.set_ylabel("north (km)")
    ax.set_title(f"flight path   (+/-{tspan:.1f} km)\n"
                 f"time marks every {step:.0f} s", fontsize=10)
    ax.grid(True, lw=0.4, alpha=0.3)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=8, frameon=False, loc="upper left")

    # ---- right: fix error relative to truth ------------------------------
    ex, ey, ea = [], [], []
    for x in r:
        if float(x["lat_fix"]) == 0:
            continue
        la, lo = float(x["lat_true"]), float(x["lon_true"])
        k2 = math.cos(math.radians(la)) * math.pi / 180 * R
        ex.append((float(x["lon_fix"]) - lo) * k2)
        ey.append((float(x["lat_fix"]) - la) * ky)
        ea.append(float(x["t"]))
    ax2 = axs[1]
    if ex:
        sc2 = ax2.scatter(ex, ey, c=ea, cmap="viridis", s=12, alpha=0.7, lw=0)
        cb = fig.colorbar(sc2, ax=ax2, shrink=0.82)
        cb.set_label("time (s)", fontsize=9)
    ax2.plot(0, 0, "k+", ms=14, mew=1.6, zorder=6)
    lim = max(1.0, sorted(math.hypot(u, v) for u, v in zip(ex, ey))
              [int(len(ex) * 0.95)] * 1.15) if ex else 1.0
    ax2.set_xlim(-lim, lim); ax2.set_ylim(-lim, lim)
    ax2.set_aspect("equal")
    ax2.set_xlabel("east error (km)"); ax2.set_ylabel("north error (km)")
    out95 = sum(1 for u, v in zip(ex, ey) if math.hypot(u, v) > lim)
    ax2.set_title(f"fix error relative to truth   (+/-{lim:.1f} km)\n"
                  f"{out95} beyond the axis", fontsize=10)
    ax2.grid(True, lw=0.4, alpha=0.3)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)

    fig.savefig(a.traj, dpi=150)
    print(f"wrote {a.traj}", file=sys.stderr)

    # ---- error over time ------------------------------------------------
    et = [t[i] for i in range(len(r)) if float(r[i]["frame_err_m"]) >= 0]
    ee = [float(r[i]["frame_err_m"]) / 1000 for i in range(len(r))
          if float(r[i]["frame_err_m"]) >= 0]

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    if ee:
        ax.plot(et, ee, ".", ms=3.5, color="#7F77DD", alpha=0.45, zorder=2,
                label="per-frame fix")
        ax.plot(et, rolling(et, ee, a.window), lw=2.0, color="#4C3FBF",
                zorder=4, label=f"mean over {a.window:.0f} s")
    orb = [(float(x["t"]), float(x["orbit_err_m"]) / 1000) for x in r
           if x.get("orbit_err_m") and float(x["orbit_err_m"]) >= 0]
    if orb:
        ax.plot([p[0] for p in orb], [p[1] for p in orb], "o", color="#E24B4A",
                ms=5, zorder=6, label="orbit fix")
    ax.set_yscale("log")
    ax.set_xlabel("time (s)"); ax.set_ylabel("position error (km)")
    ax.set_title("celestial position error", fontsize=11)
    ax.grid(True, which="both", lw=0.4, alpha=0.3)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=8.5, frameon=False)
    fig.savefig(a.err, dpi=150)
    print(f"wrote {a.err}", file=sys.stderr)

    # ---- summary --------------------------------------------------------
    imu = [math.hypot(float(x["tilt_imu_fwd_am"]),
                      float(x["tilt_imu_lat_am"])) for x in r]
    fus = [math.hypot(float(x["tilt_fused_fwd_am"]),
                      float(x["tilt_fused_lat_am"])) for x in r]
    mi, mf = sorted(imu)[len(imu)//2], sorted(fus)[len(fus)//2]
    print(f"  vertical: IMU {mi:.2f}' -> fused {mf:.2f}'"
          f"  ({mi*111.2/60:.1f} km -> {mf*111.2/60:.2f} km)", file=sys.stderr)
    if ee:
        e = sorted(ee)
        print(f"  per-frame fixes: {len(ee)}  median {e[len(e)//2]:.2f} km"
              f"  best {e[0]:.2f} km", file=sys.stderr)
    for k, (tt, e) in enumerate(orb):
        print(f"  orbit fix {k+1} at t={tt:.0f}s: {e:.2f} km", file=sys.stderr)
    if "det" in r[0]:
        d = [int(x["det"]) for x in r if x.get("det")]
        sm = [int(x["stars"]) for x in r if x.get("det")]
        if d:
            print(f"  detected mean {sum(d)/len(d):.1f}, matched mean "
                  f"{sum(sm)/len(sm):.1f} ({100*sum(sm)/max(1,sum(d)):.0f}%)",
                  file=sys.stderr)
    if "dr_err_m" in r[0]:
        de = [float(x["dr_err_m"]) for x in r
              if x.get("dr_err_m") not in (None, "") and float(x["dr_err_m"]) > 0]
        if de:
            d = sorted(de)
            print(f"  DR+fix error: median {d[len(d)//2]/1000:.2f} km"
                  f"  max {d[-1]/1000:.2f} km", file=sys.stderr)
    if "dr_open_err_m" in r[0]:
        oe = [float(x["dr_open_err_m"]) for x in r
              if x.get("dr_open_err_m") not in (None, "")
              and float(x["dr_open_err_m"]) > 0]
        if oe:
            o = sorted(oe)
            print(f"  DR ALONE (no fixes): median {o[len(o)//2]/1000:.2f} km"
                  f"  final {oe[-1]/1000:.2f} km", file=sys.stderr)
    if "boresight_deg" in r[0]:
        b = [(float(x["t"]), float(x["boresight_deg"])) for x in r
             if x.get("boresight_deg") not in (None, "")]
        if b:
            print(f"  boresight: {b[0][1]:.4f} deg -> {b[-1][1]:.4f} deg"
                  f"  ({b[-1][1]*111.2:.1f} km equivalent)", file=sys.stderr)
            if b[-1][1] * 111.2 > 3.0:
                print("    still converging: recalibration runs AT the "
                      "estimated position, so each\n    orbit improves it. "
                      "Fly longer -- 5+ orbit fixes.", file=sys.stderr)
    if "calibrated" in r[0]:
        c = [float(x["t"]) for x in r if x.get("calibrated") == "1"]
        if c:
            print(f"  mounting calibrated from t={c[0]:.0f}s"
                  f"  ({100*len(c)/len(r):.0f}% of the run)", file=sys.stderr)
    if "exposure_ms" in r[0]:
        ex = sorted(float(x["exposure_ms"]) for x in r if x.get("exposure_ms"))
        if ex:
            print(f"  exposure ms: median {ex[len(ex)//2]:.0f}"
                  f"  range {ex[0]:.0f}-{ex[-1]:.0f}", file=sys.stderr)


if __name__ == "__main__":
    main()
