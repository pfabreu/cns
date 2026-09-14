#!/usr/bin/env python3
"""Live two-pane view of a celestial navigation run.

    ./build/celestial_node --port 14556 --no-inject \
        --csv live.csv --frames-dir frames --frames-every 20 &
    ./build/fake_sitl --port 14556 --speed 25          # or ArduPilot SITL
    python3 tools/plot/live_view.py --csv live.csv --frames frames

LEFT PANE  -- ground truth against the estimate. Three things are drawn and the
   distinction matters: the TRUTH track, the DEAD RECKONED track (airspeed,
   heading and estimated wind, which drifts without limit), and the CELESTIAL
   FIXES (which do not). The point of the system is that the second is bounded
   by the third, so the interesting thing to watch is the DR track pulling back
   toward truth at each fix rather than the fixes themselves being accurate.

RIGHT PANE -- what the camera saw and what the matcher made of it. Green is a
   detection the matcher identified against the catalogue, red is one it
   rejected. Rejections are not failures: the ambiguity guard drops anything
   with a second plausible candidate nearby, and that is why a noisy detector
   costs little downstream.

Works live -- both files are polled and the newest frame is picked up -- or
after the fact on a finished run. It reads only what the node already writes,
so nothing here is a separate simulation.

Needs matplotlib. Nothing else in this project does, which is why it is a
script rather than a build target.
"""

import argparse
import csv
import glob
import math
import os

import sys

import numpy as np

# --save renders without a display, and the backend has to be chosen before
# pyplot is imported -- hence the argv peek rather than a tidier flag check.
if "--save" in sys.argv:
    import matplotlib
    matplotlib.use("Agg")

try:
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
except ImportError:
    raise SystemExit("needs matplotlib: pip install matplotlib")


def read_pgm(path):
    """Binary PGM, 8- or 16-bit. Same reader as scripts/train_unet.py."""
    with open(path, "rb") as f:
        if f.readline().strip() != b"P5":
            return None
        dims = f.readline().split()
        w, h = int(dims[0]), int(dims[1])
        maxval = int(f.readline())
        n = w * h
        if maxval > 255:
            data = np.frombuffer(f.read(2 * n), dtype=">u2").astype(np.float32)
        else:
            data = np.frombuffer(f.read(n), dtype=np.uint8).astype(np.float32)
    if data.size != n:
        return None
    return data.reshape(h, w)


def _agg(v):
    """median and RMS of a list, for the title.

    The LAST value is a trap: dr_err_m is written only on fix rows and the
    aided error is a sawtooth, so the final sample is wherever the flight
    happened to stop in the cycle. On the Sagres->Porto Santo crossing the two
    arms ended on 3.67 km and 28.11 km -- a 7.7x gap that is pure sampling;
    by median they are 6.97 and 7.79. Quote the distribution in the title."""
    if not v:
        return 0.0, 0.0
    v = sorted(v)
    med = v[len(v) // 2] if len(v) % 2 else 0.5 * (v[len(v)//2 - 1] + v[len(v)//2])
    return med, math.sqrt(sum(x * x for x in v) / len(v))


def read_track(path):
    """Truth, dead reckoning and fixes, in metres from the first truth point."""
    try:
        rows = list(csv.DictReader(open(path)))
    except OSError:
        return None
    if not rows:
        return None

    lat0 = float(rows[0]["lat_true"])
    lon0 = float(rows[0]["lon_true"])
    k = math.cos(math.radians(lat0))

    def ne(lat, lon):
        return ((lon - lon0) * 111320.0 * k, (lat - lat0) * 111320.0)

    truth, dr, dr_open, fixes = [], [], [], []
    for r in rows:
        try:
            truth.append(ne(float(r["lat_true"]), float(r["lon_true"])))
            if float(r.get("dr_lat", 0)):
                dr.append(ne(float(r["dr_lat"]), float(r["dr_lon"])))
            if float(r.get("dr_open_lat", 0)):
                dr_open.append(ne(float(r["dr_open_lat"]),
                                  float(r["dr_open_lon"])))
            if float(r.get("orbit_lat", 0)):
                x, y = ne(float(r["orbit_lat"]), float(r["orbit_lon"]))
                fixes.append((x, y, float(r.get("t", 0)),
                              float(r.get("orbit_err_m", -1))))
        except (ValueError, KeyError):
            continue

    # Error columns are written only on FIX rows and are -1 otherwise, so take
    # the last row that actually carries them rather than the last row.
    def last_valid(key, default=-1.0):
        for r in reversed(rows):
            try:
                v = float(r.get(key, default))
            except ValueError:
                continue
            if v >= 0:
                return v
        return default

    def all_valid(key):
        out = []
        for r in rows:
            try:
                v = float(r.get(key, -1))
            except ValueError:
                continue
            if v >= 0:
                out.append(v)
        return out

    dr_med, dr_rms = _agg(all_valid("dr_err_m"))
    stats = {
        "t": float(rows[-1].get("t", 0)),
        "dr_err": last_valid("dr_err_m"),      # kept, but NOT for the title
        "dr_med": dr_med,
        "dr_rms": dr_rms,
        "dr_open": last_valid("dr_open_err_m"),
        "bore": last_valid("boresight_deg", 0.0),
        "fixes": len(fixes),
    }
    return truth, dr, dr_open, fixes, stats


def newest_frame(d):
    """Newest complete image/detections pair. The node renames into place, so a
    pair is either fully written or absent -- never half read."""
    pgms = sorted(glob.glob(os.path.join(d, "f*.pgm")))
    for p in reversed(pgms):
        c = p[:-4] + ".csv"
        if os.path.exists(c):
            return p, c
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="live.csv")
    ap.add_argument("--frames", default="frames")
    ap.add_argument("--compare", default=None,
                    help="second CSV to overlay, e.g. the horizon run")
    ap.add_argument("--label", default="run A")
    ap.add_argument("--compare-label", default="run B")
    ap.add_argument("--interval", type=float, default=1.0, help="refresh, s")
    ap.add_argument("--legend-loc", default="upper left",
                    help="matplotlib legend location. Default upper left: "
                         "these missions run from north-east to south-west, so "
                         "the START of the track sits in the upper right and a "
                         "legend there covers it.")
    ap.add_argument("--fix-labels", type=int, default=1,
                    help="annotate every Nth fix with its time and error. "
                         "1 = all (default, useful live), 0 = none. An 883 km "
                         "crossing produces 100+ fixes and the labels then sit "
                         "on top of the tracks; the error-vs-time plot carries "
                         "the same timing information without the clutter.")
    ap.add_argument("--no-dr-open", action="store_true",
                    help="hide the unaided dead-reckoning track. It drifts "
                         "hundreds of km on a long crossing and sets the axis "
                         "limits, squashing the aided tracks -- and the aided "
                         "tracks are what a two-run comparison is about.")
    ap.add_argument("--save", default=None,
                    help="render ONE frame to this PNG and exit, no window. "
                         "Trajectory only -- the camera pane needs the node to "
                         "have been run with --frames-dir.")
    a = ap.parse_args()

    # Saving draws the trajectory alone. Keeping the two-pane layout would put
    # a 'waiting for frames' placeholder next to it in every saved image.
    if a.save:
        fig, ax_t = plt.subplots(1, 1, figsize=(11, 9))
        ax_i = None
    else:
        fig, (ax_t, ax_i) = plt.subplots(1, 2, figsize=(15, 7))
        fig.canvas.manager.set_window_title("celestial navigation -- live")
    shown = {"path": None}

    def update(_):
        # ---- left: trajectories -------------------------------------------
        tr = read_track(a.csv)
        ax_t.clear()
        if tr:
            truth, dr, dr_open, fixes, st = tr
            if truth:
                x, y = zip(*truth)
                ax_t.plot(x, y, "-", color="0.35", lw=1.4, label="truth")
                ax_t.plot(x[-1], y[-1], "o", color="0.2", ms=7)
            if dr:
                x, y = zip(*dr)
                ax_t.plot(x, y, "-", color="tab:blue", lw=1.2,
                          label=f"{a.label}: DR + fixes")
                ax_t.plot(x[-1], y[-1], "o", color="tab:blue", ms=7)
            # The open-loop control: same filter, no fixes applied. The gap
            # between this and the aided track IS the result -- one grows
            # without limit, the other does not.
            if dr_open and not a.no_dr_open:
                x, y = zip(*dr_open)
                ax_t.plot(x, y, "--", color="tab:red", lw=1.1, alpha=0.75,
                          label="dead reckoning ONLY (no fixes)")
                ax_t.plot(x[-1], y[-1], "s", color="tab:red", ms=6)
            if fixes:
                x = [p[0] for p in fixes]
                y = [p[1] for p in fixes]
                # Fixes share the colour of the trajectory they feed, so a
                # two-run overlay reads without a legend lookup.
                ax_t.plot(x, y, "*", color="tab:blue", ms=13, ls="none",
                          mec="k", mew=0.4,
                          label=f"{a.label}: fixes ({len(fixes)})")
                # Label each fix with its time and error. Which fixes are bad
                # and WHEN matters -- early ones are taken on an uncalibrated
                # mounting and are expected to be worse, so a plot without
                # times hides the convergence.
                labelled = (fixes[::a.fix_labels] if a.fix_labels > 0 else [])
                for px, py, pt, pe in labelled:
                    txt = f"{pt/60:.0f}'" if pe < 0 else f"{pt/60:.0f}' {pe/1000:.0f}km"
                    ax_t.annotate(txt, (px, py), fontsize=7,
                                  color="0.35", alpha=0.8,
                                  xytext=(4, 4), textcoords="offset points")
            ax_t.set_title(
                f"t = {st['t']:.0f} s    DR error median "
                f"{st['dr_med']/1000:.2f} km  RMS {st['dr_rms']/1000:.2f}"
                f"    unaided {st['dr_open']/1000:.2f} km"
                f"    boresight {st['bore']:.3f} deg")
            # ---- optional second run, same axes, matched colours ----------
            if a.compare:
                tr2 = read_track(a.compare)
                if tr2:
                    _, dr2, _, fixes2, st2 = tr2
                    if dr2:
                        x, y = zip(*dr2)
                        ax_t.plot(x, y, "-", color="tab:green", lw=1.2,
                                  label=f"{a.compare_label}: DR + fixes")
                        ax_t.plot(x[-1], y[-1], "o", color="tab:green", ms=7)
                    if fixes2:
                        ax_t.plot([p[0] for p in fixes2], [p[1] for p in fixes2],
                                  "*", color="tab:green", ms=13, ls="none",
                                  mec="k", mew=0.4,
                                  label=f"{a.compare_label}: fixes "
                                        f"({len(fixes2)})")
                    tail = ("" if a.no_dr_open
                            else f"     unaided {st['dr_open']/1000:.2f} km")
                    ax_t.set_title(
                        f"t = {st['t']:.0f} s     "
                        f"{a.label} median {st['dr_med']/1000:.2f} "
                        f"RMS {st['dr_rms']/1000:.2f} km     "
                        f"{a.compare_label} median {st2['dr_med']/1000:.2f} "
                        f"RMS {st2['dr_rms']/1000:.2f} km" + tail)
            ax_t.legend(loc=a.legend_loc, fontsize=9)
        else:
            ax_t.set_title(f"waiting for {a.csv}")
        ax_t.set_xlabel("east (m)")
        ax_t.set_ylabel("north (m)")
        ax_t.set_aspect("equal", "datalim")
        ax_t.grid(alpha=0.3)

        # ---- right: camera and matching -----------------------------------
        if ax_i is None:            # --save: trajectory only
            fig.tight_layout()
            return
        pgm, dcsv = newest_frame(a.frames)
        if pgm and pgm != shown["path"]:
            img = read_pgm(pgm)
            if img is not None:
                ax_i.clear()
                # Stars are faint against the background, so stretch hard:
                # median to the 99.9th percentile, which is roughly what the
                # detector's median+MAD threshold sees.
                lo = np.median(img)
                hi = np.percentile(img, 99.9)
                ax_i.imshow(img, cmap="gray", vmin=lo, vmax=max(hi, lo + 1),
                            origin="upper")
                nm = nu = 0
                try:
                    for r in csv.DictReader(open(dcsv)):
                        u, v = float(r["u"]), float(r["v"])
                        if int(r["matched"]):
                            ax_i.plot(u, v, "o", mfc="none", mec="lime",
                                      ms=13, mew=1.4)
                            nm += 1
                        else:
                            ax_i.plot(u, v, "o", mfc="none", mec="red",
                                      ms=9, mew=1.0)
                            nu += 1
                except OSError:
                    pass
                ax_i.set_title(f"{os.path.basename(pgm)}    "
                               f"matched {nm} (green)    rejected {nu} (red)")
                ax_i.set_xticks([])
                ax_i.set_yticks([])
                shown["path"] = pgm
        elif not pgm:
            ax_i.clear()
            ax_i.set_title(f"waiting for frames in {a.frames}/  "
                           f"(run the node with --frames-dir)")
            ax_i.set_xticks([])
            ax_i.set_yticks([])

        fig.tight_layout()

    if a.save:
        update(0)
        fig.savefig(a.save, dpi=130, bbox_inches="tight")
        print(f"wrote {a.save}")
        return

    ani = FuncAnimation(fig, update, interval=int(a.interval * 1000),
                        cache_frame_data=False)
    plt.show()


if __name__ == "__main__":
    main()
