#!/usr/bin/env python3
"""Zoom on a fix orbit: the manoeuvre, and what it does to the fix cloud.

    python3 tools/plot/plot_orbit_zoom.py runs/<run>/live.csv out.png
    python3 tools/plot/plot_orbit_zoom.py runs/<run>/live.csv out.png --compare 2 4
    python3 tools/plot/plot_orbit_zoom.py runs/<run>/live.csv out.png --list

WHY TWO PANELS AND NOT ONE. The loiter is 250 m across. The per-frame fixes it
averages are scattered over tens of kilometres. There is no axis scale that
shows both, and the ratio IS the result: a single frame carries the whole AHRS
tilt error and lands ~30 km out, and only the 360 degree heading sweep makes
those errors cancel. So the left panel is the aircraft's actual path at orbit
scale, and the right panel is the same seconds at fix scale.

Deliberately no satellite basemap. At 250 m the imagery is meaningless and it
would only obscure the track.

Per-frame fix colour is heading, because that is the variable doing the work:
the body-fixed part of the tilt error rotates with the aircraft, so fixes from
opposite headings land on opposite sides of truth and average out. A fix orbit
that does not close the heading circle does not cancel anything -- 90 deg of
sweep gives 28 km, 180 gives 12, 270 gives 7, 360 gives 6.2.
"""
import argparse
import csv
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M_PER_DEG = 111320.0


def load(path):
    rows = list(csv.DictReader(open(path)))
    g = lambda k: np.array([float(r[k]) for r in rows])
    return dict(t=g("t"), lat=g("lat_true"), lon=g("lon_true"),
                flat=g("lat_fix"), flon=g("lon_fix"), ferr=g("frame_err_m"),
                olat=g("orbit_lat"), olon=g("orbit_lon"), oerr=g("orbit_err_m"),
                yaw=g("yaw_deg"),
                tfwd=g("tilt_imu_fwd_am"), tlat=g("tilt_imu_lat_am"))


def headingSpan(yaw_deg):
    """Degrees of heading actually covered: 360 minus the largest empty gap."""
    d = np.sort(np.mod(yaw_deg, 360.0))
    if len(d) < 2:
        return 0.0
    gaps = np.append(np.diff(d), 360.0 - d[-1] + d[0])
    return 360.0 - float(gaps.max())


def orbitStats(d, t0, t1):
    m = (d["t"] >= t0) & (d["t"] <= t1)
    yaw = np.unwrap(np.radians(d["yaw"][m]))
    fm = m & (d["ferr"] > 0)
    om = m & (d["oerr"] > 0)
    return dict(
        m=m, fm=fm, om=om,
        revs=abs(yaw[-1] - yaw[0]) / (2 * math.pi),
        span=headingSpan(d["yaw"][m]),
        # The component that survives a full heading sweep is the MEAN tilt over
        # the sweep, not its magnitude: the body-fixed part rotates with the
        # aircraft and cancels, the nav-frame part does not.
        ned=math.hypot(float(np.mean(d["tfwd"][m])), float(np.mean(d["tlat"][m]))),
        medframe=float(np.median(d["ferr"][fm])) / 1000.0 if fm.any() else float("nan"),
        best=float(d["oerr"][om].min()) / 1000.0 if om.any() else float("nan"))


def loiters(d, min_rate=1.5, min_len=60.0):
    """Windows where the aircraft is turning steadily. Returns (t0, t1) list."""
    rate = np.abs(np.gradient(np.unwrap(np.radians(d["yaw"])), d["t"])) * 180 / math.pi
    turning = rate > min_rate
    out, i = [], 0
    while i < len(turning):
        if not turning[i]:
            i += 1
            continue
        j = i
        while j < len(turning) and (turning[j] or
                                    (j + 25 < len(turning) and turning[j:j + 25].any())):
            j += 1
        if d["t"][j - 1] - d["t"][i] >= min_len:
            out.append((d["t"][i], d["t"][j - 1]))
        i = j + 1
    return out


def compare(d, ls, a):
    """Two loiters side by side at fix scale, on one shared axis scale.

    The geometry is the same every time -- these orbits all cover ~355 deg in
    two revolutions -- so what separates a good fix from a bad one is not the
    manoeuvre. It is the tilt error that survives the sweep.
    """
    ia, ib = a.compare
    for i in (ia, ib):
        if not 0 <= i < len(ls):
            raise SystemExit(f"loiter {i} out of range 0..{len(ls)-1}")
    stats = [orbitStats(d, *ls[i]) for i in (ia, ib)]

    lat0 = float(np.mean(d["lat"][stats[0]["m"]]))
    lon0 = float(np.mean(d["lon"][stats[0]["m"]]))
    kx = M_PER_DEG * math.cos(math.radians(lat0))

    # one scale for both panels, or the eye is being lied to
    r = 0.0
    for st in stats:
        c_lon = float(np.mean(d["lon"][st["m"]]))
        c_lat = float(np.mean(d["lat"][st["m"]]))
        fx = (d["flon"][st["fm"]] - c_lon) * kx / 1000.0
        fy = (d["flat"][st["fm"]] - c_lat) * M_PER_DEG / 1000.0
        if len(fx) > 10:
            r = max(r, np.percentile(np.abs(fx), a.clip),
                    np.percentile(np.abs(fy), a.clip))
    r *= 1.2

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.2))
    for ax, i, st in zip(axes, (ia, ib), stats):
        c_lon = float(np.mean(d["lon"][st["m"]]))
        c_lat = float(np.mean(d["lat"][st["m"]]))
        ex = lambda lo: (lo - c_lon) * kx / 1000.0
        no = lambda la: (la - c_lat) * M_PER_DEG / 1000.0
        sc = ax.scatter(ex(d["flon"][st["fm"]]), no(d["flat"][st["fm"]]),
                        c=np.mod(d["yaw"][st["fm"]], 360.0), s=8, cmap="hsv",
                        alpha=.6, vmin=0, vmax=360, zorder=2,
                        label=f"per-frame fixes ({int(st['fm'].sum())})")
        ax.plot([0], [0], "X", color="k", ms=16, markeredgecolor="w",
                markeredgewidth=1.2, zorder=6, label="true position")
        if st["om"].any():
            ax.plot(ex(d["olon"][st["om"]]), no(d["olat"][st["om"]]), "*",
                    color="#2166ac", ms=21, markeredgecolor="w",
                    markeredgewidth=.9, zorder=7,
                    label=f"orbit fix ({int(st['om'].sum())}), best "
                          f"{st['best']:.1f} km")
        ax.set_xlim(-r, r); ax.set_ylim(-r, r)
        ax.set_aspect("equal"); ax.grid(alpha=.3)
        ax.set_xlabel("east of loiter centre, km")
        if ax is axes[0]:
            ax.set_ylabel("north of loiter centre, km")
        ax.set_title(f"orbit {i} at t = {ls[i][0]:.0f} s — fix {st['best']:.1f} km")
        ax.legend(loc="upper left", fontsize=8)
        ax.text(.5, -.155,
                f"{st['revs']:.1f} revolutions, {st['span']:.0f}° of heading  |  "
                f"NED-mean tilt {st['ned']:.1f}′",
                transform=ax.transAxes, ha="center", fontsize=9)
    cb = fig.colorbar(sc, ax=axes, fraction=.03, pad=.02)
    cb.set_label("aircraft heading at that frame, deg")
    if a.title:
        fig.suptitle(a.title)
    fig.savefig(a.out, dpi=130, bbox_inches="tight")
    print(f"  wrote {a.out}")
    for i, st in zip((ia, ib), stats):
        print(f"    orbit {i}: {st['revs']:.1f} rev, {st['span']:.0f} deg span, "
              f"NED-mean tilt {st['ned']:.1f}', median frame {st['medframe']:.0f} km, "
              f"best fix {st['best']:.1f} km")



def together(d, ls, a):
    """Two fix orbits and the leg between them: the route, then the result.

    Three scales are in play and no single axis holds them. The loiter is 250 m
    across, the leg between two fix orbits is ~12 km, and the per-frame fixes
    those orbits average are scattered over 100 km. So: left panel is the route
    at track scale, where the two loiters read as small loops at either end of
    a straight leg; right panel is what those same orbits produced, at fix
    scale, on the same centre.

    Pick ADJACENT orbits. Two orbits further apart put a leg between them long
    enough that the loops shrink to dots.
    """
    ia, ib = a.together
    for i in (ia, ib):
        if not 0 <= i < len(ls):
            raise SystemExit(f"loiter {i} out of range 0..{len(ls)-1}")
    sa, sb = orbitStats(d, *ls[ia]), orbitStats(d, *ls[ib])

    lat0 = float(np.mean([np.mean(d["lat"][sa["m"]]), np.mean(d["lat"][sb["m"]])]))
    lon0 = float(np.mean([np.mean(d["lon"][sa["m"]]), np.mean(d["lon"][sb["m"]])]))
    kx = M_PER_DEG * math.cos(math.radians(lat0))
    ex = lambda lo: (np.asarray(lo) - lon0) * kx / 1000.0
    no = lambda la: (np.asarray(la) - lat0) * M_PER_DEG / 1000.0

    COL = ["#1b9e77", "#d95f02"]
    leg = (d["t"] >= ls[ia][0]) & (d["t"] <= ls[ib][1])
    flown = float(np.sum(np.hypot(np.diff(ex(d["lon"][leg])),
                                  np.diff(no(d["lat"][leg])))))

    # The route panel is tall and narrow (a straight leg with equal aspect),
    # so give it less width and hang its legend underneath rather than let
    # it cover a loiter.
    fig, (axT, axF) = plt.subplots(1, 2, figsize=(13.4, 7.0),
                                   gridspec_kw={"width_ratios": [1, 1.55]})

    # ---- left: the route, at track scale --------------------------------
    axT.plot(ex(d["lon"][leg]), no(d["lat"][leg]), "-", color="0.45", lw=1.3,
             zorder=2, label=f"flown track, {flown:.0f} km in "
                             f"{ls[ib][1]-ls[ia][0]:.0f} s")
    for k, (i, st) in enumerate(((ia, sa), (ib, sb))):
        axT.plot(ex(d["lon"][st["m"]]), no(d["lat"][st["m"]]), "-",
                 color=COL[k], lw=2.6, zorder=4, solid_capstyle="round",
                 label=f"orbit {i}: {st['revs']:.1f} rev, {st['span']:.0f}°, "
                       f"tilt {st['ned']:.1f}'")
        cx, cy = ex(np.mean(d["lon"][st["m"]])), no(np.mean(d["lat"][st["m"]]))
        axT.annotate(f"orbit {i}", (cx, cy), textcoords="offset points",
                     xytext=(16, 12), fontsize=9, color=COL[k], weight="bold",
                     arrowprops=dict(arrowstyle="-", color=COL[k], lw=.8))
    axT.set_aspect("equal"); axT.grid(alpha=.3)
    axT.set_xlabel("east, km"); axT.set_ylabel("north, km")
    axT.set_title("the route — two 250 m loiters, one leg between them")
    axT.legend(loc="upper center", bbox_to_anchor=(0.5, -0.09),
               fontsize=8.5, frameon=False)

    # ---- right: what they produced, at fix scale ------------------------
    for k, (i, st) in enumerate(((ia, sa), (ib, sb))):
        axF.scatter(ex(d["flon"][st["fm"]]), no(d["flat"][st["fm"]]), s=9,
                    color=COL[k], alpha=.45, zorder=2 + k,
                    label=f"orbit {i} per-frame fixes ({int(st['fm'].sum())})")
        if st["om"].any():
            axF.plot(ex(d["olon"][st["om"]]), no(d["olat"][st["om"]]), "*",
                     color=COL[k], ms=21, markeredgecolor="k",
                     markeredgewidth=.9, zorder=8,
                     label=f"orbit {i} fix, best {st['best']:.1f} km")
    axF.plot(ex(d["lon"][leg]), no(d["lat"][leg]), "-", color="k", lw=1.8,
             zorder=7, label="the same flown track")
    allx = np.concatenate([ex(d["flon"][sa["fm"]]), ex(d["flon"][sb["fm"]])])
    ally = np.concatenate([no(d["flat"][sa["fm"]]), no(d["flat"][sb["fm"]])])
    r = max(np.percentile(np.abs(allx), a.clip),
            np.percentile(np.abs(ally), a.clip)) * 1.15
    nout = int(((np.abs(allx) > r) | (np.abs(ally) > r)).sum())
    axF.set_xlim(-r, r); axF.set_ylim(-r, r)
    axF.set_aspect("equal"); axF.grid(alpha=.3)
    axF.set_xlabel("east, km")
    tspan = max(np.ptp(ex(d["lon"][leg])), np.ptp(no(d["lat"][leg])))
    axF.set_title(f"what they produced — this axis is {2*r/tspan:.0f}x wider")
    axF.legend(loc="upper left", fontsize=8.5, framealpha=.95)
    if nout:
        axF.text(.99, .01, f"{nout} fixes outside the frame", ha="right",
                 transform=axF.transAxes, fontsize=7.5, color="0.4")

    if a.title:
        fig.suptitle(a.title)
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print(f"  wrote {a.out}  ({flown:.0f} km flown between orbit {ia} and {ib})")
    for i, st in ((ia, sa), (ib, sb)):
        print(f"    orbit {i}: {st['revs']:.1f} rev, {st['span']:.0f} deg, "
              f"tilt {st['ned']:.1f}', median frame {st['medframe']:.0f} km, "
              f"best fix {st['best']:.1f} km")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("out")
    ap.add_argument("--t0", type=float, default=None)
    ap.add_argument("--t1", type=float, default=None)
    ap.add_argument("--orbit", type=int, default=None,
                    help="index into the detected loiters (default: the one "
                         "with the most published orbit fixes)")
    ap.add_argument("--pad", type=float, default=8.0, help="seconds either side")
    ap.add_argument("--clip", type=float, default=96.0,
                    help="percentile of the fix cloud to frame the right "
                         "panel on. A handful of fixes land 80 km out and "
                         "would otherwise squash everything into a dot.")
    ap.add_argument("--title", default="")
    ap.add_argument("--list", action="store_true", help="list loiters and exit")
    ap.add_argument("--compare", nargs=2, type=int, metavar=("A", "B"),
                    default=None,
                    help="two loiter indices, drawn side by side on a SHARED "
                         "axis scale so the comparison is visually fair")
    ap.add_argument("--together", nargs=2, type=int, metavar=("A", "B"),
                    default=None,
                    help="two loiter indices on ONE axis, with the flown track "
                         "between them")
    a = ap.parse_args()

    d = load(a.csv)
    ls = loiters(d)
    if a.together:
        together(d, ls, a)
        return
    if a.compare:
        compare(d, ls, a)
        return
    if a.list:
        print(f"  {'#':>3}{'t0':>8}{'t1':>8}{'revs':>6}{'hdg span':>10}"
              f"{'NED-mean tilt':>15}{'med frame':>11}{'best fix':>10}")
        for i, (t0, t1) in enumerate(ls):
            st = orbitStats(d, t0, t1)
            print(f"  {i:3d}{t0:8.0f}{t1:8.0f}{st['revs']:6.1f}{st['span']:9.0f}d"
                  f"{st['ned']:14.1f}'{st['medframe']:10.1f}k{st['best']:10.1f}k")
        return

    if a.t0 is None:
        if a.orbit is None:
            counts = [int(((d["t"] >= t0) & (d["t"] <= t1) & (d["oerr"] > 0)).sum())
                      for t0, t1 in ls]
            a.orbit = int(np.argmax(counts))
        t0, t1 = ls[a.orbit]
        t0, t1 = t0 - a.pad, t1 + a.pad
    else:
        t0, t1 = a.t0, (a.t1 if a.t1 is not None else a.t0 + 300)

    m = (d["t"] >= t0) & (d["t"] <= t1)
    if m.sum() < 20:
        raise SystemExit(f"only {m.sum()} rows in {t0:.0f}..{t1:.0f} s")
    lat0 = float(np.mean(d["lat"][m]))
    lon0 = float(np.mean(d["lon"][m]))
    kx = M_PER_DEG * math.cos(math.radians(lat0))
    east = lambda lo: (lo - lon0) * kx
    north = lambda la: (la - lat0) * M_PER_DEG

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.4, 5.9))

    # ---- left: the manoeuvre, at orbit scale ----------------------------
    x, y, tt = east(d["lon"][m]), north(d["lat"][m]), d["t"][m]
    sc = axL.scatter(x, y, c=tt - tt[0], s=5, cmap="viridis", zorder=3)
    axL.plot(x, y, "-", color="0.75", lw=.7, zorder=2)
    axL.plot(x[0], y[0], "o", color="#1b7837", ms=9, zorder=4, label="entry")
    axL.plot(x[-1], y[-1], "s", color="#b2182b", ms=8, zorder=4, label="exit")
    turns = abs(np.unwrap(np.radians(d["yaw"][m]))[-1] -
                np.unwrap(np.radians(d["yaw"][m]))[0]) / (2 * math.pi)
    axL.set_aspect("equal")
    axL.grid(alpha=.3)
    axL.set_xlabel("east of loiter centre, m")
    axL.set_ylabel("north of loiter centre, m")
    axL.set_title(f"the manoeuvre — {turns:.1f} revolutions in {tt[-1]-tt[0]:.0f} s")
    axL.legend(loc="upper right", fontsize=8)
    cb = fig.colorbar(sc, ax=axL, fraction=.046, pad=.04)
    cb.set_label("seconds into the window")

    # ---- right: what those same seconds produce, at fix scale ------------
    fm = m & (d["ferr"] > 0)
    fx = east(d["flon"][fm]) / 1000.0
    fy = north(d["flat"][fm]) / 1000.0
    hdg = np.mod(d["yaw"][fm], 360.0)
    s2 = axR.scatter(fx, fy, c=hdg, s=7, cmap="hsv", alpha=.55, vmin=0, vmax=360,
                     zorder=2, label=f"per-frame fixes ({fm.sum()})")
    # At this scale the 250 m loiter is a point, so mark it rather than draw it.
    axR.plot([east(np.mean(d["lon"][m])) / 1000.0],
             [north(np.mean(d["lat"][m])) / 1000.0], "X", color="k", ms=15,
             markeredgecolor="w", markeredgewidth=1.2, zorder=7,
             label="true position")
    om = m & (d["oerr"] > 0)
    if om.any():
        axR.plot(east(d["olon"][om]) / 1000.0, north(d["olat"][om]) / 1000.0, "*",
                 color="#2166ac", ms=20, markeredgecolor="w", markeredgewidth=.8,
                 zorder=6,
                 label=f"orbit-averaged fix ({om.sum()}, median "
                       f"{np.median(d['oerr'][om])/1000:.1f} km)")
    med = np.median(d["ferr"][fm]) / 1000.0 if fm.any() else float("nan")
    if fm.sum() > 10:
        q = a.clip
        r = max(np.percentile(np.abs(fx), q), np.percentile(np.abs(fy), q)) * 1.25
        cx0 = east(np.mean(d["lon"][m])) / 1000.0
        cy0 = north(np.mean(d["lat"][m])) / 1000.0
        axR.set_xlim(cx0 - r, cx0 + r)
        axR.set_ylim(cy0 - r, cy0 + r)
        nout = int((np.abs(fx - cx0) > r).sum() + (np.abs(fy - cy0) > r).sum())
        if nout:
            axR.text(.98, .02, f"{nout} fixes outside the frame", ha="right",
                     transform=axR.transAxes, fontsize=7, color="0.4")
    axR.set_aspect("equal")
    axR.grid(alpha=.3)
    axR.set_xlabel("east of loiter centre, km")
    axR.set_ylabel("north of loiter centre, km")
    axR.set_title(f"what they average — single frames land {med:.0f} km out")
    axR.legend(loc="upper left", fontsize=8)
    cb2 = fig.colorbar(s2, ax=axR, fraction=.046, pad=.04)
    cb2.set_label("aircraft heading at that frame, deg")

    if a.title:
        fig.suptitle(a.title)
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print(f"  wrote {a.out}  (t {t0:.0f}..{t1:.0f} s, {turns:.1f} revolutions, "
          f"{fm.sum()} per-frame fixes, {om.sum()} orbit fixes)")


if __name__ == "__main__":
    main()
