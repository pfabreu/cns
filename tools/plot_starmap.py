#!/usr/bin/env python3
"""The camera's view of the sky at one celestial fix.

    python3 tools/plot_starmap.py live.csv out.png --utc 2024-12-15T19:00:00

Reconstructs, from the SAME Yale catalogue the matcher uses, which stars fall on
the sensor at the moment of a fix: zenith-pointing Alvium 1800 U-240, 1936x1216,
53.5 deg horizontal field, image +y along the aircraft heading.

This is computed from the catalogue and the logged pose, not read back from the
detector -- the node only logs star COUNTS, not identities, unless it was run
with --frames-dir. So it shows the sky that was there to be matched, and the
logged count of what actually was.
"""
import argparse, csv, datetime as dt, math, re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W, H, HFOV = 1936, 1216, math.radians(53.5)
FOCAL = (W / 2.0) / math.tan(HFOV / 2.0)
MAG_LIMIT, MATCH_MAG = 5.0, 4.0          # pipeline.hpp mag_limit / one brighter


def catalogue(path):
    """HR, RA deg J2000, Dec deg, dRA/dt "/yr, dDec/dt "/yr, Vmag."""
    out = []
    for m in re.finditer(r"\{(\d+),([-\d.]+),([-\d.]+),([-\d.]+),([-\d.]+),([-\d.]+)\}",
                         open(path).read()):
        g = m.groups()
        out.append((int(g[0]), float(g[1]), float(g[2]),
                    float(g[3]), float(g[4]), float(g[5])))
    return out


def gmst_deg(when):
    jd = 2451545.0 + (when - dt.datetime(2000, 1, 1, 12, 0, 0)).total_seconds() / 86400.0
    d = jd - 2451545.0
    t = d / 36525.0
    g = 280.46061837 + 360.98564736629 * d + 0.000387933 * t * t - t ** 3 / 38710000.0
    return g % 360.0


ap = argparse.ArgumentParser()
ap.add_argument("csv"); ap.add_argument("out")
ap.add_argument("--utc", default="2024-12-15T19:00:00",
                help="epoch the node was given; sim t is added to it")
ap.add_argument("--which", default="last",
                help="'last' fix, or a sim time in seconds to pick the nearest")
ap.add_argument("--catalog", default="src/star_catalog_data.cpp")
ap.add_argument("--place", default="")
ap.add_argument("--names", default="tools/star_names.tsv",
                help="HR -> designation/proper-name table, from BSC5")
ap.add_argument("--label-mag", type=float, default=MATCH_MAG,
                help="label stars at least this bright (default: the matcher's "
                     "own limit, so every labelled star is one it could use)")
a = ap.parse_args()

rows = [r for r in csv.DictReader(open(a.csv)) if float(r["orbit_err_m"]) >= 0]
if a.which == "last":
    fix = rows[-1]
else:
    tw = float(a.which)
    fix = min(rows, key=lambda r: abs(float(r["t"]) - tw))

t = float(fix["t"])
when = dt.datetime.fromisoformat(a.utc) + dt.timedelta(seconds=t)
lat, lon = float(fix["lat_true"]), float(fix["lon_true"])
err_km = float(fix["orbit_err_m"]) / 1000.0
yaw = math.radians(float(fix["yaw_deg"]))
n_logged = int(fix["stars"])

lst = (gmst_deg(when) + lon) % 360.0
years = (when - dt.datetime(2000, 1, 1, 12, 0, 0)).total_seconds() / (365.25 * 86400)
la = math.radians(lat)

# camera axes in ENU: boresight up, +y along heading, +x to starboard
fwd = (math.sin(yaw), math.cos(yaw), 0.0)
stb = (math.cos(yaw), -math.sin(yaw), 0.0)

pts = []
for hr, ra0, dec0, pmra, pmdec, vmag in catalogue(a.catalog):
    ra = ra0 + pmra * years / 3600.0
    dec = dec0 + pmdec * years / 3600.0
    ha = math.radians((lst - ra) % 360.0)
    dr = math.radians(dec)
    sin_alt = math.sin(dr) * math.sin(la) + math.cos(dr) * math.cos(la) * math.cos(ha)
    sin_alt = max(-1.0, min(1.0, sin_alt))
    alt = math.asin(sin_alt)
    if alt <= 0:
        continue
    az = math.atan2(-math.cos(dr) * math.sin(ha),
                    math.sin(dr) * math.cos(la) - math.cos(dr) * math.sin(la) * math.cos(ha))
    e, n, u = (math.cos(alt) * math.sin(az), math.cos(alt) * math.cos(az), math.sin(alt))
    zc = u
    if zc <= 0.2:
        continue
    xc = e * stb[0] + n * stb[1]
    yc = e * fwd[0] + n * fwd[1]
    uu = FOCAL * xc / zc + W / 2.0
    vv = FOCAL * yc / zc + H / 2.0
    if 0 <= uu < W and 0 <= vv < H:
        pts.append((uu, vv, vmag, hr))

names = {}
try:
    for ln in open(a.names):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 2:
            names[int(f[0])] = (f[2] if len(f) > 2 and f[2] else f[1])
except OSError:
    print(f"  no {a.names}, labels will be HR numbers")

fig, ax = plt.subplots(figsize=(13, 8.5))
ax.set_facecolor("#05070f")
fig.patch.set_facecolor("#05070f")
for uu, vv, vmag, hr in pts:
    # Linear in magnitude, not in flux: a flux scale makes Vega 250x the area
    # of a 6th-magnitude star and the chart becomes a few white discs.
    size = max(1.6, 1.5 + 2.1 * (6.5 - vmag))
    if vmag <= MATCH_MAG:
        ax.plot(uu, vv, "o", color="white", ms=size, zorder=4)
        ax.plot(uu, vv, "o", mfc="none", mec="lime", ms=size + 5, mew=1.0, zorder=5)
    elif vmag <= MAG_LIMIT:
        ax.plot(uu, vv, "o", color="#cfd8ff", ms=size, zorder=3)
    else:
        ax.plot(uu, vv, "o", color="#6b7699", ms=size, zorder=2)

# Label the matchable stars. Placed with a small offset and skipped where two
# would collide, since the field is dense near the galactic plane.
placed = []
for uu, vv, vmag, hr in sorted(pts, key=lambda p: p[2]):
    if vmag > a.label_mag:
        continue
    if any(abs(uu - px) < 150 and abs(vv - py) < 26 for px, py in placed):
        continue
    ax.annotate(names.get(hr, f"HR {hr}"), (uu, vv), color="#9dffc0",
                fontsize=8.5, xytext=(9, -4), textcoords="offset points",
                zorder=6)
    placed.append((uu, vv))

ax.set_xlim(0, W); ax.set_ylim(H, 0)
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values():
    s.set_color("#3a4260")
place = f" over {a.place}" if a.place else ""
ax.set_title(f"Final celestial fix{place} — error {err_km:.2f} km\n"
             f"{when:%Y-%m-%d %H:%M:%S} UTC   {lat:.3f}, {lon:.3f}   "
             f"heading {math.degrees(yaw):.0f}°   "
             f"{len(pts)} catalogue stars on sensor, {n_logged} matched in flight",
             color="white", fontsize=11)
ax.text(0.5, -0.045,
        "green rings = bright enough for the matcher (V ≤ 4.0)   ·   "
        "pale = detectable (V ≤ 5.0)   ·   grey = catalogue only (V ≤ 6.0)   ·   "
        f"Alvium 1800 U-240, 53.5° × {math.degrees(2*math.atan((H/2)/FOCAL)):.1f}°, zenith-pointing",
        transform=ax.transAxes, ha="center", va="top", color="#8b95b5", fontsize=8)
fig.tight_layout()
fig.savefig(a.out, dpi=140, facecolor=fig.get_facecolor(), bbox_inches="tight")
print(f"wrote {a.out}  ({len(pts)} stars on sensor, fix error {err_km:.2f} km)")
