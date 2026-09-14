#!/usr/bin/env python3
"""Trajectory over a satellite basemap.

    python3 tools/plot/plot_map.py a.csv out.png --compare b.csv \
        --labels "no horizon" "horizon" --mark 33.0734,-16.3500,"Porto Santo"

WHY A SEPARATE TOOL. live_view.py plots metres north/east from the first truth
point, which is the right frame for reading error and the wrong one for putting
a map underneath. Tiles are Web Mercator, so everything here works in EPSG:3857
metres and the tracks are projected into it.

Imagery: Esri World Imagery. Tiles are cached under --cache so a re-run does not
re-fetch, and the attribution Esri requires is drawn on the figure.
"""
import argparse, csv, io, math, os, urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

R = 6378137.0
TILE = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"


def merc(lat, lon):
    """lat/lon degrees -> EPSG:3857 metres."""
    x = R * math.radians(lon)
    y = R * math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
    return x, y


def tile_xy(lat, lon, z):
    n = 2 ** z
    x = (lon + 180.0) / 360.0 * n
    lr = math.radians(lat)
    y = (1.0 - math.log(math.tan(lr) + 1 / math.cos(lr)) / math.pi) / 2.0 * n
    return x, y


def tile_bounds_merc(x, y, z):
    """Web Mercator extent of one tile."""
    span = 2 * math.pi * R / (2 ** z)
    x0 = -math.pi * R + x * span
    y1 = math.pi * R - y * span
    return x0, x0 + span, y1 - span, y1


def fetch(z, x, y, cache):
    p = os.path.join(cache, f"{z}_{x}_{y}.jpg")
    if not os.path.exists(p):
        url = TILE.format(z=z, x=x, y=y)
        req = urllib.request.Request(url, headers={"User-Agent": "cns-plot/1.0"})
        with urllib.request.urlopen(req, timeout=20) as r:
            data = r.read()
        with open(p, "wb") as f:
            f.write(data)
    return Image.open(p).convert("RGB")


def basemap(ax, lat0, lat1, lon0, lon1, z, cache):
    x0f, y1f = tile_xy(lat1, lon0, z)      # north-west
    x1f, y0f = tile_xy(lat0, lon1, z)      # south-east
    xs = range(int(math.floor(x0f)), int(math.floor(x1f)) + 1)
    ys = range(int(math.floor(y1f)), int(math.floor(y0f)) + 1)
    n = 0
    for tx in xs:
        for ty in ys:
            try:
                img = fetch(z, tx, ty, cache)
            except Exception as e:
                print(f"  tile {z}/{tx}/{ty} failed: {e}")
                continue
            l, r_, b, t = tile_bounds_merc(tx, ty, z)
            ax.imshow(img, extent=(l, r_, b, t), origin="upper",
                      interpolation="bilinear", zorder=0)
            n += 1
    print(f"  {n} tiles at zoom {z}")


def track(path):
    truth, dr, fixes, dro = [], [], [], []
    for r in csv.DictReader(open(path)):
        try:
            truth.append(merc(float(r["lat_true"]), float(r["lon_true"])))
        except (ValueError, KeyError):
            pass
        try:
            if float(r["dr_lat"]):
                dr.append(merc(float(r["dr_lat"]), float(r["dr_lon"])))
        except (ValueError, KeyError):
            pass
        try:
            if float(r["dr_open_lat"]):
                dro.append(merc(float(r["dr_open_lat"]), float(r["dr_open_lon"])))
        except (ValueError, KeyError):
            pass
        try:
            if float(r.get("orbit_lat", 0)):
                fixes.append(merc(float(r["orbit_lat"]), float(r["orbit_lon"])))
        except (ValueError, KeyError):
            pass
    return truth, dr, fixes, dro


ap = argparse.ArgumentParser()
ap.add_argument("csv"); ap.add_argument("out")
ap.add_argument("--compare", default=None)
ap.add_argument("--labels", nargs=2, default=["run A", "run B"])
ap.add_argument("--zoom", type=int, default=8)
ap.add_argument("--mark", action="append", default=[],
                help="lat,lon,label  (repeatable)")
ap.add_argument("--title", default="")
ap.add_argument("--show-dr", action="store_true",
                help="draw unaided dead reckoning too. It ends hundreds of km "
                     "away, so it sets the map extent and the crossing itself "
                     "shrinks -- which is the point of the picture.")
ap.add_argument("--legend-loc", default="upper left")
ap.add_argument("--cache", default="/tmp/cns_tiles")
a = ap.parse_args()
os.makedirs(a.cache, exist_ok=True)

t1, d1, f1, o1 = track(a.csv)
allpts = list(t1) + list(d1)
if a.show_dr:
    allpts += list(o1)
t2 = d2 = f2 = None
if a.compare:
    t2, d2, f2, _o2 = track(a.compare)
    allpts += list(d2)

xs = [p[0] for p in allpts]; ys = [p[1] for p in allpts]
padx = 0.06 * (max(xs) - min(xs)); pady = 0.06 * (max(ys) - min(ys))
xmin, xmax = min(xs) - padx, max(xs) + padx
ymin, ymax = min(ys) - pady, max(ys) + pady

# invert-merc the corners to pick tiles
def inv(x, y):
    lon = math.degrees(x / R)
    lat = math.degrees(2 * math.atan(math.exp(y / R)) - math.pi / 2)
    return lat, lon
# Figure shape follows the route. A 152 km north-south corridor in a 12x10
# frame becomes a thin strip with the legend sitting on the track.
aspect = (xmax - xmin) / (ymax - ymin)
if aspect >= 1.0:
    figsize = (13.0, max(4.5, min(13.0, 13.0 / aspect)))
else:
    figsize = (max(5.0, min(13.0, 13.0 * aspect)), 13.0)
fig, ax = plt.subplots(figsize=figsize)

la0, lo0 = inv(xmin, ymin); la1, lo1 = inv(xmax, ymax)
basemap(ax, la0, la1, lo0, lo1, a.zoom, a.cache)

ax.plot(*zip(*t1), "-", color="white", lw=2.4, zorder=3, label="truth", alpha=0.9)
ax.plot(*zip(*d1), "-", color="tab:blue", lw=1.6, zorder=4, label=f"{a.labels[0]}: DR + fixes")
if f1:
    ax.plot([p[0] for p in f1], [p[1] for p in f1], "*", color="tab:blue",
            ms=10, ls="none", mec="k", mew=0.4, zorder=5,
            label=f"{a.labels[0]}: fixes ({len(f1)})")
if d2:
    ax.plot(*zip(*d2), "-", color="tab:green", lw=1.6, zorder=4, label=f"{a.labels[1]}: DR + fixes")
    if f2:
        ax.plot([p[0] for p in f2], [p[1] for p in f2], "*", color="tab:green",
                ms=10, ls="none", mec="k", mew=0.4, zorder=5,
                label=f"{a.labels[1]}: fixes ({len(f2)})")

if a.show_dr and o1:
    ax.plot(*zip(*o1), "--", color="#ff5b5b", lw=1.8, zorder=3.5, alpha=0.95,
            label="dead reckoning ONLY (no fixes)")
    ax.plot(o1[-1][0], o1[-1][1], "s", color="#ff5b5b", ms=9, zorder=6,
            mec="k", mew=0.5)

for m in a.mark:
    plat, plon, lbl = m.split(",", 2)
    mx, my = merc(float(plat), float(plon))
    ax.plot(mx, my, "o", mfc="none", mec="yellow", ms=18, mew=2.0, zorder=6)
    ax.annotate(lbl, (mx, my), color="yellow", fontsize=11, weight="bold",
                xytext=(14, 8), textcoords="offset points", zorder=6)

ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
ax.set_xticks([]); ax.set_yticks([])
if a.title:
    ax.set_title(a.title)
ax.legend(loc=a.legend_loc, fontsize=9, framealpha=0.85)
ax.text(0.995, 0.006, "Imagery © Esri, Maxar, Earthstar Geographics",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
        color="white", alpha=0.85, zorder=7)
fig.tight_layout()
fig.savefig(a.out, dpi=140, bbox_inches="tight")
print("wrote", a.out)
