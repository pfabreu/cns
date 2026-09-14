#!/usr/bin/env python3
"""Project the Yale catalogue onto LenghuSky-8 frames, to LABEL REAL STARS.

    python3 scripts/fetch_lenghusky8.py
    scripts/.venv/bin/python scripts/lenghu_astrometry.py --solve
    scripts/.venv/bin/python scripts/lenghu_astrometry.py --overlay 2018-07-13-01-37-01

WHY. Option (a) in scripts/README.md. Pasting synthetic stars onto these frames
would be actively harmful, because the frames ALREADY contain real unlabelled
stars and every one of them would become a labelled negative. But the site is
published and so is the fisheye calibration, so the real stars can be labelled
instead of invented: project BSC5 through their polynomial and the labels are
DERIVED, not fabricated. It is self-validating -- if the projection is wrong
the stars do not land on the detections and you know at once.

WHAT HAD TO BE REVERSE ENGINEERED. The calibration JSON gives the site, an
odd-order radial polynomial, an optical centre at 4096, and an azimuth offset.
Three things it does NOT give had to be recovered from the pixels:

  1. TIMESTAMPS ARE BEIJING TIME (UTC+8), not UTC. Determined by a rotation
     search: at UTC the catalogue never locks on, at UTC-8 it locks hard.
     Sanity: Lenghu is at 93.3 E, so local solar time is UTC+6.2 -- the site
     runs on Beijing time like everywhere else in China.
  2. THE IMAGE Y AXIS IS FLIPPED relative to a right-handed az convention:
     y = v0 - r sin(psi), not + . Every un-flipped combination scores at the
     noise floor.
  3. THE 4096 -> 512 MAPPING IS A PLAIN RESIZE, scale exactly 0.125, no crop.
     The dataset card says "center-cropped"; for these frames the crop is the
     whole sensor. Confirmed two ways -- the horizon circle radius r(90 deg)
     lands on the measured disc radius at 0.125, and a free fit of scale and
     centre offset moved neither.

AND ONE THING THAT IS STILL NOT UNDERSTOOD. `az_offset_rad` does not enter the
way this code assumes: after subtracting it, the residual rotation still varies
per epoch (96, 203, 97, 27 deg) rather than being constant. So the rotation is
SOLVED EMPIRICALLY PER CALIBRATION EPOCH from clear dark frames, and
`az_offset_rad` is folded in only as a starting point. This works and is
honest, but it means the convention is fitted rather than understood; if you
find the real one, `EPOCH_ROT` becomes unnecessary.

HOW WELL IT WORKS. Best clear frame (2018-07-13-01-37-01), catalogue V < 3.5,
80 strongest detections, match within 2.5 px:

    at the solved rotation      31 matched
    over all 720 rotations      median 1, p99 6
    -> 31x the null, and the global maximum

Per epoch, on clear dark frames only, matched detections against the random
rate:

    2018-05-01   17 frames   11.0x
    2018-09-27   25 frames    5.0x
    2019-04-24    2 frames   11.1x
    2019-07-05   25 frames    1.1x   <- UNSOLVED, do not use
    2023-09-27   25 frames    4.0x

Four epochs of five. The 2019-07-05 epoch does not lock on at any rotation and
its frames must be excluded until someone works out why.

CLEAR FRAMES ONLY, for SOLVING. On a cloudy frame the high-pass is cloud
texture and the "detections" are cloud edges, so the match rate collapses to
chance. Selection uses the annotators' own polygons: `sky` present, `cloud`
absent.

THE LABELS ARE GOOD, AND THE PROOF IS THE BIMODALITY. Per-frame completeness
(fraction of V<3.0 stars above 30 deg that have a detection within 3 px) over
all 211 night frames in the solved epochs is strongly bimodal: 38 frames sit at
90-100% and 53 sit at 0-10%, with the middle sparsely filled. A wrong plate
solution cannot produce that -- it would smear everything toward chance. The
low mode is simply overcast sky with no stars to find. By annotator label:

    clear  frames  n= 99   median completeness 61%
    cloudy frames  n=112   median completeness 22%

That gap is the whole experiment for scripts/unet_star_test.py.

THE ASTROMETRY IS VERIFIED INDEPENDENTLY, because a plate solve that fails
tells you nothing about which half is wrong. Polaris comes out at altitude =
latitude with the right 0.74 deg wobble, and the Sun culminates at 05:48 UTC
against 05:46 predicted from the longitude. See --selftest.

WHAT THIS IS STILL NOT. Their camera, not ours: 1624 arcsec/px against the
flight sensor's 107, 8-bit JPEG rather than radiance, multi-second exposures
with no motion blur. A detector trained on these labels is a transfer check and
the first honest measurement of what cloud costs a real detector. It is not a
flight model. See scripts/README.md.
"""
import argparse
import base64
import datetime
import glob
import json
import math
import os
import pathlib
import re
import sys

import cv2
import numpy as np

HERE = pathlib.Path(__file__).parent
LAT, LON, ALT_M = 38.9586, 93.2681, 4200.0
TZ_HOURS = 8            # filenames are Beijing time
IMG_SCALE = 0.125       # 4096 -> 512, plain resize
DISC_R_PX = 216.0       # measured horizon radius in the 512 frame

# Solved per calibration epoch: (rotation deg, scale, dx px, dy px). The
# rotation stands in for an azimuth convention that has not been worked out;
# scale and the small centre offsets absorb whatever else the published fit
# does not capture. 2019-07-05 is absent because it never locked on.
#
# Refined by coordinate descent on the median residual of V<3.5 stars matched
# to detections on clear dark frames. Median residual after refinement:
#
#     2018-05-01   0.99 px   (n=1071)
#     2018-09-27   2.31 px   (n= 802)   <- the worst; something unmodelled
#     2019-04-24   0.74 px   (n= 100)
#     2023-09-27   1.40 px   (n= 758)
#
# Note scale stayed at exactly 0.125 in every epoch even though it was free.
EPOCH_FIT = {
    "2018-05-01-00-02-44": (280.30, 0.12500, -0.85,  0.00),
    "2018-09-27-19-19-49": (129.76, 0.12500,  0.55, -0.55),
    "2019-04-24-15-39-36": (283.97, 0.12500, -0.55, -0.75),
    "2023-09-27-18-09-48": (144.00, 0.12500, -0.65, -1.40),
}
EPOCH_ROT = {k: v[0] for k, v in EPOCH_FIT.items()}


# --------------------------------------------------------------------------
# astronomy
# --------------------------------------------------------------------------
def julianDay(dt):
    y, m = dt.year, dt.month
    if m <= 2:
        y -= 1
        m += 12
    a = y // 100
    b = 2 - a + a // 4
    d = dt.day + (dt.hour + dt.minute / 60 + dt.second / 3600) / 24
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + b - 1524.5


def altAz(ra_deg, dec_deg, dt_utc):
    """Topocentric altitude/azimuth in radians. Azimuth from north, eastward."""
    n = julianDay(dt_utc) - 2451545.0
    lst = math.radians((((18.697375 + 24.065709824279 * n) % 24) * 15 + LON) % 360)
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    p = math.radians(LAT)
    ha = lst - ra
    alt = np.arcsin(np.clip(np.sin(dec) * math.sin(p) +
                            np.cos(dec) * math.cos(p) * np.cos(ha), -1, 1))
    az = np.arctan2(-np.sin(ha) * np.cos(dec),
                    np.sin(dec) * math.cos(p) - np.cos(dec) * math.sin(p) * np.cos(ha))
    return alt, np.mod(az, 2 * np.pi)


def sunAltDeg(dt_utc):
    """Low-precision solar altitude, good to ~0.01 deg. Enough for day/night."""
    n = julianDay(dt_utc) - 2451545.0
    L = math.radians((280.460 + 0.9856474 * n) % 360)
    g = math.radians((357.528 + 0.9856003 * n) % 360)
    lam = L + math.radians(1.915) * math.sin(g) + math.radians(0.020) * math.sin(2 * g)
    eps = math.radians(23.439)
    ra = math.degrees(math.atan2(math.cos(eps) * math.sin(lam), math.cos(lam)))
    dec = math.degrees(math.asin(math.sin(eps) * math.sin(lam)))
    return math.degrees(altAz(ra, dec, dt_utc)[0])


def loadCatalog(src=None):
    """(ra, dec, vmag, hr) from the generated C++ table. One source of truth."""
    src = src or (HERE.parent / "src/star_catalog_data.cpp")
    rows = re.findall(r"\{(\d+),([-\d.]+),([-\d.]+),([-\d.]+),([-\d.]+),([-\d.]+)\}",
                      pathlib.Path(src).read_text())
    if not rows:
        sys.exit(f"no catalogue rows parsed from {src}")
    a = np.array([[float(r[1]), float(r[2]), float(r[5]), float(r[0])] for r in rows])
    return a


# --------------------------------------------------------------------------
# dataset
# --------------------------------------------------------------------------
def calibrations(cal_dir):
    idx = json.load(open(os.path.join(cal_dir, "calibration_index.json")))
    out = []
    for c in idx["calibrations"]:
        out.append((datetime.datetime.strptime(c["start_file"], "%Y-%m-%d-%H-%M-%S"),
                    c["start_file"],
                    json.load(open(os.path.join(cal_dir, c["file"])))))
    return sorted(out, key=lambda t: t[0])


def slotFor(t, slots):
    j = 0
    for i, (tt, _, _) in enumerate(slots):
        if t >= tt:
            j = i
    return slots[j]


def frameTime(stem):
    return datetime.datetime.strptime(stem, "%Y-%m-%d-%H-%M-%S")


def utcOf(stem):
    return frameTime(stem) - datetime.timedelta(hours=TZ_HOURS)


def loadFrame(data_dir, stem):
    d = json.load(open(os.path.join(data_dir, stem + ".json")))
    im = cv2.imdecode(np.frombuffer(base64.b64decode(d["imageData"]), np.uint8),
                      cv2.IMREAD_GRAYSCALE)
    labels = [s["label"] for s in d["shapes"]]
    return im, labels


def detectSources(im, nmax=80, k=7.0):
    """Local maxima well above a robust noise floor. Not a science detector."""
    a = im.astype(np.float32)
    hp = a - cv2.GaussianBlur(a, (0, 0), 2.5)
    s = 1.4826 * np.median(np.abs(hp - np.median(hp)))
    if s <= 0:
        return np.empty((0, 2)), hp, 1.0
    mx = cv2.dilate(hp, np.ones((5, 5), np.float32))
    ys, xs = np.nonzero((hp >= mx) & (hp > k * s))
    o = np.argsort(-hp[ys, xs])[:nmax]
    return np.stack([xs[o], ys[o]], 1).astype(float), hp, s


def project(cat, dt_utc, fit, epoch_fit, vlim=5.0, alt_min_deg=20.0):
    """Catalogue -> 512x512 pixel coordinates. Returns x, y, vmag.

    `epoch_fit` is an EPOCH_FIT tuple, or a bare rotation in degrees for the
    older call style.
    """
    if np.isscalar(epoch_fit):
        rot_deg, scale, dx, dy = float(epoch_fit), IMG_SCALE, 0.0, 0.0
    else:
        rot_deg, scale, dx, dy = epoch_fit
    b = cat[cat[:, 2] < vlim]
    alt, az = altAz(b[:, 0], b[:, 1], dt_utc)
    k = alt > math.radians(alt_min_deg)
    th = np.pi / 2 - alt[k]
    r = np.zeros_like(th)
    for i, c in enumerate(fit["a_coeffs"]):
        r = r + c * th ** (2 * i + 1)
    psi = az[k] + fit["az_offset_rad"] + math.radians(rot_deg)
    x = (fit["u0"] + r * np.cos(psi)) * scale + dx
    y = (fit["v0"] - r * np.sin(psi)) * scale + dy     # y flipped, see docstring
    return x, y, b[k, 2]


def clearDarkFrames(data_dir):
    """Annotator says sky and not cloud, and the Sun is 18 deg down."""
    out = []
    for p in sorted(glob.glob(os.path.join(data_dir, "*.json"))):
        stem = os.path.basename(p)[:-5]
        labels = [s["label"] for s in json.load(open(p))["shapes"]]
        if "cloud" in labels or "sky" not in labels:
            continue
        if sunAltDeg(utcOf(stem)) < -18:
            out.append(stem)
    return out


def matchCount(P, det, tol=2.5):
    if len(P) == 0 or len(det) == 0:
        return 0
    d = np.linalg.norm(P[:, None, :] - det[None, :, :], axis=2)
    return int((d.min(0) < tol).sum())      # detections explained by a star


# --------------------------------------------------------------------------
def cmdSolve(a, cat, slots):
    frames = clearDarkFrames(a.data)
    by = {}
    for s in frames:
        by.setdefault(slotFor(frameTime(s), slots)[1], []).append(s)
    print(f"  {len(frames)} clear dark frames over {len(by)} calibration epochs\n")
    print(f"  {'epoch':<22}{'n':>4}{'rot':>7}{'matched':>9}{'null med':>10}"
          f"{'null p99':>10}{'ratio':>8}")
    for key in sorted(by):
        grp = by[key][:a.max_frames]
        fit = [c for _, k, c in slots if k == key][0]["fit"]
        dets = []
        for s in grp:
            im, _ = loadFrame(a.data, s)
            d, _, _ = detectSources(im, a.ndet)
            if len(d) >= 10:
                dets.append((d, utcOf(s)))
        if not dets:
            continue
        scores = []
        for rot in np.arange(0, 360, 0.5):
            tot = 0
            for d, t in dets:
                x, y, _ = project(cat, t, fit, rot, a.vlim)
                ok = (x > 2) & (x < 509) & (y > 2) & (y < 509)
                tot += matchCount(np.stack([x[ok], y[ok]], 1), d)
            scores.append(tot)
        scores = np.array(scores)
        i = int(scores.argmax())
        med, p99 = np.median(scores), np.percentile(scores, 99)
        print(f"  {key:<22}{len(dets):>4}{np.arange(0,360,0.5)[i]:>7.1f}"
              f"{scores[i]:>9}{med:>10.1f}{p99:>10.1f}"
              f"{scores[i]/max(med,0.5):>7.1f}x")
    print("\n  Ratios near 1 mean the epoch did not lock on; exclude it.")


def cmdOverlay(a, cat, slots):
    stem = a.overlay
    im, labels = loadFrame(a.data, stem)
    t, key = utcOf(stem), slotFor(frameTime(stem), slots)[1]
    if key not in EPOCH_FIT:
        sys.exit(f"epoch {key} has no solved plate (see EPOCH_FIT)")
    fit = [c for _, k, c in slots if k == key][0]["fit"]
    det, hp, s = detectSources(im, a.ndet)
    x, y, v = project(cat, t, fit, EPOCH_FIT[key], a.vlim)
    ok = (x > 2) & (x < 509) & (y > 2) & (y < 509)
    n = matchCount(np.stack([x[ok], y[ok]], 1), det)
    print(f"  {stem}  epoch {key}  rot {EPOCH_ROT[key]} deg  labels {labels}")
    print(f"  {n}/{len(det)} detections explained by V<{a.vlim} stars "
          f"({int(ok.sum())} in the disc)")
    # Crop to the disc before upscaling. The black surround is half the frame
    # and none of the evidence, and high-passed speckle compresses badly enough
    # that a full 1024 frame costs 1.1 MB of PNG for no extra information.
    cx, cy, rad = 254, 259, 220
    x0, y0 = max(0, cx - rad), max(0, cy - rad)
    x1, y1 = min(512, cx + rad), min(512, cy + rad)
    bg = cv2.cvtColor((np.clip(hp / (a.stretch * s), 0, 1) * 255).astype(np.uint8),
                      cv2.COLOR_GRAY2BGR)[y0:y1, x0:x1]
    sc = 880.0 / bg.shape[1]
    bg = cv2.resize(bg, None, fx=sc, fy=sc, interpolation=cv2.INTER_NEAREST)
    xs, ys, vs = project(cat, t, fit, EPOCH_FIT[key], a.label_mag)
    o = (xs > 2) & (xs < 509) & (ys > 2) & (ys < 509)
    for xi, yi, vi in zip(xs[o], ys[o], vs[o]):
        cv2.circle(bg, (int((xi - x0) * sc), int((yi - y0) * sc)),
                   max(5, int(9 * sc - 1.6 * sc * vi)), (0, 235, 0), 2, cv2.LINE_AA)
    cv2.putText(bg, f"{stem} Beijing = {t.strftime('%Y-%m-%d %H:%M')} UTC",
                (12, 28), cv2.FONT_HERSHEY_SIMPLEX, 0.62, (0, 255, 255), 2, cv2.LINE_AA)
    cv2.putText(bg, f"green = BSC5 V<{a.label_mag} via their fisheye polynomial   "
                    f"{n}/{len(det)} detections matched",
                (12, 54), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 255), 1, cv2.LINE_AA)
    cv2.imwrite(a.out, bg, [cv2.IMWRITE_PNG_COMPRESSION, 9])
    print(f"  wrote {a.out}")


def cmdSelftest():
    print("  Polaris should sit at altitude = latitude, wobbling 0.74 deg:")
    t = datetime.datetime(2019, 1, 1, 18)
    for h in (0, 6, 12, 18):
        al, az = altAz(37.9545, 89.2641, t + datetime.timedelta(hours=h))
        print(f"    +{h:2d}h  alt {math.degrees(al):6.2f} (lat {LAT:.2f})  "
              f"az {math.degrees(az):7.2f}")
    best = max(((sunAltDeg(datetime.datetime(2019, 6, 21) + datetime.timedelta(minutes=m)), m)
                for m in range(0, 1440, 2)))
    print(f"\n  Sun culminates {best[1]//60:02d}:{best[1]%60:02d} UTC at "
          f"alt {best[0]:.2f}")
    print(f"  expected       {int((12-LON/15)):02d}:{int(((12-LON/15)%1)*60):02d} UTC at "
          f"alt {90-LAT+23.44:.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=str(HERE / "data/lenghusky8_seg/data"))
    ap.add_argument("--calib", default=str(HERE / "data/calibrations"))
    ap.add_argument("--solve", action="store_true",
                    help="re-fit the per-epoch rotation from clear dark frames")
    ap.add_argument("--overlay", metavar="STEM",
                    help="render the catalogue over one frame")
    ap.add_argument("--selftest", action="store_true",
                    help="check the astrometry against Polaris and the Sun")
    ap.add_argument("--vlim", type=float, default=5.0, help="magnitude for matching")
    ap.add_argument("--label-mag", type=float, default=4.0, help="magnitude to draw")
    ap.add_argument("--ndet", type=int, default=80)
    ap.add_argument("--max-frames", type=int, default=25)
    ap.add_argument("--out", default="overlay.png")
    ap.add_argument("--stretch", type=float, default=12.0,
                    help="display white point, in sigma of the high-pass. "
                         "Low values drown the stars in JPEG speckle.")
    a = ap.parse_args()

    if a.selftest:
        cmdSelftest()
        return
    if not os.path.isdir(a.calib):
        sys.exit(f"calibration dir {a.calib} missing. Fetch and unzip:\n"
                 f"  curl -L -o calibration.zip https://huggingface.co/datasets/"
                 f"ruiyicheng/LenghuSky-8/resolve/main/data/calibration.zip\n"
                 f"  unzip calibration.zip -d {HERE}/data")
    cat = loadCatalog()
    slots = calibrations(a.calib)
    if a.solve:
        cmdSolve(a, cat, slots)
    elif a.overlay:
        cmdOverlay(a, cat, slots)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
