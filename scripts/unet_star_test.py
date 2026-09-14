#!/usr/bin/env python3
"""UNET TEST: can a UNet find stars through cloud that a threshold misses?

Named a test, not a trainer, because the previous name (train_star_detector)
read like it produced the navigation system's star detector. It does not,
and it never will -- see the caveats below.

    scripts/.venv/bin/python scripts/unet_star_test.py --epochs 120

PROOF OF CONCEPT ONLY. This never goes near the navigation system. It is a
21 arcmin/px all-sky camera on a fixed mount at an observatory, with no motion
smear; the flight sensor is 107 arcsec/px, strapdown, smeared. A detector
trained here would never have seen the streaked star ours exists to find. What
it can establish is narrower and worth having: that a learned detector recovers
stars through real cloud that a threshold on a matched background does not.

WHERE THE LABELS COME FROM. Not from a human, and not invented. The Yale BSC5
catalogue is projected through LenghuSky-8's own published fisheye calibration
for the site and time -- see lenghu_astrometry.py, which had to reverse engineer
the Beijing timestamps, the y flip and the 4096->512 scale before this worked.
Label positions are good to 0.7-2.3 px depending on epoch.

The labels do NOT depend on whether the star is visible. That is the point. A
star behind thick cloud is still labelled at its true position, so recall
measured against these labels is a direct measurement of how much cloud costs a
detector -- and both detectors are scored against the SAME labels, so whatever
the cloud takes away, it takes from both.

WHY THE LABELS ARE TRUSTWORTHY. Per-frame completeness of the classical
detector against them is strongly bimodal -- 38 frames at 90-100%, 53 at 0-10%.
A wrong plate solution smears toward chance and cannot produce two modes. The
low mode is overcast sky.

SPLIT BY NIGHT. Frames hours apart on the same night share cloud state and
nearly the same sky. Holding out whole dates is the same discipline
train_cloud_unet.py uses, for the same reason.

THE COMPARISON. Both detectors produce a score map; peaks are extracted and
matched to catalogue positions within 3 px. Sweeping each detector's threshold
traces a precision-recall curve, so the two are compared at matched precision
rather than at whatever operating point flatters one of them. The classical
baseline is the same high-pass-over-robust-noise detector used to solve the
plate in the first place, which is the honest incumbent here.

EXCLUDED. The 2019-07-05 calibration epoch, which never locked on. Daylight and
twilight frames (Sun above -15 deg). Anything outside the fisheye disc, which
is dome and sky-less and would otherwise dominate both the loss and the
precision.
"""
import argparse
import glob
import json
import math
import os
import pathlib
import random
import sys

import cv2
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

HERE = pathlib.Path(__file__).parent
sys.path.insert(0, str(HERE))
from lenghu_astrometry import (EPOCH_FIT, calibrations, detectSources, frameTime,
                               loadCatalog, loadFrame, project, slotFor, sunAltDeg,
                               utcOf)
from train_cloud_unet import UNet

DISC_C = (254.1, 258.8)
DISC_R = 214.0          # the optical horizon, measured from a median frame

# THE EVALUATION REGION IS NOT THE DISC, and getting this wrong cost a whole
# round of results. Labels are cut at `--alt-min` (25 deg by default), which
# through the fisheye polynomial lands at r = 162 px in every epoch. The disc
# runs to 214. So the annulus 162..214 -- 42% of the disc area, and precisely
# where the dome rim and horizon glow live -- CAN CONTAIN NO LABEL, and every
# peak there is an automatic false positive.
#
# That penalised the classical detector far more than the UNet: the high-pass
# fires hard on the rim, while the network had learned from the training labels
# that the rim is never a star. Scoring inside the labelled region only, the
# gap narrows considerably. Margin is 4 px for the 3 px match tolerance.
EVAL_R = 158.0


# --------------------------------------------------------------------------
def buildIndex(data_dir, calib_dir, vlim, alt_min, sun_max):
    """One entry per usable night frame: image, star pixel positions, labels."""
    slots = calibrations(calib_dir)
    cat = loadCatalog()
    out = []
    for p in sorted(glob.glob(os.path.join(data_dir, "*.json"))):
        stem = os.path.basename(p)[:-5]
        if sunAltDeg(utcOf(stem)) >= sun_max:
            continue
        key = slotFor(frameTime(stem), slots)[1]
        if key not in EPOCH_FIT:
            continue
        fit = [c for _, k, c in slots if k == key][0]["fit"]
        im, labels = loadFrame(data_dir, stem)
        x, y, v = project(cat, utcOf(stem), fit, EPOCH_FIT[key], vlim, alt_min)
        r = np.hypot(x - DISC_C[0], y - DISC_C[1])
        ok = r < EVAL_R - 4
        if ok.sum() < 10:
            continue
        out.append(dict(stem=stem, date=stem[:10], epoch=key, im=im,
                        xy=np.stack([x[ok], y[ok]], 1), vmag=v[ok],
                        cloud="cloud" in labels))
    return out


def discMask(h=512, w=512, radius=DISC_R):
    yy, xx = np.mgrid[0:h, 0:w]
    return (np.hypot(xx - DISC_C[0], yy - DISC_C[1]) < radius).astype(np.float32)


def evalMask(h=512, w=512):
    """Where a label can exist. See EVAL_R -- this is not the disc."""
    return discMask(h, w, EVAL_R)


def heatmap(xy, sigma=1.2, h=512, w=512):
    """Gaussian at each star. Peak 1.0, so BCE reads as 'star here'."""
    m = np.zeros((h, w), np.float32)
    rad = int(math.ceil(3 * sigma))
    for x, y in xy:
        xi, yi = int(round(x)), int(round(y))
        x0, x1 = max(0, xi - rad), min(w, xi + rad + 1)
        y0, y1 = max(0, yi - rad), min(h, yi + rad + 1)
        if x0 >= x1 or y0 >= y1:
            continue
        gy, gx = np.mgrid[y0:y1, x0:x1]
        g = np.exp(-((gx - x) ** 2 + (gy - y) ** 2) / (2 * sigma ** 2))
        m[y0:y1, x0:x1] = np.maximum(m[y0:y1, x0:x1], g)
    return m


class StarSet(Dataset):
    def __init__(self, items, augment, sigma):
        self.items, self.augment, self.sigma = items, augment, sigma
        # The LOSS mask is the evaluation region, not the disc. Training the
        # network on the rim, where no label can exist, teaches it that the
        # rim is never a star -- knowledge the classical baseline has no way
        # to acquire, and which then shows up as a spurious win.
        self.mask = evalMask()

    def __len__(self):
        return len(self.items)

    def __getitem__(self, i):
        it = self.items[i]
        im = it["im"].astype(np.float32) / 255.0
        tg = heatmap(it["xy"], self.sigma)
        mk = self.mask
        if self.augment:
            # The disc is very nearly centred, so quarter turns and flips are
            # valid skies. Without them 150-odd frames is not enough.
            k = random.randint(0, 3)
            if k:
                im, tg, mk = (np.rot90(im, k), np.rot90(tg, k), np.rot90(mk, k))
            if random.random() < 0.5:
                im, tg, mk = im[:, ::-1], tg[:, ::-1], mk[:, ::-1]
        f = lambda a: torch.from_numpy(np.ascontiguousarray(a))[None]
        return f(im), f(tg), f(mk)


# --------------------------------------------------------------------------
def peaks(score, thresh, mask, nmax=600):
    """Local maxima above `thresh`, inside the disc."""
    mx = cv2.dilate(score, np.ones((5, 5), np.float32))
    sel = (score >= mx) & (score > thresh) & (mask > 0)
    ys, xs = np.nonzero(sel)
    if len(xs) > nmax:
        o = np.argsort(-score[ys, xs])[:nmax]
        ys, xs = ys[o], xs[o]
    return np.stack([xs, ys], 1).astype(float)


def prCurve(scores, items, mask, thresholds, tol=3.0):
    """Precision/recall against the catalogue, pooled over frames."""
    out = []
    for t in thresholds:
        tp = fp = fn = 0
        for sc, it in zip(scores, items):
            det = peaks(sc, t, mask)
            gt = it["xy"]
            if len(gt) == 0:
                continue
            if len(det) == 0:
                fn += len(gt)
                continue
            d = np.linalg.norm(gt[:, None, :] - det[None, :, :], axis=2)
            hit = (d.min(1) < tol)
            tp += int(hit.sum())
            fn += int((~hit).sum())
            fp += int((d.min(0) >= tol).sum())
        p = tp / max(tp + fp, 1)
        r = tp / max(tp + fn, 1)
        out.append((t, p, r, 2 * p * r / max(p + r, 1e-9)))
    return out


def classicalScore(im):
    """The incumbent: high-pass over a robust local noise estimate, in sigma."""
    a = im.astype(np.float32)
    hp = a - cv2.GaussianBlur(a, (0, 0), 2.5)
    s = 1.4826 * np.median(np.abs(hp - np.median(hp)))
    return hp / max(s, 1e-6)


def bestAtPrecision(curve, pmin):
    ok = [c for c in curve if c[1] >= pmin]
    return max(ok, key=lambda c: c[2]) if ok else None


# --------------------------------------------------------------------------
def recallAtPrecision(curve, pmin):
    """Best recall among operating points that reach `pmin` precision."""
    ok = [c for c in curve if c[1] >= pmin]
    return max(ok, key=lambda c: c[2]) if ok else None


def report(model, te, mask, dev, a):
    """Head to head on held-out nights, UNet against the classical detector."""
    amp = dev.type == "cuda"
    unet_scores = []
    with torch.no_grad():
        for it in te:
            x = torch.from_numpy(it["im"].astype(np.float32) / 255.0)[None, None].to(dev)
            unet_scores.append(torch.sigmoid(model(x))[0, 0].float().cpu().numpy())
    cls_scores = [classicalScore(it["im"]) for it in te]
    ut = np.concatenate([np.linspace(0.02, 0.95, 32), [0.97, 0.98, 0.99]])
    # The classical detector needs a WIDE sweep or it never reaches the
    # precision the UNet does, and the comparison silently flatters the UNet.
    # Its precision genuinely CAPS near 0.45 -- past ~30 sigma the survivors
    # are saturated blobs rather than the brightest stars, so precision falls.
    ct = np.concatenate([np.linspace(2.0, 30.0, 29), np.linspace(35.0, 200.0, 20)])

    curves = {}
    print(f"\n  held-out: {len(te)} frames on {len({it['date'] for it in te})} nights")
    for name, group in (("CLEAR", [i for i, it in enumerate(te) if not it["cloud"]]),
                        ("CLOUDY", [i for i, it in enumerate(te) if it["cloud"]]),
                        ("ALL", list(range(len(te))))):
        if not group:
            continue
        gi = [te[i] for i in group]
        cu = prCurve([unet_scores[i] for i in group], gi, mask, ut)
        cc = prCurve([cls_scores[i] for i in group], gi, mask, ct)
        curves[name] = (cu, cc)
        pmax_c = max(c[1] for c in cc)
        print(f"\n  --- {name}: {len(group)} frames, "
              f"{sum(len(it['xy']) for it in gi)} catalogue stars ---")
        print(f"      classical precision tops out at {pmax_c:.3f}; "
              f"UNet reaches {max(c[1] for c in cu):.3f}")
        print(f"      {'at precision':<14}{'UNet recall':>12}{'classical':>11}{'gain':>8}")
        for pmin in (0.15, 0.20, 0.25, 0.30, 0.40):
            if pmin > pmax_c:
                continue
            bu, bc = recallAtPrecision(cu, pmin), recallAtPrecision(cc, pmin)
            if bu and bc and bc[2] > 0:
                print(f"      >= {pmin:<11.2f}{bu[2]:>12.3f}{bc[2]:>11.3f}"
                      f"{bu[2]/bc[2]:>7.1f}x")
        fu, fc = max(c[3] for c in cu), max(c[3] for c in cc)
        print(f"      best F1        {fu:>12.3f}{fc:>11.3f}{fu/max(fc,1e-9):>7.2f}x")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
        for ax, name in zip(axes, ("CLEAR", "CLOUDY")):
            if name not in curves:
                continue
            cu, cc = curves[name]
            ax.plot([c[2] for c in cu], [c[1] for c in cu], "-o", ms=3,
                    label="UNet", color="#1f77b4")
            ax.plot([c[2] for c in cc], [c[1] for c in cc], "-s", ms=3,
                    label="classical high-pass", color="#d62728")
            ax.set_xlabel("recall against BSC5"); ax.set_ylabel("precision")
            ax.set_title(f"{name} nights"); ax.grid(alpha=.3)
            ax.set_xlim(0, .75); ax.set_ylim(0, 1.02); ax.legend()
        fig.suptitle("Star detection on real LenghuSky-8 night sky, held-out nights "
                     "(labels: Yale BSC5 through their fisheye calibration)")
        fig.tight_layout()
        fig.savefig(a.plot, dpi=130)
        print(f"\n  wrote {a.plot}")



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=str(HERE / "data/lenghusky8_seg/data"))
    ap.add_argument("--calib", default=str(HERE / "data/calibrations"))
    ap.add_argument("--vlim", type=float, default=4.5)
    ap.add_argument("--alt-min", type=float, default=25.0)
    ap.add_argument("--sun-max", type=float, default=-15.0)
    ap.add_argument("--sigma", type=float, default=1.2)
    ap.add_argument("--base", type=int, default=16)
    ap.add_argument("--epochs", type=int, default=120)
    ap.add_argument("--batch", type=int, default=4)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--patience", type=int, default=25)
    ap.add_argument("--pos-weight", type=float, default=40.0)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--cpu", action="store_true")
    ap.add_argument("--out", default=str(HERE / "runs"))
    ap.add_argument("--eval-only", metavar="CKPT",
                    help="skip training, score this checkpoint")
    ap.add_argument("--plot", metavar="PNG", default=None,
                    help="write the precision-recall figure here")
    a = ap.parse_args()

    if not a.cpu and not torch.cuda.is_available():
        sys.exit("CUDA not available; pass --cpu if you really mean it")
    dev = torch.device("cpu" if a.cpu else "cuda")
    amp = dev.type == "cuda"
    random.seed(a.seed); np.random.seed(a.seed); torch.manual_seed(a.seed)

    items = buildIndex(a.data, a.calib, a.vlim, a.alt_min, a.sun_max)
    if not items:
        sys.exit("no usable frames; run fetch_lenghusky8.py and unzip calibration.zip")
    dates = sorted({it["date"] for it in items})
    random.Random(a.seed).shuffle(dates)
    ntest = max(1, int(0.25 * len(dates)))
    test_d, val_d = set(dates[:ntest]), set(dates[ntest:ntest + max(1, ntest // 2)])
    tr = [it for it in items if it["date"] not in test_d | val_d]
    va = [it for it in items if it["date"] in val_d]
    te = [it for it in items if it["date"] in test_d]
    nc = lambda g: sum(it["cloud"] for it in g)
    print(f"  device {dev}" + (f"  {torch.cuda.get_device_name(0)}" if amp else ""))
    print(f"  {len(items)} night frames / {len(dates)} nights, "
          f"catalogue V<{a.vlim}, alt>{a.alt_min:.0f} deg")
    print(f"  train {len(tr):3d} ({nc(tr)} cloudy)  val {len(va):3d} ({nc(va)} cloudy)  "
          f"test {len(te):3d} ({nc(te)} cloudy)")
    print(f"  stars per frame: median {int(np.median([len(it['xy']) for it in items]))}")

    mask = discMask()
    if a.eval_only:
        ck = torch.load(a.eval_only, map_location="cpu", weights_only=False)
        model = UNet(1, 1, ck["base"]).to(dev)
        model.load_state_dict(ck["model"])
        model.eval()
        report(model, te, evalMask(), dev, a)
        return
    ltr = DataLoader(StarSet(tr, True, a.sigma), batch_size=a.batch, shuffle=True,
                     num_workers=4, pin_memory=amp, drop_last=len(tr) > a.batch,
                     persistent_workers=True)
    lva = DataLoader(StarSet(va, False, a.sigma), batch_size=a.batch, num_workers=2)

    model = UNet(1, 1, a.base).to(dev)
    print(f"  UNet base={a.base} 1ch->1ch  "
          f"{sum(p.numel() for p in model.parameters())/1e6:.2f} M params")
    opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=a.epochs)
    scaler = torch.amp.GradScaler("cuda", enabled=amp)
    pw = torch.tensor([a.pos_weight], device=dev)

    def loss_fn(logit, tg, mk):
        # Stars are ~0.05% of pixels, so plain BCE predicts empty sky forever.
        bce = F.binary_cross_entropy_with_logits(logit, tg, pos_weight=pw,
                                                 reduction="none")
        bce = (bce * mk).sum() / mk.sum().clamp(min=1)
        p = torch.sigmoid(logit) * mk
        t = tg * mk
        dice = 1 - (2 * (p * t).sum() + 1) / ((p * p).sum() + (t * t).sum() + 1)
        return bce + dice

    best, best_state, bad = float("inf"), None, 0
    print(f"\n  {'epoch':>5}{'train':>10}{'val':>10}")
    for ep in range(a.epochs):
        model.train(); tot = n = 0
        for x, t, m in ltr:
            x, t, m = x.to(dev), t.to(dev), m.to(dev)
            opt.zero_grad(set_to_none=True)
            with torch.autocast("cuda", torch.float16, enabled=amp):
                l = loss_fn(model(x), t, m)
            scaler.scale(l).backward(); scaler.step(opt); scaler.update()
            tot += float(l.detach()); n += 1
        sched.step()
        model.eval(); vt = vn = 0
        with torch.no_grad():
            for x, t, m in lva:
                x, t, m = x.to(dev), t.to(dev), m.to(dev)
                with torch.autocast("cuda", torch.float16, enabled=amp):
                    vt += float(loss_fn(model(x), t, m)); vn += 1
        vl = vt / max(vn, 1)
        if (ep + 1) % 5 == 0 or ep == 0:
            print(f"  {ep+1:>5}{tot/max(n,1):>10.4f}{vl:>10.4f}")
        if vl < best - 1e-5:
            best, bad = vl, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
            if bad >= a.patience:
                print(f"  early stop at epoch {ep+1}"); break
    model.load_state_dict(best_state); model.eval()

    report(model, te, evalMask(), dev, a)

    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    torch.save({"model": model.state_dict(), "args": vars(a), "base": a.base},
               out / "star_detector.pt")
    print(f"\n  wrote {out}/star_detector.pt")


if __name__ == "__main__":
    main()
