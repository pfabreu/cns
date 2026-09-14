#!/usr/bin/env python3
"""Train a small UNet on LenghuSky-8 to segment sky / cloud / contamination.

    python3 scripts/fetch_lenghusky8.py
    scripts/.venv/bin/python scripts/train_cloud_unet.py --split group
    scripts/.venv/bin/python scripts/train_cloud_unet.py --split random   # their protocol

WHAT THIS IS. Real night sky, real cloud, hand-labelled. Everything else this
project knows about cloud is synthetic and therefore circular (scripts/README.md).
This is the first non-circular thing in that column, and it is still not the
flight sensor -- see "the honest gap" at the bottom of this docstring.

THE LABELS. labelme JSON, one per image, with the 512x512 JPEG embedded as
base64 in `imageData`, so the 22 MB zip is self-contained and the 20 GB of
image tars are not needed.

    sky 0, cloud 1, contamination 2, UNLABELLED 3 = ignore

`contamination` is the moon, satellite trails, the dome, frost on the window.

THE POLYGONS ARE NOT EXHAUSTIVE. This is the one thing to get right. Most
images carry two or three polygons and leave the rest of the frame unlabelled,
so an unlabelled pixel means "the annotator did not say", NOT "sky". Filling
the background with class 0 would train the model to call every ambiguous
region sky and would score it against a label nobody wrote. Upstream's
`train_segmentation_Unet.py` initialises the mask to 3 and fills polygons over
it; CrossEntropyLoss(ignore_index=3) then drops those pixels from loss AND from
the metrics. Reproduced exactly here.

SPLIT BY NIGHT, NOT BY FILE -- and this is where we depart from upstream.
Their `bootstrap_segmentation.py` shuffles the 1,111 files and cuts 80/10/10.
But the images are timestamped and they cluster: 249 distinct dates, a median
of 3 images per date, and 226 consecutive pairs less than an hour apart. Cloud
does not reorganise itself in an hour, so a random cut puts near-duplicate
frames on both sides of the wall and reports a number inflated by leakage.
`--split group` holds out whole dates. `--split random` reproduces theirs, so
the gap between the two is measurable rather than asserted -- and that gap is
the actual result of running this.

This is the same lesson already written into `train_unet.py` for the synthetic
data ("SPLIT BY TRAJECTORY, NOT BY FRAME"), which is some comfort that it is a
real effect and not special pleading.

THE MODEL. 4 down, 4 up, 16 base channels, 1.96 M parameters -- the depth and
width named in NOTES-private.md. Upstream's baseline starts at 64 channels, so
this is ~16x fewer per layer and ~16x fewer parameters. The point is not to win their benchmark; it is to find out what a
model that fits inside a 100 ms frame budget can actually do on real cloud.

THE HONEST GAP, restated because it is easy to lose. LenghuSky-8 is an all-sky
fisheye on a fixed mount taking multi-second exposures at a premier
astronomical site. The flight sensor is a 53 deg camera on a moving aircraft at
0.1 s with motion blur. Different field, different PSF, different noise, no
smear. A good number here is evidence that the cloud-appearance half of the
problem is learnable from real data. It is not a flight result and must not be
quoted as one.
"""
import argparse
import base64
import collections
import datetime
import json
import pathlib
import random
import sys
import time

import cv2
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

CLASSES = ["sky", "cloud", "contamination"]
IGNORE = 3

# The annotators were inconsistent. 7 of 2,846 polygons carry a misspelling;
# dropping them would silently delete labelled area, so map them.
LABEL2ID = {
    "sky": 0,
    "cloud": 1, "clode": 1,
    "contamination": 2, "contination": 2, "contimination": 2,
}


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def loadJsonImageMask(path, gray):
    """Decode one labelme file to (image HWC uint8, mask HW uint8).

    Rasterisation matches upstream exactly: mask starts at IGNORE, polygons are
    filled in file order, so a later polygon overwrites an earlier one where
    they overlap. That is the annotator's own layering and we keep it.
    """
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    buf = np.frombuffer(base64.b64decode(d["imageData"]), np.uint8)
    img = cv2.imdecode(buf, cv2.IMREAD_GRAYSCALE if gray else cv2.IMREAD_COLOR)
    if not gray:
        img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    if img.ndim == 2:
        img = img[:, :, None]

    h, w = d["imageHeight"], d["imageWidth"]
    mask = IGNORE * np.ones((h, w), np.uint8)
    for s in d["shapes"]:
        cid = LABEL2ID.get(s["label"])
        if cid is None:
            continue
        pts = np.round(np.array(s["points"], np.float32)).astype(np.int32)
        pts[:, 0] = np.clip(pts[:, 0], 0, w - 1)
        pts[:, 1] = np.clip(pts[:, 1], 0, h - 1)
        # 4 shapes are tagged `linestrip` rather than `polygon`; upstream fills
        # anything with >=3 points and so do we.
        if len(pts) >= 3:
            cv2.fillPoly(mask, [pts], cid)
    return img, mask


class CloudSet(Dataset):
    """Decodes once into RAM: 1,111 x 512 x 512 is ~1.1 GB at 3 channels."""

    def __init__(self, files, gray, augment):
        self.augment = augment
        self.imgs, self.masks = [], []
        for p in files:
            im, mk = loadJsonImageMask(p, gray)
            self.imgs.append(im)
            self.masks.append(mk)

    def __len__(self):
        return len(self.imgs)

    def __getitem__(self, i):
        im, mk = self.imgs[i], self.masks[i]
        if self.augment:
            # An all-sky frame is zenith-centred, so flips and quarter turns
            # are all physically plausible skies. With 889 training images this
            # is not optional.
            k = random.randint(0, 3)
            if k:
                im, mk = np.rot90(im, k, (0, 1)), np.rot90(mk, k, (0, 1))
            if random.random() < 0.5:
                im, mk = im[:, ::-1], mk[:, ::-1]
        x = torch.from_numpy(np.ascontiguousarray(im.transpose(2, 0, 1)))
        y = torch.from_numpy(np.ascontiguousarray(mk))
        return x.float().div_(255.0), y.long()


def splitFiles(files, mode, seed, val_ratio=0.1, test_ratio=0.1):
    """Return (train, val, test).

    mode 'random' shuffles files, which is upstream's protocol.
    mode 'group'  shuffles DATES and keeps every image of a date together.
    """
    rng = random.Random(seed)
    if mode == "random":
        fs = sorted(files)
        rng.shuffle(fs)
        n = len(fs)
        ntest = int(n * test_ratio)
        nval = int(n * val_ratio)
        return fs[ntest + nval:], fs[ntest:ntest + nval], fs[:ntest]

    by_date = collections.defaultdict(list)
    for p in sorted(files):
        by_date[p.stem[:10]].append(p)          # yyyy-mm-dd
    dates = sorted(by_date)
    rng.shuffle(dates)
    n = len(dates)
    ntest = max(1, int(n * test_ratio))
    nval = max(1, int(n * val_ratio))
    pick = lambda ds: [p for d in ds for p in by_date[d]]
    return (pick(dates[ntest + nval:]), pick(dates[ntest:ntest + nval]),
            pick(dates[:ntest]))


# ---------------------------------------------------------------------------
# model -- 4 down, 4 up, `base` channels at full resolution
# ---------------------------------------------------------------------------
class DoubleConv(nn.Module):
    def __init__(self, cin, cout):
        super().__init__()
        self.f = nn.Sequential(
            nn.Conv2d(cin, cout, 3, padding=1, bias=False),
            nn.BatchNorm2d(cout), nn.ReLU(inplace=True),
            nn.Conv2d(cout, cout, 3, padding=1, bias=False),
            nn.BatchNorm2d(cout), nn.ReLU(inplace=True))

    def forward(self, x):
        return self.f(x)


class UNet(nn.Module):
    def __init__(self, cin=3, ncls=3, base=16):
        super().__init__()
        b = base
        self.inc = DoubleConv(cin, b)
        self.d1, self.d2 = DoubleConv(b, b * 2), DoubleConv(b * 2, b * 4)
        self.d3, self.d4 = DoubleConv(b * 4, b * 8), DoubleConv(b * 8, b * 16)
        self.u4 = DoubleConv(b * 16 + b * 8, b * 8)
        self.u3 = DoubleConv(b * 8 + b * 4, b * 4)
        self.u2 = DoubleConv(b * 4 + b * 2, b * 2)
        self.u1 = DoubleConv(b * 2 + b, b)
        self.out = nn.Conv2d(b, ncls, 1)
        self.pool = nn.MaxPool2d(2)

    def _up(self, x, skip):
        # `size=` rather than `scale_factor=` so odd spatial dims survive the
        # round trip. The flight frame decimated by 4 is 484x304, not a power
        # of two, and a scale_factor graph silently mismatches by a pixel.
        x = F.interpolate(x, size=skip.shape[-2:], mode="bilinear",
                          align_corners=False)
        return torch.cat([skip, x], 1)

    def forward(self, x):
        x1 = self.inc(x)
        x2 = self.d1(self.pool(x1))
        x3 = self.d2(self.pool(x2))
        x4 = self.d3(self.pool(x3))
        x5 = self.d4(self.pool(x4))
        y = self.u4(self._up(x5, x4))
        y = self.u3(self._up(y, x3))
        y = self.u2(self._up(y, x2))
        y = self.u1(self._up(y, x1))
        return self.out(y)


# ---------------------------------------------------------------------------
# metrics -- confusion matrix over VALID pixels only
# ---------------------------------------------------------------------------
def confusion(model, loader, device, amp):
    model.eval()
    cm = np.zeros((3, 3), np.int64)
    loss_sum, nb = 0.0, 0
    crit = nn.CrossEntropyLoss(ignore_index=IGNORE)
    with torch.no_grad():
        for x, y in loader:
            x, y = x.to(device, non_blocking=True), y.to(device, non_blocking=True)
            with torch.autocast("cuda", torch.float16, enabled=amp):
                logits = model(x)
                loss_sum += float(crit(logits, y))
            nb += 1
            p = logits.argmax(1).view(-1)
            t = y.view(-1)
            keep = t != IGNORE
            if keep.any():
                idx = (t[keep] * 3 + p[keep]).cpu().numpy()
                cm += np.bincount(idx, minlength=9).reshape(3, 3)
    return cm, loss_sum / max(1, nb)


def reportCm(cm, title):
    tp = np.diag(cm).astype(float)
    sup = cm.sum(1).astype(float)          # true per class
    pred = cm.sum(0).astype(float)
    prec = np.divide(tp, pred, out=np.zeros(3), where=pred > 0)
    rec = np.divide(tp, sup, out=np.zeros(3), where=sup > 0)
    den = prec + rec
    f1 = np.divide(2 * prec * rec, den, out=np.zeros(3), where=den > 0)
    union = sup + pred - tp
    iou = np.divide(tp, union, out=np.zeros(3), where=union > 0)
    acc = tp.sum() / max(1.0, cm.sum())
    print(f"\n  {title}")
    print(f"    {'class':<14}{'prec':>8}{'recall':>8}{'F1':>8}{'IoU':>8}"
          f"{'pixels':>12}")
    for i, c in enumerate(CLASSES):
        print(f"    {c:<14}{prec[i]:8.3f}{rec[i]:8.3f}{f1[i]:8.3f}{iou[i]:8.3f}"
              f"{int(sup[i]):12d}")
    print(f"    {'':<14}{'':>8}{'':>8}{f1.mean():8.3f}{iou.mean():8.3f}"
          f"   <- macro")
    print(f"    overall pixel accuracy {acc:.4f}   ({int(cm.sum())} labelled px)")
    return acc, f1, iou


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    here = pathlib.Path(__file__).parent
    ap.add_argument("--data", default=str(here / "data/lenghusky8_seg/data"))
    ap.add_argument("--split", choices=["group", "random"], default="group",
                    help="group = hold out whole dates (honest); "
                         "random = upstream's file-level shuffle (leaky)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--base", type=int, default=16, help="UNet base channels")
    ap.add_argument("--gray", action="store_true",
                    help="1-channel input, so the export is drop-in for "
                         "PriorNet, which feeds a mono frame")
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--patience", type=int, default=20)
    ap.add_argument("--no-augment", action="store_true")
    ap.add_argument("--cpu", action="store_true", help="force CPU (slow)")
    ap.add_argument("--out", default=str(here / "runs"))
    ap.add_argument("--tag", default="")
    a = ap.parse_args()

    if not a.cpu and not torch.cuda.is_available():
        sys.exit("CUDA not available. Re-run with --cpu if you really mean it.")
    device = torch.device("cpu" if a.cpu else "cuda")
    amp = device.type == "cuda"

    random.seed(a.seed)
    np.random.seed(a.seed)
    torch.manual_seed(a.seed)

    files = sorted(pathlib.Path(a.data).glob("*.json"))
    if not files:
        sys.exit(f"no json in {a.data}; run scripts/fetch_lenghusky8.py first")
    tr, va, te = splitFiles(files, a.split, a.seed)

    print(f"  device  {device}" + (f"  {torch.cuda.get_device_name(0)}"
                                   if device.type == "cuda" else ""))
    print(f"  split   {a.split}  seed {a.seed}")
    ndate = lambda fs: len({p.stem[:10] for p in fs})
    print(f"  train {len(tr):4d} img / {ndate(tr):3d} dates")
    print(f"  val   {len(va):4d} img / {ndate(va):3d} dates")
    print(f"  test  {len(te):4d} img / {ndate(te):3d} dates")
    if a.split == "group":
        overlap = ({p.stem[:10] for p in tr} &
                   ({p.stem[:10] for p in va} | {p.stem[:10] for p in te}))
        assert not overlap, f"date leaked across the split: {overlap}"
        print("  no date appears on both sides of the split")

    t0 = time.time()
    dtr = CloudSet(tr, a.gray, not a.no_augment)
    dva = CloudSet(va, a.gray, False)
    dte = CloudSet(te, a.gray, False)
    print(f"  decoded {len(files)} images in {time.time()-t0:.1f} s")

    ltr = DataLoader(dtr, batch_size=a.batch, shuffle=True, num_workers=4,
                     pin_memory=amp, drop_last=True, persistent_workers=True)
    lva = DataLoader(dva, batch_size=a.batch, num_workers=2, pin_memory=amp)
    lte = DataLoader(dte, batch_size=a.batch, num_workers=2, pin_memory=amp)

    cin = 1 if a.gray else 3
    model = UNet(cin, 3, a.base).to(device).to(memory_format=torch.channels_last)
    npar = sum(p.numel() for p in model.parameters())
    print(f"  model   UNet base={a.base} in={cin}ch  {npar/1e6:.2f} M params")

    crit = nn.CrossEntropyLoss(ignore_index=IGNORE)
    opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=a.epochs)
    scaler = torch.amp.GradScaler("cuda", enabled=amp)

    best, best_state, bad, best_ep = float("inf"), None, 0, -1
    print(f"\n  {'epoch':>5}{'train':>9}{'val':>9}{'val acc':>9}"
          f"{'mIoU':>8}{'s':>7}")
    for ep in range(a.epochs):
        te0 = time.time()
        model.train()
        tot, nb = 0.0, 0
        for x, y in ltr:
            x = x.to(device, non_blocking=True).to(memory_format=torch.channels_last)
            y = y.to(device, non_blocking=True)
            opt.zero_grad(set_to_none=True)
            with torch.autocast("cuda", torch.float16, enabled=amp):
                loss = crit(model(x), y)
            scaler.scale(loss).backward()
            scaler.step(opt)
            scaler.update()
            tot += float(loss.detach())
            nb += 1
        sched.step()
        cm, vloss = confusion(model, lva, device, amp)
        tp = np.diag(cm).astype(float)
        vacc = tp.sum() / max(1.0, cm.sum())
        union = cm.sum(1) + cm.sum(0) - tp
        miou = float(np.divide(tp, union, out=np.zeros(3),
                               where=union > 0).mean())
        print(f"  {ep+1:5d}{tot/max(1,nb):9.4f}{vloss:9.4f}{vacc:9.4f}"
              f"{miou:8.3f}{time.time()-te0:7.1f}")
        if vloss < best - 1e-5:
            best, best_state, bad, best_ep = vloss, {
                k: v.detach().clone() for k, v in model.state_dict().items()}, 0, ep
        else:
            bad += 1
            if bad >= a.patience:
                print(f"  early stop (no val improvement for {a.patience})")
                break

    model.load_state_dict(best_state)
    print(f"\n  best epoch {best_ep+1}, val loss {best:.4f}")
    cmv, _ = confusion(model, lva, device, amp)
    reportCm(cmv, f"VAL  ({a.split} split)")
    cmt, _ = confusion(model, lte, device, amp)
    acc, f1, iou = reportCm(cmt, f"TEST ({a.split} split)  <- the number")

    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = a.tag or f"{a.split}_b{a.base}_{'gray' if a.gray else 'rgb'}_s{a.seed}"
    torch.save({"model": model.state_dict(), "args": vars(a),
                "in_channels": cin, "base": a.base,
                "test_acc": float(acc), "test_f1": f1.tolist(),
                "test_iou": iou.tolist()}, out / f"cloud_unet_{tag}.pt")
    json.dump({"split": a.split, "seed": a.seed, "base": a.base,
               "in_channels": cin, "params": npar, "best_epoch": best_ep + 1,
               "test_accuracy": float(acc),
               "test_f1": dict(zip(CLASSES, f1.tolist())),
               "test_iou": dict(zip(CLASSES, iou.tolist())),
               "confusion": cmt.tolist(),
               "n_train": len(tr), "n_val": len(va), "n_test": len(te)},
              open(out / f"cloud_unet_{tag}.json", "w"), indent=2)
    print(f"\n  wrote {out}/cloud_unet_{tag}.pt and .json")


if __name__ == "__main__":
    main()
