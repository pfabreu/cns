#!/usr/bin/env python3
"""Train a small UNet to produce a DETECTION PRIOR, and export ONNX.

    for s in 1 2 3 4 5 6 7 8; do ./build/gen_dataset data/train --frames 500 --seed $s; done
    ./build/gen_dataset data/val --frames 200 --seed 99
    python scripts/train_unet.py --train data/train --val data/val --epochs 20
    # -> star_unet.onnx, for cv::dnn

WHAT THIS IS FOR, AND WHAT IT IS NOT ALLOWED TO DO
--------------------------------------------------
The output is NOT a decision. It is fed to `DetectorConfig::prior`, where it
lowers the detection threshold between `threshold_k` and `threshold_k_low`:

    k(x,y) = k_high - (k_high - k_low) * p(x,y)

The detection still has to clear a real significance floor computed from actual
photons in the original image, and the centroid is always measured on the
original image. So:

  * the network CANNOT hallucinate a star -- a confident prediction over empty
    sky produces nothing, because there is no flux to exceed even k_low;
  * the network CANNOT delete a star -- a prior of zero leaves the nominal
    threshold untouched.

`testPriorIsBounded` asserts the exact equivalence: prior=1 everywhere is
identical to hand-setting k_low, prior=0 is identical to k_high. The worst an
adversarial model can do is move you to an operating point you could have picked
yourself. That is a bound on the damage, and it holds with no assumption that
the model is any good.

This is why a segmentation network is defensible here and a network that
directly outputs detections is not.

WHAT THE MEASUREMENTS SAY BEFORE YOU START
------------------------------------------
Matched stars per frame, this simulator:

  clear 38.3 | cloud 25.8 | flare 19.3 | cloud+flare 13.8

Lens flare is the dominant contaminant and it is purely additive, and a 60-line
classical mesh background already recovers it (1.76x, `bg_mesh_px`). What is
LEFT for a network is cloud, whose damage is attenuation plus a raised shot
noise floor -- photons that never arrived, and Poisson noise that cannot be
subtracted. So this model has to beat the mesh, not beat nothing.

The honest case for it is narrower than "networks are better at cloud":

  * sharp cloud EDGES, where a median mesh smooths across the discontinuity.
    Note the simulator's cloud is a sum of smooth sinusoids and has no sharp
    edges, so it cannot currently show this;
  * faint-star recovery under attenuation, where context can justify accepting
    a low-SNR detection that a fixed threshold rejects. This is the one thing
    genuinely outside the classical pipeline's reach, and it is exactly what
    the prior mechanism above is shaped to exploit safely.

AND THE TRAP
------------
A network trained on cloud we invented will beat a classical filter on cloud we
invented. Teague & Chahl's contribution is not that a UNet segments synthetic
stars; it is that theirs transferred to REAL FLIGHT DATA at 85.9%
identification. Until there is real night-sky footage from the actual airframe,
this is pipeline plumbing, not evidence.

SPLIT BY TRAJECTORY, NOT BY FRAME. Consecutive frames are nearly identical, so a
random frame-level split leaks validation into training. Each --seed of
gen_dataset is one trajectory; the commands above hold whole seeds out.

See scripts/README.md and "Learned detector" in NOTES-private.md.
"""

import argparse
import glob
import os

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset


def read_pgm(path):
    """Minimal binary PGM reader, 8- or 16-bit. Avoids a Pillow dependency."""
    with open(path, "rb") as f:
        assert f.readline().strip() == b"P5", f"not a binary PGM: {path}"
        dims = f.readline().split()
        w, h = int(dims[0]), int(dims[1])
        maxval = int(f.readline())
        n = w * h
        if maxval > 255:
            data = np.frombuffer(f.read(2 * n), dtype=">u2").astype(np.float32)
        else:
            data = np.frombuffer(f.read(n), dtype=np.uint8).astype(np.float32)
    return data.reshape(h, w)


class StarFrames(Dataset):
    """Image/mask pairs, cropped to tiles.

    Full frames are 1936x1216 and mostly empty sky. Training on random CROPS is
    both cheaper and better: it multiplies the effective dataset and stops the
    network learning a fixed vignetting pattern as a positional prior.
    """

    def __init__(self, root, tile=256, tiles_per_frame=4, train=True):
        self.imgs = sorted(glob.glob(os.path.join(root, "*_img.pgm")))
        if not self.imgs:
            raise SystemExit(f"no *_img.pgm under {root} -- run gen_dataset first")
        self.tile, self.tpf, self.train = tile, tiles_per_frame, train

    def __len__(self):
        return len(self.imgs) * self.tpf

    def __getitem__(self, i):
        path = self.imgs[i // self.tpf]
        img = read_pgm(path)
        mask = read_pgm(path.replace("_img.pgm", "_mask.pgm"))

        t = self.tile
        h, w = img.shape
        if self.train:
            y, x = np.random.randint(0, h - t), np.random.randint(0, w - t)
        else:
            k = i % self.tpf
            y = min((k // 2) * (h - t), h - t)
            x = min((k % 2) * (w - t), w - t)

        img = img[y : y + t, x : x + t]
        mask = (mask[y : y + t, x : x + t] > 127).astype(np.float32)

        # PER-TILE normalisation, not global. The background level moves with
        # cloud, flare and exposure; a network handed raw ADU would spend its
        # capacity learning those offsets instead of what a star looks like.
        # This mirrors the median+MAD the C++ detector already computes.
        med = np.median(img)
        mad = np.median(np.abs(img - med)) + 1e-6
        img = (img - med) / (1.4826 * mad)

        if self.train and np.random.rand() < 0.5:
            img, mask = img[:, ::-1].copy(), mask[:, ::-1].copy()
        if self.train and np.random.rand() < 0.5:
            img, mask = img[::-1].copy(), mask[::-1].copy()

        return torch.from_numpy(img)[None], torch.from_numpy(mask)[None]


def block(cin, cout):
    return nn.Sequential(
        nn.Conv2d(cin, cout, 3, padding=1, bias=False),
        nn.BatchNorm2d(cout), nn.ReLU(inplace=True),
        nn.Conv2d(cout, cout, 3, padding=1, bias=False),
        nn.BatchNorm2d(cout), nn.ReLU(inplace=True),
    )


class UNet(nn.Module):
    """Deliberately small: 4 levels, 16 base channels, ~0.5M parameters.

    Two reasons not to scale this up without measuring. It has to run at 10 Hz
    on a Raspberry Pi alongside everything else -- check ./build/bench, the
    detector is 3.8 ms of a 100 ms budget. And UNet is a poor
    architectural fit for tiny objects: a 2 px star is sub-pixel by the third
    downsample, so the deep path contributes little and the skip connections do
    the work. If this underperforms, a dilated stack with no downsampling is
    the more sensible next try, not a bigger UNet.
    """

    def __init__(self, base=16):
        super().__init__()
        b = base
        self.d1, self.d2 = block(1, b), block(b, b * 2)
        self.d3, self.d4 = block(b * 2, b * 4), block(b * 4, b * 8)
        self.pool = nn.MaxPool2d(2)
        self.u3 = nn.ConvTranspose2d(b * 8, b * 4, 2, stride=2)
        self.c3 = block(b * 8, b * 4)
        self.u2 = nn.ConvTranspose2d(b * 4, b * 2, 2, stride=2)
        self.c2 = block(b * 4, b * 2)
        self.u1 = nn.ConvTranspose2d(b * 2, b, 2, stride=2)
        self.c1 = block(b * 2, b)
        self.out = nn.Conv2d(b, 1, 1)

    def forward(self, x):
        s1 = self.d1(x)
        s2 = self.d2(self.pool(s1))
        s3 = self.d3(self.pool(s2))
        x = self.d4(self.pool(s3))
        x = self.c3(torch.cat([self.u3(x), s3], 1))
        x = self.c2(torch.cat([self.u2(x), s2], 1))
        x = self.c1(torch.cat([self.u1(x), s1], 1))
        return self.out(x)  # logits; the C++ side applies the sigmoid


def dice_bce(logits, target, pos_weight=50.0):
    """Dice + weighted BCE.

    The positive class is ~0.3% of pixels (measured on generated data), so plain
    BCE is minimised by predicting all-background: it will report 99.7% accuracy
    while detecting nothing. Dice is scale-free in class frequency; the weighted
    BCE term keeps gradients alive early.
    """
    bce = F.binary_cross_entropy_with_logits(
        logits, target, pos_weight=torch.tensor(pos_weight, device=logits.device)
    )
    p = torch.sigmoid(logits)
    num = 2 * (p * target).sum() + 1.0
    den = p.sum() + target.sum() + 1.0
    return bce + (1 - num / den)


@torch.no_grad()
def f1_score(logits, target, thr=0.5):
    """Pixel F1 -- comparable with Teague & Chahl's reported numbers, and NOT
    the figure of merit. See the note printed at the end of training."""
    p = (torch.sigmoid(logits) > thr).float()
    tp = (p * target).sum()
    fp = (p * (1 - target)).sum()
    fn = ((1 - p) * target).sum()
    return (2 * tp / (2 * tp + fp + fn + 1e-6)).item()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--val", required=True)
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--tile", type=int, default=256)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--onnx", default="star_unet.onnx")
    a = ap.parse_args()

    dev = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {dev}")

    tr = DataLoader(StarFrames(a.train, a.tile, 8, True), batch_size=a.batch,
                    shuffle=True, num_workers=4, drop_last=True)
    va = DataLoader(StarFrames(a.val, a.tile, 4, False), batch_size=a.batch,
                    num_workers=2)

    net = UNet().to(dev)
    print(f"parameters: {sum(p.numel() for p in net.parameters())/1e6:.2f}M")
    opt = torch.optim.AdamW(net.parameters(), lr=a.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, a.epochs)

    best = -1.0
    for ep in range(a.epochs):
        net.train()
        tot = 0.0
        for x, y in tr:
            loss = dice_bce(net(x.to(dev)), y.to(dev))
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += loss.item()
        sched.step()

        net.eval()
        f1s = []
        with torch.no_grad():
            for x, y in va:
                f1s.append(f1_score(net(x.to(dev)), y.to(dev)))
        f1 = float(np.mean(f1s)) if f1s else 0.0
        print(f"epoch {ep+1:3d}  loss {tot/max(1,len(tr)):.4f}  val F1 {f1:.4f}")
        if f1 > best:
            best, _ = f1, torch.save(net.state_dict(), "star_unet.pt")

    print(f"best val F1 {best:.4f}")

    # ONNX for OpenCV 5's cv::dnn -- operator coverage went from ~22% to over
    # 80%, so a plain UNet loads without adding a second runtime. Dynamic axes
    # so the C++ side can feed whatever tile size fits its budget.
    net.load_state_dict(torch.load("star_unet.pt"))
    net.eval()
    torch.onnx.export(
        net, torch.randn(1, 1, a.tile, a.tile, device=dev), a.onnx,
        input_names=["image"], output_names=["logits"],
        dynamic_axes={"image": {0: "n", 2: "h", 3: "w"},
                      "logits": {0: "n", 2: "h", 3: "w"}},
        opset_version=13,
    )
    print(f"wrote {a.onnx}")
    print(
        "\nNEXT, in order:\n"
        "  1. The pixel F1 above is NOT the figure of merit. Feed sigmoid(logits)\n"
        "     to DetectorConfig::prior and compare MATCHED STARS PER FRAME\n"
        "     against the mesh background (testMeshBackground) -- the mesh gets\n"
        "     1.76x under flare, so that is the number to beat.\n"
        "  2. ./build/bench with the network in the loop. 100 ms budget.\n"
        "  3. transit_scenario --imaging for end-to-end km.\n"
        "  4. Do not believe any of it until it has been checked against real\n"
        "     night-sky footage."
    )


if __name__ == "__main__":
    main()
