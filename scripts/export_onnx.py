#!/usr/bin/env python3
"""Export a trained cloud UNet to ONNX and check it loads under cv::dnn.

    scripts/.venv/bin/python scripts/export_onnx.py scripts/runs/cloud_unet_group_b16_rgb_s42.pt
    python3 scripts/export_onnx.py --check-only scripts/runs/cloud_unet_...onnx   # system cv2

WHY ONNX AT ALL. Only so the existing C++ can consume it. `PriorNet`
(include/celestial/prior_net.hpp) already loads ONNX through `cv::dnn` and is
built and tested; ONNX means a learned model needs no second runtime and no new
dependency. Nothing about the ONNX format is load bearing for the result.

OPSET 12 AND dynamo=False. The repo builds against whatever OpenCV is
installed -- 4.6.0 here, with OpenCV 5 the stated target. The TorchScript
exporter at opset 12 emits a graph both read; torch 2.9+ defaults to the dynamo
exporter, whose output the 4.x importer often will not take. Pinned on purpose.

SHAPE, AND THE ONE PLACE OpenCV 4.6 AND 5.0 ACTUALLY DIVERGE. `PriorNet::infer`
decimates the flight frame by 4 and runs on whatever comes out -- 1936x1216 / 4
= 484x304, neither 512 nor a power of two. The UNet upsamples with `size=` taken
from the skip connection rather than `scale_factor=`, so odd dimensions survive
the round trip. That emits a `Shape` node, and there the two importers part:

    dynamic axes (default)   OpenCV 5.0 OK      4.6 REFUSES
                             "parseShape: !isDynamicShape assertion failed"
    --static W H             OpenCV 5.0 OK      4.6 OK

With a fixed input size the exporter constant-folds `Shape` away and 4.6 is
happy. So export dynamic if you build against OpenCV 5, and `--static 484 304`
if you build against 4.x -- which is what is installed here, so that is the
variant the C++ can actually load today. Measured, not assumed: this script
loads both under whichever cv2 it is run with. Batch stays 1 either way,
because cv::dnn feeds one frame.

THREE THINGS THAT DO NOT LINE UP WITH PriorNet AS IT STANDS, and each one is a
deliberate non-decision rather than an oversight:

  1. CHANNELS. Trained on 3-channel RGB all-sky; `PriorNet::infer` wraps a
     1-channel 16-bit mono frame. Train with --gray for a drop-in, at whatever
     accuracy colour was worth.
  2. SCALE. This model wants float in [0,1]. `PriorNet` calls
     `blobFromImage(small)` with the default scalefactor 1.0 on a CV_16U frame,
     so it would feed 0..65535. Whoever wires this up must set the scale.
  3. SEMANTICS. `DetectorConfig::prior` is P(star here), and lowering the
     threshold where the prior is high. This model outputs P(sky), P(cloud),
     P(contamination). Which of those becomes the prior is a real design
     question -- P(sky) says "be permissive where it is clear", the opposite
     says "recover faint stars under thin cloud", and scripts/README.md argues
     the second is the only thing outside the classical pipeline's reach.
     Not decided here, because it should be decided by measurement.
"""
import argparse
import pathlib
import sys
import time

import numpy as np


def buildModel(ckpt_path):
    import torch
    from train_cloud_unet import UNet
    ck = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    m = UNet(ck["in_channels"], 3, ck["base"])
    m.load_state_dict(ck["model"])
    m.eval()
    return m, ck


def checkCv2(onnx_path, cin, hw=(304, 484)):
    """Load under cv::dnn and run one frame. Same importer the C++ uses."""
    import cv2
    print(f"  cv2 {cv2.__version__}")
    net = cv2.dnn.readNetFromONNX(str(onnx_path))
    blob = np.random.rand(1, cin, hw[0], hw[1]).astype(np.float32)
    net.setInput(blob)
    out = net.forward()
    t0 = time.time()
    n = 10
    for _ in range(n):
        net.setInput(blob)
        net.forward()
    dt = (time.time() - t0) / n
    print(f"    input  {blob.shape}")
    print(f"    output {out.shape}  range {out.min():.3f}..{out.max():.3f}")
    print(f"    {dt*1000:.1f} ms/frame at {hw[1]}x{hw[0]} "
          f"(the flight frame decimated by 4)")
    assert out.shape[:2] == (1, 3) and out.shape[2:] == hw, out.shape
    print("    shape and dynamic spatial axes OK")
    return dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("checkpoint", help=".pt from train_cloud_unet.py, or a "
                                       ".onnx with --check-only")
    ap.add_argument("--out", default=None)
    ap.add_argument("--opset", type=int, default=12)
    ap.add_argument("--check-only", action="store_true",
                    help="skip export, just load the .onnx under cv::dnn. Use "
                         "with the SYSTEM python so it exercises the same "
                         "OpenCV the C++ links.")
    ap.add_argument("--channels", type=int, default=3,
                    help="only for --check-only, which cannot read the .pt")
    ap.add_argument("--static", nargs=2, type=int, metavar=("W", "H"),
                    default=None,
                    help="export at a FIXED size instead of dynamic spatial "
                         "axes. Needed for OpenCV 4.x. Use the decimated "
                         "flight frame: --static 484 304")
    ap.add_argument("--check-size", nargs=2, type=int, metavar=("W", "H"),
                    default=[484, 304])
    a = ap.parse_args()

    if a.check_only:
        cw, ch = a.check_size
        checkCv2(pathlib.Path(a.checkpoint), a.channels, (ch, cw))
        return

    import torch
    sys.path.insert(0, str(pathlib.Path(__file__).parent))
    model, ck = buildModel(a.checkpoint)
    cin = ck["in_channels"]
    out = pathlib.Path(a.out or str(pathlib.Path(a.checkpoint).with_suffix(".onnx")))

    print(f"  {a.checkpoint}")
    print(f"    in_channels {cin}  base {ck['base']}  "
          f"test acc {ck.get('test_acc', float('nan')):.4f}")

    if a.static:
        w, h = a.static
        dummy = torch.randn(1, cin, h, w)
        dyn = None
        print(f"    exporting STATIC {w}x{h} (OpenCV 4.x compatible)")
    else:
        dummy = torch.randn(1, cin, 512, 512)
        dyn = {"image": {2: "h", 3: "w"}, "logits": {2: "h", 3: "w"}}
        print(f"    exporting DYNAMIC h,w (needs OpenCV 5)")
    torch.onnx.export(
        model, dummy, str(out),
        opset_version=a.opset, dynamo=False, do_constant_folding=True,
        input_names=["image"], output_names=["logits"], dynamic_axes=dyn)
    print(f"  wrote {out}  ({out.stat().st_size/1e6:.2f} MB)")

    import onnx
    onnx.checker.check_model(onnx.load(str(out)))
    print("  onnx.checker OK")

    # Numerical agreement with torch. For a dynamic export this deliberately
    # uses a size that is NOT the export size -- that is the whole point.
    cw, ch = (a.static if a.static else a.check_size)
    with torch.no_grad():
        x = torch.randn(1, cin, ch, cw)
        ref = model(x).numpy()
    import cv2
    net = cv2.dnn.readNetFromONNX(str(out))
    net.setInput(x.numpy())
    got = net.forward()
    err = float(np.abs(ref - got).max())
    print(f"  max |torch - cv::dnn| = {err:.2e} at {cw}x{ch} "
          f"({'OK' if err < 1e-3 else 'MISMATCH'})")

    print()
    checkCv2(out, cin, (ch, cw))
    print(f"\n  now check the OpenCV the C++ actually links:")
    print(f"    python3 scripts/export_onnx.py --check-only "
          f"--channels {cin} {out}")


if __name__ == "__main__":
    main()
