#!/usr/bin/env python3
"""Fetch the LenghuSky-8 labelled cloud-segmentation set (~22 MB).

    python3 scripts/fetch_lenghusky8.py            # -> scripts/data/lenghusky8_seg/

WHY ONLY THIS FILE. The full dataset is 84 GB: 100 image tars (~20 GB), 100
DINOv3 logit tars (~40 GB), plus background masks and astrometric calibration.
None of it is needed to train a segmenter, because `data/segmentation.zip`
embeds each annotated image as base64 inside its own labelme JSON. 22 MB gets
you all 1,111 images AND all 1,111 label sets. Do not download the tars.

The GitHub README says "252 images used for benchmarking"; that is a subset.
The zip holds 1,111, which is the number the paper quotes for its balanced
labelled set, so it is the whole thing.

SOURCE. https://huggingface.co/datasets/ruiyicheng/LenghuSky-8 (Apache-2.0),
paper arXiv:2603.16429, code https://github.com/ruiyicheng/LenghuSky-8.

WHAT THIS IS NOT. All-sky fisheye, fixed mount, multi-second exposures at an
astronomical site. Not a 53 deg star camera on a moving aircraft: no motion
blur, a different PSF, a far wider field. It is evidence about CLOUD
APPEARANCE on real sky and a transfer check, and it is not training data for
the flight sensor. See scripts/README.md.
"""
import argparse
import pathlib
import sys
import urllib.request
import zipfile

URL = ("https://huggingface.co/datasets/ruiyicheng/LenghuSky-8/"
       "resolve/main/data/segmentation.zip")
EXPECT_BYTES = 22_104_744
EXPECT_FILES = 1111

ap = argparse.ArgumentParser()
ap.add_argument("--out", default=str(pathlib.Path(__file__).parent / "data"),
                help="directory to extract into (default scripts/data)")
ap.add_argument("--force", action="store_true", help="re-download if present")
a = ap.parse_args()

out = pathlib.Path(a.out)
dest = out / "lenghusky8_seg"
if dest.exists() and not a.force:
    n = len(list((dest / "data").glob("*.json")))
    print(f"already present: {dest} ({n} json)")
    sys.exit(0)

out.mkdir(parents=True, exist_ok=True)
zp = out / "segmentation.zip"
print(f"downloading {URL}")


def _progress(blocks, bs, total):
    got = blocks * bs
    pct = 100.0 * got / total if total > 0 else 0.0
    print(f"\r  {got/1e6:6.1f} / {total/1e6:.1f} MB  {pct:5.1f}%", end="")


urllib.request.urlretrieve(URL, zp, _progress)
print()
got = zp.stat().st_size
if got != EXPECT_BYTES:
    print(f"  WARNING: got {got} bytes, expected {EXPECT_BYTES}. "
          f"Upstream may have changed; check the dataset card.")

with zipfile.ZipFile(zp) as z:
    z.extractall(dest)
zp.unlink()

n = len(list((dest / "data").glob("*.json")))
print(f"extracted {n} labelme json to {dest}/data")
if n != EXPECT_FILES:
    print(f"  WARNING: expected {EXPECT_FILES}")
