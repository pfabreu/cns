# scripts/

Off-repo experiments. Nothing here is on the flight path, and nothing here is
built by CMake.

## Should there be a network here at all?

A UNet stub used to live in this directory, ported from Teague & Chahl (2026)
on the strength of their reported F1 0.77 under cloud. It was removed after the
question was actually measured. The measurement is below; it changed the answer
twice.

### What cloud and flare actually cost

Matched stars per frame, six frames per condition (`/tmp` probe, reproduced in
`testMeshBackground`):

| condition | matched | vs clear |
|---|---|---|
| clear | 38.3 | 100% |
| cloud: glow + attenuation | 25.8 | 67% |
| cloud: GLOW only | 31.5 | 82% |
| cloud: ATTENUATION only | 33.0 | 86% |
| **flare only** | **19.3** | **50%** |

Two things fall out. **Lens flare is the dominant contaminant**, worse than
cloud, and it is purely additive -- a smooth gradient with no attenuation.
Cloud's damage splits roughly evenly between glow and attenuation.

### A 60-line classical mesh takes the flare

A SExtractor-style background -- sigma-clipped estimate per 64 px cell, median
filtered across nodes, bilinearly interpolated, subtracted:

| condition | plain | mesh | gain |
|---|---|---|---|
| clear | 38.3 | 36.8 | 0.96x |
| cloud 0.6 | 25.8 | 26.2 | **1.01x** |
| flare 0.5 | 19.3 | 34.0 | **1.76x** |
| cloud + flare | 13.8 | 24.2 | 1.75x |

It recovers nearly all of the flare loss and **does nothing whatever for
cloud**. That is now in the library as `DetectorConfig::bg_mesh_px`.

### Why the neural background estimator is not obviously worth building

The plan was: predict the background field with a small network at 1/8
resolution, subtract, hand the residual to the classical detector. Attractive
because a 1/8-resolution output structurally CANNOT delete a 2 px star, and it
is ~64x less compute than full-resolution segmentation.

The measurement undercuts it. The additive component -- the part any background
estimator can address -- is already taken by the classical mesh. What is left is
cloud, and **cloud damage is not a background problem**. It attenuates
starlight, which is photons that never arrived, and its glow raises the SHOT
NOISE floor (sigma^2 ~ B). Subtracting a background cannot un-add Poisson noise.
That was the mechanism I had wrong: I assumed cloud hurt mainly through
threshold bias, and it does not.

So a learned background estimator would have to beat the mesh's 1.76x, not beat
nothing. On this evidence it has very little room.

### What would still justify a network

**Sharp cloud edges.** A median mesh smooths across a discontinuity and leaves
residuals on both sides, where a network sees the edge as an edge. The cloud
model here is a sum of smooth sinusoids and has no sharp edges at all, so this
simulator is structurally incapable of showing that failure -- which is exactly
the limitation that makes the result above provisional.

**Faint-star recovery under attenuation.** A segmenter can in principle use
context to accept a low-SNR detection a threshold would reject. That is the one
thing genuinely outside the classical pipeline's reach.

**Their real-data result.** Teague & Chahl report 85.9% identification on real
flight data, +23.3% over the baseline. That is a measurement on real sky and it
stands regardless of anything here.

### The stub is back, with a safety property

`train_unet.py` returns, but its output is a **prior**, not a decision. It feeds
`DetectorConfig::prior`, which lowers the detection threshold between
`threshold_k` and `threshold_k_low`:

    k(x,y) = k_high - (k_high - k_low) * p(x,y)

The detection still has to clear a real significance floor computed from actual
photons, and the centroid is always measured on the original image. So the
network **cannot hallucinate a star** -- a confident prediction over empty sky
produces nothing, because there is no flux to exceed even k_low -- and it
**cannot delete one**, because a zero prior leaves the nominal threshold
untouched.

`testPriorIsBounded` asserts the exact equivalence:

| prior | equivalent to |
|---|---|
| 1 everywhere | hand-setting `threshold_k = threshold_k_low` |
| 0 everywhere | the ordinary detector, unchanged |

Measured: 25 detections / 21 matched, and 7 / 7, matching the hand-set
thresholds exactly in both directions. **The worst an adversarial model can do
is move you to an operating point you could have chosen yourself.** That is a
bound on the damage, and it holds with no assumption that the model is good.

This is the reason a segmentation network is defensible here while a network
that outputs detections directly is not. It also happens to be the mechanism
that suits the one job genuinely outside the classical pipeline's reach --
accepting a low-SNR detection under attenuation that a fixed threshold would
reject, which is precisely "lower the bar where there is reason to".

### Running it: PriorNet

`include/celestial/prior_net.hpp` loads the exported ONNX through `cv::dnn` and
produces the prior. Verified end to end against a hand-built ONNX:

```
PriorNet::available() = true
load(missing file)    = false          <- graceful
infer without model   = false, prior empty
load(toy model)       = true
first call after load = 17.3 ms
inference (steady)    =  7.3 ms        <- at 1/4 resolution
prior                 = 2354176 values, range 0.0004-0.9990
matched: no prior 10, with prior 17
```

Inference runs on a DECIMATED image and the prior is upsampled. That is not
only for speed: a coarse prior cannot target an individual pixel, so the model
is structurally incapable of fine-grained mischief even inside its permitted
range. 7.3 ms sits comfortably inside the 100 ms frame budget alongside the
24.8 ms matched filter -- check `./build/bench` before scaling the model up.

**Warm up at load, not in flight.** `cv::dnn` initialises layers, allocates
buffers and selects implementations lazily on the first forward pass, so an
un-warmed model costs **165 ms on frame one** against 7 ms thereafter. On an
aircraft that lands on whichever frame happens to be first. `PriorNet::load`
runs a dummy pass to pay it up front. (An earlier version of this note quoted
17.8 ms as the inference cost; that was a first-call measurement and wrong in
both directions.)

**OpenCV 5 is the target, 4.6 is what was tested.** Only the stable subset is
used -- `readNetFromONNX`, `blobFromImage`, `forward`, `resize` -- which is
unchanged between 4.x and 5.x, so both work. One deliberate detail: PriorNet
does NOT call `setPreferableBackend` or `setPreferableTarget`. OpenCV 5 ships
three DNN engines behind one API and defaults to `ENGINE_AUTO`, trying the new
graph engine first and falling back to the classic one; the new engine is where
the 80%+ ONNX coverage, shape inference and operator fusion live. Pinning a
backend forces the classic engine, which is only correct if you want CUDA or
OpenVINO and is a straight loss on CPU. Leaving both unset gets the best
available engine on 5 and identical behaviour on 4.

**OpenCV is optional.** `find_package(OpenCV QUIET)`; without it PriorNet
compiles to a stub that reports unavailable and every caller falls through to
the ordinary detector. Both builds are tested (`testPriorNetFallback`). It is
used for exactly one thing -- ONNX inference, so a learned model does not drag
in a second runtime -- and NOT for thresholding, connected components,
centroiding, the matched filter or the mesh background, all of which are
hand-written because each is measured and documented and a library call would
lose the reasoning.

### The order of work

1. Enable `bg_mesh_px` where flare is expected. Free, already done, 1.76x.
2. Temporal median across frames -- stars move with attitude, flare and hot
   pixels do not. Cheap, classical, and targets the same contaminant from a
   different angle.
3. Get real night-sky footage. Everything below this line is unfalsifiable
   without it.
4. Add sharp-edged cloud to the renderer, since the current model cannot
   produce the failure a network would fix.
5. Only then the learned prior, benchmarked against the mesh -- 1.76x under
   flare is the number to beat, not zero.

`tools/gen_dataset` writes the image/mask pairs, and the labels (`StarLabel`:
centroid, both streak endpoints, truncation flag) are exact.

## Real night-sky data

Synthetic training with real validation is the only version of this experiment
that means anything. Nothing below is a drop-in: they are all-sky fisheye
cameras on fixed mounts taking multi-second exposures, so they have no motion
blur, a completely different PSF and a much wider field. They are useful for
CLOUD APPEARANCE and as a transfer check, not as training data for a 53 deg
camera on a moving aircraft.

**LenghuSky-8** (arXiv:2603.16429, March 2026) is the best fit found.
429,620 frames at 512x512 over eight years (2018-2025) from a premier
astronomical site, 81.2% night-time coverage, and -- the reason it matters --
**star-aware cloud masks**, background masks, and per-pixel altitude-azimuth
calibration. Manually labelled balanced subset of 1,111 images. Astrometric
calibration to ~0.37 deg at zenith.

**Dev et al., nighttime sky/cloud segmentation** (ICIP 2017),
github.com/Soumyabrata/nighttime-imaging. Much smaller, Singapore, segmentation
ground truth, code included. Good for a quick sanity check.

**LSST all-sky camera archive** (SMTN-005). 30 s exposures through every night,
~30 GB/night, Canon RGB Bayer, with SExtractor photometry of ~2000 bright stars
alongside sky brightness. Their own note is directly relevant to us: *in cloudy
conditions the coordinate solution fails and the stars become mis-identified* --
which is the failure mode this whole experiment is about. Access is via
`lsst-dev.ncsa.illinois.edu`, so not an anonymous download.

**What none of these give you** is night-sky footage from your airframe, with
your optics, at your frame rate, with real motion blur. That is the missing
input for validating the matched filter, for measuring the true star
identification rate, and for making any UNet result meaningful. It is the
highest-value thing to collect on the first flight.
