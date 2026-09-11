# Working notes — private

Not for publication. This is the lab notebook: what was tried, what failed, and
what had to be retracted. The README carries the durable technical content; this
carries the history, which is mostly a record of being wrong in instructive ways.

Kept because nearly every dead end here looked correct beforehand, and several
were only caught when real flight data contradicted the simulator. A future
version of me will otherwise retry them.

## Star-aided attitude, tested at last -- and it does not help the FIX

Implemented in `analyse_log --star-aided`: run `StarAidedAttitude` on real gyro
from the log plus the absolute star attitude from each frame, then recompute the
orbit fix with the corrected attitude and compare.

| orbit | attitude error, EKF3 -> star-aided | fix (hdg-wtd) | fix (circle) |
|---|---|---|---|
| 1 | 0.944 -> **0.396 deg** | 5.79 -> 5.78 km | 7.49 -> 5.81 km |
| 2 | 0.952 -> **0.395 deg** | 2.94 -> 2.94 km | 3.17 -> 2.96 km |
| 3 | 1.002 -> **0.395 deg** | 11.53 -> 11.48 km | 12.23 -> 11.54 km |

**Attitude improves 2.4x. The fix does not move.** Same degeneracy as every
other attempt: converting the star attitude from ECEF to NED needs a POSITION,
and an error there enters as a constant NED-frame rotation -- which is not
body-fixed, so the orbit does not average it out and it passes straight into the
fix. Replacing the attitude wholesale just reproduces the reference position it
was given.

The residual 0.395 deg is the unknown 0.4 deg MOUNTING error, which is now the
dominant term in the star-derived attitude. That one IS body-fixed, which is
why the circle fit improves a little (7.49 -> 5.81) and heading-weighted does
not.

FIRST ATTEMPT WAS CIRCULAR and reported 5.79 -> 0.01 km. It converted ECEF to
NED using `seg.centre`, the TRUE orbit centre, so truth was feeding the
correction. Worth remembering what a circular test looks like: it produces a
number too good to be true, and 10 m of position from a 0.4 deg uncalibrated
mounting is exactly that.

### The de-drift variant, implemented and measured

Correct only the TIME-VARYING part of EKF3's attitude, keeping its absolute
tilt. The construction uses RELATIVE rotations, which is what makes the unknown
mounting cancel exactly:

    dS_i = S_i S_0^T   = T_i T_0^T                <- mounting cancels
    X_i  = B_i^T B_0   = C_ekf,i^T (P^T dS_i P) C_ekf,0

so `C_ekf,i X_i = C_true,i B_0` -- the residual is frozen at its frame-0 value,
which is constant and body-fixed, exactly what the orbit average removes.

An earlier version conjugated into the body frame at frame ZERO and applied the
result at frame i. Those are different frames once the aircraft has turned; the
expression above applies the correction in the body frame at i.

**Measured, and the mechanism is CORRECT:**

| reference position | residual scatter | fix |
|---|---|---|
| estimated (honest) | 0.0521 deg | 5.79 km |
| TRUE (oracle, `CNS_ORACLE_P=1`) | **0.0011 deg** | **0.00 km** |

Given a known position it removes essentially all the time-varying error and the
fix collapses to zero. Given an estimated one it reproduces that estimate.

**Why, and the step I got wrong.** `P^T (S_i S_0^T) P` conjugates the relative
star rotation -- which is LARGE, up to 360 deg round the orbit -- by the
position error. Conjugating a large rotation by a small angle perturbs it to
FIRST order; it is second order only for SMALL rotations, and I asserted second
order without checking which case applied. 0.0521 deg x 111 km = 5.8 km is
precisely the reference error that went in.

Doing it incrementally frame-to-frame does not escape it either: each step's
error is small but they share the same epsilon, so they accumulate
systematically to the same total.

**This is the fourth independent confirmation of the same degeneracy.** Tilt and
position are one observable to a celestial fix, and every route that turns an
absolute or relative attitude into a position correction reintroduces the
dependency it was trying to escape. The oracle row is the cleanest statement of
it in the project: the estimator is right, the information is absent.

## The same mistake, made twice

`transit_scenario` dropped its timed fix trigger early on, for a reason recorded
at the time: without a vertical reference a straight-leg fix carries the whole
body-fixed error and is worth ~40 km. **`Pipeline` kept its timed trigger for
another dozen sessions.** The lesson was learned in one file and not applied to
the other, and it surfaced only when a real SITL transit produced fixes
consistently displaced behind the aircraft -- which reads as a frame or sign bug
and is nothing of the kind.

The tell was in the data the whole time: fixes with 12 of 12 heading bins came
out 3.6-23.5 km, fixes with 1-2 bins came out 20.9-39.7 km. Nobody looked at the
bin count against the error until the symptom forced it.

Worth generalising: when the same rule applies to two code paths, the second one
does not inherit it, and the failure will look like a new bug rather than an old
one.

## The horizon sensor breaks the celestial compass

Restored the horizon module, added Lepton 2.5 / Boson 320 presets, and measured
it on the DETERMINISTIC transit -- `fake_sitl` cannot A/B anything, its runs
differ because the node samples a UDP stream and frame timing varies.

Six seeds:

| | compass ON | compass OFF |
|---|---|---|
| no horizon | **7.21 km** | 10.66 km |
| 1x Lepton 2.5 | 11.00 km | 8.91 km |

The sensor works. Tilt error halves (23.96 -> 12.22 arcmin) and individual fixes
go from 5-7 km to 0.5-2 km, exactly as the isolated probe predicted (6.16 ->
2.46 km). What it breaks is the COMPASS: heading residual against a 3 deg
magnetometer bias goes 0.096 -> 0.480 deg.

Eliminated:
* not the mount calibration -- the effect persists with `--boresight 0`
* not fix quality -- the fixes improve substantially
* not a fake_sitl artefact -- this is the seeded deterministic scenario

Suspect: the alias in `orbit.hpp`. An AHRS yaw bias appears as a rotation about
the LOCAL VERTICAL expressed in body axes, and `mountTiltOnly` must remove
exactly that component for the compass to observe anything. The horizon changes
the estimated attitude and therefore the estimated vertical, so the axis being
removed is no longer the right one. Not confirmed.

**The lesson worth keeping**: an isolated probe said 2.5x and the end-to-end
system said 1.5x WORSE. Two subsystems that each improve attitude were competing
for the same information, and only the full pipeline showed it. The probe was
not wrong -- it answered a narrower question than the one that mattered.

## Scoreboard of things that did not work

| attempt | why it failed |
|---|---|
| Post-hoc tilt estimation, three variants | A celestial fix cannot separate tilt from position; they are one observable |
| Drift-aware orbit estimator | Drift at the orbit frequency is exactly degenerate with position |
| Star-tracker de-drift over an orbit | Comparing an ECEF attitude to NED reintroduces the position dependence |
| Star-aided attitude replacing EKF3 wholesale | Attitude 2.4x better, fix unchanged -- inherits the reference position |
| Mid-transit loiters | 2-3x worse; each recalibration writes the instantaneous bias into the mounting |
| Splitting the tilt filter states | Bias and drift load identically on the measurement; over-parameterised |
| GNC-TLS robust solver | Lost to RANSAC by an order of magnitude under gross outliers |
| Polar-coordinate matched filter | Rotation centre lands outside the frame; 6x SLOWER than the naive version |
| Neural background estimator | The classical mesh already takes the recoverable (additive) part |
| Two circle-fit divergence guards | First was mis-scaled; second fired one step too late |

## Things I got wrong and had to retract

* **"The attitude error is body-fixed."** The paper's premise, inherited without
  question. SITL shows a clean one-per-revolution component that does not
  average away.
* **"The NED-mean predicts the residual to 3%."** True GPS-aided (ratio 1.03),
  false denied (0.56). An aided-case result generalised without warrant.
* **"Do not fix in the first few minutes."** Drawn from EKF3 settling with GPS.
  Denied there is nothing to settle. Would have been published as a rule.
* **`ahrs_turn_coupling = 0.02`.** Invented. Measured 0.0009 aided, 0.0041
  denied -- 5x to 20x too pessimistic, and it drove the orbit-radius decision.
* **"The maneuver term will make the numbers worse."** It improved them, because
  the horizon camera absorbed it.
* **"Don't reach for a network."** True for centroiding a Gaussian PSF, and it
  says nothing about detection under cloud, which is the harder problem.
* **"17.8 ms inference."** A first-call measurement. Steady state is 7.3 ms and
  the un-warmed first frame is 165 ms.
* **Pinning the OpenCV DNN backend.** Would have forced the classic engine on
  OpenCV 5 and thrown away the reason to want 5.
* **Doubting the GNSS denial.** EKF3's own status bits said `dead_reckoning`
  set, `using_gps` cleared. I had an anomaly I could not explain and blamed the
  setup instead of my check. Twice in one session.

## Mount calibration absorbs the yaw bias

The scenario used to report `residual 2.934 deg` -- essentially no correction.
An earlier revision of this document blamed the simulator, claiming
`simulateObservationsWithAttitude` rendered the observations from the estimated
attitude. **That diagnosis was wrong.** That function renders from
`eulerToDcm(ft.roll, ft.pitch, ft.yaw)` -- the TRUTH -- and the supplied
estimate is only stored as `FrameData::C_l_b_est` / `yaw_est`. It reads neither
`ahrs_bias` nor `ahrs_noise`, so the `e2.ahrs_bias.setZero()` in the scenario is
a no-op on that path and is probably what made the setup look self-consistent.
The observations always carried a real, recoverable heading error.

The actual cause is `recalibrateMount`. It solves

```
C_b_c = C_l_b_est^T * R
```

and `C_l_b_est` is wrong by the magnetometer bias. Writing the estimate as
`E = Rz(b) * T` and the truth as `R = T * M`:

```
E^T R  =  T^T Rz(-b) T * M
```

The leading factor is a rotation by `b` about `T^T * z_ned`, the local down
direction in body axes. Level, that IS the boresight -- so a magnetometer bias
and a rotated camera bracket are the same object, and the calibration returns
the bracket, because a mounting is all it can return. Then `celestialHeading`
computes `C_l_b = R * C_b_c^T`, which unwinds a bias-contaminated mounting
straight back to `C_l_b_est` and reports no error at all.

This is a SILENT failure. The compass returns a small, stable, converged-looking
number no matter how large the bias. Instrumented, the scenario showed it
plainly -- the first fix, taken during the loiter before calibration, recovered
2.842 of the 3 deg (residual 0.158); the calibration then booked 2.939 deg into
the mounting's clocking and every later fix returned about -0.06.

### Why iterating does not fix it

The obvious repair -- estimate heading, correct the attitude, recalibrate, repeat
-- does not work, and fails deceptively. Substituting `T = Rz(-b) E`, the
per-frame mounting implied by a candidate bias `b` is

```
M_i(b) = [rotation by (b - b_true) about E_i^T z_ned] * M
```

`E_i^T z_ned` depends on roll and pitch and **not on yaw**. So in a constant-bank
loiter that axis is pinned no matter how far heading sweeps, `M_i(b)` is equally
consistent across frames for EVERY `b`, and the bias is unidentifiable. The
iteration has a null direction and converges to its initialisation rather than
to the truth. Tried directly, it fixed heading and silently walked the boresight
from 0.389 to 0.836 deg to pay for it.

Bias is identifiable only to the extent that roll and pitch VARY across the
data -- bank diversity, not heading sweep. At 9 deg of bank the leverage is
`sin(bank)`, about 0.16, so heading to 0.1 deg would need mounting consistency
near 0.016 deg. That remains the principled route and is not implemented.

### What is implemented: tilt / clocking split

The two halves of the mounting are not equally compromised. Deviation from
nominal, in the CAMERA frame, splits into x/y (tilt of the boresight) and z
(clocking about it). Only z is aliased with yaw. `mountTiltOnly` keeps the tilt
and drops the clocking, so the scenario calibrates as before and then discards
the contaminated component:

```cpp
C_b_c = mountTiltOnly(recalibrateMount(window, r.position),
                      meanVerticalBody(window));
```

The alias axis is the LOCAL VERTICAL in body axes, not body z. They coincide
level and differ by the bank angle in a turn, and using body z leaves
`bias * sin(bank)` of the yaw behind AS TILT, which lands in the boresight --
0.47 deg at 9 deg of bank and a 3 deg bias, against a calibration that otherwise
reaches 0.02 deg. Both axes are asserted in `testCompassSurvivesCalibration`,
because choosing wrong shows up nowhere else.

This costs nothing the position solver needs. `OrbitResult` scores the boresight
DIRECTION, which does not depend on clocking; measured, it moves from 0.0161 to
0.0246 deg.

**Measured, transit scenario** (`make scenario`, 3 deg bias injected):

| heading source | residual | DR + fixes, mean | DR without fixes, peak |
|---|---|---|---|
| magnetometer | 3.00 deg uncorrected | 3.65 km | 63.06 km |
| celestial compass | **-0.011 deg** | **2.75 km** | 57.13 km |

Regression coverage is `testCompassSurvivesCalibration` in `test/test_orbit.cpp`,
which asserts BOTH halves: that the full mounting absorbs the bias (residual
3.000 deg, so the failure mode stays documented and cannot return unnoticed) and
that the tilt-only mounting recovers it (0.033 deg), with the boresight
direction unchanged to 0.0002 deg.

### What this assumes

Clocking is no longer estimated in flight; it is assumed nominal. That is an
assumption about the airframe, not a result, and the heading numbers above are
conditional on it. It is a reasonable one -- clocking is the mounting DOF that a
bench measurement pins most easily, and it does not drift -- but it must be
stated wherever the heading figure is quoted. It is also flattered here: the
simulator injects a pure TILT (`boresight_axis = (1,0,0)`) and no clocking error
at all, so "nominal clocking" is exactly right rather than approximately right.
A clocking error would appear as heading error roughly one-for-one.

The magnetometer is also NOT out of the loop; see "How the two headings combine"
below.

## Robust solver: RANSAC, not GNC

`singleFrameFix` used Graduated Non-Convexity with a truncated least-squares
cost (Yang et al., RA-L 2020), chosen over the paper's RANSAC on the grounds
that it is deterministic, uses every measurement at every iteration, and needs
no initial guess. All three are true and the choice was still wrong.

Measured on a realistic near-zenith star field -- a 53 deg field pointed up sees
nothing below ~64 deg elevation, and sampling the whole sky flatters GNC --
20 stars, 60 trials per cell, mean fix error:

| misidentified | plain LS | GNC | RANSAC |
|---|---|---|---|
| 0% | 0.34 km | 0.34 km | **0.34 km** |
| 10% | 25.48 km | 6.22 km | **0.34 km** |
| 25% | 45.06 km | 31.69 km | **0.42 km** |

RANSAC is an order of magnitude better under outliers and identical without
them; GNC barely improves on plain least squares. Two likely causes: the
initial least-squares seed is badly corrupted by gross outliers and the
graduated schedule converges from it to a poor local minimum, and a
misidentified star here is GROSS (0.3-3 deg), which is the regime a minimal-seed
hypothesis test handles best.

There is also a symptom worth chasing: GNC IMPROVES as its noise bound is
loosened, 6.65 -> 2.43 km going from 0.02 to 0.5 deg, which is backwards. That
points at the residual metric, which divides the plane residual by
`sin(zenith)` clamped at 0.05 -- inflating residuals by up to 20x for exactly
the near-zenith stars this camera sees, so legitimate measurements get
truncated as outliers. Not confirmed.

It also reset every weight to 1 on entry, silently discarding the refraction
weighting `singleFrameFix` applies two lines earlier. RANSAC does not touch
weights, so the switch fixed that too. GNC has since been REMOVED rather than
left in place, since a robust solver nobody calls is a liability.

## Star-aided attitude: the mechanism, implemented for testing

**This is the one remaining avenue against the dominant error term, and unlike
the three that failed, it has not been shown not to work.** It is implemented
(`star_aided.hpp`, `StarAidedAttitude`) so the ArduPilot experiment has
something to run. It is NOT validated on real hardware, and no number here
should be quoted as a system result.

### Why it is not the thing that failed

Three attempts to remove the AHRS tilt error from the FIX failed against one
wall: a celestial fix cannot separate a tilt error from a position error,
because they are a single observable. Post-hoc estimation therefore has nothing
to work with, whatever the estimator.

This estimates a different quantity. Not absolute tilt -- GYRO BIAS, which is
observable from the SEQUENCE of attitudes rather than any one of them:

* star directions in ECEF do not depend on the observer's position, so
  `kabsch(v_cam, ecef_dirs)` is an absolute attitude with **no position in it**;
* two of them give the true rotation between their epochs, exactly;
* the gyro's integral over the same interval differs from that by the bias.

Drift-as-a-RATE is observable where tilt-as-an-OFFSET is not. The two results do
not conflict. And it matters because an aided gyro's attitude error stops
growing and becomes very nearly CONSTANT over an orbit -- and a constant
body-fixed error is exactly what the orbit average removes.

### What was built

`StarAidedAttitude` is a 6-state multiplicative EKF: attitude error and gyro
bias. `predict` takes a body rate and removes Earth rate (0.6 deg over a 150 s
orbit -- far too large to drop); `update` fuses an absolute body-to-ECEF star
attitude. `ErrorModel` gained a gyro model (`gyro_bias_sigma`, `gyro_bias_tau`,
`gyro_arw`) and `simulateGyro` produces samples from the truth attitude
sequence, because **the simulator previously had no gyro at all** -- the AHRS
attitude was synthesised directly, which is precisely why this could not be
tested before.

`testStarAidedAttitude`, 302 s of flight, star fix every 2 s at 60 arcsec:

| | mean | peak |
|---|---|---|
| gyro alone | 0.0792 deg | 0.1455 deg |
| star-aided | **0.0167 deg** | **0.0423 deg** |

with the gyro bias recovered at 0.20 deg/hr. The filter is correct. That is all
this shows.

### What it does NOT show

It does not show the transit improves. The gyro model is one this project
invented, and the real question -- does aiding a Cube Orange's attitude with a
hobby-grade IMU reduce the tilt error enough to matter, given what EKF3 already
does internally -- can only be answered against the real thing.

### Integrating with ArduPilot

EKF3 has no general-purpose attitude measurement input. It fuses external-nav
YAW (`EK3_SRC*_YAW = 6`) and position/velocity, but not roll and pitch. So:

1. **Emit yaw as external-nav yaw.** Easy, and nearly pointless here --
   `celestialHeading` already recovers heading to ~0.1 deg.
2. **Run this filter OUTSIDE EKF3**, on `RAW_IMU` / `SCALED_IMU` gyro plus the
   star attitude, and use its output in place of EKF3's attitude in the fix
   path. No autopilot changes. This is what the class is shaped for, and it is
   the recommended experiment.
3. **Add a tilt fusion path inside EKF3.** The correct answer, a much larger
   job, and it needs care: EKF3's attitude feeds the star matcher, so the loop
   has to be closed without the matcher depending on its own output.

### Prior art

Sodern's Astradia does exactly this in hardware -- see "Related work" below.
Their residual budget after attitude drift is removed is dominated by
accelerometer bias and star-tracker-to-IMU alignment, which are the same two
terms this project reached independently. That is the strongest external
evidence that the approach is sound; what is unproven is doing it cheaply.

## Learned detector: a roadmap

Teague & Chahl (2026) hold F1 0.77 under CLOUD where every classical baseline
collapses below 0.26. Cloud is the one contaminant the matched filter does not
address -- it is structured, non-stationary, and looks like signal at some
scales. This is the plan to find out whether a network is worth it here,
scoped to stay small.

**Stage 1 only.** Segmentation. Skip their inter-frame keypoint detector: it is
the harder half and its main extra benefit -- precise camera-to-AHRS temporal
alignment -- is a separate problem worth solving separately.

**SUPERSEDED BY MEASUREMENT.** The renderer now has cloud and flare
(`SensorModel::cloud_amount` / `flare_amount`) and `tools/gen_dataset` emits
image/mask pairs, but the training script was removed and the roadmap below is
no longer the recommendation. Measured: lens flare, not cloud, is the dominant
contaminant (-50% matched stars), it is purely additive, and a 60-line
SExtractor-style mesh background recovers it (1.76x, now
`DetectorConfig::bg_mesh_px`). What remains is cloud, whose damage is
attenuation plus a raised SHOT NOISE floor -- neither recoverable by
subtracting anything, and neither a background problem. See
`scripts/README.md` for the numbers and for what would still justify a network.

### Step 1: give the simulator weather (the only real prerequisite)

`renderFrame` models motion blur, PSF, shot and read noise, and hot pixels. It
does NOT model cloud or lens flare, which is exactly the regime where a network
is supposed to win. Nothing can be evaluated until it does.

Keep it crude on purpose:

* **cloud** -- two or three octaves of value noise, smoothly advected between
  frames, used to ATTENUATE stars multiplicatively and RAISE the background
  additively. Cloud is not an overlay; it dims what is behind it and glows.
* **flare** -- one or two broad radial gradients anchored to a bright off-axis
  source, fixed in the CAMERA frame so it does not move with the stars. This is
  what kills Bernsen thresholding (their F1 0.02) and what a temporal median
  removes cheaply.

Add as `SensorModel::cloud_amount` / `flare_amount`, defaulting to zero so every
existing number is unchanged.

### Step 2: emit training pairs

The labels already exist. `StarLabel` carries the true centroid (`u,v`), both
streak endpoints (`u0,v0` / `u1,v1`) and a truncation flag, so a segmentation
mask is the `(u0,v0)->(u1,v1)` segment rasterised and dilated by the PSF sigma.
A `tools/gen_dataset` writing PNG image/mask pairs over randomised attitude,
rate, exposure, cloud and flare is a few hundred lines and reuses everything.

Aim for ~20k frames. Hold out by TRAJECTORY, not by frame: consecutive frames
are nearly identical and a random split will report a number that means nothing.

### Step 3: train, off-repo

A small UNet -- 4 down, 4 up, 16 base channels -- in PyTorch, in `scripts/`,
not in the C++ tree. Dice + BCE, since the positive class is well under 1% of
pixels. This should train in hours on one GPU. Export ONNX.

### Step 4: inference in C++

OpenCV 5's DNN module is the natural path: ONNX operator coverage went from
~22% to over 80%, it benchmarks competitively with ONNX Runtime on CPU, and it
keeps the dependency count at one. A `CentroidMethod`-style switch in
`DetectorConfig` selecting `Segmentation` keeps it behind the existing API.

Budget check first. The matched filter already misses the frame budget by 2.5x
(see above); a UNet forward pass on 2.35 Mpx will not be free. Measure with
`./build/bench` before believing any of it.

### Step 5: evaluate against what exists

The harness is already there and the metric is already chosen: MATCHED STARS
PER FRAME via `testMatchedFilter`, and end-to-end km via
`transit_scenario --imaging`. Compare four ways -- plain, matched filter,
UNet, UNet + matched filter -- across clear, cloud and flare.

### The trap to avoid

**A network trained on cloud we invented will beat a classical filter on cloud
we invented.** That result is circular and worth nothing. The comparison is only
meaningful on REAL night-sky footage from the actual airframe and optics, which
is the same thing needed to validate the matched filter and to measure the star
identification rate honestly. Teague & Chahl's contribution is not that a UNet
segments synthetic stars; it is that theirs transferred to real flight data at
85.9% identification.

So: do not start this before there is real footage to test against. Until then
the honest ordering is top-hat, temporal median, predicted-position windowing,
and making the matched filter fast enough to fly.

## MEASURED IN ARDUPILOT SITL

Two of the three numbers this project invented have now been measured against a
real EKF3. One was wrong by a factor of twenty, and the central assumption of
the paper's method turns out to be false as stated.

**Read the caveat first.** Both runs had SITL's GPS ENABLED and EKF3 fusing it.
That is what the paper's own flight did, and it is what `cns_sitl.parm` sets up,
but it means EKF3 had good velocity for its centripetal correction and good
aiding for attitude. GNSS-denied it has only airspeed and a wind estimate.
**Everything below is an OPTIMISTIC BOUND**, not the GNSS-denied number.

### The attitude error is NOT body-fixed

The whole method assumes it. `analyse_log` decomposes it, seven full
revolutions, 150-173 m loiter. Body-frame error against heading, orbit 4:

| heading | roll | pitch | yaw |
|---|---|---|---|
| 0 | -0.049 | +0.429 | +1.005 |
| 90 | +0.016 | +0.154 | +0.759 |
| 180 | +0.174 | +0.017 | +0.779 |
| 270 | +0.130 | +0.287 | +1.149 |
| **peak-to-peak** | **0.223** | **0.466** | **0.477** |

Flat would mean body-fixed. It is not flat: there is a clean one-cycle-per-
revolution pattern, same phase across orbits 2-5. So a real EKF3 error has a
substantial component that does NOT average away over a heading sweep.

**But the decomposition predicts the position error, which is the useful
result.** The NED-mean horizontal -- the part that survives averaging -- against
the measured fix error:

| orbit | predicted | measured (hdg-wtd) | ratio |
|---|---|---|---|
| 1 | 11.65 km | 15.15 | 1.30 |
| 2 | 9.15 | 10.46 | 1.14 |
| 3 | 9.31 | 10.02 | 1.08 |
| 4 | 9.52 | 9.95 | 1.05 |
| 5 | 10.01 | 10.22 | 1.02 |
| 6 | 2.91 | 3.24 | 1.11 |
| 7 | 4.14 | 3.26 | 0.79 |

Mean ratio **1.03** excluding the first orbit, which is a settling transient.
So the claim to make is not "the error is body-fixed and the orbit removes it"
but **"the error decomposes, the orbit removes the body-fixed part, and the
surviving NED component is measurable in flight and predicts the residual."**
That is a better claim, because it is checkable rather than assumed.

### EKF3 settling, aided -- an observation, NOT an operational rule

With GPS, the NED-mean horizontal fell from 0.086 deg (orbits 2-5) to 0.026 deg
(orbit 6), yaw error going 0.887 -> 0.440 -> -0.092, and position error followed:
~10 km early, 3.2-3.3 km once settled.

**This was briefly written up here as a rule -- "do not fix in the first few
minutes" -- and that was wrong.** It does not survive GNSS denial: denied, the
error does not settle at all, it wanders and changes sign (see below). And it is
in any case an observation about THIS SITL configuration, whose IMU noise and
bias are whatever `cns_sitl.parm` sets rather than a Cube Orange's. The
mechanism -- no absolute attitude reference means attitude drifts -- is general;
the magnitude and the settling behaviour are not.

### The heading-weighted mean, on real EKF3 data

First-iteration error, the paper's Eq. (23) against heading weighting:

| orbit | naive | heading-weighted | gain |
|---|---|---|---|
| 2 | 19.64 km | 9.78 | 2.0x |
| 4 | 17.99 | 9.95 | 1.8x |
| 6 | 17.95 | 1.67 | **10.7x** |
| 7 | 17.38 | 4.54 | 3.8x |

Consistent and large, on real flight data rather than simulation. This is now
the best-evidenced contribution in the project.

### Turn coupling: 0.0009 aided, 0.0041 denied, against 0.02 invented

Weighted fits through the origin over the banks actually flown: +5..+45 deg
across 100/200/400 m loiters aided (~30k samples), +5..+25 deg at 250 m denied
(~8.5k samples). **Denial costs 4.6x**, which is the mechanism behaving exactly
as predicted -- centripetal compensation needs velocity, GPS supplies it, and
airspeed plus a wind estimate does not. The denied value is now the default.
The invented 0.02 was still 5x too pessimistic.

Two qualifications. It is NOT linear -- implied coupling is 0.0002 at 15-20 deg
rising to 0.0018 at 45 deg, so compensation degrades with bank rather than
leaking a fixed fraction; a linear model is conservative low and optimistic
high. And the 50-65 deg bins reverse sign with n = 19-63, which are transients
rather than steady turns and are discarded, as are all negative-bank bins
(n = 176 total).

**This moved a design decision.** The orbit-radius trade-off was dominated by
the calibration term, which is several times weaker than modelled. Re-measured
at the DENIED coupling of 0.004: 150 m gives 8.78 km, **250 m gives 7.27 km**,
400 m gives 8.76 km, 600 m gives 10.23 km. The calibration loiter moved from
400 m to 250 m, and 250 m is the optimum at both the aided and denied coupling,
so the choice is not sensitive to which value is right.

### The circle fit could diverge

Orbit 1 of the SITL log: 48.88 km at iteration 1, **73.90 km after six
iterations**. `iterateOrbit` assumed that recalibrating on a better fix yields a
better fix, which fails when the averaging is poorly conditioned -- a partial
revolution, an unsettled EKF3, a near-degenerate arc. More iterations was making
it worse. Now guarded: the loop stops when the estimate moves further than the
previous fit's own spread can justify.

### GNSS DENIED: the numbers that actually apply

Same airframe, same flight, GPS denied at t = 137.8 s (`SIM_GPS1_ENABLE 0`,
confirmed in the log by GPS Status 6 -> 1 and NSats 10 -> 3), then four full
revolutions at 250 m.

| | GPS aided | GNSS denied |
|---|---|---|
| NED-mean horizontal, orbits 2-5 | 0.085 deg | **0.158 deg** |
| ...once settled | 0.032 deg | does not settle |
| position error range | 3.2 - 10.2 km | 4.0 - 15.9 km |

**Denial roughly doubles the irreducible error** -- 1.9x the aided value, 5x the
aided-and-settled one. The 3.2 km measured with GPS was flattering by about a
factor of five.

**EKF3 does not converge denied; it wanders.** Body-frame roll runs +0.26, +0.16,
-0.25, -0.28 across the four orbits and yaw +0.93, +0.76, -0.39, -0.51 -- it
changes SIGN. Again: specific to this SITL setup, instructive rather than a
general result.

**The NED-mean predictor does NOT survive denial.** Measured over predicted:

| | ratios | mean |
|---|---|---|
| aided | 1.14 1.08 1.05 1.02 1.11 0.79 | **1.03** |
| denied | 0.78 0.82 0.28 0.37 | **0.56** |

It over-predicts by up to 3.6x. The earlier claim that the decomposition
predicts the residual to 3% was an AIDED-CASE result and should not have been
generalised. The likely reason is visible in the data: aided, the error is
quasi-static so its orbit mean is meaningful; denied it drifts within each orbit
and flips between them, so a mean over the revolution does not represent what
the position solve sees. It errs conservative, which is the safe direction, but
it is not a validated predictor.

### Which estimator, and the divergence guard

First-iteration error, denied, mean over four orbits:

| estimator | iteration 1 | converged |
|---|---|---|
| naive, Eq. (23) | 22.74 km | 10.37 km |
| heading-weighted | 11.78 km | 10.36 km |
| **circle fit** | **6.93 km** | 9.14 km |

The two corrections fix DIFFERENT errors, which is why the ordering changes with
conditions. Heading weighting removes a SAMPLING bias; the circle fit removes a
MISALIGNMENT bias by modelling the circle a body-fixed error traces and taking
its axis. Denied, the body-fixed error is large (body-frame mean 0.19-0.32 deg),
the circle is big, and the fit has plenty to grip -- so it wins. Aided and
settled, the body-frame mean was 0.042-0.066 deg, the circle is small, and
heading weighting won instead. Both beat naive at first iteration in every
orbit.

### The circle fit must not be iterated

`iterateOrbit` recalibrates the mounting between passes. That is right for the
averaging methods and WRONG for the circle fit, because **the fitted radius IS
the misalignment** -- recalibrating from the fit's own output applies the same
correction twice, and the fixed-point iteration can walk to a wrong answer.

Eight real orbits, four aided and four denied:

| | mean |
|---|---|
| circle fit, iteration 1 | **6.64 km** |
| circle fit, "converged" | 9.37 km |

Iterating helps in **2 of 8**. Worst cases 5.13 -> 15.19 km denied, and
31.68 -> 67.53 km aided. `iterateOrbit` now returns a single pass for
CircleFit; the averaging methods still iterate, since for them recalibration is
doing real work. `testCircleFitIsNotIterated` pins it.

CONFIRMED on the same SITL logs after the change. Circle fit, four denied
orbits: 15.19 -> 5.13, 10.08 -> 4.50, 4.25 -> 8.97, 7.03 -> 9.12; mean
**9.14 -> 6.93 km**. The aided transient orbit went 67.53 -> 31.68 km. Note it
is not uniform -- orbits 3 and 4 were being HELPED by iteration and got worse.
The wins are simply larger than the losses.

After the change, on GNSS-denied data:

| estimator | mean | worst |
|---|---|---|
| naive, Eq. (23) | 10.37 km | 15.21 km |
| heading-weighted | 10.36 km | 15.90 km |
| **circle fit, 1 pass** | **6.93 km** | **9.12 km** |

### The simulator and SITL disagree about which estimator to use

Switching the transit scenario to CircleFit throughout, six seeds:

| | simulator | SITL, denied |
|---|---|---|
| heading-weighted | **7.36 km** | 10.36 km |
| circle fit | 7.59 km | **6.93 km** |

**They disagree, and the real data should win.** The scenario default is
therefore LEFT ALONE rather than changed on the strength of a simulator that
cannot reproduce the effect in question.

The gap is identifiable. `ErrorModel` generates a drift that is essentially
body-frame, while the SITL logs show a pronounced ONE-PER-REVOLUTION NAV-FRAME
component -- pitch swinging 0.47 deg peak-to-peak with consistent phase across
orbits 2-5. The circle fit grips exactly that structure; the simulator does not
generate it, so it cannot see the advantage.

**This is the next thing to fix in the model**: add a nav-frame error term to
`ErrorModel` so the simulator reproduces the 1/rev signature, then re-run this
comparison. Until then the simulator cannot arbitrate estimator choice, and any
estimator recommendation should cite the SITL numbers rather than the scenario.

**Two guards were tried before this one, and both were wrong.** The first
bounded the step against the fitted circle radius -- useless, because the radius
is the misalignment expressed as a distance, tens of km for a fraction of a
degree of mounting error, so the bound sat far above any real step and never
fired. The second was a contraction test, stopping when a step grew rather than
shrank; that is scale-free and correct in principle, but it fires one step too
late. The damage is done by the SECOND pass, which then sits at a wrong fixed
point with tiny subsequent movement, so nothing downstream looks wrong. Neither
guard was as good as not iterating.

**Keeping the first estimate is worth more than the choice between
estimators** -- 6.64 vs 9.37 km, against 10.36 vs 10.37 for heading-weighted
against naive on the same orbits.

### What is still unmeasured

`ahrs_drift_sigma` / `ahrs_drift_tau` still come from the paper's Figure 5, not
from this airframe. And every number above needs remeasuring with GPS denied,
which is the configuration the project is actually about.

## Related work, as of September 2026

### The paper's own sequel

Teague & Chahl, "Star detection for low-altitude atmospheric star trackers",
J. Opt. Soc. Am. A 43(9) 1502-1517, published 1 September 2026
(doi:10.1364/JOSAA.600488). Same authors, now Adelaide University and DST Group.
It replaces the detection front end of the 2024 Drones paper this project is
built on.

Two-stage: a UNet segmenting stars, then a multi-task INTER-FRAME keypoint
detector that maps the geometric bisection across the temporal boundary of two
consecutive frames -- which resolves the directional ambiguity and endpoint
fading that defeat single-frame thresholding when a star is a motion-blurred
streak. Trained on synthetic imagery with blur, cloud and flare, validated on
real flight data.

| | reported |
|---|---|
| segmentation F1, clear | 0.85 vs 0.40 best baseline (Niblack) |
| segmentation F1, cloud / flare | 0.77 / 0.72 vs 0.26 / 0.31 |
| keypoint error, dynamic | 0.51 px mean |
| reprojection error, real flight | 0.60 px vs 0.92 px CoG, 1.28 px endpoint |
| **star ID success, 1000 real frames** | **85.9%, +23.3% over baseline** |

**Three things this project should take from it.**

1. **Our simulator is far too kind to the matcher.** `--imaging` reports 0.3% of
   frames unusable. They measure 85.9% star ID success on real flight data with
   a BETTER detector than ours, implying roughly 14% frame loss at best and
   ~37% with the baseline method. Our +7% imaging penalty is therefore a floor,
   not an estimate. The gap is cloud, flare and real motion blur, none of which
   `renderFrame` models.

2. **Our centroid noise assumption is optimistic.** We use 10-20 arcsec (0.2-0.4
   px); their real-data reprojection error is 0.60 px even with the learned
   detector. Since the fix is attitude-dominated this probably does not move the
   headline, but the assumption should be labelled.

3. **Their inter-frame keypoint incidentally solves our timestamp problem.**
   The bisection point is defined AT the boundary between two frames, so it
   gives precise temporal alignment between the star sensor and the AHRS. Camera
   -to-AHRS skew is listed as an untested risk here and is worth 1.4 arcmin per
   10 ms at 2.4 deg/s.

### Sodern Astradia, and independent confirmation of our error budget

Sodern (France) launched Astradia in 2025, a DAYTIME endoatmospheric star
tracker for aircraft, coupled to an IMU in a strapdown configuration with no
pointing mechanism. Reported: 1.2 Nm (2.2 km) at 1 h and 2.4 Nm (4.4 km) at
10 h, 95%, from a civil-aviation-grade INS.

This matters here for two reasons.

**It is the architecture recommended as the open experiment above.** The star
tracker supplies an absolute attitude reference that kills inertial attitude
drift, rather than the drift being estimated out afterwards -- exactly the
"feed the star attitude into EKF3" item. It is a productised system, so the
approach is sound; what is unproven is doing it with a hobby-grade IMU.

**Their residual error budget matches ours.** With attitude drift removed, they
report the dominant terms become ACCELEROMETER BIASES and STAR-TRACKER-TO-IMU
ALIGNMENT. Those are precisely the two this project converged on independently:
maneuver coupling (an accelerometer-derived vertical problem) and the mounting
calibration. Arriving at the same two terms from a different direction is the
strongest external check available on the analysis here.

The caveat is cost class: a daytime tracker plus an aviation INS is a different
budget from a Raspberry Pi and a Cube Orange, and their 2.2 km at 1 h is not
comparable to our 8.8 km over 194 km on that hardware.

### Also worth a look

- Chasan & Bingul, "Visible and short-wave infrared star tracker objective
  designs for UAV navigation", JOSA A 43(9) 1472-1479 (2026) -- same issue,
  optical design for the same application.
- Sensors 26(15) 4790 (2026), a 2020-onwards review of automated celestial
  navigation, covering star trackers, solar tracking and horizon detection.

### Star pipeline, wired through

`--imaging` now runs the real pipeline in the transit -- render, detect,
centroid, match -- instead of ideal star vectors. It was previously a flag that
was parsed, printed and never used.

Measured on a short transit, four seeds: **6.38 km ideal, 6.81 km with plain
detection (+7%), 6.43 km with the matched filter (+1%)**. Per frame, plain
detection finds ~11 stars and matches ~9 with 0.3% of frames unusable; with the
filter it finds ~124 and matches ~12, and no frame is unusable.

So real detection and matching is NOT a meaningful error source, which is what
the attitude budget predicted: at 1 deg of attitude = 100 km, a 10 arcsec
centroid sits three orders of magnitude below the dominant term. It is off by
default because it is ~50x slower and this scenario is also a regression test.

The residual risk from the star pipeline is availability, not accuracy -- a
matcher that fails entirely yields no fix at all, and no fix means dead
reckoning alone at 36 km. That is a binary risk and the 0.3% unusable-frame
figure is the number to watch, not the 7%.

### Matched filter: use the gyro instead of learning the streak

Teague & Chahl (2026) attack blurred-star detection with a UNet, and their
baselines -- adaptive Gaussian, Niblack, Bernsen -- are all single-frame spatial
thresholds that know nothing about the streak they are hunting. That is
deliberate: their method "obviates the need for angular rate sensors
altogether", which is the right goal for a standalone star tracker.

**We have a gyro.** The AHRS gives the body rate, so the smear direction and
length are predictable in every pixel, and a star becomes a signal of KNOWN
shape in approximately Gaussian noise. The optimal linear detector for that is a
matched filter: integrate along the streak and recover the SNR that smearing
spread out. Twenty lines, no training data, no inference cost. Much of what
their network must learn, we can simply know.

`DetectorConfig::omega_cam` / `exposure_s` / `focal_px` enable it. MATCHED stars
per frame -- the figure that feeds the fix, not raw detections:

| exposure | smear | plain | matched filter |
|---|---|---|---|
| 20 ms | 4.8 px | 4.2 | **11.8** |
| 50 ms | 12.1 px | 9.2 | **46.8** |
| 100 ms | 24.2 px | 13.8 | **53.8** |
| 200 ms | 48.3 px | 28.4 | **76.0** |

**2.7 to 5x more identified stars.** `testMatchedFilter` asserts at least 2x.

Three things went wrong on the way, all worth keeping:

**The smear is not uniform.** In an orbit the dominant rate is about the
BORESIGHT, which ROTATES the field rather than translating it. A star at radius
r from the principal point smears tangentially by `r*omega*t`, so direction and
length both vary across the frame. A single global kernel is wrong; the
rotational optical flow is evaluated per pixel.

**Per-TILE kernels flood the detector.** Reusing one kernel over a 64 px tile
leaves a step in the filtered map at every tile edge, and every step is a local
maximum: 73-2603 false detections per frame against 0.1-3.4 for a single global
kernel. Evaluating the flow per pixel costs a few flops against a tap loop that
is already there, so tiling bought nothing.

**Connected components are the wrong detector for a filtered map.** The filter
smooths, so its noise is spatially correlated and a median+MAD threshold floods.
The signal in a matched-filter output is a PEAK, so detect peaks, with a
suppression radius covering the streak length. And where there is NO smear --
at the principal point, under boresight rotation -- the kernel collapses to one
tap and the "filtered" map is the raw image, where peak detection at 5 sigma
fires on every noise spike. That region belongs to the ordinary detector, which
now runs alongside; the two are merged with a duplicate check.

The filter still produces many false peaks. That is fine and expected: they do
not land near a predicted star, so the matcher's ambiguity guard drops them.
**Judge this by matched stars, never by detection count** -- the raw counts look
alarming and mean nothing.

#### It is too slow to fly, and that is fixable

`./build/bench`, on a desktop x86 core, 2.35 Mpx frames, 10 Hz camera:

| stage | per frame | max rate |
|---|---|---|
| detect, plain | 4.95 ms | 202 Hz |
| **detect, matched filter** | **247 ms** | **4.0 Hz** |
| match to catalogue | 0.31 ms | 3224 Hz |
| per-frame fix (RANSAC) | 1.00 ms | 1001 Hz |

**50x the plain detector, and 2.5x over the 100 ms frame budget** -- on a
desktop. A Raspberry Pi 5 is several times slower again. As written this cannot
keep up with the camera, and every accuracy number above assumes it does.

The implementation is naive: a per-pixel loop over up to 129 taps across
2.35 Mpx, roughly 300 M memory-bound operations per frame, no SIMD. Nothing
else in the pipeline is close to the budget.

#### The polar idea, which was wrong

The intended fix was a coordinate change. To first order the smear field is a
rotation about the principal point plus a uniform translation, and a translation
can be absorbed into the centre, so the whole field is a PURE ROTATION about the
instantaneous centre `c = (w_x, w_y) * f / w_z`. About that point every streak
subtends the same angle at every radius, so in polar coordinates the kernel is a
fixed width and a running sum makes it O(1) regardless of streak length.

Implemented and measured: **1585 ms**, six times SLOWER than the naive filter it
replaced.

The centre is not near the image. In a coordinated turn the body rates are
roughly `q = W sin(bank)` and `r = W cos(bank)`, so at 23 deg of bank the
translation terms are ~40% of the rotation term and the centre lands ~810 px
outside the frame. The polar grid then has to span radii to ~1950 px, and
resolving one pixel of arc at that radius needs ~12000 angular bins -- the
warped image comes out LARGER than the original. Polar is the right move when
the rotation centre is inside or near the frame. For a banked orbit it is not.

#### What worked: adaptive decimation

Detection does not need full resolution. A star smeared over 24 px is still 6 px
at 4x decimation, and the centroid is computed at FULL resolution afterwards
regardless, so localisation is untouched.

The saving is close to cubic: decimating by `d` cuts pixels by `d^2` and kernel
length by `d`, so cost goes as `N*s/d^3`.

Decimation must be ADAPTIVE. A flat 4x costs accuracy at short smear -- a 12 px
streak becomes 3 px and the filter stops helping (40.8 -> 10.6 matched stars).
The factor is chosen per frame as the largest that keeps the decimated smear at
or above `mf_target_smear_px` (8 px). That also bounds the cost naturally,
because the expensive frames are the long-smear ones that tolerate the most
decimation.

| | before | after |
|---|---|---|
| detect, matched filter | 312 ms | **24.8 ms** |
| vs plain detector | 50x | 4x |
| vs 100 ms frame budget | 3.1x OVER | **4x under** |

Matched stars per frame are unchanged (40.8 / 38.4 / 54.4 at 50 / 100 / 200 ms),
so this is a 12x speedup for no measurable accuracy cost, and the filter now
fits the frame budget with room for a Raspberry Pi to be several times slower.

### Drift-aware orbit estimator: tried, does not work, and why

Averaging removes a body-fixed error because it rotates in the navigation frame
as heading sweeps. It does not remove AHRS tilt DRIFT, and that residual is the
entire gap: one orbit at 150 m gives 4.14 km with drift modelled and 0.66 km
without. So a joint estimator that models the drift instead of averaging over it
looked like the largest remaining prize in the project.

It was implemented. Small angles, body tilt as a constant plus a smooth
time-varying part, fitted jointly with position:

    y_i = c + M(psi_i) * (m + sum_k a_k P_k(t_i))

`y_i` is the per-frame fix offset in NED, `c` the position offset wanted, `m`
the constant body tilt, `P_k` Legendre polynomials on the window excluding the
constant term, and `M(psi)` the measured body-tilt-to-NED map (`d_along =
+R*pitch`, `d_cross = -R*roll`). With `K = 0` it reproduces the naive mean
exactly, which confirms the model is right.

**Every drift basis makes it worse**, and at one revolution catastrophically
(K=4 gives 28 km against 2.88 naive). The reason is a hard identifiability
limit, not a tuning problem.

Over one revolution, time maps linearly to heading, so a drift term `P_k(t)`
is also `P_k(psi)`. Modulated by `M(psi)`, a linear drift produces frequency
components at `1 +/- 1` -- including **DC**. A body tilt drifting at the orbit
frequency displaces the fix by a CONSTANT amount in NED, which is exactly what a
position error looks like. The two are not separable by any estimator.

And the component that survives averaging is precisely the one near the orbit
frequency. Sweeping the drift correlation time against a 38 s revolution:

| drift tau | naive | best drift-aware |
|---|---|---|
| 15 s | 4.71 km | 5.08 (worse) |
| 38 s | 3.15 km | 3.22 (worse) |
| 120 s | 1.99 km | 1.82 |
| 600 s | 1.08 km | 0.82 |

Fast drift is white over the orbit and averages out. Slow drift is nearly
body-fixed and averages out. The worst case is `tau` comparable to the
revolution period -- which is exactly where the real system sits, 60 s against
38 s. In the only regime where the estimator helps, averaging was already
working.

**So the 0.66 km figure is not reachable by better estimation.** It was quoted
earlier in this document as the prize for a drift-aware estimator; that was
wrong and the estimator is not in the tree.

The mechanism does explain the orbit-radius finding from first principles: error
scales with the ratio of revolution period to drift correlation time, which is
why 150 m (38 s) beats 400 m (100 s). But the lever is weak. Reaching 1 km at
`tau = 60 s` would need a ~4 s revolution, a radius of a couple of metres. What
is left is `1/sqrt(N)` in revolutions -- 1.96 km at four revolutions -- which is
a mission-time trade, not an algorithm.

### Is MAGSAC worth applying?

No, and the reason is worth recording because it also rules out the next three
ideas of the same shape.

MAGSAC's contribution is marginalising over an unknown noise scale, so no inlier
threshold has to be chosen, plus sigma-consensus weighted refinement instead of
a hard inlier cut. Both are real improvements to RANSAC in general. Neither buys
anything here, because **the robust stage is not losing anything to recover.**
At 0% outliers RANSAC and plain least squares both give 0.336 km -- that is the
measurement noise floor, and RANSAC is already sitting on it. Under outliers it
stays within 20% of that floor. There is no gap for a better robust method to
close.

Three cheap improvements were implemented and measured, 80 trials per cell,
20 stars:

| outliers | current | + exhaustive sampling | + normalised residual | + count-first scoring |
|---|---|---|---|---|
| 0% | 0.336 km | 0.336 | 0.336 | 0.336 |
| 10% | 0.348 | 0.348 | 0.346 | 0.348 |
| 25% | 0.405 | 0.405 | 0.406 | 0.422 |
| 40% | 0.368 | 0.368 | **0.305** | 0.403 |

Exhaustive enumeration of all C(n,3) minimal sets gives results IDENTICAL to 100
random draws -- 100 samples already finds the same model, so the sampling is not
the limitation. Count-first scoring is slightly worse. Only the normalised
residual helps, and only at implausible outlier rates; it is adopted anyway
because it makes `tolerance` mean the angle it claims to.

Note also that `RansacConfig::seed` defaults to 1 and the fix path never varies
it, so the solver is ALREADY deterministic run to run. The determinism argument
that motivated GNC was answered by a default value the whole time.

The general lesson: before adopting a better robust estimator, check what the
current one is losing. Here it is losing nothing, and the real uncertainty is
upstream -- the simulator generates no misidentifications at all, so the outlier
columns above are entirely synthetic. Wiring `--imaging` through would tell us
the actual outlier rate, and that number decides whether any of this matters.

**None of this moves any number in this document**, because the simulator
produces no misidentifications at all -- the outlier-free row is the only one
the scenario exercises. That is exactly why it went unnoticed, and it is the
second reason to finish item 4 above. `testRobustSolverChoice` pins the
comparison, including a guard that fires if GNC is ever fixed.

## How the two headings combine

They are used TOGETHER, as a complementary pair. The magnetometer is not
removed; it is corrected.

- `yaw_est` -- magnetometer-derived, via EKF3 -- supplies heading every frame at
  2 Hz. Dead reckoning integrates at frame rate and needs it there.
- `celestialHeading` supplies an absolute bias correction, but only at fixes
  (every 500 s here) and only when enough frames carry 3+ identified stars.

The scenario adds them at every propagation step:

```cpp
dr.predict(dt, 25.0, fd.yaw_est + hdg_corr);
```

so the magnetometer provides SHORT-TERM relative heading and the stars provide
the LONG-TERM absolute reference -- the same division of labour as a gyro/mag
complementary filter, one level up. What the compass removes is the
magnetometer's slowly varying error: hard/soft iron, and the declination lookup
at a drifting last-known position. What it does NOT remove is dependence on the
magnetometer between fixes. A spoof faster than the fix interval passes through
essentially uncorrected, and going fully magnetometer-free would mean
propagating heading on gyro alone between fixes.

Two caveats in the current implementation. `hdg_corr` is overwritten wholesale
at each fix with no smoothing or gating on `HeadingResult::spread`, so one bad
fix moves heading directly. And the `dead reckoning ALONE` row also uses
`yaw_est + hdg_corr`, so it is "no position fixes" rather than "no celestial
input" -- which is why that row moved when the compass started working. Both are
worth tightening before the row is quoted as a baseline.

## What you can test, and how

Three levels. The first two run anywhere; only the third needs ArduPilot.

### 1. The test suite and the deterministic scenario

```bash
make test        # 6 tests, ~2 minutes
make scenario    # the headline result, ~5 seconds
```

`make scenario` should print exactly this, every time, on any machine:

```
vertical reference    17.53' ->  1.54'   (32.5 km -> 2.85 km)
celestial fixes      18 over 183 km
dead reckoning ALONE mean  28.93 km, peak  57.15 km  (31.2% of distance)
dead reckoning + FIXES mean   2.15 km                    BOUNDED
```

If those numbers differ, something changed. It is a regression test as much as
a demonstration.

Worth sweeping, since each takes seconds:

```bash
./build/transit_scenario --dr-quality poor        # 15-30% drift
./build/transit_scenario --fix-revs 4             # longer fix orbits
./build/transit_scenario --rate 10                # the paper's frame rate
./build/transit_scenario --seed 3                 # drift makes this matter
./build/transit_scenario --boresight 0.0          # calibrated camera
./build/transit_scenario --imaging                # full render chain, slower
./build/transit_scenario --km 400                 # longer transit
```

### 2. The live node, without ArduPilot

```bash
make fake
```

Runs `celestial_node` against `fake_sitl`, which replays a generated trajectory
as real MAVLink over UDP. This exercises the transport -- stream requests,
message pairing, timing skew, the stale-pose guard -- none of which the
scenario touches.

```bash
./build/fake_sitl --pattern orbit --revs 4      # loiter: calibrates the mounting
./build/fake_sitl --pattern straight            # no maneuver at all
./build/fake_sitl --pattern transit --revs 6    # legs and loiters
```

Use `--pattern orbit` when testing calibration. The transit pattern opens with
a straight leg, so the coverage gate correctly refuses to calibrate and the
boresight never converges -- that is the gate working, not a failure.

### 3. ArduPilot SITL -- the part only you can run

Everything above uses simulated or replayed attitude. SITL is the only thing
that puts a REAL EKF3 in the loop, with real timing, real sensor error and real
mode changes. Every timing bug in this project was found there and none were
reachable from levels 1 and 2.

Three experiments, in order of value:

**a. Is EKF3's attitude error body-fixed?** The assumption everything rests on,
and still unverified.

```bash
python3 tools/dump_log.py logs/00000001.BIN > flight.csv
./build/analyse_log flight.csv
```

Validate the analyser first on a known input:

```bash
./build/analyse_log --synth 0.10 0.03 > fake.csv && ./build/analyse_log fake.csv

## Log analysis and the attitude-error decomposition


The goal is **not** to inject anything into ArduPilot yet. It is to answer one
question that everything else rests on:

> Milestone 1 injected the AHRS error as a constant **body-frame** vector. The
> whole orbit-averaging mechanism depends on that. Is it true of a real EKF3?

A body-fixed error rotates with heading, traces the circle of Figure 7, and
averages away. An earth-fixed error does not — the same distinction that made
refraction survive the orbit. This was baked in by construction in Milestone 1
and never verified, and the paper does not verify it either: their Figure 5
shows 0.2–0.3° of wander over 89 s but never decomposes it by heading.

That question needs no MAVLink code at all. Log-based, offline, done in an
afternoon.

---

### Phase 3a — log analysis (start here)

#### 1. Fly the orbits

```bash
sim_vehicle.py -v ArduPlane -f plane --console --map
```

**Make SITL lie to you first.** The default IMU is nearly ideal, so EKF3 will
look far better than a Cube Orange ever does and you will get a beautiful,
meaningless result:

```
param set SIM_ACC1_BIAS_X 0.05      # m/s/s -- this is your tilt error
param set SIM_ACC1_BIAS_Y 0.05
param set SIM_GYR1_RND 0.5
param set SIM_ACC1_RND 0.5
param set SIM_GYR1_SCALE 0.01
param set SIM_VIB_FREQ 30
param set SIM_WIND_SPD 8
param set SIM_WIND_TURB 3
param set SIM_MAG1_OFS_X 20
```

Accelerometer bias matters most: it maps straight to the tilt error that
becomes 111 km/deg. Sweep it and check the position error scales as predicted —
that is a clean closure test between Milestones 1 and 3.

Then fly loiters of a few different radii and both directions, ideally 2+
revolutions each. Both directions matters: see the CW/CCW note below.

#### 2. Dump the log

```bash
pip install pymavlink
python3 tools/dump_log.py logs/00000001.BIN > flight.csv
```

Uses `ATT` (EKF3 attitude — what the algorithm gets), `SIM` (**true** attitude
and position), `GPS` (absolute time), `ARSP` (airspeed). `SIM` only exists in
SITL logs; on a real flight there is no attitude truth, so section [1] below is
unavailable and you can only score the position estimate.

#### 3. Analyse

```bash
./build/analyse_log flight.csv            # ideal star vectors
./build/analyse_log flight.csv --imaging  # plus the full imaging chain
```

Three sections per detected orbit:

1. **EKF3 attitude error decomposition** — the question above.
2. **Position with ideal star vectors** (Milestone 1 path).
3. **Position through the imaging chain** (Milestone 2 path), with `--imaging`.

Keeping 2 and 3 separate is deliberate: it keeps the attitude cost and the
imaging cost independently attributable.

---

### Reading the output

#### The decomposition

```
body-frame mean  roll -0.0000  pitch +0.1000  yaw +0.0000
NED mean         N    +0.0299  E     -0.0000  D   +0.0106

body-frame mean magnitude : 0.1000 deg  (11.1 km if it did NOT average)
NED-mean horizontal       : 0.0299 deg  -> PREDICTED RESIDUAL 3.33 km
```

- **body-frame mean** is the part the orbit removes. Large is fine.
- **NED mean horizontal** is the part it cannot. This is the number that
  matters. A body-fixed error averages to ~0 in NED over a full revolution, so
  whatever survives there is earth-fixed.

#### The heading table

```
heading      n      roll     pitch       yaw
0          126   +0.0287   +0.0924   +0.0008
90         127   -0.0076   +0.0715   +0.0030
180        126   -0.0287   +0.1076   -0.0008
270        126   +0.0076   +0.1285   -0.0030
range            +0.0573   +0.0570   +0.0061
```

**Flat means body-fixed and the paper's method is sound.** A sinusoidal
variation, as above, is the signature of an earth-fixed component: it is
constant in NED, so it appears in body axes modulated at the orbit frequency.

#### Verifying the analyser itself

```bash
./build/analyse_log --synth 0.10 0.03 > fake.csv    # 0.10 deg body-fixed, 0.03 deg earth-fixed
./build/analyse_log fake.csv
```

On that input the tool reports:

| | injected | recovered |
|---|---|---|
| body-fixed pitch | 0.100° | **0.1000°** |
| earth-fixed north tilt | 0.030° | **0.0299°** |
| predicted residual | — | 3.33 km |
| **measured position error** | — | **3.34 km** |

The separation is exact and the NED mean predicts the un-averagable residual to
0.3%. Note also that the 0.10° body-fixed part — worth 11.1 km if it did not
average — contributes nothing to the final error. Run this before trusting the
tool on real data.

With `--imaging` the same case gives 3.42 km against 3.34 km, so the imaging
chain costs ~80 m on top of the attitude error. Consistent with Milestone 2.

---

### Specifically worth looking for

**Centrifugal compensation.** EKF3 corrects sensed acceleration for the turn.
An error in that correction is a tilt that is roughly constant in the body frame
during a steady turn — so it averages. But a **turn-direction-dependent**
component would not cancel between CW and CCW. The paper's Table 2 has both
(orbits 5–9): CCW 600 m gave 1.73 km against CW 600 m at 2.54–3.48 km. Might be
noise. Might be real. Fly both directions and compare the NED means.

**Bank-angle dependence.** Fly several radii. If the NED mean scales with bank
angle, the culprit is in the acceleration correction rather than the gyros.

**Wind.** `SIM_WIND_SPD` with a GPS-guided loiter reproduces the Section 4.4
mechanism, where the autopilot varies bank to hold a ground radius. Check
whether the heading-weighted average still recovers it on real EKF3 attitude
the way it did in simulation (14.95 km → 0.42 km).

---

### Phase 3b — injection (only after 3a)

```
AHRS_EKF_TYPE 3
GPS1_TYPE 14          # MAV -- accept GPS_INPUT
EK3_SRC1_POSXY 3
EK3_GPS_CHECK 0
ARMING_CHECK ...      # relax as needed
```

Send `GPS_INPUT` (msg 232) with `h_acc ≈ 4000` m, one fix per orbit. Expect
friction: EKF3's innovation gating, `EK3_POSNE_M_NSE`, and the arming/failsafe
logic all assume GPS-like rates and noise. A 4 km fix every three minutes is
nothing like a GPS stream. Budget real time for parameter fighting, and keep
SITL's own GPS enabled for failsafes exactly as the paper's flight did.

Version-check the disable parameter — it has been `SIM_GPS_DISABLE` and
`SIM_GPS1_ENABLE` depending on release.

Highest value-per-effort item here is actually **yaw**, not position: set the
yaw source to external and inject a heading. It costs almost nothing and
protects pitch and roll from magnetometer-driven innovation garbage.

### Phase 3c — actually GPS-denied

The thing nobody mentions: **you cannot fly a `LOITER` without a position
estimate.** Section 4.4's fixed-attitude orbit is not a loiter at all — it is a
commanded bank angle with a heading reference, i.e. FBWA with a held roll
input. So 3c needs RC override in FBWA or a small custom mode. Worth knowing
now so it does not surprise you.

---

## Findings, condensed

The blow-by-blow of SITL bring-up lives in git history. What survives is the
set of things that would cost a day each to rediscover.

### Timing — five bugs, one cause

Every one was two quantities that must refer to the same instant, and did not.
None were reachable from `fake_sitl`; four were invisible on the ground because
the error scales with turn rate.

| symptom | cause |
|---|---|
| `stars 0`, apparent 6 deg tilt | SIMSTATE arriving slower than ATTITUDE |
| both streams pinned at 4.0 Hz | MAVProxy re-sending REQUEST_DATA_STREAM over the per-message interval |
| every frame rejected on skew | SIMSTATE has no timestamp, so post-hoc arrival comparison always shows one message period of difference -- **pair at arrival instead** |
| fine on the ground, dead in a turn | image rendered for the exposure START, navigation using the END |
| a 10 km error floor that ignored more stars | flux-weighted centroid sits at exposure MIDDLE, navigation used the start |

The last one is the general lesson: **if a floor does not respond to more
stars, better estimators or more averaging, it is a bias.** Look for something
systematic in the timing chain before reaching for a better fit.

### Detection and matching

* `det` and `matched` must be logged SEPARATELY. "0 stars" is ambiguous -- a
  dark frame and a frame full of stars whose predictions all miss look
  identical and need opposite fixes.
* **Render deep, match shallow.** Mutual nearest-neighbour breaks when
  predictions outnumber detections: at V=6.0 there are ~235 predictions against
  ~15 detections, and almost nothing matches.
* Exposure has an optimum set by turn rate, and **both ends fail the same way**
  -- too few stars. Too short loses photons; too long spreads them below
  threshold.
* Turbulence (`SIM_WIND_TURB`) is real attitude motion at ~1 Hz and smears the
  image. `fake_sitl` cannot reproduce it, which is why the loopback reported
  20 detections where SITL gave 3.

### The calibration loop

The mounting and the position are coupled: recalibration runs AT the estimated
position, so a fix N km off injects N/111.2 deg of mounting error. That
converges when fixes are good and DIVERGES when they are not.

Three guards, each added after the corresponding failure:

* **Reject implausible fixes** before recalibrating. A bad fix is recoverable;
  recalibrating on one is not.
* **Require coverage.** Heading bins counted from the frames that SURVIVED, not
  from the accumulator -- a hard turn drops the star count and leaves the
  survivors clustered while the accumulator still claims a revolution.
* **Smooth, do not replace.** The true mounting is static, so successive
  estimates are repeated measurements of one fixed quantity. Replacing outright
  let it oscillate between 0.009 and 0.090 deg.

### Estimator selection

The circle fit is a **bootstrap**, not a better mean. It recovers position from
~110 deg of arc while the mounting is unknown, and must be switched off
afterwards: with no circle left, fitting one to a noisy shallow arc is
ill-conditioned. The trap is that a large reported radius is the SYMPTOM of
that failure, not evidence against it.


---

# SITL crossing campaign, 2026-09-10

Three GNSS-denied crossings flown end to end. The durable numbers are in
RESULTS.md; this is what it took to get them and what had to be withdrawn.

## The unattended SITL test could not start at all

Four silent failures, each of which reported something other than its cause.

* **`--out` is a MAVProxy option.** `sim_vehicle.py` consumes `opts.out` only
  inside `start_mavproxy()`, and `--no-mavproxy` returns before that, so BOTH
  `--out` flags were discarded and nothing arrived on either port. The symptom
  was "no heartbeat", i.e. a dead vehicle, not a dropped argument. Fixed with
  `-A --serial1=udpclient:...`, which reaches the ArduPlane binary, where
  `--serial1` is a real MAVLink port. `udpclient` and not `udpin` because
  `celestial_node` BINDS its port (`mavlink_source.cpp`), so SITL must send.
* **SITL blocks at startup** until a client connects to SERIAL0's TCP port.
  MAVProxy is normally that client. With `--no-mavproxy` nothing connected, so
  the simulator never left its wait loop and NO link produced a byte --
  including the node's, which made it look like a node fault.
* **SITL persists parameters in `./eeprom.bin`** (its cwd is the repo root). The
  GNSS denial performed by one run was therefore still in force at the next
  boot: "TIMEOUT waiting for GPS 3D fix" and a refused arm, with nothing to say
  why, because the denial was applied by a process that had already exited. The
  three denial parameters appear in `cns_demo.parm` only as COMMENTS, so
  `--add-param-file` does not restore them. Launch with `-w`.
* **ArduPlane runs in its own process group** via `run_in_terminal_window.sh`,
  so killing `sim_vehicle`'s group does not reliably stop it. A leaked instance
  holds 5760 and the next run fails to start.

Also: `--location -46.6,...` fails argparse, because a value beginning with "-"
is indistinguishable from a flag unless it parses as a plain number, and
"lat,lon,alt,heading" does not. Every southern-hemisphere departure hits this;
the northern ones never did.

## The horizon sensor comparison that was not one

`HorizonConfig` default-constructs with ONE camera -- a FLIR Boson 640 -- not an
empty vector, and nothing cleared it. `pipeline.hpp` says "Empty `cameras`
disables it, which is the default"; that described the intent, not the code.

Consequences, all discovered after several hours of comparisons:

* Every `--horizon-compare` run was **Boson 640 vs Lepton 2.5**, not horizon
  versus none.
* **`--horizon` was a downgrade**: 0.06 arcmin tilt floor to 5.13. The
  per-frame fix went 10.04 -> 16.46 km, exactly as the resolution predicts.
* `pipeline.cpp`'s `cameras.empty() ? 6000.0 : 2700.0` fix sigma **never took
  the 6000 branch**, so both arms were priced identically.

An entire experiment was built on top of that: noticing that horizon fixes were
better (9.65 vs 11.05 km median) while the filtered track was worse, I
attributed it to the filter over-trusting them at 2700 against 6000, and
"re-priced" one arm to 5240. Both arms were at 2700. The ~1 km deltas it
produced have no mechanism behind them and sit inside the run-to-run spread.
`transit_scenario` is not affected -- it calls `hc.cameras.clear()` first --
which is why its measured table showed the Lepton helping all along.

## The statistic that fooled me three times in one session

`dr_err_m` is written only on fix rows, and the aided error is a SAWTOOTH: it
grows along a leg and collapses at each fix orbit. The final value is therefore
wherever the flight happened to stop in that cycle, and it swings by an order of
magnitude.

| quoted | what it was | what the aggregate said |
|---|---|---|
| "horizon 4.00 vs 10.81 km" | last sample | mean 11.62 vs 10.17 -- horizon WORSE |
| "re-priced is good at 3.67 km" | last sample, off a plot title | median 6.97 vs 7.79 |
| "Snares: 1.99 vs 16.81 km" | median of a 6-9 sample final window | endpoints 12.35 vs 20.04 -- both MISS |

Fixed in `report()` and in `live_view.py`'s titles, which both quote median and
RMS now -- and then walked into it a third time in prose, after fixing it twice
in code. The lesson is not "fix the tool", it is that a number which looks like
a triumph should be recomputed a second way before it is said out loud.

## Retracted this session

* "The horizon sensor wins on the filtered result" -- from the last-sample trap.
* "Three runs show the same pattern" -- two of the three were the same artifact,
  and the data to recheck one had been overwritten.
* "Sensor quality washes out after filtering" -- true for Canberra (-11%) and
  Porto Santo (-21%), false for the Snares (-41%). The pattern is that the
  horizon matters most where the dead reckoning it corrects is worst.
* "Both configurations find the island, regardless" -- true for Porto Santo,
  false for a 3.5 km unlit target.

## Method notes worth keeping

* **Two nodes, one flight.** `--horizon-compare` runs a second `celestial_node`
  on `--serial2`, so both arms consume identical telemetry from one aircraft.
  That removes the run-to-run variance the repo already warns about (5.80 vs
  7.53 km on an unchanged configuration) -- but NOT variance in the fix quality
  being fed in, which still differed by 2.3 km between two runs of one config.
* **The autopilot is GNSS-denied too.** With `--no-inject` the fixes never reach
  ArduPilot, so the aircraft flies to where EKF3 believes the waypoint is. Over
  13.7 h that was 18.6 km from the real one. An arrival test keyed on proximity
  to the final waypoint therefore never fires; key it on closest approach.
* **A watcher that stops a run needs an "is it flying yet" guard.** A stagnation
  test fired on a parked aircraft during the pre-arm EKF wait and killed a run
  before it armed.
* **The sky model does not brighten.** `sky_mag_per_arcsec2 = 21.5` is fixed, so
  any flight longer than the night it starts in gets stars it should not have.
  The Porto Santo crossing ran 2.7 h past astronomical dawn.
