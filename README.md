# Celestial Navigation System

A GNSS-denied celestial navigation system for fixed-wing aircraft, inspired by
the SR-71's Astro-Inertial Navigation System and proven out in ArduPilot SITL.
It is built for flying over terrain with no visual features to fix against —
ocean, desert, ice.

A zenith-pointing strapdown star camera and air-data dead reckoning, fused so
that absolute position error stays **bounded indefinitely** while the aircraft
emits nothing. One camera, no vertical reference, no extra sensors.

This is a **navigation and simulation library**; hardware is out of scope. It
implements the estimator, the star pipeline and the sensor fusion end to end,
and exercises them in a repeatable simulation, against a MAVLink replayer, and
against ArduPilot SITL.

![Sagres to Porto Santo, 883 km](docs/images/portosanto-map-dr.png)

**883 km across open Atlantic, GNSS denied.** Red is dead reckoning with no
celestial fixes — it ends **346 km** out. The aided track holds a 6.9 km median
and finds Porto Santo, an 11 km island, with margin.

![Star field at one celestial fix](docs/images/starmap.png)

**What the camera actually matches.** The zenith-pointing Alvium's
53.5° × 35.1° field, reconstructed from the same Yale catalogue the matcher
uses. Green rings are bright enough to identify (V ≤ 4.0), pale dots are
detectable but too faint. 36 were matched in flight and that fix came out
6.74 km from truth.

---

## Contents

- [How it works](#how-it-works)
- [Build and test](#build-and-test)
- [Running it](#running-it)
- [Results](#results)
- [Known limitations](#known-limitations)
- [References](#references)

---

## How it works

```
                        +------------- ORBIT (360 deg of heading) ----------+
                        |  body-fixed error traces a circle and averages out |
                        +------------------------+---------------------------+
                                                 |
 star camera -> detect -> centroid -> match --+-> per-frame fix -> AVERAGE --+
                                              |   (~43 km)                   |
 IMU / EKF3 attitude -------------------------+                              |
   (the ONLY vertical reference)              |                              v
                                              |                      absolute fix
                                              |                        (~4 km)
                                              |                              |
                                              +-> Kabsch --+-> mounting -----+
                                                           |  (calibrate     |
                                                           |   ONCE, low bank)
                                                           |                 v
                                                           +-> heading -> dead reckoning -> trajectory
                                                              (celestial        ^    |
                                                               compass)         |    |
 airspeed + wind ---------------------------------------------------------------+    |
                                                           +---- dr_ne ------------- +
                                                    (displacement only, not a dependency)
```

Three things to read off it.

**The orbit is the estimator, not a flight-planning detail.** With no vertical
reference, the mounting error and the AHRS tilt bias are both *body-fixed*, and
the only thing that removes them is averaging over a heading sweep. On a
straight leg there is no sweep, so a fix carries the whole error and is worth
about 43 km. Fly a circle and it is worth about 4 km. **The aircraft must
manoeuvre in order to navigate**, which is why the mission is legs separated by
fix orbits: the orbits make fixes, the legs are where dead reckoning drifts.

**The same Kabsch rotation yields three products** — the position fix, the
mounting calibration and the heading — which is why the celestial compass costs
almost nothing to add.

**Mounting and heading split off before the average**, because they are
per-frame quantities while the fix is not.

### Code layout

| module | job |
|---|---|
| `types` | `Geodetic`, `Epoch`, units, haversine. **Frame conventions live here** |
| `star_catalog` | Yale BSC5 subset, proper motion, magnitudes, common names |
| `sky_model` | RA/Dec to ECEF/NED, refraction, sub-stellar points, **and the position solver** (least squares + RANSAC) |
| `attitude` | DCM/Euler, `kabsch`, `averageRotation`, and `StarAidedAttitude`, a gyro-bias MEKF |
| `imaging` | render, detect, centroid, match; matched filter; mesh background |
| `orbit` | the estimator: orbit fix, mount calibration, compass, simulators |
| `deadreckon` | 4-state air-data filter (N, E, wind_N, wind_E) |
| `prior_net` | optional: ONNX detection prior via `cv::dnn` |

Tools live in `tools/`: `pipeline` (nav logic, no MAVLink), `celestial_node`
(live UDP node), `fake_sitl` (replayer), `analyse_log`, `visualise`, `bench`,
`gen_dataset`. See [tools/README.md](tools/README.md).

**The celestial core is independent of dead reckoning.** `deadreckon` consumes
fixes; nothing in the core includes it. Delete it and the fix still works.

### The algorithms

**Position solve** — plane-intersection weighted least squares, SVD-solved:
each star gives a plane through Earth's centre and the zenith is where they
intersect. Wrapped in RANSAC with a 3-star minimal set, residuals normalised by
`sin(zenith)` so the tolerance means an angle, deterministic via a fixed seed.

**Attitude** — Kabsch/Wahba SVD for the camera-to-NED rotation from matched
pairs; chordal rotation averaging for combining per-frame mount estimates.

**Orbit averaging**, three variants — naive mean of zenith vectors,
heading-weighted mean, and a small-circle fit whose axis is the position and
whose angular radius is the misalignment.

**Star pipeline** — median + MAD robust thresholding, connected-component
labelling, Gaussian-fit centroiding (5-parameter Gauss-Newton, verified against
the Cramér-Rao bound), nearest-neighbour matching with a 2x ambiguity guard.

**Sky model** — ERFA precession/nutation/sidereal time, proper motion, Bennett
and Saemundsson refraction with elevation-dependent weighting.

**Dead reckoning** — 4-state Kalman filter, air-data driven, with a wind random
walk that absorbs systematic airspeed and heading error. Fixes enter as position
measurements and are gated on a Mahalanobis distance against the filter's own
covariance, so a bad fix cannot overwrite a good prior.

#### Matched filter

A star under motion blur is a streak of **known** shape — the AHRS gives the
body rate, so direction and length are predictable per pixel, and the optimal
linear detector for a known signal in Gaussian noise is a matched filter.
Matched stars per frame:

| exposure | smear | plain | matched filter |
|---|---|---|---|
| 20 ms | 4.8 px | 4.2 | **11.8** |
| 50 ms | 12.1 px | 9.2 | **46.8** |
| 100 ms | 24.2 px | 13.8 | **53.8** |
| 200 ms | 48.3 px | 28.4 | **76.0** |

The smear is **not uniform**: in an orbit the dominant rate is about the
boresight, which rotates the field rather than translating it, so the rotational
optical flow is evaluated per pixel. The filter runs on a **decimated** image,
the factor chosen per frame to keep the decimated smear above 8 px — 312 ms down
to 24.8 ms with no measurable accuracy cost. Centroiding is always at full
resolution.

#### Mesh background

A SExtractor-style background: sigma-clipped estimate per 64 px cell,
median-filtered across nodes, bilinearly interpolated, subtracted
(`DetectorConfig::bg_mesh_px`, off by default).

| condition | plain | mesh | gain |
|---|---|---|---|
| clear | 38.3 | 36.8 | 0.96x |
| cloud 0.6 | 25.8 | 26.2 | 1.01x |
| flare 0.5 | 19.3 | 34.0 | **1.76x** |

Note *which* contaminant it fixes. Lens flare is a steep additive gradient and
the mesh recovers nearly all of it. Cloud it does not help with at all, because
cloud's damage is not threshold bias: it attenuates starlight and its glow
raises the **shot noise** floor. Subtracting a background cannot un-add Poisson
noise or recover photons that never arrived.

### Sensor models

Everything the simulator injects, and where the numbers come from.

| term | value | source |
|---|---|---|
| camera | 1936x1216, 53.5 deg FOV, 10 Hz | Alvium 1800 U-240 |
| boresight error | 0.4 deg, uncalibrated at departure | assumed |
| AHRS tilt bias | 0.25 / 0.15 deg | assumed |
| AHRS tilt drift | 0.15 deg, tau = 60 s | the paper's Figure 5 |
| manoeuvre coupling | 0.0041 | **measured in SITL, GNSS-denied** |
| centroid noise | 10 arcsec | assumed |
| wind | 5 m/s | the paper's flight conditions |

**Manoeuvre coupling** deserves a note. An accelerometer senses specific force,
so in a coordinated turn its vertical is pulled toward the aircraft's own down
axis and the AHRS under-reads the bank. `ErrorModel::ahrs_turn_coupling` is the
fraction that survives EKF3's compensation. Measured: **0.0009 GPS-aided,
0.0041 denied** — denial costs 4.6x, because centripetal compensation needs
velocity and airspeed is a poorer substitute than GPS.

The **orbit radius** is a genuine trade-off. A fix wants a short period, because
drift only averages out if the orbit is quick against its correlation time; a
calibration wants low bank, because manoeuvre coupling puts a bank-proportional
tilt into the AHRS. At the measured coupling: 150 m gives 8.78 km, **250 m gives
7.27 km**, 400 m gives 8.76 km, 600 m gives 10.23 km. 250 m is the optimum at
both the aided and denied coupling, so the choice is not sensitive to which is
right.

### Optional extras

All off by default, all tested in their absent state.

**Horizon sensor** — an LWIR camera looking forward gives a *non-inertial*
vertical reference: it observes the AHRS tilt error directly and, unlike an
accelerometer, is not confused by the acceleration of a turn, which is exactly
when the orbit needs it most. Worth about 30% on fix quality; see
[Results](#a-horizon-sensor-is-worth-about-30). A short moving average on the
tilt measurement is worth a further 37%, but the window must stay under ~3 s.

**Star-aided attitude** (`analyse_log --star-aided`) — a 6-state multiplicative
EKF estimating gyro bias from absolute star attitude. Star directions in ECEF do
not depend on position, so the *sequence* of star attitudes observes gyro bias
even though a single fix cannot separate tilt from position. On a real log it
improves attitude by 2.4x (0.95 → 0.40 deg) and leaves the fix essentially
unchanged.

**Learned detection prior** (`prior_net.hpp`) — an ONNX model run through
`cv::dnn`, producing a **prior** rather than a decision. It lowers the detection
threshold between `threshold_k` and `threshold_k_low`; the detection still has
to clear a significance floor computed from actual photons, and the centroid is
always measured on the original image. So it **cannot hallucinate** a star and
**cannot delete** one. `testPriorIsBounded` asserts the exact equivalence: the
worst an adversarial model can do is move you to an operating point you could
have chosen yourself. 7.3 ms at 1/4 resolution. See
[scripts/README.md](scripts/README.md).

---

## Build and test

```bash
sudo apt install cmake libeigen3-dev liberfa-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
```

Two optional dependencies, both auto-detected:

```bash
# MAVLink -- needed for celestial_node and fake_sitl
git clone --depth 1 https://github.com/mavlink/c_library_v2 ~/c_library_v2
cmake -S . -B build -DMAVLINK_DIR=~/c_library_v2

# OpenCV -- only for running an ONNX detection prior through cv::dnn
sudo apt install libopencv-dev
```

Without MAVLink those two targets are skipped with a message that is easy to
miss — check `ls build/`. Without OpenCV, `PriorNet` compiles to a stub and
nothing else changes. Both build paths are tested.

| binary | covers |
|---|---|
| `test_roundtrip` | geometry, sub-stellar points, position solver, RANSAC |
| `test_orbit` | orbit averaging, mount calibration, compass, paper replication |
| `test_imaging` | render, detect, centroid, match, matched filter, mesh background |
| `test_horizon` | horizon sensor model and tilt fusion |
| `transit_scenario` | end-to-end smoke test, bounded by CI |

---

## Running it

### Deterministic scenario

Seeded and reproducible. This is the right place to compare options.

```bash
./build/transit_scenario                  # default
./build/transit_scenario --seed 3         # drift makes this matter
./build/transit_scenario --imaging        # real star pipeline
./build/transit_scenario --horizon lepton # with a horizon sensor
./build/transit_scenario --magnetic-heading   # disable the celestial compass
```

Quote a seed sweep, not a run: the spread on an unchanged configuration is
6.1 to 8.5 km.

### Against a MAVLink replayer

```bash
./build/celestial_node --port 14556 --no-inject --csv live.csv &
./build/fake_sitl --port 14556 --speed 25
```

Expect ~8 fixes and a boresight self-calibrating from 0.40 deg. This exercises
the transport; it is not reproducible run to run, so do not A/B with it.

### Against ArduPilot SITL

```bash
sim_vehicle.py -v ArduPlane \
  --add-param-file=$PWD/ardupilot/params/cns_sitl.parm \
  --console --map --out=udp:127.0.0.1:14556
```

`cns_sitl.parm` matters: SITL's default IMU is nearly ideal, EKF3 will hold
attitude far better than a Cube Orange, and the run will prove nothing. It
injects accelerometer bias, gyro noise, vibration and wind.

To deny GNSS, **after** takeoff and once established in a loiter:

```
param set SIM_GPS1_ENABLE 0
param set EK3_SRC1_POSXY 0
param set EK3_SRC1_VELXY 0
```

Set from boot they leave EKF3 with no horizontal source, so it never aligns and
pre-arm fails. Confirm denial from EKF3's own status bits — `using_gps` and
`horiz_pos_abs` clear, `dead_reckoning` sets — not from GPS fix status, which
keeps logging rows with a zero fix.

### The whole SITL test, unattended

```bash
pip install pymavlink
python3 tools/sitl/make_mission.py > ardupilot/missions/demo.waypoints
python3 tools/sitl/run_sitl_test.py --ardupilot ~/ardupilot
```

Starts SITL and `celestial_node`, uploads the mission, arms, flies AUTO, denies
GNSS at `--deny-after` seconds, runs to the end or `--max-minutes`, kills
everything and prints a summary.

| flag | what it does |
|---|---|
| `--run-name NAME` | write CSVs and logs to `runs/<timestamp>-NAME/` so one run does not overwrite the last |
| `--location=LAT,LON,ALT,HDG` | move SITL's home. Use the `=` form; a negative latitude looks like a flag otherwise |
| `--utc "YYYY-MM-DDTHH:MM:SS"` | the observation epoch. The camera looks at the ZENITH, so latitude and this decide which sky is overhead |
| `--horizon-compare` | run TWO nodes on one flight, one with a horizon sensor and one without — the only sound way to A/B in SITL |
| `--max-minutes N` | wall-clock cap, default 40; a long crossing needs it raised |

The mission is a 3-turn calibration loiter, then legs each ending in a fix
orbit at 250 m. The departure loiter gets an extra turn because the mounting is
estimated **once** there and everything downstream inherits it.

### Analysis and plotting

```bash
python3 tools/sitl/dump_log.py logs/00000001.BIN > flight.csv
./build/analyse_log flight.csv              # the real analysis, SITL-only
./build/bench                               # timing against the frame budget

python3 tools/plot/live_view.py --csv live.csv --best 3
python3 tools/plot/plot_map.py live.csv out.png --show-dr
python3 tools/plot/plot_starmap.py live.csv out.png --simple
```

`analyse_log` decomposes the EKF3 attitude error into the part a heading sweep
removes and the part it cannot, and predicts the resulting fix error. It needs
`SIM` messages, so it is SITL-only.

---

## Results

### 1. Does the paper's method reproduce?

```bash
./build/test_orbit          # testPaperReplication
```

Teague & Chahl report **4 km** from one orbit through 360° of compass heading,
on real flight data. Reproduced here with the same class of error — 0.4°
uncalibrated mounting, AHRS tilt bias, drift at the rate their own Figure 5
measures, manoeuvre coupling, 5 m/s wind, six seeds:

| orbit | period | mean error |
|---|---|---|
| 150 m radius, 1 rev | 38 s | **4.14 km** |
| 400 m radius, 1 rev | 100 s | 8.80 km |
| 150 m radius, 4 rev | 151 s | **1.96 km** |

A single frame under the same conditions is 43 km, so the orbit buys a factor
of ten.

**The radius is the finding.** A body-fixed error averages away over a heading
sweep; AHRS drift is not body-fixed, so it only averages as `sqrt(2*tau/T)` —
which makes the orbit **period** a design variable. The paper reports 4 km and
never reports a radius.

### 2. Can it navigate a transit?

```bash
./build/transit_scenario
```

192 km over water, GNSS-denied: a calibrating loiter at departure, then legs
separated by fix orbits.

| quantity | result |
|---|---|
| camera boresight, self-calibrated in flight | 0.40 deg → **0.11 deg** |
| heading, against a 3 deg magnetometer bias | residual **0.10 deg** |
| celestial fixes | 10 over 192 km |
| dead reckoning ALONE | mean 29 km, peak **58 km** (30% of distance) |
| dead reckoning + fixes | mean **7.3 km**, bounded |

The transit number is dominated by dead reckoning between orbits, not by fix
quality: orbit more often and it improves.

### 3. Can it cross an ocean and find an island?

Full SITL, GNSS denied 240 s after takeoff, fixes **not** fed back to the
autopilot — the navigation estimate is a passive observer while ArduPilot flies
the mission on its own dead reckoning.

| crossing | flown | unaided DR | aided, median | aided, p90 |
|---|---:|---:|---:|---:|
| Canberra demo | 154 km | 56 km | 6.1 km | 15.2 km |
| Bluff → the Snares | 262 km | 98 km | 8.2 km | 17.4 km |
| Sagres → Porto Santo | 883 km | 346 km | 6.9 km | 24.7 km |

**Bounded against unbounded is the whole claim**, and it holds at every scale
tested: unaided error grows to 30–48% of distance flown, the aided estimate does
not grow at all.

![Canberra demo, 154 km](docs/images/canberra-map-dr.png)

**Canberra, 154 km.** Unaided dead reckoning ends 56 km out; the aided estimate
holds a 6.1 km median.

![Bluff to the Snares, 262 km](docs/images/snares-map-dr.png)

**Bluff → the Snares, 262 km.** Unaided DR ends 98 km out, aided median 8.2 km.
The Snares are 3.5 km, uninhabited and unlit, and this one is **missed** by
12 km. This is where the system as configured runs out.

Whether the bound is *enough* depends on the target. Porto Santo is 11 km long
with a town on it and was found with margin; the Snares are 3.5 km, unlit, and
were missed. Somewhere between those is the limit, and it has not been mapped.

Per-crossing detail is in **[RESULTS.md](RESULTS.md)**.

### What every fix looks like

![Canberra trajectory, dead reckoning and every fix](docs/images/canberra-trajectory-fixes.png)

Every fix the system produced on one flight. Grey is truth, red dashed is dead
reckoning with no fixes applied, blue is the same filter with fixes fed in, and
each blue star is one of the **37 fixes** it was given.

The fixes are scattered, and that is the honest picture: individual ones land
tens of kilometres out. The estimate is not carried by any single fix being
accurate — it is carried by the fixes being *unbiased*, so the filter is pulled
back toward truth every time one arrives. Median **6.05 km**, RMS 8.77 km,
against 57 km unaided. The three best are gold; the best is **1.2 km**.

The blue track sawtooths — it drifts between orbits and is yanked back at each
fix — which is why every number here is a median or a p90 over a whole flight
rather than an instantaneous value.

### Fix quality is heading coverage

Because the dominant errors are body-fixed, fix quality is a function of how
much heading the orbit sweeps, and it is steep. 40 seeds:

| revs | sweep | no horizon | 1x Lepton |
|---|---|---|---|
| 0.50 | 180 deg | 12.24 km | 8.65 km |
| 0.75 | 270 deg | 6.99 km | 3.98 km |
| **1.00** | **360 deg** | **6.18 km** | **2.57 km** |
| 1.25 | 450 deg | 10.16 km | 4.73 km |
| **2.00** | **720 deg** | **3.94 km** | **1.66 km** |
| **3.00** | **1080 deg** | **3.32 km** | **1.38 km** |

**Whole revolutions are local minima, and 1.25 revolutions is 64% worse than
1.00 despite covering more sky.** A body-fixed error only cancels if every
heading is sampled equally, and a fractional sweep leaves the over-sampled
sector weighted. `Pipeline` therefore trims the solve window to the most recent
whole number of revolutions, which is free — it discards frames that were
actively hurting.

Below one revolution a fix is still emitted, with an honest per-fix `sigma_m`
derived from the coverage so the filter can weight it; below `min_heading_bins`
it is refused outright, because there the error is a bias rather than noise and
no sigma describes it honestly.

### A horizon sensor is worth about 30%

An LWIR camera looking forward, measured over 8 seeds at 250 m:

| configuration | tilt error | orbit fix |
|---|---|---|
| none | 0.261 deg | 6.12 km |
| **1x FLIR Lepton 2.5** | **0.078 deg** | **1.34 km** |
| 2x Lepton 2.5, fore/aft | 0.058 deg | 1.23 km |

Confirmed end to end in SITL with `--horizon-compare`, so both arms see
identical telemetry:

| mission | orbit fix | filtered median |
|---|---:|---:|
| Canberra, 152 km | −29.8% | −10.8% |
| Bluff → the Snares, 205 km | −12.9% | **−41.3%** |
| Sagres → Porto Santo, 883 km | −33.0% | −20.5% |

Fix quality improves ~30% consistently. The effect on the *filtered* track
varies with how bad the dead reckoning being corrected is.

### Validated against ArduPilot SITL

Not just simulation. `analyse_log` on SITL flights, GPS-aided and then denied
mid-flight.

**The paper's central assumption is false as stated.** It assumes the AHRS
attitude error is body-fixed and therefore averages away over a heading sweep.
Measured, it is not: there is a clean one-per-revolution component, peak-to-peak
0.47 deg in pitch, with consistent phase across revolutions. Part of the error
survives.

**But the surviving part is measurable.** The NED-mean horizontal component
against the measured fix error, six aided orbits: mean ratio **1.03**. So the
claim becomes checkable in flight rather than assumed — the error decomposes,
the orbit removes the body-fixed part, and what survives predicts the residual.

**Denial roughly doubles the irreducible error**: NED-mean horizontal 0.085 →
0.158 deg, position error 3.2–10.2 km aided against 4.0–15.9 km denied.

| estimator, GNSS-denied, 4 orbits | mean | worst |
|---|---|---|
| naive | 10.37 km | 15.21 km |
| heading-weighted | 10.36 km | 15.90 km |
| **circle fit, 1 pass** | **6.93 km** | **9.12 km** |

### Performance

| stage | per frame | max rate |
|---|---|---|
| detect, plain | 5.9 ms | 169 Hz |
| detect, matched filter | 24.8 ms | 40 Hz |
| match to catalogue | 0.10 ms | 10070 Hz |
| per-frame fix (RANSAC) | 1.24 ms | 806 Hz |

The whole flight path fits a 10 Hz camera with 4x margin. Real detection and
matching cost about 1% end to end — **the star pipeline is not the constraint.
Attitude is**, and one degree of it is 111 km of position.

---

## Known limitations

* **`ahrs_drift_sigma` and `ahrs_drift_tau` are not measured for this airframe.**
  They come from the paper's Figure 5.
* **The simulator does not reproduce the one-per-revolution nav-frame error**
  the real EKF3 shows, so the simulator and SITL disagree about which estimator
  is better. The real data is taken as authoritative.
* **Timestamp skew between camera and AHRS is not modelled.** At 2.4 deg/s,
  10 ms is 1.4 arcmin — larger than several terms already budgeted.
* **No real night-sky footage has been tested against.** Everything about
  detection under cloud is synthetic and therefore circular.
* **EKF3 behaviour under denial is specific to this SITL configuration.** The
  mechanism generalises; the magnitudes do not.
* **The sky is modelled as permanently dark**, with no sun and no twilight, so a
  flight longer than the night it departs in is given stars it could not really
  see. The 883 km crossing overruns its night by ~2.7 h; its night-window
  figures are quoted separately in RESULTS.md.
* **Small unlit targets are beyond it.** 3.5 km missed, 11 km found; the
  boundary has not been mapped.
* **Fixes are not fed back to the autopilot.** Every result here is a passive
  estimate. Closing that loop (`--inject`) is implemented but was not flown.

### Next

1. Add the one-per-revolution nav-frame error to the simulator.
2. Measure `ahrs_drift_sigma` and `ahrs_drift_tau` on a real airframe —
   `tools/analysis/measure_drift_tau.py` does it from a `.BIN`.
3. Model camera-to-AHRS timestamp skew.
4. Get real night-sky footage from the airframe.
5. Close the loop with `--inject`.
6. Map the small-target limit.
7. Solve mounting tilt and clocking jointly.

---

## References

Teague, S. & Chahl, J. *An Algorithm for Affordable Vision-Based GNSS-Denied
Strapdown Celestial Navigation.* Drones 8(11) 652, 2024.

Teague, S. & Chahl, J. *Star detection for low-altitude atmospheric star
trackers.* J. Opt. Soc. Am. A 43(9) 1502–1517, 2026. The successor: a UNet plus
an inter-frame keypoint detector, reporting 85.9% star identification on real
flight data.

Delabie, T. et al., on Gaussian-fit centroiding — cited by the paper, which then
uses a 3x3 centre-of-gravity window instead.
