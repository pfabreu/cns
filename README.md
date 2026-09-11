# Celestial navigation for a GNSS-denied fixed-wing UAV

A strapdown star camera and air-data dead reckoning, fused to bound absolute
position error indefinitely without emitting anything. One camera, no vertical
reference, no extra sensors.

Built on Teague & Chahl, *"An Algorithm for Affordable Vision-Based GNSS-Denied
Strapdown Celestial Navigation"*, Drones 2024, 8(11) 652 — reproduced, then
corrected and extended in four measurable ways, and validated against a real
ArduPilot EKF3 in SITL rather than only in simulation.

This is a **navigation and simulation library**. It implements the estimator,
the star pipeline and the sensor fusion end to end, and exercises them in a
repeatable simulation, against a MAVLink replayer, and against ArduPilot SITL.
Hardware, daytime operation, sun/moon sighting and the polarisation compass are
deliberately out of scope.

---

## Contents

- [Results](#results)
- [Build and test](#build-and-test)
- [Running it](#running-it)
- [How it works](#how-it-works)
- [The algorithms](#the-algorithms)
- [Sensor models](#sensor-models)
- [Validation against ArduPilot SITL](#validation-against-ardupilot-sitl)
- [What this adds to the paper](#what-this-adds-to-the-paper)
- [Optional extras](#optional-extras)
- [Known limitations](#known-limitations)
- [References](#references)

---

## Results

Two numbers matter, and they answer different questions.

### 1. Does the paper's method reproduce?

```bash
./build/test_orbit          # testPaperReplication
```

Teague & Chahl report **4 km** from one orbit through 360° of compass heading,
on real flight data with a Cube Orange AHRS. Reproduced here with the same class
of error — 0.4° uncalibrated mounting, AHRS tilt bias, drift at the rate their
own Figure 5 measures, maneuver coupling, 5 m/s wind, six seeds:

| orbit | period | mean error |
|---|---|---|
| 150 m radius, 1 rev | 38 s | **4.14 km** |
| 400 m radius, 1 rev | 100 s | 8.80 km |
| 150 m radius, 4 rev | 151 s | **1.96 km** |

A single frame under the same conditions is 43 km, so the orbit buys a factor of
ten.

**The radius is the finding.** A body-fixed error averages away over a heading
sweep; AHRS drift is not body-fixed, so it only averages as `sqrt(2*tau/T)` —
which makes the orbit **period** a design variable. The paper reports 4 km and
never reports a radius; on this model that number needs an orbit short against
the drift correlation time.

### 2. Can it navigate a transit?

```bash
./build/transit_scenario
```

192 km over water, GNSS-denied: a calibrating loiter at departure, then legs
separated by fix orbits. The paper stops at a single fix; this adds dead
reckoning between fixes, a celestial compass for heading, and the estimator
corrections below.

| quantity | result |
|---|---|
| camera boresight, self-calibrated in flight | 0.40 deg -> **0.11 deg** |
| heading, against a 3 deg magnetometer bias | residual **0.10 deg** |
| celestial fixes | 10 over 192 km |
| dead reckoning ALONE | mean 29 km, peak **58 km** (30% of distance) |
| dead reckoning + fixes | mean **7.3 km**, bounded |

Seed spread is 6.1–8.5 km; quote a sweep, not a run. `--seed N` is provided.

**A fix is an orbit.** With no vertical reference a straight-leg fix carries the
whole AHRS tilt error and is worth about 40 km, so the aircraft must maneuver to
navigate. That is the paper's method and it is the operating mode here. The
transit number is dominated by dead reckoning between orbits, not fix quality:
orbit more often and it improves.

### 3. Can it cross an ocean and find an island?

Full SITL, GNSS denied 240 s after takeoff, celestial fixes NOT fed back to the
autopilot — the navigation estimate is a passive observer while ArduPilot flies
the mission on its own dead reckoning.

![Sagres to Porto Santo](docs/images/portosanto-map-dr.png)

883 km, Sagres to Porto Santo, one night in December. Unaided dead reckoning
(red) ends **346 km** out, in open Atlantic. The aided tracks stay on truth.

| crossing | flown | unaided DR | aided, median | aided, p90 |
|---|---:|---:|---:|---:|
| Canberra demo | 154 km | 56 km | 6.1 km | 15.2 km |
| Bluff → the Snares | 262 km | 98 km | 8.2 km | 17.4 km |
| Sagres → Porto Santo | 883 km | 346 km | 6.9 km | 24.7 km |

**Bounded against unbounded is the whole claim**, and it holds at every scale
tested: unaided error grows to 30–48 % of distance flown, the aided estimate
does not grow at all.

Whether that is *enough* depends on the target:

* **Porto Santo** — 11 km long, a town, lights visible ~50 km. Found with
  margin, with or without a horizon sensor.
* **The Snares** — 3.5 km, uninhabited, unlit. **Missed** by both
  configurations (12 km and 20 km at the final fix). This is where the system
  as configured runs out.

Per-crossing numbers, the exposure sweep behind `min_exposure_s`, and the
horizon-sensor comparison are in **[RESULTS.md](RESULTS.md)**.

### Performance

```bash
./build/bench
```

| stage | per frame | max rate |
|---|---|---|
| detect, plain | 5.9 ms | 169 Hz |
| detect, matched filter | 24.8 ms | 40 Hz |
| match to catalogue | 0.10 ms | 10070 Hz |
| per-frame fix (RANSAC) | 1.24 ms | 806 Hz |

The whole flight path fits a 10 Hz camera with 4x margin.

---

## Build and test

```bash
sudo apt install cmake libeigen3-dev liberfa-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build            # 4 suites
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
miss. Without OpenCV, `PriorNet` compiles to a stub and nothing else changes.
Both build paths are tested.

### The test suites

| binary | covers |
|---|---|
| `test_roundtrip` | geometry, sub-stellar points, position solver, RANSAC |
| `test_orbit` | orbit averaging, mount calibration, compass, paper replication |
| `test_imaging` | render, detect, centroid, match, matched filter, mesh background |
| `transit_scenario` | end-to-end smoke test, bounded by CI |

---

## Running it

### Deterministic scenario

```bash
./build/transit_scenario                  # default
./build/transit_scenario --seed 3         # drift makes this matter
./build/transit_scenario --imaging        # real star pipeline, ~50x slower
./build/transit_scenario --fix-revs 4     # longer fix orbits
./build/transit_scenario --rate 10        # the paper's frame rate
./build/transit_scenario --magnetic-heading   # disable the celestial compass
```

### Against a MAVLink replayer

```bash
./build/celestial_node --port 14556 --no-inject --csv live.csv &
./build/fake_sitl --port 14556 --speed 25
```

Expect ~8 fixes, heading-weighted beating naive on each, boresight
self-calibrating from 0.40 deg, and a 25-column CSV. If this does not work, do
not move on to SITL.

**`fake_sitl` runs are NOT reproducible and must not be used for A/B
comparisons.** The node samples a UDP stream, so frame timing varies run to run
and different frames land in each window; two runs differ even with identical
settings (measured: unaided DR 5.80 km against 7.53 km on the same
configuration). It exercises the transport, which the in-process scenario cannot.
For comparing options, use a seeded in-process experiment.

### Against ArduPilot SITL

```bash
sim_vehicle.py -v ArduPlane \
  --add-param-file=$PWD/ardupilot/params/cns_sitl.parm \
  --console --map --out=udp:127.0.0.1:14556
```

`cns_sitl.parm` matters: SITL's default IMU is nearly ideal, EKF3 will hold
attitude far better than a Cube Orange, and the run will prove nothing. It
injects accelerometer bias, gyro noise, vibration and wind.

To deny GNSS mid-flight, once established in a loiter:

```
param set SIM_GPS1_ENABLE 0
param set EK3_SRC1_POSXY 0
param set EK3_SRC1_VELXY 0
```

Set those **after** takeoff. From boot they leave EKF3 with no horizontal
source, so it never aligns and pre-arm fails with "AHRS: not using configured
AHRS type".

**Confirm denial from EKF3's own status bits**, not from GPS fix status:

```bash
python3 -c "
from pymavlink import mavutil
m = mavutil.mavlink_connection('logs/00000001.BIN')
bits=['attitude','horiz_vel','vert_vel','horiz_pos_rel','horiz_pos_abs','vert_pos',
      'terrain_alt','const_pos_mode','pred_horiz_pos_rel','pred_horiz_pos_abs',
      'takeoff_detected','takeoff','touchdown','using_gps','gps_glitching',
      'gps_quality_good','initialized','rejecting_airspeed','dead_reckoning']
seen=set()
while True:
    msg=m.recv_match(type=['XKF4'])
    if msg is None: break
    if msg.SS not in seen:
        seen.add(msg.SS)
        on=[n for i,n in enumerate(bits) if (msg.SS>>i)&1]
        print(msg.SS, '->', ', '.join(on))
"
```

Denial has taken when `using_gps` and `horiz_pos_abs` clear and
`dead_reckoning` sets. GPS fix status and message counts both prove less than
they appear to — ArduPilot keeps logging `GPS` rows with a zero fix.

### Running the whole SITL test unattended

```bash
pip install pymavlink
python3 tools/make_mission.py > ardupilot/missions/demo.waypoints
python3 tools/run_sitl_test.py --ardupilot ~/ardupilot
python3 tools/run_sitl_test.py --ardupilot ~/ardupilot --horizon --speedup 20
```

Starts SITL and `celestial_node`, uploads the mission, arms, flies AUTO, denies
GNSS at `--deny-after` seconds, runs to the end or `--max-minutes`, kills
everything and prints a summary. Then:

```bash
python3 tools/live_view.py --csv live.csv          # trajectories
python3 tools/dump_log.py logs/00000001.BIN > flight.csv
./build/analyse_log flight.csv                     # the real analysis
```

Two things it does that hand-driving does not. **The denial time is a
parameter**, so the analysis can rely on it instead of on whatever moment you
managed to type at. And **everything is killed on exit** — including on Ctrl-C
and on failure — with each child in its own process group, because
`sim_vehicle.py` spawns several processes and killing only the parent leaves
ArduPlane holding the port so the next run fails.

Useful flags:

| flag | what it does |
|---|---|
| `--run-name NAME` | write CSVs and node logs to `runs/<timestamp>-NAME/` instead of the repo root, so one run does not overwrite the last |
| `--location=LAT,LON,ALT,HDG` | move SITL's home. Use the `=` form; a negative latitude looks like a flag otherwise |
| `--utc "YYYY-MM-DDTHH:MM:SS"` | the observation epoch. The camera looks at the ZENITH, so latitude and this decide which sky is overhead — the default is mid-morning at Atlantic longitudes and would give no stars |
| `--horizon-compare` | run TWO nodes on the same flight, one with a horizon sensor and one without. The only sound way to A/B in SITL: both see identical telemetry |
| `--max-minutes N` | wall-clock cap, default 40 — a long crossing needs it raised |

Logs land in `/tmp/sitl.log` and `/tmp/node.log`, or in the run directory under
`--run-name`. If arming fails, the pre-arm message is in the first.

### A trajectory that demonstrates the claim

```bash
python3 tools/make_mission.py > ardupilot/missions/demo.waypoints
```

96 km: a 3-turn calibration loiter, then eight 12 km legs each ending in a
2-turn fix orbit at 250 m. In MAVProxy:

```
param set SIM_SPEEDUP 10          # 84 min of flight in ~8 min
param set BATT_MONITOR 0          # no low-battery failsafe
param set WP_LOITER_RAD 250
wp load ardupilot/missions/demo.waypoints
mode auto
arm throttle
# once established, after the calibration loiter:
param set SIM_GPS1_ENABLE 0
param set EK3_SRC1_POSXY 0
param set EK3_SRC1_VELXY 0
```

**Why this shape.** The claim is not that celestial navigation is accurate — over
a few minutes plain dead reckoning beats it. The claim is that DR error grows
*without limit* while the fix *bounds* it, so the trajectory has to be long
enough for the two curves to separate visibly. Unaided DR runs at roughly 30% of
distance travelled, so at ~7 km of bounded error the crossover is near 25 km and
you want several times that. Over 96 km, unaided DR should reach ~29 km against
a bounded ~7 km — a 4x separation that is obvious on the plot.

**Why leg-then-orbit.** With no vertical reference a straight-leg fix carries the
whole AHRS tilt error and is worth ~40 km. The aircraft has to maneuver to
navigate, so fixes happen in the orbits and the legs are where DR drifts.

The orbit numbers are the measured optima: 250 m is best at both the GPS-aided
and GNSS-denied maneuver coupling, so the choice does not depend on which value
is right, and the departure loiter gets an extra turn because the mounting is
estimated **once** there and everything downstream inherits it. The legs fan by
4° each so heading varies across the flight rather than repeating.

### Log analysis

```bash
python3 tools/dump_log.py logs/00000001.BIN > flight.csv
./build/analyse_log flight.csv
./build/analyse_log flight.csv --imaging     # add the real star pipeline
```

`analyse_log` decomposes the EKF3 attitude error into the part a heading sweep
removes and the part it cannot, and predicts the resulting fix error. It needs
`SIM` messages, so it is SITL-only.

### Optional: a horizon sensor

An LWIR camera looking forward gives a **non-inertial vertical reference** — it
observes the AHRS tilt error directly, and unlike an accelerometer it is not
confused by the acceleration of a turn, which is exactly when the orbit needs it
most.

```bash
./build/celestial_node --port 14556 --no-inject --horizon        # 1x Lepton 2.5
./build/celestial_node --port 14556 --no-inject --horizon-pair   # fore/aft pair

# A/B it deterministically -- run the scenario twice and overlay the tracks
./build/transit_scenario --seed 1 --csv plain.csv
./build/transit_scenario --seed 1 --horizon lepton --csv horizon.csv
python3 tools/live_view.py --csv plain.csv --compare horizon.csv \
    --label "no horizon" --compare-label "1x Lepton"
```

Compare in the SCENARIO, not the live node: `fake_sitl` samples UDP and its runs
are not reproducible, so an A/B there measures timing jitter. The scenario is
seeded and deterministic, and the viewer draws each run's fixes in the same
colour as its trajectory.

Measured, one revolution at 250 m, mount calibrated, 8 seeds:

| configuration | tilt error | orbit fix |
|---|---|---|
| none | 0.261 deg | 6.12 km |
| **1x FLIR Lepton 2.5** | **0.078 deg** | **1.34 km** |
| 2x Lepton 2.5, fore/aft | 0.058 deg | 1.23 km |

Confirmed end to end in SITL, one flight with `--horizon-compare` so both arms
see identical telemetry:

| mission | orbit fix | filtered median |
|---|---:|---:|
| Canberra, 152 km | −29.8 % | −10.8 % |
| Bluff → the Snares, 205 km | −12.9 % | **−41.3 %** |
| Sagres → Porto Santo, 883 km | −33.0 % | −20.5 % |

Fix quality improves ~30 % consistently. The effect on the *filtered* track
varies with how bad the dead reckoning being corrected is — largest on the
205 km leg, where unaided DR reached 48 % of distance flown.

**A short moving average on the tilt measurement is worth 37%** and is the
single cheapest improvement here (2.23 -> 1.34 km). The original code argued no
filter was needed because "the orbit machinery downstream already performs the
time-averaging" — true for a research-grade core at 0.5 px of edge noise, and
false for a Lepton at 1.2 px, which is NOISE-limited: halving its edge noise is
worth 2.2x on the fix.

The window must stay SHORT. Measured, 30 seeds:

| averaging | tilt error | orbit fix |
|---|---|---|
| none | 0.200 deg | 2.31 km |
| 1.0 s | 0.104 deg | 1.48 km |
| **2.5 s** | 0.084 deg | **1.38 km** |
| 4.0 s | 0.075 deg | 1.93 km |

Past ~3 s the tilt error keeps falling while the fix gets WORSE. The average
lags, and during an orbit a lag is a heading-correlated error — exactly the kind
the orbit average cannot remove. **Optimising the tilt number past this point
makes the navigation worse**, which is a useful reminder that tilt is not the
objective.

Note also the fore/aft pair now adds almost nothing (1.34 -> 1.23 km): with
averaging the residual is no longer dominated by the noise a second camera would
help with. The pair earned its keep against anomalous dip, which averaging does
not touch but which is now a smaller share of the total.

#### Choosing a thermal camera

The requirement is set by the ATMOSPHERE, not the sensor. Anomalous dip is about
3 arcmin (0.05 deg) and irreducible, so a camera much better than that is
wasted. And the horizon spans the whole frame, so a line fit across N columns
averages the per-column edge error down by sqrt(N) — which is what lets a very
coarse array work at all.

Estimated tilt precision, assuming 0.15 px sub-pixel edge fit per column:

| sensor | pixels | vFOV | deg/px | est. tilt | verdict |
|---|---|---|---|---|---|
| Panasonic AMG8833 (Grid-EYE) | 8x8 | 60 | 7.50 | 0.40 deg | too coarse |
| Melexis MLX90641 | 16x12 | 35 | 2.92 | 0.109 deg | marginal |
| **Melexis MLX90640** | 32x24 | 35 | 1.46 | **0.039 deg** | at the atmospheric floor |
| **FLIR Lepton 2.5** | 80x60 | 38 | 0.63 | **0.011 deg** | comfortable |
| FLIR Lepton 3.5 | 160x120 | 43 | 0.36 | 0.004 deg | overkill |
| FLIR Boson 320 | 320x256 | 27 | 0.11 | 0.001 deg | research grade |

The **MLX90640** is the interesting one: 32x24 thermopile array, NETD 0.1 K,
I2C, under 23 mA, tens of dollars. Against a sea/sky contrast of tens of kelvin
that sensitivity is enormous margin — you are detecting an EDGE, not measuring
temperature. The **Lepton 2.5** is the safe choice, a couple of hundred dollars
for 4x the margin, and it is what `HorizonCamera::lepton25()` models.

Two caveats on the table. The 0.15 px sub-pixel edge fit is an ASSUMPTION and
drives every row; cheap silicon optics with real distortion and pixel
non-uniformity could halve it, which would push the MLX90640 into "marginal".
And whether sea/sky thermal contrast survives 101 km of atmosphere at 800 m
altitude is untested here — marine thermal imagers see the horizon routinely,
but not at that slant range.

**A single cheap Lepton is worth 4.6x on the orbit fix** in isolation. But see
the next section before buying one: on the full transit it still makes things
worse, because it degrades the celestial compass by more than it improves the
fix.

Note also what it does NOT do: at 16 km a straight-leg fix is still worse than
simply orbiting, so the heading-coverage gate stays and the flight plan is
unchanged. Only a Boson-class sensor makes straight-leg fixes viable — that is
where the paper-era 3.3 km transit came from.

#### A real bug this exposed: fuseTilt leaked into yaw

`fuseTilt` built its correction as `(dtheta_fwd, dtheta_lat, 0)` — zeroing the
BODY-Z component — and its comment claimed "yaw is untouched". **That is false
whenever the aircraft is banked.** Yaw is rotation about the LOCAL VERTICAL, not
about body z, and the two differ by the bank angle, so a correction about a
horizontal body axis carries a component along the nav vertical and moves the
heading.

Measured with a PERFECT sensor and NO atmosphere, one orbit at 14.3 deg mean
bank: **-0.176 deg of yaw change per frame**, consistently signed. It was
identical with the atmosphere off and with a perfect sensor, which is what
proved it was structural rather than noise or anomalous dip.

Fixed by projecting the correction orthogonal to the local vertical in body
axes. Yaw change is now +0.0012 deg. This is the mirror of the trap
`mountTiltOnly` exists to avoid, and that function's own documentation states it
from the other side: *"Level it is body z; banked it tilts by the bank angle."*
The same mistake was made twice, in opposite directions, in two files.

Effect on the celestial compass, recovering a 3 deg magnetometer bias, 20 seeds:

| | before | after |
|---|---|---|
| no horizon | 0.089 deg | 0.089 deg |
| 1x Lepton | **0.310 deg** | **0.098 deg** |

#### What is still wrong: the mounting calibration

With the yaw leak fixed, everything the horizon touches directly is now better —
and the transit is STILL worse.

Deterministic transit, seed 1, with and without a single Lepton:

| | no horizon | 1x Lepton |
|---|---|---|
| AHRS tilt error | 23.96' | **8.27'** |
| heading residual (3 deg bias) | 0.096 deg | **0.006 deg** |
| per-fix error | 6.72 / 13.38 / 6.32 / 5.37 km | **1.83 / 6.98 / 1.21 / 1.73 km** |
| **camera boresight after calibration** | **0.157 deg** | **0.634 deg** |
| transit, 6 seeds | **6.1 - 8.5 km** | 7.1 - 9.7 km |

Tilt a third of what it was, compass sixteen times better, fixes three to five
times better — and the whole flight worse. The one row going the wrong way is
the MOUNTING CALIBRATION: 0.634 deg is worse than the 0.400 deg it started from.

Isolated by removing the mounting error entirely (`--boresight 0`, nothing to
calibrate): the horizon becomes **neutral**, 8.84 km against 8.85 km over four
seeds. So the loss is specifically `recalibrateMount` operating on
horizon-corrected attitudes, not anything about the fix or the compass.

**Diagnosed.** `recalibrateMount` is not the problem — with a horizon it is
SEVEN TIMES BETTER (0.056 deg against 0.905 with a 3 deg yaw bias present). The
damage comes from `mountTiltOnly`, which runs afterwards.

`mountTiltOnly` is a REPAIR, not a refinement. A yaw bias gets absorbed into the
mounting as clocking, hiding it from the compass, and stripping the vertical
component undoes that. But the absorption only happens when the AHRS attitude
carries a tilt error — with a good vertical reference the calibration comes out
clean, and the repair then removes real mounting tilt instead. Measured, 20
seeds, 3 deg yaw bias, boresight error after calibration:

| | `recalibrateMount` | `+ mountTiltOnly` |
|---|---|---|
| no horizon | 0.905 deg | **0.286 deg** (repaired) |
| 1x Lepton | **0.056 deg** | 0.678 deg (damaged) |

It cannot simply be switched off, because that kills the compass — the heading
residual against a 3 deg bias goes 0.006 -> **2.875 deg**, the mounting having
swallowed the bias. So this is a genuine trade-off:

| | compass | boresight |
|---|---|---|
| mountTiltOnly ON | 0.006 deg | 0.678 deg |
| mountTiltOnly OFF | 2.875 deg | 0.056 deg |

Two hypotheses tested and both wrong: it is NOT anomalous dip (the damage is
identical with the atmosphere switched off) and it is NOT the choice of removal
axis (stripping about the camera boresight rather than the mean local vertical
is a no-op, 0.056 vs 0.056). What is actually needed is a calibration that
solves for tilt and clocking JOINTLY, rather than solving for both and
subtracting one afterwards.

**Do not buy the sensor until this is resolved.** The sensor works, and one real
bug in the way was found and fixed; the software still does not know what to do
with a good vertical reference.

Note also that with the mount UNCALIBRATED all three configurations give ~46 km
on a straight leg, because there the 0.4 deg mounting dominates and the horizon
does not touch mounting. The calibration orbit is required either way.

The Lepton's `edge_sigma_px = 1.2` is an **estimate**, not a measurement, and it
drives the whole table. It is the first thing to check against real hardware.

### Fix quality is heading coverage

Without a vertical reference the mounting error and the AHRS tilt bias are both
**body-fixed**, and only averaging over a heading sweep removes them. So fix
quality is a function of coverage, and it is steep. Measured, mount calibrated,
8 seeds:

| sweep | no horizon | 1x Lepton 2.5 |
|---|---|---|
| 36 deg | 27.1 km | 16.2 km |
| 90 deg | 23.1 km | 14.4 km |
| 180 deg | 9.6 km | 9.3 km |
| 270 deg | 4.8 km | **2.6 km** |
| 360 deg | 6.2 km | **2.7 km** |

The full curve, 40 seeds, is more interesting than that fragment. **The first
half revolution buys almost nothing** (18 deg gives 35.2 km, 90 deg gives 28.3),
the error only collapses between 180 and 360 deg, and **whole revolutions are
local minima**:

| revs | sweep | no horizon | 1x Lepton |
|---|---|---|---|
| 0.50 | 180 deg | 12.24 km | 8.65 km |
| 0.75 | 270 deg | 6.99 km | 3.98 km |
| **1.00** | **360 deg** | **6.18 km** | **2.57 km** |
| 1.25 | 450 deg | 10.16 km | 4.73 km |
| 1.50 | 540 deg | 6.79 km | 4.34 km |
| **2.00** | **720 deg** | **3.94 km** | **1.66 km** |
| 2.25 | 810 deg | 6.39 km | 3.23 km |
| **3.00** | **1080 deg** | **3.32 km** | **1.38 km** |

**1.25 revolutions is 64% worse than 1.00 despite covering MORE sky.** A
body-fixed error only cancels if every heading is sampled equally, and a
fractional sweep leaves the over-sampled sector weighted. That is the paper's
Eq. (23) defect arriving from a different direction, and it shows heading
weighting only partly compensates.

`Pipeline` therefore TRIMS the solve window to the most recent whole number of
revolutions whenever it holds at least one. This is free — it discards frames
that were actively hurting. Note the horizon's benefit is also largest at whole
revolutions (2.4x) and smallest on partial arcs (1.4x), because on a partial arc
the residual is dominated by uncancelled body-fixed error, which the horizon
does not touch.

Below one full revolution `Pipeline` still emits a fix — a partial arc is worth
having when dead reckoning has drifted far enough — but it reports an honest
per-fix `sigma_m` from the coverage
(`sigma_full * (12/bins)^0.8`, which fits both columns to ~20%) and lets the
dead-reckoning filter weight it. A 16 km fix is useless against a DR estimate
that is 5 km uncertain and valuable against one that is 40 km uncertain, and the
Kalman update already knows which — it just has to be told the truth about
sigma. Below `min_heading_bins` (4 of 12) the fix is still refused outright,
because there the error is a body-fixed BIAS rather than noise and no sigma
describes it honestly.

### A fix is an orbit

Without a vertical reference the mounting error and the AHRS tilt bias are both
**body-fixed**, and the only thing that removes them is averaging over a heading
sweep. On a straight leg there is no sweep, so the fix carries the whole error —
and because the error is body-fixed while the heading is constant, it lands in a
fixed direction *relative to the aircraft*, which looks like a broken solver.

Measured on a SITL transit:

| fix window | heading bins covered | error |
|---|---|---|
| in an orbit | 12 of 12 | 3.6 – 23.5 km |
| mid-leg | 1 – 2 of 12 | 20.9 – 39.7 km |

The mid-leg fixes were displaced consistently *behind* the aircraft — an
uncalibrated 0.32 deg mounting (35 km) projected through a southward heading.

`Pipeline` therefore emits a fix only when the window covers at least
`min_heading_bins` of 12 (default 8, i.e. 240 degrees). Long straight legs
produce no fixes at all, which is correct: dead reckoning carries you between
orbits, and that is what the orbits are for.

### Fix rejection

A bad fix is not merely useless, it is **corrosive**: it overwrites the position
prior that the matcher and the next solve both depend on, and it can be
recalibrated upon. Left ungated, that runs away — measured on a 65-minute SITL
flight, fixes ranged 3 to 194 km, the mounting was never successfully calibrated
(boresight stuck at 0.400 deg for the entire flight), and dead reckoning **with**
fixes finished worse than unaided: 64.8 km against 48.5 km.

Dead reckoning is the independent check. Over one fix interval it cannot be far
wrong, and the filter's covariance says how far, so a fix that disagrees with DR
by much more than DR's own uncertainty is the thing that is wrong.

The test is **Mahalanobis, not a circle**. Dead-reckoning uncertainty is
genuinely anisotropic — along-track error comes from airspeed scale, cross-track
from heading — so a circular gate either admits bad fixes across the narrow axis
or rejects good ones along the wide one. `DeadReckoner::mahalanobis` measures the
distance in the shape the filter actually believes, with the eigenvalues floored
at 1.5 km so a freshly-updated filter cannot become so confident it rejects
everything. `Pipeline` rejects beyond 4 sigma, and a rejected fix leaves the
prior and the mounting untouched, which is the whole point.

**Note what this can and cannot do.** Gating removes bad fixes; it cannot make
good ones. If the fix scatter is 17 km because of surviving attitude error, a
tighter gate discards the tail without narrowing what remains — and a gate tight
enough to force agreement with DR would simply be reporting DR back, since DR is
itself anchored to earlier fixes.

Before the first fix there is no DR to check against, so the gate falls back to
the departure prior — seeded from the known takeoff position, which is available
even in the GNSS-denied case because you know where you launched.

For the same reason **dead reckoning starts at the departure position, not at
the first fix.** Denial happens in flight, not on the ramp, so starting the DR
at the first celestial fix throws away a position you already know and inherits
the fix's error instead. Measured: a first fix 27 km off put the DR track 27 km
from truth before it had flown anywhere.

### Live view

Two panes, updating during a run: ground truth against the estimate on the
left, and what the camera saw with the matcher's verdict on the right.

```bash
./build/celestial_node --port 14556 --no-inject \
    --csv live.csv --frames-dir frames --frames-every 20 &
./build/fake_sitl --port 14556 --speed 25        # or ArduPilot SITL
python3 tools/live_view.py --csv live.csv --frames frames
```

The **left** pane draws four things and the distinction is the point: truth, the
corrected track (dead reckoning with fixes applied), the **open-loop control**
(the same filter with no fixes at all, dashed red) and the fixes themselves.

The gap between the orange and red tracks IS the result. One grows without
limit, the other does not, and watching orange get pulled back at each fix says
more than any individual fix being accurate. Note that over a SHORT flight red
will beat orange — dead reckoning is good over minutes and only fails over
hours, so the two curves need distance to separate.

The **right** pane shows the rendered frame with green for detections the
matcher identified against the catalogue and red for those it rejected.
Rejections are not failures: the ambiguity guard drops anything with a second
plausible candidate nearby, which is why a noisy detector costs little
downstream.

`--frames-dir` writes a PGM and a detections CSV every Nth frame, renamed into
place so the viewer never reads a half-written file. Off by default; at 10 Hz,
dumping every frame is 47 MB/s. Needs matplotlib, which nothing else here does.

### Other tools

```bash
./build/bench                                  # timing against the frame budget
./build/visualise match out/                   # match diagrams as SVG
./build/visualise frames out/ 120              # star frames as PGM
./build/gen_dataset data/train --frames 500 --seed 1
```

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

**The orbit is the estimator, not a flight-planning detail.** Without a heading
sweep the AHRS tilt error passes straight through and a fix is worth ~43 km.

**The same Kabsch rotation yields three products** — the position fix, the
mounting calibration and the heading — which is why the celestial compass costs
almost nothing to add.

**Mounting and heading split off before the average**, because they are
per-frame quantities while the fix is not.

### Code layout

| module | job |
|---|---|
| `types` | `Geodetic`, `Epoch`, units, haversine. **Frame conventions live here** |
| `star_catalog` | Hipparcos subset, proper motion, magnitudes, common names |
| `sky_model` | RA/Dec to ECEF/NED, refraction, sub-stellar points, **and the position solver** (least squares + RANSAC) — the solver operates on the `StarSight` values this module defines |
| `attitude` | DCM/Euler, `kabsch`, `averageRotation`, **and `StarAidedAttitude`**, the gyro-bias MEKF built from those helpers |
| `imaging` | render, detect, centroid, match; matched filter; mesh background |
| `orbit` | the estimator: orbit fix, mount calibration, compass, simulators |
| `deadreckon` | 4-state air-data filter (N, E, wind_N, wind_E) |
| `prior_net` | optional: ONNX detection prior via `cv::dnn` |

Tools: `pipeline` (nav logic, no MAVLink), `celestial_node` (live UDP node),
`fake_sitl` (replayer), `analyse_log`, `visualise`, `bench`, `gen_dataset`.

**The celestial core is independent of dead reckoning.** `deadreckon` consumes
fixes; nothing in the core includes it. Delete it and the fix still works.

---

## The algorithms

**Position solve** — plane-intersection weighted least squares, SVD-solved (the
paper's Eq. 7–10): each star gives a plane through Earth's centre and the zenith
is where they intersect. Wrapped in RANSAC with a 3-star minimal set (the
paper's Algorithm 1), residuals normalised by `sin(zenith)` so the tolerance
means an angle, and deterministic via a fixed seed.

**Attitude** — Kabsch/Wahba SVD for the camera-to-NED rotation from matched
pairs; chordal rotation averaging for combining per-frame mount estimates.

**Orbit averaging**, three variants — naive mean of zenith vectors (Eq. 23),
heading-weighted mean, and small-circle fit whose axis is the position and whose
angular radius is the misalignment.

**Star pipeline** — median + MAD robust thresholding, connected-component
labelling, Gaussian-fit centroiding (5-parameter Gauss-Newton, verified against
the Cramer-Rao bound), nearest-neighbour matching with a 2x ambiguity guard.
Optional matched filter and SExtractor-style mesh background, below.

**Sky model** — ERFA precession/nutation/sidereal time, proper motion, Bennett
and Saemundsson refraction with elevation-dependent weighting.

**Dead reckoning** — 4-state Kalman filter, air-data driven, with a wind random
walk that absorbs systematic airspeed and heading error.

### Matched filter: use the gyro instead of learning the streak

A star under motion blur is a streak of **known** shape — the AHRS gives the
body rate, so direction and length are predictable per pixel, and the optimal
linear detector for a known signal in Gaussian noise is a matched filter.
Matched stars per frame, which is the figure that feeds the fix:

| exposure | smear | plain | matched filter |
|---|---|---|---|
| 20 ms | 4.8 px | 4.2 | **11.8** |
| 50 ms | 12.1 px | 9.2 | **46.8** |
| 100 ms | 24.2 px | 13.8 | **53.8** |
| 200 ms | 48.3 px | 28.4 | **76.0** |

Two subtleties. The smear is **not uniform**: in an orbit the dominant rate is
about the boresight, which rotates the field rather than translating it, so the
rotational optical flow is evaluated per pixel. And the filter runs on a
**decimated** image, with the factor chosen per frame to keep the decimated
smear above 8 px — that took it from 312 ms to 24.8 ms with no measurable
accuracy cost. Centroiding is always at full resolution.

### Mesh background

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

---

## Sensor models

Everything the simulator injects, and where the numbers come from.

| term | value | source |
|---|---|---|
| camera | 1936x1216, 53.5 deg FOV, 10 Hz | the paper's Alvium 1800 U-240 |
| boresight error | 0.4 deg, uncalibrated at departure | assumed |
| AHRS tilt bias | 0.25 / 0.15 deg | assumed |
| AHRS tilt drift | 0.15 deg, tau = 60 s | the paper's Figure 5 |
| maneuver coupling | 0.0041 | **measured in SITL, GNSS-denied** |
| centroid noise | 10 arcsec | assumed |
| wind | 5 m/s | the paper's flight conditions |

**Maneuver coupling** deserves a note. An accelerometer senses specific force,
so in a coordinated turn its vertical is pulled toward the aircraft's own down
axis and the AHRS under-reads the bank. `ErrorModel::ahrs_turn_coupling` is the
fraction that survives EKF3's compensation. Measured: **0.0009 GPS-aided,
0.0041 denied** — denial costs 4.6x, because centripetal compensation needs
velocity and airspeed is a poorer substitute than GPS. It is not linear
(~0.002 at 10–15 deg of bank, ~0.007 at 25 deg), the denied fit only reached
25 deg, and the top bin dominates it — treat 0.004 as the middle of a
0.0015–0.007 range.

The **orbit radius** is a genuine trade-off. A fix wants a short period, because
drift only averages out if the orbit is quick against its correlation time; a
calibration wants low bank, because maneuver coupling puts a bank-proportional
tilt into the AHRS that `recalibrateMount` books into the mounting. At the
measured coupling: 150 m gives 8.78 km, **250 m gives 7.27 km**, 400 m gives
8.76 km, 600 m gives 10.23 km. 250 m is the optimum at both the aided and denied
coupling, so the choice is not sensitive to which is right.

---

## Validation against ArduPilot SITL

Not just simulation. `analyse_log` on SITL flights, GPS-aided and then denied
mid-flight (confirmed by EKF3 reporting `dead_reckoning` and clearing
`using_gps` and `horiz_pos_abs`).

**The paper's central assumption is false as stated.** It assumes the AHRS
attitude error is body-fixed and therefore averages away over a heading sweep.
Measured, it is not: there is a clean one-per-revolution component, peak-to-peak
0.47 deg in pitch, with consistent phase across revolutions. Part of the error
survives.

**But the surviving part is measurable.** The NED-mean horizontal component
against the measured fix error, six aided orbits: ratios 1.14 / 1.08 / 1.05 /
1.02 / 1.11 / 0.79, **mean 1.03**. So the claim becomes checkable in flight
rather than assumed — *the error decomposes, the orbit removes the body-fixed
part, and what survives predicts the residual.*

**Denial roughly doubles the irreducible error**: NED-mean horizontal 0.085 ->
0.158 deg, and position error 3.2–10.2 km aided against 4.0–15.9 km denied.

| estimator, GNSS-denied, 4 orbits | mean | worst |
|---|---|---|
| naive, Eq. (23) | 10.37 km | 15.21 km |
| heading-weighted | 10.36 km | 15.90 km |
| **circle fit, 1 pass** | **6.93 km** | **9.12 km** |

**The circle fit must not be iterated.** `iterateOrbit` recalibrates the
mounting between passes, which is right for the averaging methods and wrong
here: the fitted radius **is** the misalignment, so recalibrating from the fit's
own output applies the same correction twice. Across eight real orbits,
iteration 1 gives 6.64 km and "converged" gives 9.37 km, and iterating helps in
2 of 8. Fixed, and confirmed on the logs: 9.14 -> 6.93 km denied, 67.53 -> 31.68
aided.

---

## What this adds to the paper

**1. Orbit period is a design variable.** The paper reports 4 km and never a
radius. The averaging only removes body-fixed error, and AHRS drift is not
body-fixed, so it averages as `sqrt(2*tau/T)`. 150 m gives 4.14 km, 400 m gives
8.80 km, everything else equal. Their own Figure 5 supplies the drift rate, so
this follows from their data.

**2. Their conclusion about GPS is an estimator artefact.** They observe 18.27
km for a GPS-guided orbit against 2.29 km flying fixed attitude and conclude GPS
guidance is harmful. It is Eq. (23): an arithmetic mean assumes uniform sampling
in heading, and a GPS-guided track in wind does not provide it. Weighting by
heading increment gives 4.03 -> **0.11 km** with no change to the flight plan.
Caveat: that win is masked once realistic AHRS drift is modelled (6.94 -> 5.86
km) — it fixes a sampling bias, and drift dominates. Both bounds are asserted in
`testHeadingWeightedMean`.

**3. Navigation, not just a fix.** The paper stops at one orbit. This adds dead
reckoning between fixes, a celestial compass recovering a 3 deg magnetometer
bias to 0.10 deg from the same Kabsch rotation, and the fusion that bounds error
over a 192 km transit at 7.3 km against 29 km unaided.

**4. Calibrate once, at low bank.** Maneuver coupling puts a bank-proportional
tilt into the AHRS which `recalibrateMount` books into the mounting.
Recalibrating at every fix orbit made the boresight wander between 0.11 and
0.46 deg. The mounting is a bracket; estimate it once, where the bank is lowest.

Also: the star pipeline uses Gaussian-fit centroiding rather than the paper's
3x3 centre-of-gravity window, which is 7–10x better and reaches the Cramer-Rao
bound; and the detector adds a gyro-informed matched filter, where the paper's
successor reaches for a neural network instead.

---

## Optional extras

Both off by default, both tested in their absent state.

**Star-aided attitude** (`attitude.hpp`, and `analyse_log --star-aided`) — a 6-state multiplicative EKF
estimating gyro bias from absolute star attitude. Star directions in ECEF do not
depend on position, so the *sequence* of star attitudes observes gyro bias even
though a single fix cannot separate tilt from position. Verified in simulation
only: unaided attitude drifts to 0.146 deg peak over 302 s, aided holds 0.042
deg.

Now also runnable on a real log — `dump_log.py` carries the gyro and
`analyse_log --star-aided` recomputes the fix with the corrected attitude.
Measured, it improves the ATTITUDE by 2.4x (0.95 -> 0.40 deg) and leaves the FIX
essentially unchanged, because converting the star attitude from ECEF to NED
needs a position and an error there enters as a constant NED rotation that the
orbit cannot average away. See NOTES-private.md.

**Learned detection prior** (`prior_net.hpp`) — an ONNX model run through
`cv::dnn`, producing a **prior** rather than a decision. It lowers the detection
threshold between `threshold_k` and `threshold_k_low`; the detection still has
to clear a real significance floor from actual photons, and the centroid is
always measured on the original image. So it **cannot hallucinate** a star
(confidence over empty sky yields nothing) and **cannot delete** one (a zero
prior leaves the nominal threshold untouched). `testPriorIsBounded` asserts the
exact equivalence: prior = 1 everywhere is identical to hand-setting the low
threshold, prior = 0 is identical to the ordinary detector. The worst an
adversarial model can do is move you to an operating point you could have chosen
yourself.

7.3 ms at 1/4 resolution, warmed up at load so the first frame does not pay the
165 ms lazy-init cost. `scripts/` has the training stub and the dataset
generator.

---

## Known limitations

* **`ahrs_drift_sigma` and `ahrs_drift_tau` are not measured for this airframe.**
  They come from the paper's Figure 5. Of the three parameters this project
  invented, two have since been measured and one was 20x off, so treat the third
  with suspicion.
* **The simulator does not reproduce the one-per-revolution nav-frame error**
  the real EKF3 shows. Consequence: the simulator and SITL **disagree about
  which estimator is better** — the simulator prefers heading-weighted (7.36 vs
  7.59 km), the real data prefers the circle fit (6.93 vs 10.36 km). The real
  data is taken as authoritative and the scenario default left unchanged. Adding
  a nav-frame error term is the next thing to fix in the model.
* **Timestamp skew between camera and AHRS is not modelled.** At 2.4 deg/s,
  10 ms is 1.4 arcmin — larger than several terms already budgeted.
* **The transit assumes perfect star identification by default.** `--imaging`
  runs the real pipeline and costs +1% with the matched filter, but the
  simulator has no cloud sharp enough, and no real night-sky footage has been
  tested against.
* **EKF3 behaviour under denial is specific to this SITL configuration**, whose
  IMU noise and bias are whatever `ardupilot/params/cns_sitl.parm` sets rather
  than a Cube Orange's. The mechanism generalises; the magnitudes do not.
* **The sky is modelled as permanently dark.** `sky_mag_per_arcsec2 = 21.5` is a
  constant, with no sun and no twilight, so a flight longer than the night it
  departs in will be given stars it could not really see. The 883 km crossing
  runs 13.7 h and overruns its night by ~2.7 h; its night-window figures are
  quoted separately in RESULTS.md for that reason.
* **Small unlit targets are beyond it.** A 3.5 km uninhabited island was missed
  by 12 km with a horizon sensor and 20 km without, after 205 km of flight. An
  11 km island with a town on it was found with margin. Somewhere between those
  is the limit, and it has not been mapped.
* **Fixes are not fed back to the autopilot.** Every result here is a passive
  estimate; ArduPilot flew on its own GNSS-denied dead reckoning throughout, and
  over 13.7 h its own navigation drifted 18.6 km. Closing that loop
  (`--inject`) is implemented but was not part of this campaign.

---

## References

Teague, S. & Chahl, J. *An Algorithm for Affordable Vision-Based GNSS-Denied
Strapdown Celestial Navigation.* Drones 8(11) 652, 2024.

Teague, S. & Chahl, J. *Star detection for low-altitude atmospheric star
trackers.* J. Opt. Soc. Am. A 43(9) 1502–1517, 2026. The successor: a UNet plus
an inter-frame keypoint detector, reporting 85.9% star identification on real
flight data.

Yang, H. et al. *Graduated Non-Convexity for Robust Spatial Perception.* RA-L
2020. Implemented, measured, and **not** used — RANSAC beats it by an order of
magnitude on gross outliers here.

Bertin, E. & Arnouts, S. *SExtractor: Software for source extraction.* A&AS 117,
1996. The mesh background follows its approach.

Delabie, T. et al., on Gaussian-fit centroiding — cited by the paper, which then
uses a 3x3 centre-of-gravity window instead.
