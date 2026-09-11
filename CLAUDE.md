# CLAUDE.md

Working notes for an agent operating in this repository. The README describes
what the system does; this describes how to work on it without repeating
mistakes that have already been made and measured.

## What this is

GNSS-denied celestial navigation for a fixed-wing UAV: a strapdown star camera
plus air-data dead reckoning, bounding absolute position error indefinitely. It
reproduces and extends Teague & Chahl, Drones 2024, 8(11) 652.

**The one fact that explains everything else:** 1 degree of attitude error is
about 111 km of position error. Every other subsystem — the star detector, the
centroider, the position solver, the robust estimator, the mounting calibration
— has been measured and found NOT to be the constraint. Attitude is. If a
change does not reduce attitude error, or does not change how attitude error
averages out over a heading sweep, it will not move the headline number.

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build        # 5 suites, all must pass
```

Two optional dependencies, both auto-detected:

```bash
git clone --depth 1 https://github.com/mavlink/c_library_v2 ~/c_library_v2
cmake -S . -B build -DMAVLINK_DIR=~/c_library_v2   # celestial_node, fake_sitl
sudo apt install libopencv-dev                     # PriorNet only
```

Without MAVLink those two targets are skipped with a message that is easy to
miss. Check `ls build/` before concluding a target is broken.

## Running things

```bash
./build/transit_scenario --seed 3          # end-to-end, DETERMINISTIC
./build/test_orbit                         # includes testPaperReplication
./build/bench                              # timing against the frame budget
./build/analyse_log flight.csv             # real ArduPilot logs
python3 tools/run_sitl_test.py --ardupilot ~/ardupilot   # unattended SITL
```

## How to measure things here

These rules were each learned by getting them wrong first.

**Never A/B on `fake_sitl`.** It samples a UDP stream, so frame timing varies
run to run and different frames land in each window. Two runs of the SAME
configuration differed by 5.80 vs 7.53 km. Use `transit_scenario`, which is
seeded and deterministic, or `analyse_log` on a recorded `.BIN`, which
reprocesses deterministically. `fake_sitl` is for exercising the transport
(stream requests, message pairing, timing skew), not for comparing options.

**Quote a seed sweep, not a run.** Six seeds minimum for the transit; the spread
is 6.1 to 8.5 km on an unchanged configuration. Twenty to forty seeds for a
subtle effect. A 28% difference on 8 seeds turned out to be noise sitting next
to a real effect that only appeared at 40.

**A result that is too good is usually circular.** A star-aided attitude
experiment reported 5.79 km improving to 0.01 km. It was converting frames using
the TRUE orbit centre. Ten metres of position from a 0.4 degree uncalibrated
mounting is not a triumph, it is a bug. When a number looks wonderful, find the
line where truth leaked in before reporting it.

**Relink probes after rebuilding the library.** Scratch programs link
`build/libcns.a` statically. Rebuilding the library does not update them, and
you will get identical output and conclude your change did nothing. This
happened twice.

**Isolated probes can invert end to end.** An isolated measurement said a
horizon sensor was worth 2.5x. The full transit said 1.5x WORSE, because two
subsystems were competing for the same information. The probe was not wrong; it
answered a narrower question. Always confirm a win end to end.

**Do not optimise a proxy.** Averaging the horizon tilt measurement keeps
improving tilt error past 3 seconds while making the position fix worse, because
the lag becomes a heading-correlated error. Tilt is not the objective.

## What is settled — do not re-litigate

Each of these was measured. Re-deriving them wastes a session.

* **Tilt and position are ONE observable to a celestial fix** (1 arcmin = 1 nm).
  Four separate attempts to estimate tilt out of the fix have failed against
  this: post-hoc tilt estimation, a drift-aware orbit estimator, star-tracker
  de-drift, and star-aided attitude replacing EKF3. The cleanest statement is in
  `analyse_log`'s `CNS_ORACLE_P` diagnostic: with the TRUE position the
  correction works perfectly (0.00 km); with an estimated one it reproduces the
  estimate. The estimator is right, the information is absent.
* **A fix requires a heading sweep.** Body-fixed errors only cancel over 360
  degrees. Fix quality against sweep: 90 deg gives 28 km, 180 gives 12, 270
  gives 7, 360 gives 6.2. And WHOLE revolutions are local minima — 1.25 revs is
  64% worse than 1.00 despite covering more sky.
* **RANSAC beats GNC-TLS** by an order of magnitude on gross outliers. GNC was
  implemented, measured and removed.
* **The star pipeline is not the constraint.** Real detection and matching costs
  +1% end to end. Centroid noise from 2 to 60 arcsec is invisible.
* **Mounting calibration does not improve an orbit fix.** The orbit already
  removes the body-fixed error that calibration corrects: 6.34 km at every
  boresight from 0.40 deg down to 0.00.
* **The paper's assumption is false as stated.** Real EKF3 attitude error is not
  body-fixed; there is a one-per-revolution component. What survives averaging
  is the NED-mean horizontal, measured at 0.158 deg GNSS-denied.

## Open problems worth working on

1. **`mountTiltOnly` and a horizon sensor conflict.** With a horizon the
   calibration comes out clean (0.056 deg) and `mountTiltOnly` then damages it
   (0.678 deg); switching it off kills the compass (heading residual 0.006 ->
   2.875 deg). Two hypotheses tested and both wrong — it is not anomalous dip,
   and it is not the choice of removal axis. Needs a calibration that solves for
   tilt and clocking JOINTLY rather than solving both and subtracting one.
2. **`ahrs_drift_sigma` and `ahrs_drift_tau` are unmeasured** for any real
   airframe; they come from the paper's Figure 5. Of the three parameters this
   project invented, two have since been measured and one was 20x off.
3. **The simulator lacks the one-per-revolution nav-frame error** the real EKF3
   shows. Consequence: the simulator and SITL DISAGREE about which estimator is
   better. The real data is authoritative.
4. **No real night-sky footage** has ever been tested against. Everything about
   detection under cloud is synthetic and therefore circular.

## Dead ends — do not retry without new information

Full reasoning in `NOTES-private.md`.

| attempt | why it failed |
|---|---|
| Post-hoc tilt estimation (3 variants) | tilt and position are one observable |
| Drift-aware orbit estimator | drift at the orbit frequency is degenerate with position |
| Star-tracker de-drift | comparing ECEF attitude to NED reintroduces position |
| Mid-transit calibration loiters | 2-3x worse; each recalibration writes the instantaneous bias into the mounting |
| GNC-TLS robust solver | lost to RANSAC by 10x |
| Polar-coordinate matched filter | rotation centre lands outside the frame; 6x SLOWER |
| Neural background estimator | the classical mesh already takes the recoverable part |
| Moon parallax for position | needs 0.09 px astrometry on an 18 px disc; phase error is 400 km |

## Style

Comments explain WHY, especially where a number is surprising or a simpler
approach was tried and failed. Several constants in this codebase are load
bearing and non-obvious; if you change one, say what you measured.

State uncertainty explicitly. `ahrs_turn_coupling` was carried for a long time
marked "ASSUMED, NOT MEASURED", which is exactly why it was the first thing
checked against real logs — and it was 20x off.

When a test fails, check the test before the code. Several apparent bugs in this
session were faulty probes: units confused (metres for km), a stale binary, a
`pgrep` pattern matching its own parent.
