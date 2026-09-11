# An error-state filter for the celestial fix — proposal, not a result

Prompted by IEEE Xplore document 11124185,
<https://ieeexplore.ieee.org/document/11124185>. Only the abstract is
accessible; everything below is inference from it and has not been checked
against the paper's method.

> Modern aviation often relies on global navigation satellite systems (GNSS) for
> positioning accuracy. However, today's electronic warfare technology presents
> an integrity risk for such systems. [...] This paper describes a new way to
> solve the aircraft positioning problem through an error-state Kalman filter
> (ESKF) for autonomous astronomical-based navigation by adapting an observation
> model developed for planetary rovers to the context of aviation. Its
> efficiency to reduce inertial navigation system errors is demonstrated through
> Monte Carlo simulations of randomly generated aircraft trajectories.

## Why this is architecturally different from what is here

This project computes a **position** from attitude and star directions, then
feeds that position to a four-state Kalman filter — north, east, wind_n, wind_e
(`deadreckon.hpp`). There is no attitude state anywhere in the estimator; the
AHRS is an input, and the celestial system is a passive observer of it.

An ESKF inverts that: the state is the **INS error**, attitude error included,
and star observations are measurements on it directly. Attitude is estimated
rather than trusted.

## The hypothesis worth testing

`CLAUDE.md` records that tilt and position are one observable, and four attempts
to separate them failed. That is true **per fix** and it is geometry — no filter
changes it. But a filter does not have to separate them within one fix. It can
exploit the fact that the two errors have **different dynamics**:

* position error is driven by the dead-reckoning model — wind, airspeed scale,
  heading bias
* attitude error is driven by gyro drift, a Gauss-Markov process with
  correlation time `ahrs_drift_tau`

Two processes with different spectra can be separable over a trajectory even
when each individual measurement is rank-deficient. This is the standard
argument for why an INS/GNSS filter observes accelerometer bias only during
manoeuvres.

The recorded failure of the drift-aware orbit estimator — *"drift at the orbit
frequency is degenerate with position"* — is consistent with this rather than
contrary to it: its only excitation **was** at the orbit frequency. A filter
running across a whole transit sees orbits separated by straight legs.

## The measurement model, and what the orbit actually buys

A fix's error decomposes as

    fix_error = δp + R_earth · R(ψ) · δθ_body

with ψ the heading. The tilt term rotates with the aircraft; the position term
does not. Stacking fixes gives `H = [I₂ , R_earth · R(ψ)]` per row, so:

| trajectory | rank of the stack | tilt states |
|---|---|---|
| straight leg, ψ constant | 2 of 4 | unobservable — a fixed tilt is a fixed position offset |
| full orbit, ψ sweeps 360° | 4 | separable |

That is the same mechanism the existing circle fit exploits, which matters for
what follows.

## The parameter that decides it

The circle fit already uses the **within-orbit** rotation. An ESKF's marginal
contribution is linking tilt **across** orbits and through the straight legs —
and that requires δθ_body to persist from one orbit to the next.

    one full heading sweep, 250 m loiter at 25 m/s   62.8 s
    ahrs_drift_tau (default)                         60.0 s

**The attitude error decorrelates in about the time it takes to fly one orbit.**
At that ratio there is nothing left to carry between orbits, and the ESKF should
reduce to roughly what the circle fit already does. The prediction is therefore
specific and falsifiable:

* `tau ≪ orbit period` — δθ is white between fixes, unobservable, and no
  architecture helps.
* `tau ≈ orbit period` — the present case. Marginal. Consistent with four
  failed attempts.
* `tau ≫ orbit period` — δθ is effectively a constant bias over many orbits,
  observable given heading variation, and an ESKF should beat the circle fit.

`ahrs_drift_tau` is one of the two parameters this project has never measured;
it is taken from the source paper's Figure 5. **Which regime the real airframe
is in is unknown**, so the question of whether this architecture can work is
currently gated on a number nobody has checked. Two open problems are one.

## Step 1 result: the geometry works, the memory does not

Run on the Canberra flight, 13 760 frames, nine orbits separated by straight
legs (`tools/eskf_observability.py`). Posterior 1-sigma on tilt after a 120 s
window, against a 9 arcmin prior:

| tau | straight leg | full orbit |
|---|---|---|
| 6 s | 2.75′ (3.3x) | 2.47′ (3.6x) |
| 60 s | 6.20′ (1.5x) | **1.74′ (5.2x)** |
| 600 s | 8.10′ (1.1x) | 0.69′ (13x) |
| infinite | 7.73′ (1.2x) | 0.05′ (194x) |

Two things fall straight out.

**Straight legs teach the filter essentially nothing** — 1.1x to 1.5x at the
taus that matter. Exactly as the rank argument predicts: a fixed tilt on a fixed
heading is a fixed position offset, and no amount of integration separates them.

**Orbits do work, at every tau.** Even at 60 s an orbit sharpens tilt from 9′ to
1.74′. So the information is there, and the circle fit is already collecting it.

Which raises the real question — not "can an orbit see tilt" but "can anything
carry that knowledge to the NEXT orbit", because carrying it is the only thing
an ESKF would add over the circle fit. Measured on the same flight:

    gap between consecutive orbits, median          950 s
    ahrs_drift_tau (assumed)                         60 s
    -> 15.8 correlation times

    tau       gap/tau     tilt memory surviving to the next orbit
      6 s       158           0.0000 %
     60 s        15.8         0.0000 %
    600 s         1.6        20.5    %
   6000 s         0.2        85.4    %

**At the assumed tau, nothing survives.** Whatever an orbit learns about tilt is
gone — to fifteen decimal places — before the next orbit begins. An ESKF would
therefore reduce to per-orbit tilt estimation, which is what the circle fit
already does, and it would add complexity for no information.

That is a quantitative explanation for four failed attempts, and it is a better
reason than "we tried and it did not work". It also says precisely what would
change the answer: at tau = 600 s a fifth of the estimate survives the leg, and
cross-orbit linking starts to be worth something.

So the question is not architectural. **It is: what is `ahrs_drift_tau`?**

## Step 2 result: tau was the wrong question, because the model shape is wrong

No real airframe is needed to answer this. `ahrs_drift_tau` is not a gyro
property — it is the correlation time of the EKF3 TILT ERROR, a closed-loop
property of the filter. A SITL `.BIN` carries both `ATT` (what EKF3 believed)
and `SIM` (truth), so the tilt error is directly recoverable and its
autocorrelation can simply be read. Measured on the Canberra flight, 215 105
samples at 25 Hz over 2.4 h (`tools/measure_drift_tau.py`):

| lag | 1 s | 10 s | 60 s | 120 s | 300 s | 630 s | 950 s |
|---|---|---|---|---|---|---|---|
| roll (lateral) | +0.640 | +0.571 | +0.390 | +0.123 | **−0.518** | **+0.444** | −0.379 |
| pitch (forward) | +0.284 | +0.026 | +0.018 | +0.010 | −0.002 | +0.009 | −0.004 |

**Neither axis is a Gauss-Markov process, and neither has tau = 60 s.**

* **Pitch is effectively white beyond ~10 s** (0.026 at 10 s). Far faster than
  the assumed 60 s. Nothing to carry anywhere.
* **Roll oscillates rather than decays** — negative at 300 s, positive again at
  630 s. That is not drift, it is structure locked to the mission: 8 legs of
  480 s plus their loiters is a ~606 s cycle, against an observed ~630 s. It
  survives with the turn samples removed and with them included, so it is not an
  artifact of gap interpolation. This is the one-per-revolution nav-frame error
  `CLAUDE.md` records, seen in the time domain.
* **The magnitudes are 2–5x the assumed value.** `ahrs_drift_sigma` is 0.15 deg
  (9 arcmin); measured is 19 arcmin roll and 45 arcmin pitch over the whole
  flight (16 and 23 on straight legs alone). Of the two never-measured
  parameters, this one is understated.

For the ESKF the conclusion is unchanged and now rests on measurement rather
than assumption: **at the 950 s inter-orbit gap the autocorrelation is ≤ 0 on
both axes.** Nothing survives, so there is nothing to carry between orbits, so
an error-state filter reduces to the per-orbit estimate the circle fit already
produces.

It is stronger than the tau=60 s prediction, in fact: the earlier analysis
assumed a decaying exponential and asked whether it decayed too fast. The
measurement says the process is not exponential at all — pitch is white, roll is
periodic — and a Gauss-Markov tilt model is the wrong shape for either.

**Caveat, and it is not small.** This is SITL's IMU model from
`ardupilot/params/cns_sitl.parm` with the real EKF3 algorithm on top. It is a
measurement of the filter, not of an airframe. A real Cube Orange in a vibrating
fixed-wing may well be slower and larger. What it does establish is that the
number in use is not defensible as a modelling assumption, and that the
apparatus to replace it exists and is cheap.

## Remaining plan, if tau turns out to be long

1. ~~Observability, analytically.~~ Done above.
2. **If it survives that, sweep tau in `transit_scenario`.** It already has
   `--drift-tau`, is seeded and deterministic, and the short config
   (`--km 25 --legs 1`) runs in ~3 minutes, so a 20-seed sweep across
   tau ∈ {6, 60, 600} s is about an hour of CPU. This measures whether the
   *existing* estimator is tau-limited before any ESKF is written.
3. **Only then implement.** Extend the four-state filter with two tilt-error
   states (forward/lateral) under Gauss-Markov dynamics at `tau`, fed by
   per-frame fixes rather than the orbit-averaged one. Compare against the
   circle fit in `transit_scenario` over a seed sweep — never in SITL, whose
   runs differ by 5.80 vs 7.53 km on an unchanged configuration.

## Two other things the abstract suggests

**The rover observation model is the hard part to transfer.** A stationary rover
gets local vertical for free: with no acceleration, the accelerometer measures
gravity. An aircraft cannot, because specific force includes acceleration — which
is precisely why a fix here needs a heading sweep. Whatever the adaptation is,
it is the crux of that paper, and the horizon sensor in this project is an
alternative answer to the same gap.

**Monte Carlo over random trajectories averages over the dominant variable.**
Measured here, fix quality against heading coverage: 90° gives 28 km, 180° gives
12, 270° gives 7, 360° gives 6.2 — and whole revolutions are local minima, with
1.25 revs 64 % worse than 1.00 despite covering more sky. A mean over randomly
generated trajectories marginalises over exactly that. It is a reasonable way to
report robustness and a poor way to find the design variable.

## What would make this project wrong

If an ESKF at the measured `ahrs_drift_tau` beats the circle fit on a seed sweep
in `transit_scenario`, then *"tilt and position are one observable"* is too
strong as written: true per fix, false across a trajectory. That sentence in
`CLAUDE.md` should then be narrowed rather than deleted, because the per-fix
statement would still be correct and is what the `CNS_ORACLE_P` diagnostic
actually demonstrates.
