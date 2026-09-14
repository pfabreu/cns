#!/usr/bin/env python3
"""Are the tilt states observable alongside position? Analytic, no simulation.

    python3 tools/analysis/eskf_observability.py runs/<run>/live.csv

WHAT THIS ANSWERS. NOTES-eskf.md proposes replacing the four-state position
filter with an error-state filter that also carries tilt. Before writing one,
ask whether tilt is observable at all from this trajectory. It is a property of
the heading profile and the drift correlation time, so it needs a logged flight
and no new one.

THE MODEL. A celestial fix's error decomposes as

    fix_error_NED = dp + k * R(psi) * dtheta_body

with psi the heading. The tilt term rotates with the aircraft, the position term
does not, so each fix contributes a 2x4 row block

    H_i = [ I2 | k * R(psi_i) ]

Units are chosen to make the blocks comparable, or the condition number means
nothing: position in km, tilt in arcmin, k = 1.852 km/arcmin (1 arcmin = 1 nm).

THE DRIFT, AND WHY A GRAMIAN IS THE WRONG TOOL FOR IT. dtheta is Gauss-Markov
with correlation time tau. A deterministic observability Gramian discounts the
tilt columns by exp(-Dt/tau) and stops there -- which makes SHORT tau look
observable, because a decaying state is easy to tell apart from a constant
position offset by its decay alone. That is an artifact: the real process does
not decay to zero, it is re-randomised by process noise and stays stationary.
Omitting that noise measures the wrong thing, and measures it backwards.

So run the actual covariance recursion instead:

    predict:  P <- F P F' + Q,  F = diag(1, 1, e^-dt/tau, e^-dt/tau)
              Q_tilt = sigma_theta^2 (1 - e^-2dt/tau)   <- the part that matters
    update:   standard Kalman with H above and R = (meas noise)^2

and report the STEADY-STATE posterior on tilt. Long tau -> little process noise
per step, tilt persists across the heading sweep, observable. Short tau -> the
process noise floods it faster than the sweep can resolve it, and no amount of
geometry recovers it.

READING THE OUTPUT. sigma_tilt is the posterior 1-sigma on tilt, in arcmin, for
1 km of per-fix measurement noise: sqrt of the mean of the tilt diagonal of
W^-1. Small means observable. Compare it against the AHRS tilt error the filter
would be trying to correct -- around 16 arcmin in these runs. An estimator whose
posterior is worse than the prior is not learning anything.
"""
import argparse, csv, math

import numpy as np

K_KM_PER_ARCMIN = 1.852


def load(path):
    """(t, heading_rad) for every frame carrying a valid per-frame fix."""
    out = []
    for r in csv.DictReader(open(path)):
        try:
            if float(r["frame_err_m"]) < 0:
                continue
            out.append((float(r["t"]), math.radians(float(r["yaw_deg"]))))
        except (ValueError, KeyError):
            continue
    return out


def heading_span(psis):
    """Circular range covered, degrees. 360 minus the largest empty gap."""
    d = sorted((math.degrees(p) % 360.0) for p in psis)
    if len(d) < 2:
        return 0.0
    gaps = [d[i + 1] - d[i] for i in range(len(d) - 1)] + [360.0 - d[-1] + d[0]]
    return 360.0 - max(gaps)


def run_filter(win, tau, sigma_theta, meas_km, pos_walk_km_s):
    """Covariance recursion over one window. Returns posterior tilt sigma."""
    P = np.diag([50.0 ** 2, 50.0 ** 2, sigma_theta ** 2, sigma_theta ** 2])
    R = np.eye(2) * meas_km ** 2
    prev_t = win[0][0]
    for t, psi in win:
        dt = max(1e-3, t - prev_t)
        prev_t = t
        if tau is None:
            f, q = 1.0, 0.0
        else:
            f = math.exp(-dt / tau)
            q = sigma_theta ** 2 * (1.0 - math.exp(-2.0 * dt / tau))
        F = np.diag([1.0, 1.0, f, f])
        Q = np.diag([(pos_walk_km_s * dt) ** 2, (pos_walk_km_s * dt) ** 2, q, q])
        P = F @ P @ F.T + Q
        c, s = math.cos(psi), math.sin(psi)
        H = np.array([[1.0, 0.0, K_KM_PER_ARCMIN * c, -K_KM_PER_ARCMIN * s],
                      [0.0, 1.0, K_KM_PER_ARCMIN * s,  K_KM_PER_ARCMIN * c]])
        S = H @ P @ H.T + R
        K = P @ H.T @ np.linalg.inv(S)
        P = (np.eye(4) - K @ H) @ P
    return math.sqrt(max(0.0, (P[2, 2] + P[3, 3]) / 2.0))


ap = argparse.ArgumentParser()
ap.add_argument("csv")
ap.add_argument("--window", type=float, default=120.0, help="seconds")
ap.add_argument("--straight-max", type=float, default=45.0,
                help="heading span below which a window counts as a straight leg")
ap.add_argument("--orbit-min", type=float, default=300.0,
                help="heading span above which a window counts as an orbit")
ap.add_argument("--taus", default="6,60,600,inf")
ap.add_argument("--sigma-theta", type=float, default=9.0,
                help="Gauss-Markov stationary tilt sigma, arcmin "
                     "(ahrs_drift_sigma 0.15 deg = 9 arcmin)")
ap.add_argument("--meas-km", type=float, default=1.0,
                help="per-frame fix noise NOT explained by tilt, km. The 30 km "
                     "per-frame error in these runs IS the tilt, i.e. signal")
ap.add_argument("--pos-walk", type=float, default=0.003,
                help="position random walk, km/s (~10 km/h of DR drift)")
a = ap.parse_args()

rows = load(a.csv)
if len(rows) < 50:
    raise SystemExit(f"only {len(rows)} usable frames in {a.csv}")
print(f"  {len(rows)} frames, {rows[0][0]:.0f}..{rows[-1][0]:.0f} s "
      f"({(rows[-1][0]-rows[0][0])/3600:.2f} h)")

# Slide a window over the flight and classify each by how much heading it covers.
wins, i = [], 0
while i < len(rows):
    t0 = rows[i][0]
    w = [r for r in rows[i:] if r[0] <= t0 + a.window]
    if len(w) >= 10:
        wins.append(w)
    j = i
    while j < len(rows) and rows[j][0] < t0 + a.window / 2:
        j += 1
    i = max(j, i + 1)

straight = [w for w in wins if heading_span([p for _, p in w]) < a.straight_max]
orbit = [w for w in wins if heading_span([p for _, p in w]) > a.orbit_min]
print(f"  {len(wins)} windows of {a.window:.0f} s: "
      f"{len(straight)} straight (<{a.straight_max:.0f} deg), "
      f"{len(orbit)} orbit (>{a.orbit_min:.0f} deg)\n")

taus = [None if x.strip() == "inf" else float(x) for x in a.taus.split(",")]
print(f"  prior tilt sigma {a.sigma_theta:.1f} arcmin, measurement noise "
      f"{a.meas_km:.1f} km, position walk {a.pos_walk*3600:.1f} km/h\n")
print(f"  {'tau':>6}  {'trajectory':<10} {'n':>4} {'posterior tilt':>15} {'vs prior':>10}")
print(f"  {'':>6}  {'':<10} {'':>4} {'arcmin':>15} {'':>10}")
for tau in taus:
    label = "inf" if tau is None else f"{tau:.0f}s"
    for name, group in (("straight", straight), ("orbit", orbit)):
        if not group:
            continue
        ss = [run_filter(w, tau, a.sigma_theta, a.meas_km, a.pos_walk)
              for w in group]
        med = float(np.median(ss))
        print(f"  {label:>6}  {name:<10} {len(group):>4} {med:15.3f} "
              f"{a.sigma_theta / med:9.1f}x")
    print()

print("  'vs prior' is how much the geometry sharpens the tilt estimate. 1x means")
print("  it learned nothing; the filter would just be repeating its own prior.")
