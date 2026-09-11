#!/usr/bin/env python3
"""Generate an ArduPilot mission that demonstrates what this system claims.

    python3 tools/make_mission.py > missions/demo.waypoints
    python3 tools/make_mission.py --legs 12 --leg-km 15 > missions/long.waypoints

WHY THIS SHAPE. The claim is not "celestial navigation is accurate" -- over a
few minutes plain dead reckoning beats it. The claim is that DR error grows
WITHOUT LIMIT while the celestial fix BOUNDS it, so the trajectory has to be
long enough for those two curves to separate visibly. Unaided DR runs at roughly
30% of distance travelled; at 7 km of bounded error the crossover is around
25 km, and you want several times that to make a convincing plot.

Hence: leg, orbit, leg, orbit. The legs let DR drift. The orbits are where the
fixes happen -- with no vertical reference a straight-leg fix carries the whole
AHRS tilt error and is worth ~40 km, so the aircraft must maneuver to navigate.

The numbers are the measured optima, not guesses:

  * 250 m orbits. A fix wants a SHORT PERIOD, because AHRS drift only averages
    out if the orbit is quick against its ~60 s correlation time; a CALIBRATION
    wants LOW BANK, because maneuver coupling puts a bank-proportional tilt into
    the AHRS. 250 m is the optimum at both the GPS-aided and GNSS-denied
    coupling, so the choice is insensitive to which is right.
  * 3 turns at departure to calibrate the mounting, 2 per fix orbit after.
  * 800 m, matching the altitude the sensor models assume.
"""

import argparse
import math

# SITL's default start (CMAC, Canberra). Change with --lat/--lon, and pass the
# same place to sim_vehicle.py with --location if you move it.
LAT0, LON0 = -35.363261, 149.165230

WPL_HEADER = "QGC WPL 110"
NAV_WAYPOINT, NAV_LOITER_TURNS, NAV_TAKEOFF, NAV_RTL = 16, 18, 22, 20
FRAME_REL = 3


def offset(lat, lon, north_m, east_m):
    dlat = north_m / 111320.0
    dlon = east_m / (111320.0 * math.cos(math.radians(lat)))
    return lat + dlat, lon + dlon


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lat", type=float, default=LAT0)
    ap.add_argument("--lon", type=float, default=LON0)
    ap.add_argument("--alt", type=float, default=800.0)
    ap.add_argument("--legs", type=int, default=8)
    ap.add_argument("--leg-km", type=float, default=12.0)
    ap.add_argument("--radius", type=float, default=250.0)
    ap.add_argument("--turns", type=int, default=2)
    ap.add_argument("--calib-turns", type=int, default=3)
    ap.add_argument("--track", type=float, default=190.0,
                    help="degrees; legs fan by --fan to vary heading")
    ap.add_argument("--fan", type=float, default=4.0,
                    help="degrees added to the track PER LEG. The default 4 is "
                         "a gentle spread over the 8-leg demo (32 deg total), "
                         "but it scales with --legs: 24 legs fan 92 deg and "
                         "the route curves away from the intended track "
                         "entirely. Use 0 for a straight point-to-point "
                         "crossing, where the heading variation comes from the "
                         "fix orbits rather than the legs.")
    a = ap.parse_args()

    out = [WPL_HEADER]
    i = 0

    def row(cmd, p1=0, p2=0, p3=0, p4=0, lat=0, lon=0, alt=0, frame=FRAME_REL):
        nonlocal i
        out.append(f"{i}\t{1 if i == 0 else 0}\t{frame}\t{cmd}\t"
                   f"{p1:.6f}\t{p2:.6f}\t{p3:.6f}\t{p4:.6f}\t"
                   f"{lat:.7f}\t{lon:.7f}\t{alt:.6f}\t1")
        i += 1

    # 0: home
    row(NAV_WAYPOINT, lat=a.lat, lon=a.lon, alt=a.alt, frame=0)
    # 1: takeoff
    row(NAV_TAKEOFF, p4=0, lat=a.lat, lon=a.lon, alt=a.alt)

    # 2: calibration loiter. Longest of the orbits, because the mounting is
    # estimated ONCE here and everything downstream inherits it.
    lat, lon = offset(a.lat, a.lon, 1500, 0)
    row(NAV_LOITER_TURNS, p1=a.calib_turns, p3=a.radius,
        lat=lat, lon=lon, alt=a.alt)

    # then leg / fix-orbit pairs
    for k in range(a.legs):
        hdg = math.radians(a.track + a.fan * k)   # fan, so heading varies
        d = a.leg_km * 1000.0
        lat, lon = offset(lat, lon, d * math.cos(hdg), d * math.sin(hdg))
        row(NAV_WAYPOINT, lat=lat, lon=lon, alt=a.alt)
        row(NAV_LOITER_TURNS, p1=a.turns, p3=a.radius,
            lat=lat, lon=lon, alt=a.alt)

    row(NAV_RTL)

    total = a.legs * a.leg_km
    print("\n".join(out))
    print(f"# {a.legs} legs x {a.leg_km:.0f} km = {total:.0f} km, "
          f"{a.legs + 1} orbits at {a.radius:.0f} m", file=__import__("sys").stderr)
    print(f"# at 25 m/s that is {total * 1000 / 25 / 60:.0f} min of legs plus "
          f"~{(a.legs * a.turns + a.calib_turns) * 2 * math.pi * a.radius / 25 / 60:.0f} min "
          f"of orbiting", file=__import__("sys").stderr)
    print(f"# unaided DR should reach ~{0.3 * total:.0f} km; the bounded "
          f"result should stay near 7 km", file=__import__("sys").stderr)


if __name__ == "__main__":
    main()
