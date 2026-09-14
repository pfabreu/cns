#!/usr/bin/env python3
"""Run the whole SITL test unattended: start, fly, deny GNSS, stop, report.

    python3 tools/sitl/run_sitl_test.py --ardupilot ~/ardupilot
    python3 tools/sitl/run_sitl_test.py --ardupilot ~/ardupilot --horizon --speedup 20

Starts SITL and celestial_node, uploads the mission, takes off, denies GNSS once
established, flies to the end, then kills everything and prints a summary. The
CSV it leaves behind can be fed to live_view.py afterwards -- the frames pane
will be empty unless --frames-dir was used, which is fine for trajectories.

WHY A SCRIPT. Doing this by hand means typing into MAVProxy at the right
moments, which is unrepeatable and impossible remotely. It also means the
denial time is whatever you managed, rather than a number the analysis can rely
on. Here it is a parameter.

EVERYTHING IS KILLED ON EXIT, including on Ctrl-C and on failure. Each child
runs in its own process group so the whole tree dies with it -- sim_vehicle.py
in particular spawns several processes and killing only the parent leaves
ArduPlane running and the next run fails on a bound port.
"""

import argparse
import os
import shlex
import signal
import subprocess
import sys
import time

try:
    from pymavlink import mavutil
except ImportError:
    raise SystemExit("needs pymavlink: pip install pymavlink")


class Child:
    """A subprocess in its own process group, killed reliably."""

    def __init__(self, name, cmd, log):
        self.name = name
        self.log_path = log
        self.log = open(log, "wb")
        print(f"  starting {name}: {' '.join(cmd[:3])} ...")
        self.p = subprocess.Popen(cmd, stdout=self.log, stderr=subprocess.STDOUT,
                                  start_new_session=True)

    def alive(self):
        return self.p.poll() is None

    def kill(self):
        if self.p.poll() is not None:
            self.log.close()
            return
        try:
            os.killpg(os.getpgid(self.p.pid), signal.SIGTERM)
            try:
                self.p.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(os.getpgid(self.p.pid), signal.SIGKILL)
                self.p.wait(timeout=5)
        except (ProcessLookupError, PermissionError):
            pass
        self.log.close()
        print(f"  stopped {self.name}")


# Options whose VALUE can legitimately start with "-": a negative latitude, or
# a string of arguments meant for a child process.
_NEG_VALUE_OPTS = {"--location", "--node-args", "--horizon-node-args"}


def joinNegativeValues(argv):
    """Rewrite `--opt -value` to `--opt=-value` for the options listed above.

    argparse cannot tell a value beginning with "-" from another option unless
    the value parses as a plain negative number, and "lat,lon,alt,heading" does
    not. So EVERY southern-hemisphere departure fails with

        error: argument --location: expected one argument

    which points at the wrong thing entirely -- the argument is right there.
    Found flying Bluff -> the Snares at -46.6; the northern-hemisphere runs
    never hit it because their latitudes are positive."""
    out, i = [], 0
    while i < len(argv):
        if (argv[i] in _NEG_VALUE_OPTS and i + 1 < len(argv)
                and argv[i + 1].startswith("-")):
            out.append(f"{argv[i]}={argv[i + 1]}")
            i += 2
            continue
        out.append(argv[i])
        i += 1
    return out


def sitl_pids():
    """PIDs of running ArduPlane SITL binaries.

    The bracket keeps pgrep from matching this script's own command line. A
    self-match would report us as a stray and the reap below would kill the
    reaper -- the same trap that has bitten probes in this repo before."""
    r = subprocess.run(["pgrep", "-f", "bin/ardupla[n]e"],
                       capture_output=True, text=True)
    return {int(x) for x in r.stdout.split() if x.isdigit()}


def reap_stray_sitl(pre_existing):
    """Kill any ArduPlane THIS run started that is somehow still alive.

    WHY THIS IS NEEDED. sim_vehicle.py does not exec the binary directly; it
    goes through run_in_terminal_window.sh, which puts ArduPlane in its OWN
    process group and session. Killing sim_vehicle's group therefore does not
    reliably take ArduPlane with it. The wrapper has a watchdog that usually
    notices within a second, but "usually" leaves 5760 still bound and the next
    run failing to start, which is exactly the failure the module docstring
    promises not to have.

    SIGTERM before SIGKILL: ArduPlane flushes its dataflash .BIN on the way
    out, and that log is the real analysis input for analyse_log.

    Only PIDs that appeared after we started are touched, so an unrelated SITL
    belonging to someone else survives."""
    for attempt in range(10):
        stray = sitl_pids() - pre_existing
        if not stray:
            return True
        for pid in stray:
            try:
                os.kill(pid, signal.SIGTERM if attempt == 0 else signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
        time.sleep(0.5)
    stray = sitl_pids() - pre_existing
    if stray:
        print(f"  WARNING: ArduPlane STILL RUNNING: {sorted(stray)} -- kill by hand")
        return False
    return True


def wait_for(m, cond, timeout, what):
    """Poll until cond(m) or timeout. Returns True on success."""
    t0 = time.time()
    while time.time() - t0 < timeout:
        m.recv_match(blocking=True, timeout=1)
        if cond(m):
            return True
    print(f"  TIMEOUT waiting for {what} after {timeout:.0f}s")
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ardupilot", default=os.path.expanduser("~/ardupilot"))
    ap.add_argument("--mission", default="ardupilot/missions/demo.waypoints")
    ap.add_argument("--csv", default="live.csv")
    ap.add_argument("--frames-dir", default=None)
    ap.add_argument("--horizon", action="store_true",
                    help="run celestial_node with a Lepton 2.5 horizon sensor")
    ap.add_argument("--horizon-compare", action="store_true",
                    help="run TWO nodes on the same flight, one with a horizon "
                         "sensor and one without, writing --csv and "
                         "--horizon-csv. The only honest way to A/B the "
                         "horizon: both see identical telemetry.")
    ap.add_argument("--horizon-csv", default="live_horizon.csv",
                    help="second node's CSV under --horizon-compare")
    ap.add_argument("--horizon-port", type=int, default=14558)
    ap.add_argument("--fix-sigma", type=float, default=0.0,
                    help="override the plain node's full-sweep fix sigma, m")
    ap.add_argument("--horizon-fix-sigma", type=float, default=0.0,
                    help="override the HORIZON node's full-sweep fix sigma, m. "
                         "Set it to the accuracy the horizon actually delivers "
                         "rather than the accuracy it is assumed to.")
    ap.add_argument("--node-args", default="",
                    help="extra arguments passed verbatim to the FIRST node, "
                         "e.g. \"--horizon\". Lets both arms of "
                         "--horizon-compare differ in one variable only")
    ap.add_argument("--horizon-node-args", default="",
                    help="extra arguments for the SECOND node")
    ap.add_argument("--run-name", default=None,
                    help="write this run's CSVs and node logs to "
                         "runs/<timestamp>-<name>/ instead of the repo root. "
                         "Without it every run overwrites live.csv, and an "
                         "830 km crossing that cost an hour is gone the moment "
                         "the next run starts.")
    ap.add_argument("--location", default=None,
                    help="SITL home as LAT,LON,ALT,HEADING (sim_vehicle -l). "
                         "Must match the --lat/--lon the mission was built "
                         "with, or the aircraft starts nowhere near it.")
    ap.add_argument("--utc", default=None,
                    help="observation epoch for the node, 'YYYY-MM-DD HH:MM:SS'. "
                         "The camera looks at the ZENITH, so latitude and this "
                         "decide which sky is overhead -- a daytime epoch at "
                         "the flight's longitude gives no stars at all.")
    ap.add_argument("--speedup", type=float, default=10.0)
    ap.add_argument("--deny-after", type=float, default=240.0,
                    help="seconds after AUTO starts before GNSS is denied")
    ap.add_argument("--max-minutes", type=float, default=40.0,
                    help="wall-clock cap; the run is killed regardless")
    ap.add_argument("--node-port", type=int, default=14556)
    ap.add_argument("--ctrl-port", type=int, default=5760,
                    help="SITL SERIAL0 TCP port. Connecting to it is also what "
                         "releases SITL from its startup wait, so it is not "
                         "optional.")
    ap.add_argument("--keep-going", action="store_true",
                    help="do not stop when the mission completes")
    a = ap.parse_args(joinNegativeValues(sys.argv[1:]))

    # repo root: tools/sitl/run_sitl_test.py -> up three, not two. This moved
    # when tools/ was split into src, sitl, plot and analysis.
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    # Per-run output directory. Resolved BEFORE the nodes start, so the CSV
    # paths handed to them already point somewhere unique.
    run_dir = None
    if a.run_name:
        stamp = time.strftime("%Y%m%d-%H%M%S")
        run_dir = os.path.join(root, "runs", f"{stamp}-{a.run_name}")
        os.makedirs(run_dir, exist_ok=True)
        a.csv = os.path.join(run_dir, os.path.basename(a.csv))
        a.horizon_csv = os.path.join(run_dir, os.path.basename(a.horizon_csv))
        print(f"  run directory: {run_dir}")
    node = os.path.join(root, "build", "celestial_node")
    if not os.path.exists(node):
        raise SystemExit(f"no {node} -- build first, with -DMAVLINK_DIR=...")
    mission = os.path.join(root, a.mission)
    if not os.path.exists(mission):
        raise SystemExit(f"no mission at {mission} -- "
                         f"python3 tools/sitl/make_mission.py > {a.mission}")

    sim = os.path.join(a.ardupilot, "Tools", "autotest", "sim_vehicle.py")
    if not os.path.exists(sim):
        raise SystemExit(f"no sim_vehicle.py under {a.ardupilot}")

    # Snapshot before starting anything, so the reap in finally can tell our
    # ArduPlane apart from one that was already running.
    pre_existing = sitl_pids()

    kids = []
    try:
        # --- SITL ------------------------------------------------------------
        # --no-mavproxy: we drive it with pymavlink instead, which is the whole
        # point. The node and this script still get one link each, but NOT via
        # --out.
        #
        # WHY NOT --out. It is a MAVProxy option. sim_vehicle.py consumes
        # opts.out inside start_mavproxy() only, and --no-mavproxy returns
        # before that, so BOTH --out flags were silently discarded and nothing
        # ever arrived on either port. The failure looked like a dead vehicle
        # ("no heartbeat") rather than a dropped argument, which is why it cost
        # a session to find. -A passes arguments to the ArduPlane binary
        # itself, where --serial1 is a real MAVLink port.
        #
        # udpclient, not udpin: celestial_node BINDS its port
        # (mavlink_source.cpp), so SITL has to be the sender.
        # -w WIPES THE EEPROM, and this run is not repeatable without it.
        # SITL persists parameters in ./eeprom.bin, so the GNSS denial this
        # script performs (SIM_GPS1_ENABLE 0, EK3_SRC1_POSXY/VELXY 0) is still
        # in force at the NEXT boot. Symptom: "TIMEOUT waiting for GPS 3D fix"
        # and a refused arm on every run after the first successful one, with
        # nothing in the log to say why -- the denial is invisible because it
        # was applied by a previous process. Those three are comments in
        # cns_demo.parm, not settings, so --add-param-file does not restore
        # them. Wiping also makes each run start from one known state, which
        # is the whole point of an unattended test.
        cmd = [sys.executable, sim, "-v", "ArduPlane", "--no-mavproxy", "-w",
               "--speedup", str(int(a.speedup)),
               f"--add-param-file={root}/ardupilot/params/cns_sitl.parm",
               f"--add-param-file={root}/ardupilot/params/cns_demo.parm",
               *(["-l", a.location] if a.location else []),
               "-A", f"--serial1=udpclient:127.0.0.1:{a.node_port}"
                     + (f" --serial2=udpclient:127.0.0.1:{a.horizon_port}"
                        if a.horizon_compare else "")]
        kids.append(Child("SITL", cmd, "/tmp/sitl.log"))

        # SERIAL0 is a TCP SERVER and SITL BLOCKS at startup until a client
        # connects to it. MAVProxy normally is that client; with --no-mavproxy
        # this script has to be, or the simulator never leaves its wait loop
        # and NO link produces a byte -- including the node's, which makes this
        # look like a node bug when it is not.
        #
        # retries= and the long timeout cover sim_vehicle.py rebuilding
        # ArduPlane first: a cold waf build measured 47 s, which on its own ate
        # half of the previous 90 s budget.
        if run_dir:
            with open(os.path.join(run_dir, "command.txt"), "w") as f:
                f.write(" ".join(sys.argv) + "\n")
                f.write(f"\nSITL: {' '.join(cmd)}\n")
        print(f"  connecting on tcp:127.0.0.1:{a.ctrl_port} ...")
        m = mavutil.mavlink_connection(f"tcp:127.0.0.1:{a.ctrl_port}",
                                       retries=60)
        if not m.wait_heartbeat(timeout=180):
            raise SystemExit("no heartbeat -- see /tmp/sitl.log")
        print(f"  heartbeat from system {m.target_system}")

        # --- navigation node --------------------------------------------------
        ncmd = [node, "--port", str(a.node_port), "--no-inject",
                "--csv", os.path.join(root, a.csv)]
        if a.utc:
            ncmd += ["--utc", a.utc]
        if a.fix_sigma > 0:
            ncmd += ["--fix-sigma", str(a.fix_sigma)]
        ncmd += shlex.split(a.node_args)
        if a.horizon:
            ncmd.append("--horizon")
        if a.frames_dir:
            ncmd += ["--frames-dir", os.path.join(root, a.frames_dir)]
        node_log = (os.path.join(run_dir, "node.log") if run_dir
                    else "/tmp/node.log")
        kids.append(Child("celestial_node", ncmd, node_log))

        # SECOND NODE, same flight. Comparing a horizon run against a separate
        # no-horizon run cannot settle anything: SITL is not deterministic, and
        # two runs of ONE configuration have differed by 5.80 vs 7.53 km. Both
        # nodes here consume the same MAVLink stream from the same aircraft, so
        # the horizon is the only difference between the two CSVs.
        if a.horizon_compare:
            hcmd = [node, "--port", str(a.horizon_port), "--no-inject",
                    "--horizon", "--csv", os.path.join(root, a.horizon_csv)]
            if a.utc:
                hcmd += ["--utc", a.utc]
            if a.horizon_fix_sigma > 0:
                hcmd += ["--fix-sigma", str(a.horizon_fix_sigma)]
            hcmd += shlex.split(a.horizon_node_args)
            hlog = (os.path.join(run_dir, "node_horizon.log") if run_dir
                    else "/tmp/node_h.log")
            kids.append(Child("celestial_node+horizon", hcmd, hlog))

        # --- mission ---------------------------------------------------------
        print("  waiting for EKF and GPS ...")
        # Arming is refused until EKF3 has a position solution. Waiting for a
        # 3D fix is the cheap proxy; the arming call below is the real test.
        wait_for(m, lambda c: (c.messages.get("GPS_RAW_INT")
                               and c.messages["GPS_RAW_INT"].fix_type >= 3),
                 120, "GPS 3D fix")
        time.sleep(10)   # let EKF3 settle before asking it to arm

        print(f"  uploading {a.mission} ...")
        upload_mission(m, mission)

        print("  arming ...")
        m.set_mode_apm("AUTO")
        m.mav.command_long_send(m.target_system, m.target_component,
                                mavutil.mavlink.MAV_CMD_COMPONENT_ARM_DISARM,
                                0, 1, 0, 0, 0, 0, 0, 0)
        if not wait_for(m, lambda c: (c.messages.get("HEARTBEAT")
                                      and c.messages["HEARTBEAT"].base_mode
                                      & mavutil.mavlink.MAV_MODE_FLAG_SAFETY_ARMED),
                        60, "arm"):
            print("  ARM FAILED -- check /tmp/sitl.log for the pre-arm message")
            return 1
        print("  armed, flying AUTO")

        t_start = time.time()
        denied_at = None
        deadline = t_start + a.max_minutes * 60

        while time.time() < deadline:
            m.recv_match(blocking=True, timeout=1)
            el = time.time() - t_start

            if denied_at is None and el >= a.deny_after:
                print(f"  DENYING GNSS at t+{el:.0f}s wall clock")
                for name, val in (("SIM_GPS1_ENABLE", 0),
                                  ("EK3_SRC1_POSXY", 0),
                                  ("EK3_SRC1_VELXY", 0)):
                    m.mav.param_set_send(m.target_system, m.target_component,
                                         name.encode(), float(val),
                                         mavutil.mavlink.MAV_PARAM_TYPE_REAL32)
                    time.sleep(0.3)
                denied_at = el

            hb = m.messages.get("HEARTBEAT")
            if hb and not (hb.base_mode
                           & mavutil.mavlink.MAV_MODE_FLAG_SAFETY_ARMED):
                if not a.keep_going and el > 60:
                    print(f"  disarmed at t+{el:.0f}s -- mission complete")
                    break

        if time.time() >= deadline:
            print(f"  hit the {a.max_minutes:.0f} minute cap")

        # Let the node flush its last rows before it is killed.
        time.sleep(3)

    except KeyboardInterrupt:
        print("\n  interrupted")
    finally:
        for c in reversed(kids):
            c.kill()
        reap_stray_sitl(pre_existing)

    report(os.path.join(root, a.csv))
    if a.horizon_compare:
        print("\n  --- second node: WITH horizon sensor ---")
        report(os.path.join(root, a.horizon_csv))
    return 0


def upload_mission(m, path):
    """Upload a QGC WPL 110 file. Written out rather than shelling to MAVProxy
    so the whole run has one dependency and one failure mode."""
    from pymavlink import mavwp
    wp = mavwp.MAVWPLoader()
    wp.load(path)
    m.waypoint_clear_all_send()
    m.waypoint_count_send(wp.count())
    for _ in range(wp.count()):
        msg = m.recv_match(type=["MISSION_REQUEST", "MISSION_REQUEST_INT"],
                           blocking=True, timeout=10)
        if msg is None:
            print("  mission upload stalled")
            return False
        m.mav.send(wp.wp(msg.seq))
    ack = m.recv_match(type="MISSION_ACK", blocking=True, timeout=10)
    print(f"  mission uploaded, {wp.count()} items, ack={ack.type if ack else '?'}")
    return True


def report(csv_path):
    """Whatever the run produced. Deliberately terse -- analyse_log on the .BIN
    is the real analysis, and this only says whether the flight was usable."""
    import csv as _csv
    import statistics as st
    if not os.path.exists(csv_path):
        print(f"\n  no {csv_path} -- the node produced nothing, see /tmp/node.log")
        return
    rows = list(_csv.DictReader(open(csv_path)))
    if not rows:
        print(f"\n  {csv_path} is empty -- see /tmp/node.log")
        return

    def f(r, k, d=0.0):
        try:
            return float(r[k])
        except (ValueError, KeyError):
            return d

    fx = [f(r, "orbit_err_m") / 1000 for r in rows if f(r, "orbit_err_m") >= 0]
    dr = [f(r, "dr_err_m") / 1000 for r in rows if f(r, "dr_err_m") >= 0]
    op = [f(r, "dr_open_err_m") / 1000 for r in rows if f(r, "dr_open_err_m") >= 0]
    bs = [f(r, "boresight_deg") for r in rows if f(r, "orbit_err_m") >= 0]

    print(f"\n  {csv_path}: {len(rows)} rows, {len(fx)} fixes")
    if fx:
        print(f"    fix error   median {st.median(fx):6.2f} km   "
              f"best {min(fx):.2f}   worst {max(fx):.2f}")
    if bs:
        print(f"    boresight   {bs[0]:.4f} -> {bs[-1]:.4f} deg"
              f"   (0.40 uncalibrated; falling is good)")
    if dr and op:
        # MEDIAN AND RMS, NOT THE LAST SAMPLE.
        #
        # dr_err_m is written only on fix rows, and the aided error is a
        # SAWTOOTH: it grows along a leg and collapses at each fix orbit. The
        # final value is therefore just wherever the flight happened to stop in
        # that cycle, and it swings by an order of magnitude. Measured on the
        # Sagres->Porto Santo crossing, night window, the two configurations
        # ended on 10.81 km and 4.00 km -- a 2.7x gap that reads as a decisive
        # win and is pure sampling. The aggregates over the SAME data:
        #
        #                  final    mean   median    RMS     p90
        #   no horizon     10.81   10.17    7.49   13.89   22.41
        #   horizon         4.00   11.62    9.11   14.79   28.58
        #
        # i.e. the arm that "won" on the last sample is equal or worse on every
        # aggregate. Quote a distribution; the final sample is kept only for
        # continuity with older logs.
        import math
        def _rms(v):
            return math.sqrt(sum(x * x for x in v) / len(v))
        def _p90(v):
            return sorted(v)[min(len(v) - 1, int(0.9 * len(v)))]
        print(f"    DR + fixes  median {st.median(dr):6.2f} km   "
              f"RMS {_rms(dr):6.2f}   p90 {_p90(dr):6.2f}   "
              f"(final {dr[-1]:.2f}, n={len(dr)})")
        print(f"    unaided     median {st.median(op):6.2f} km   "
              f"RMS {_rms(op):6.2f}   (final {op[-1]:.2f})")
        if st.median(dr) > st.median(op):
            print("    NOTE: fixes made it WORSE than unaided dead reckoning.")
            print("          Over a short flight that is expected -- DR is good")
            print("          over minutes. Fly further before concluding.")
    print(f"\n  next:  python3 tools/plot/live_view.py --csv {csv_path}")
    print(f"         python3 tools/sitl/dump_log.py <logs/NNN.BIN> > flight.csv")
    print(f"         ./build/analyse_log flight.csv")


if __name__ == "__main__":
    sys.exit(main())
