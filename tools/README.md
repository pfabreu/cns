# tools/

Four groups, split by what they talk to rather than by language.

## `src/` — C++ built by CMake

Everything here is a build target; nothing else in `tools/` is compiled.

| | |
|---|---|
| `analyse_log.cpp` | reprocess a recorded ArduPilot `.BIN` deterministically. The authority when the simulator and SITL disagree |
| `celestial_node.cpp` | the flight node: MAVLink in, celestial fixes out |
| `mavlink_source.{cpp,hpp}` | MAVLink transport, shared by the node and the replayer |
| `pipeline.{cpp,hpp}` | the node's processing chain, separated so it can be tested without a transport |
| `fake_sitl.cpp` | MAVLink replayer. For exercising the transport, NEVER for A/B -- it samples a UDP stream and two runs of one configuration differed by 5.80 vs 7.53 km |
| `gen_dataset.cpp` | image/mask pairs with exact `StarLabel` truth |
| `bench.cpp` | timing against the 100 ms frame budget |
| `visualise.cpp` | render frames for eyeballing |

`celestial_node` and `fake_sitl` are skipped without MAVLink, with a message
that is easy to miss -- check `ls build/` before concluding a target is broken.

## `sitl/` — driving and reading ArduPilot

| | |
|---|---|
| `run_sitl_test.py` | the unattended harness: launches SITL, flies the mission, reaps every process it started |
| `make_mission.py` | generate a `.waypoints` leg/orbit trajectory |
| `dump_log.py` | `.BIN` to CSV |

`run_sitl_test.py` finds the repo root by walking up **three** directories from
itself. That changed when this split happened; if you move it again, fix it.

## `plot/` — figures

`live_view.py` (trajectory and error against time, and the frame pane),
`plot_map.py` (track over Esri World Imagery), `plot_starmap.py` (the camera's
star field at one fix), `plot_orbit_zoom.py` (`--list` tabulates every fix orbit in a run with its
heading span and surviving tilt; `--compare A B` draws two side by side on one shared
scale, `--together A B` puts the route and the fixes in one figure),
`plot_trajectory.py`. `star_names.tsv` is the HR to
proper-name table and `plot_starmap.py` now resolves it relative to its own
location, so the script works from any working directory.

## `analysis/` — offline measurement and codegen

`eskf_observability.py` (is tilt observable from this trajectory?),
`measure_drift_tau.py` (the correlation time of the EKF3 tilt error, measured
rather than assumed), `make_catalog.py` (regenerates
`src/star_catalog_data.cpp` from BSC5 -- run once, not part of the build).

---

Two rules that live in `CLAUDE.md` and are worth repeating here, because both
were learned by getting them wrong: **never A/B on `fake_sitl`**, and **quote a
seed sweep, not a run** -- the transit spreads 6.1 to 8.5 km on an unchanged
configuration.
