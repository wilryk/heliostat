# STATE — beamdown (rewritten 2026-07-30 ~10:00)

## What this project is
Python package `beamdown/` drives Quadoa Optical CAD to trace a 645-heliostat
beam-down solar field for a paper comparing annualized collected energy across
geometries: 3 secondary layouts (axicon / prime_focus / cassegrain) x
(focused / flat heliostats). DNI from PVGIS TMY, monthly-mean default.

## Running RIGHT NOW
- **full8 sweep** (12 dates, 161 timesteps, axicon, occluders traced, 120k
  rays): resumed after a power-event kill at 93/161; now ~97/161, ETA ~19:00.
  Holds the ONLY Quadoa licence seat. Log: `analysis_output/full8.log`,
  lock: `analysis_output/.full8.lock`. On completion it auto-runs
  `scripts/report_energy.py` (annual MWh, sine fit, declination pairs).
- **Vetting analysis** (subagent, no licence): scalar vs traced occlusion from
  the full5/full6/full7 ladder → verdict decides whether remaining sweeps skip
  `--occluders`. Output will land in `analysis_output/vet_occlusion/` +
  `scripts/vet_occlusion_scalars.py` + README section.

## Suite status
`python -m tests.verify --no-quadoa` → 12/12. `python tests/test_gui.py
analysis_output/full7` → PASS. Axicon numbers bit-identical through all
refactors (verify stage 3b). Reference: full7 annual = 10,237.0 MWh @ monthly
DNI, optical eta 0.6040.

## Decisions made (owner)
- 120,000 rays everywhere — SNR on single-heliostat irradiance plots drives it.
- One traceRays call per heliostat (measured: 646 ms vs 698 ms for 2 calls;
  chunking is pure loss at this scale — `scripts/probe_ray_cost.py`).
- Comparison sweeps trace 7 distinct declinations, not 12 dates (5 dates are
  declination duplicates; interpolation surface identical).
- Prime focus receiver at z = 47,000 mm (H = 20 m ABOVE axicon vertex,
  mirroring receiver 20 m below). Cassegrain aims at its own F1 = 38,986 mm.
- Cassegrain design final: rim 15 m radius @ z 32,460; K = -5.8789,
  |R_v| = 32,181.3 mm, vertex z = 29,589 (`scripts/design_cassegrain.py`).

## Open decisions
- Scalar-vs-traced occlusion verdict (vetting agent, in flight). If clean:
  5 remaining sweeps run scalar at 7 declinations (~est. cheaper; measure
  no-occluder trace speed when seat frees).
- Chunk-size crossover above 120k rays (owner saw small-chunks-win last year;
  measured opposite below 120k) — extended probe queued, needs seat.

## Queued work (in order)
1. full8 finishes → read report, patch nothing (manifest already corrected),
   run `scripts/verify_prime_focus_model.py` (needs seat, minutes).
2. Launch axicon-flat sweep (no new model needed) — first of the 5.
3. Owner builds cassegrain hyperboloid in Quadoa manually (numbers above);
   prime-focus model ALREADY BUILT: `models/heliostat_field_prime_focus.optx`.
4. Prime-focus and cassegrain sweeps (x focused/flat).
5. Cross-geometry comparison report + paper figures.
6. AFTER full8 + agents go quiet: relocate the whole tree (with .git) from
   C:\gitlab to C:\gitlab\heliostats (owner wants C:\gitlab as a multi-repo
   workspace). Same-volume rename, instant; blocked only by open handles.
   Then add the owner's GitLab remote (URL not yet provided) and push.

## Traps that bite (details in CLAUDE.md / README)
- ONE licence seat. Never workers>1, never retry seat failures, never import
  quadoa while a `.lock` dir exists under analysis_output/.
- Running sweeps re-read config.toml at report time — never edit its VALUES
  while a sweep runs.
- `--rays` historically never reached workers (fixed 2026-07-30); every run
  before full8's resume traced 120k regardless of what was asked.
- `setRayDistributionCount1` is literal on sequences 0/3, grid-density on 1/2.
- Coincident heliostats 144=192 and 241=289 double-count ~0.3% of field power.

## Recently landed (this session, pre-commit)
Energy tab + DNI-mode selector; corrected time grid; pluggable secondaries;
flat-mirror option; Trace tab (GUI sweep launcher); ray-cost probe + measured
cost model; prime-focus model file; cassegrain design script.
