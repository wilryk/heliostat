# STATE — beamdown (rewritten 2026-07-31, relocation finished)

## Git / location — READ THIS FIRST
**Home is C:\gitlab\heliostats (this tree). Start all sessions here.**
Remote: https://github.com/wilryk/heliostat.git (branch main). Push from
this tree ONLY.
RELOCATION COMPLETE (2026-07-31): junction deleted, `analysis_output/`
(34 GB) and `legacy/old_output/` (2.2 GB) physically moved into this tree;
12/12 + test_gui PASS re-verified afterwards. ONE leftover: deleting the
old tree's stale code + .git at C:\gitlab was blocked by the permission
classifier — owner must run the rm themselves (everything tracked lives
here, verified identical incl. vet_occlusion_scalars.py; the only real
diffs were the intended path fixups). Until then, do NOT work in C:\gitlab.

## What this project is
Python package `beamdown/` drives Quadoa Optical CAD to trace a 645-heliostat
beam-down solar field for a paper comparing annualized collected energy across
geometries: 3 secondary layouts (axicon / prime_focus / cassegrain) x
(focused / flat heliostats). DNI from PVGIS TMY, monthly-mean default.

## Running RIGHT NOW
- Nothing. Licence seat FREE. full8 FINISHED 2026-07-31 04:49:
  **10,152.2 MWh annual, eta 0.5990 — THE axicon reference** (occluders
  traced, 120k rays, 12 dates, 161 steps). full7's 10,237.0 is RETIRED as
  reference (stays as read-path regression pin for test_gui): both runs
  are TRACED (scalar run in the vetting ladder is full5, NOT full7); the
  -0.829% decomposes EXACTLY as grid correction -1.175% + date coverage
  +0.351%. Old grid never sampled below el 8.78° and extrapolated 14.4%
  of annual hours high; corrected grid reaches el 1.75°, extrap 1.5%.
  Matched sun positions agree to -0.15% (noise floor). Resume boundary
  clean (4/5 declination pairs straddle it, ≤0.023%). Vetting verdict
  untouched (it compared full5 vs full7 aperture energy). Residual watch
  items: -0.15% matched-sun offset consistently negative (noise-level);
  slot overflow at el<5° adds ~+0.01% annual to full8.
- **Vetting RESOLVED** (scalar vs traced occlusion, full5/6/7 ladder):
  scalar path is 0.338% ± 0.004% LOW on annual aperture energy — one-sided,
  fully explained (eta_shade × eta_block double-charges overlapping losses;
  union form via shading.occlusion_efficiency cuts it to 0.114%). Secondary
  channel exact in aggregate (+0.002%). OWNER-APPROVED CONSEQUENCE: remaining
  comparison sweeps run WITHOUT --occluders, scalars in post. Scalars remain
  invalid for: low-sun instantaneous power (few % low), heavily-occluded
  single mirrors (+5.4%), and ALL through-focus/spot-shape work (traced only).
  Full verdict: analysis_output/vet_occlusion/verdict.txt; README section
  "Scalar vs traced occlusion"; rerunnable scripts/vet_occlusion_scalars.py.

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
- Product vs UNION form for scalar occlusion in the 5 remaining sweeps: union
  (shading.occlusion_efficiency) halves the known 0.34% low bias to 0.11%.
  Small change to how eff is baked in sweep.py's scalar branch — decide
  before launching axicon-flat. Recommendation: union.
- Chunk-size crossover above 120k rays (owner saw small-chunks-win last year;
  measured opposite below 120k) — extended probe queued, needs seat.

## Queued work (in order)
1. Owner deletes stale code + .git in C:\gitlab (classifier blocked agent).
2. Run `scripts/verify_prime_focus_model.py` (seat now free, minutes).
3. Decide union vs product scalar-occlusion form (see Open decisions), then
   launch axicon-flat sweep (no new model needed) — first of the 5.
4. Owner builds cassegrain hyperboloid in Quadoa manually (numbers above);
   prime-focus model ALREADY BUILT: `models/heliostat_field_prime_focus.optx`.
5. Prime-focus and cassegrain sweeps (x focused/flat).
6. Cross-geometry comparison report + paper figures.

## Delegation policy (owner, 2026-07-31)
Spawn subagents with EXPLICIT model. Haiku/Sonnet: data pulls, searches,
fact checks. Opus: code review, verification, ordinary implementation.
Lead (Fable) keeps: audits, hard architecture, complex implementation.
Lead also owns short/mid/long-term memory continuity across sessions.

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
