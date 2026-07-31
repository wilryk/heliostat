# STATE — beamdown (rewritten 2026-07-31 ~07:20, session end)

## Git / location
**Home is C:\gitlab\heliostats. Start all sessions here.** Remote:
https://github.com/wilryk/heliostat.git (main, pushed through 704aa08).
Relocation FINISHED 2026-07-31: analysis_output (34 GB) + legacy/old_output
physically here, junction gone, baseline re-verified. One leftover: OWNER
must hand-delete the stale code + .git in C:\gitlab (classifier blocks
agents). Do not work in C:\gitlab.

## What this project is
Package `beamdown/` drives Quadoa to trace a 645-heliostat beam-down field
for a paper: annualized energy across (axicon / prime_focus / cassegrain)
x (focused / flat heliostats). DNI: PVGIS TMY monthly-mean, 1,751.9 kWh/m2.

## Running RIGHT NOW (as of 2026-07-31 ~16:00)
- **prime_focus sweep** (launched 15:56 by the chain watcher): 645 x 94
  steps, focused, F1 47000, n_mirrors 1, verified .optx, union scalars —
  header overrides all confirmed. Holds the ONLY seat. Expect ~9 h (ends
  ~01:00 Aug 1); then the watcher auto-launches **prime_focus_flat**.
- **axicon_flat FINISHED 15:55 clean: 3,781.2 MWh, eta 0.2231** (37% of
  focused full8 — spillage, as designed). Worst residual 0.79%, R² 0.88,
  extrap 1.5%. Report at end of analysis_output/axicon_flat.log.
- Sweep + watcher are background children of the 2026-07-31 session. If a
  later session finds a lock but NO live python: stale lock — remove the
  lock dir, relaunch the same run script (sweeps RESUME from stored
  steps), rerun scripts/run_after_axicon_flat.sh (idempotent).

## Reference numbers
- **Axicon reference: full8 = 10,152.2 MWh, eta 0.5990** (traced occluders,
  corrected grid, 120k rays homogeneous, resume boundary verified clean).
- full7 = 10,237.0 MWh is the read-path REGRESSION PIN only: its old grid
  never sampled below el 8.78°, extrapolated 14.4% of hours, read high.
  The -0.83% delta = grid -1.175% + dates +0.351%, exact (journal 07-31).
- Vetting (full5 scalar vs full7 traced, aperture energy): product form
  0.338% low; UNION form 0.114%. Scalars invalid for low-sun instantaneous
  power, single heavily-occluded mirrors, and all spot-shape work.

## Landed 2026-07-31 (all pushed)
- UNION occlusion form (e4f0e0f): sweep scalar branch uses
  shading.occlusion_efficiency; manifest key occlusion_form (absent =
  historical product); ALL readers via store.occlusion_weight_columns.
  Bonus fix: GUI heliostat-flux double-charged occlusion on traced runs.
- Prime-focus model VERIFIED 10/10 (scripts/verify_prime_focus_model.py).
- Harness: .claude/agents fetcher(haiku)/implementer(opus)/reviewer(opus)
  — load at session START; permission allowlist .claude/settings.json.

## Suite status
`python -m tests.verify --no-quadoa` → 12/12; `python tests/test_gui.py
analysis_output/full7` → PASS (must stay so after any change). Axicon
solve() bit-identical (stage 3b).

## Queued work (in order)
1. Chain drains: axicon_flat → prime_focus → prime_focus_flat (automated,
   see Running). Then read the three report_energy outputs.
2. Owner deletes stale C:\gitlab code + .git.
3. Owner builds cassegrain hyperboloid .optx by hand (rim z 32,460 r
   15,000; F1 38,986; K -5.8789; |R_v| 32,181.3; vertex z 29,589 —
   scripts/design_cassegrain.py). Then cassegrain focused/flat sweeps:
   clone run_prime_focus*.sh with --secondary cassegrain,
   --focus-height-mm 38986, --rim-height-mm 32460, --n-mirrors 2.
4. Cross-geometry comparison report + paper figures. Figure models: 25cfg
   MYSTERY SOLVED 2026-07-31 — shipped figure_model_25cfg.optx was never
   populated (configs 1-24 all zeros; old model_edit.build_figure_model
   needed a seat for its second half, which never ran). NEW licence-free
   generator: `python scripts/build_figure_model.py --date D --hour H
   [--flat] --check` (25 configs = the 25 downselected heliostats at one
   instant; sun is a single_param, CANNOT vary per config). Once the seat
   frees: run `python scripts/verify_figure_model.py <built.optx>` (it
   refuses under any lock). Old model_edit.build_figure_model is dead
   code — removal offered as a spawned task chip.
5. Chunk-size probe >120k rays (owner's contrary prior) — needs seat.

## Delegation policy (owner, 2026-07-31)
Spawn subagents with EXPLICIT model. Haiku/Sonnet: data pulls, searches,
fact checks. Opus: review, verification, ordinary implementation. Lead
(Fable): audits, architecture, complex implementation, memory continuity.

## Traps that bite (details in CLAUDE.md / README)
- ONE licence seat: never workers>1, never retry seat failures, never
  import quadoa while any analysis_output/.*.lock exists.
- Never edit config.toml VALUES while a sweep runs (report re-reads it).
- `--suggest-dates N` ADDS to config's 12 dates (does NOT replace) — use
  explicit --dates; first axicon_flat launch died to this, 2 min in.
- `setRayDistributionCount1`: literal on sequences 0/3, grid-density 1/2.
- Coincident heliostats 144=192, 241=289: ~0.3% double-count, deliberate.
- cfg.optics.throughput applies at READ time; manifests record the truth.
