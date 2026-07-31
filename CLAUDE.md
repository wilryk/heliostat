# beamdown — session rules

## Memory protocol (do this FIRST)
Before exploring any code: read `memory/STATE.md`, then the newest file in
`memory/journal/`. That pair is the resume point — do not re-derive project
state from the codebase. At session end, rewrite STATE.md (overwrite, ≤80
lines), add one journal entry (≤40 lines), and commit. Full protocol:
`memory/README.md`.

## Invariants that protect running work
- **ONE Quadoa licence seat.** Never `--workers` > 1. Never retry a seat
  failure (each attempt pops a modal H0038 dialog). Never `import quadoa` or
  open a `beamdown.session.Session` while any `analysis_output/.*.lock`
  directory exists — that seat belongs to a multi-hour sweep.
- **Never edit config.toml VALUES while a sweep runs** — the run re-loads it
  for its end-of-run report. Comments are safe. Run options belong on the
  CLI (`python -m beamdown sweep --help`) or the GUI Trace tab, which
  propagate to workers via the override mechanism.
- **Never write into `analysis_output/<run>/` by hand** — runs are written by
  RunStore only; the lock dir is not politeness.

## Traps that produce plausible-but-wrong results
- Quadoa sequence/surface indices are 0-based (GUI is 1-based).
  `setRayDistributionCount1` is a LITERAL ray count on sequences 0/3 and a
  per-axis GRID DENSITY on 1/2 of the main model — probe per sequence.
- The trace is normalised to exactly 1000 W/m²; DNI, reflectivity and
  occlusion scalars are post-processing multipliers. Changing them never
  requires re-tracing.
- Heliostats 144=192 and 241=289 are byte-identical positions (field-file
  quirk): ~0.3% double-count, warned at load, deliberately not auto-removed.
- `cfg.optics.throughput` is applied at READ time — changing `n_mirrors`
  rescales how OLD runs read. Manifests record what was traced; trust them.

## Verification baseline
`python -m tests.verify --no-quadoa` must stay 12/12 and
`python tests/test_gui.py analysis_output/full7` must stay PASS after any
change. Axicon solve() is regression-pinned bit-identical (stage 3b);
full7 @ monthly DNI = 10,237.0 MWh is the read-path regression pin.
The axicon REFERENCE annual is full8 = 10,152.2 MWh (corrected time grid;
full7's old grid never sampled below 8.78° elevation and read high).
