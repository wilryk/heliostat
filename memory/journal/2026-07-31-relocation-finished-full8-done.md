# 2026-07-31 — relocation finished, full8 done

Attempted: complete the C:\gitlab → C:\gitlab\heliostats relocation and
absorb the full8 result.

Done (verified):
- **Old-vs-new tree audit**: 42 file diffs were line-endings only; the 11
  real diffs were exactly the intended path fixups (REPO=, ModelFolder,
  cd targets) plus newer STATE/journal here. vet_occlusion_scalars.py
  committed and identical. Nothing tracked exists only in the old tree.
- **Junction deleted; data physically moved** (same-volume rename):
  analysis_output (34 GB, 25 entries) and legacy/old_output (2.2 GB) now
  live in this tree. 12/12 verify + test_gui full7 PASS re-run afterwards.
- **full8 FINISHED 04:49**: annual 10,152.2 MWh, eta 0.5990, sine fit
  R²=0.878 peak doy 201, declination pairs ±0.0002, worst residual 1.50%.
  Flag: traced full8 is 0.83% BELOW scalar full7 (10,237.0) — opposite
  sign to vetting's "scalar 0.338% low"; full7 used the old time grid so
  not directly comparable. Owner to adjudicate the axicon reference number.

Blocked: deleting old tree's code + .git — permission classifier refused
recursive rm twice (both bash and PowerShell). Owner will delete manually.
Old tree is stale-but-present; do not work there.

Owner decisions (this session): subagent delegation ladder with explicit
model per spawn (Haiku/Sonnet fact-work, Opus review+implementation, Fable
leads high-judgment) — recorded in STATE and in coordinator private memory
(delegation-model-tiers.md). Lead is responsible for short/mid/long-term
memory continuity.

Rest of session (compressed): full8-vs-full7 ADJUDICATED — my "opposite
sign to vetting" flag was a false alarm on a wrong premise (full7 is
TRACED; the scalar rung is full5). The -0.83% = grid correction (-1.175%)
+ 12-vs-4 dates (+0.351%), exact; matched-sun eta agrees to -0.15%
(noise). full8 10,152.2 MWh = axicon reference; full7 = regression pin
(CLAUDE.md updated). UNION occlusion form landed (e4f0e0f, owner
decision): manifest occlusion_form key, all readers via
store.occlusion_weight_columns; bonus fix — GUI single-heliostat flux
double-charged occlusion on traced runs. Prime-focus model VERIFIED
10/10. axicon_flat LAUNCHED (94 steps, 7 explicit declinations; first
launch used --suggest-dates which ADDS to config's dates — 214 steps,
killed at 2 min, partial deleted, trap recorded). Chain watcher armed:
axicon_flat → prime_focus → prime_focus_flat (scripts written, chains
only on clean "done"). Harness landed: .claude/agents fetcher=haiku /
implementer=opus / reviewer=opus + permission allowlist; delegation
ladder + memory stewardship in coordinator private memory.

Late: 25cfg figure model solved (Opus agent) — shipped .optx never
populated (generator's licence half never ran; configs 1-24 zeros). New
scripts/build_figure_model.py (licence-free, exact round-trip 12/12) +
deferred verify_figure_model.py, committed c43c27e. Sun is single_param:
one instant per model, by construction.

Loose ends: chain is session-tied (recovery steps in STATE queued work
item 1); owner still to delete old C:\gitlab code+.git; cassegrain model
is owner's manual build; chunk-size probe still queued behind the seat;
verify_figure_model.py to run when the seat frees.
