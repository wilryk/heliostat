# 2026-07-30 — rays bug, power loss, new layouts ready

Attempted: keep full8 alive; make the 3x2 geometry comparison runnable.

Landed (all verified, 12/12 + test_gui PASS; initial commit 860d1a9, pushed
to github.com/wilryk/heliostat):
- **--rays override bug**: CLI ray count NEVER reached sweep workers (they
  reload config.toml); every historical run traced 120k rays regardless.
  Fixed (workers replay full override dict, recorded in manifest). Every
  prior "ray count vs speed" claim was measuring 120k vs 120k — retracted
  in README/config/run_full8.sh.
- **Ray cost measured for the first time** (`scripts/probe_ray_cost.py`,
  15-cell grid): setup 0.26 ms/heliostat, 5.29 us/ray, 94% of trace time is
  rays at 120k. Chunking pure loss (646/698/776 ms for 1/2/4 calls). Probe
  aborts if emitted != requested (grid-density trap).
- **Power event** killed full8 at 93/161 (exit 127, licence key dropout).
  Resumed with --rays 120000 --rays-per-trace 120000: homogeneous run,
  manifest self-corrected (was falsely 60000), one-call speedup free.
- **Flat-mirror option** (--flat-mirrors): zeroes heliostat Zernike c3/c4/c5,
  pointing bit-identical; ast-based guard makes bypass impossible. Verified
  zeroing == deactivating the form (base radius inf, only active form).
- **Prime-focus model** `models/heliostat_field_prime_focus.optx` built by
  guarded text surgery: detector z -> single_param pf_height=47000, seq 3
  rewired to sun->helio_surf->prime_focus. Deferred licence checks in
  `scripts/verify_prime_focus_model.py`. Base .optx is NOT well-formed XML
  (duplicate attr, pre-existing); float_ap proven non-clipping from full7 rays.

Owner decisions: 120k rays stays (single-heliostat SNR); 7 declinations not
12 dates for comparisons; prime focus at 47 m (symmetric throw), NOT at
cassegrain F1; scalar occlusion for remaining sweeps IF vetting passes.

Failed/withdrawn: my "rays are nearly free" cost model (artifact of the
override bug); 14 h full8 estimate (real: ~24 h).

Loose ends: vetting agent in flight (scalar vs traced verdict); extended
chunk probe >120k queued (owner's contrary prior experience); cassegrain
.optx is owner's manual build; full8 finishes ~19:00, auto-reports to log.

Late addition: VETTING VERDICT landed — scalar 0.338% low (double-charged
shade x block overlap; union form 0.114%), secondary scalar exact in
aggregate, traced-only for through-focus. Remaining sweeps go scalar.
Also: repo relocated to C:\gitlab\heliostats (clone + analysis_output
junction; old tree is a husk hosting full8 until done). Open: bake union
vs product form into the 5 comparison sweeps — decide before axicon-flat.
