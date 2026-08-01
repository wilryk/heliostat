# STATE — beamdown (rewritten 2026-08-01 ~13:00)

## Git / location
**Home is C:\gitlab\heliostats.** Remote https://github.com/wilryk/heliostat.git
(main). 2026-08-01 work is committed LOCALLY through 7bc0790 + this commit —
**NOT pushed** (classifier blocks agents; owner runs `git push`). Owner still
to hand-delete stale C:\gitlab code + .git.

## Seat / runs
**Seat FREE since 2026-08-01 12:23** — the axicon_flat → prime_focus →
prime_focus_flat chain drained cleanly, no intervention ever needed.
All four comparison runs traced (94 steps, 7 declinations, scalar UNION):
- axicon full8 **10,152.2 MWh** (eta 0.5990, traced occluders; reference)
- axicon_flat **3,781.2** (0.2231) — flat keeps 37%
- prime_focus **12,096.3** (0.7137, F1 47,000)
- prime_focus_flat **8,333.5** (0.4917) — flat keeps 69%
Concentration (cumulative field spot, traced): pf r90 474 mm / ~23,000 suns;
axicon 595 mm / ~11,000 suns — axicon wins ground-receiver concentration;
cassegrain is worst (best ideal estimate 636 mm, real will be larger).

## Cassegrain: SETTLED and script-built (2026-08-01)
Owner constraints: blocking comparable to axicon (→ F1 36,000 = tip+9 m;
union occ 0.9445 vs axicon 0.9431 at the avg-AOI instant) + 30 m dia cap +
dish as low as coverage allows. Design: rim r 15,000 at z 30,000,
**K −6.582109, |R_v| 31,548.867, vertex z 27,151.783** (design_cassegrain.py
--rim-height-mm 30000 --f1-height-mm 36000; coverage 645/645, F2 miss 1e-10).
Models built by scripts/build_cassegrain_model.py (--rim-z-mm/--f1-mm):
- models/heliostat_field_cassegrain.optx (sweep; run_cassegrain{,_flat}.sh)
- models/figure_model_25cfg_20260220_0927_cass30.optx (GUI review)
TRAP: in these models sec_height/rec_offset are UNREFERENCED literals —
session writes to them are no-ops; never run an axicon sweep against them.
**Licence-gated verification PENDING before first cassegrain sweep**
(self-test trace; centroid prediction still (0,0) — stigmatic relay).

## Fixed mirror figures (owner's scenarios, landed)
`sweep --fixed-shapes <csv>` freezes c3/c4/c5 per heliostat (pointing still
tracks); mutually exclusive with --flat-mirrors; manifest key fixed_shapes.
Tables: data/fixed_shapes_pf36000_{spherical,mean_cos,median}.csv (annual
c3≈0 by symmetry; mean_cos vs median differ ~40% in c5 → test needed).
4-day family READY, not launched: run_prime_focus_f36{,_sphere,_meancos,
_median}.sh — F1 36,000, dates 12-21/02-20/04-21/06-21, ~6 h each.

## Design studies (all licence-free, scripts committed)
- AOI: annual DNI-weighted 33.3°; representative instant **2026-02-20
  hour 9.45401542341586** (el 42.24); axicon aim rays cross axis tip+2.6
  to +10.0 m (mean +7.8) — owner's 9 m estimate confirmed (aoi_stats.py).
- Estimator (cos × union occ × rho^n, 94-grid, annual_energy): validated
  axicon +0.64%, prime focus +0.19% vs traced (scan_cassegrain_annual.py).
- Axicon 27,000/20° is Pareto-efficient under the owner's sagittal cap
  (7.115e-6 /mm at inner mirrors; hits at 1,943 mm). Higher tip = +energy
  but inner hits collapse (33 km tip: 658 mm, crowding 1.41 — unusable).
  Angle pinned by coverage to ~18–21° (scan_axicon_annual.py, exact).
- beamdown/design_eval.py is the single source of design math (GUI + all
  scripts import it).

## GUI (committed d326cc7, 7bc0790)
Design tab: layout radio + sliders, live cross-section from real sag
equations, plain-language readouts, sagittal-cap traffic light, .optx
export buttons (any explored geometry; never --force).
plot_style.py paper style everywhere (white, 2 pt lines, tight margins);
every tab: "Save figure…" (600 dpi PNG + vector PDF) + "Save data (CSV)".

## Suite status
verify --no-quadoa 12/12 (figure stage 21 checks); test_gui full7 PASS
(new Design-tab + export sections). Axicon solve() bit-identical pin holds.

## Queued (in order)
1. Owner: GUI review of cass30 model; owner runs `git push`.
2. Seat: verify 25cfg models (verify_figure_model.py) + cassegrain sweep
   model self-test. Minutes each; do BEFORE any cassegrain sweep.
3. 4-day fixed-figure runs (owner picks order) — then mean-vs-median call.
4. run_cassegrain.sh + run_cassegrain_flat.sh (settled design, ~9 h each).
5. Cross-geometry paper report incl. concentration column (traced r50/r90).
6. Chunk-size probe >120k rays (needs seat).

## Traps (new this session; older ones in CLAUDE.md/README)
- Heavy CPU jobs slow a running sweep even at below-normal priority
  (measured ~30%: 397→533 s/step). Nothing heavy beside a sweep.
- numpy 2 repr: format CSV cells with float(v), never repr(np scalar).
- Owner prefs: plain-language summaries; paper plot style; sagittal cap
  as a hard design rule (see private memory).
