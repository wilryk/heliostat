# 2026-08-01 — four runs done, cassegrain settled, design GUI

Attempted: absorb the finished sweep chain; answer the owner's AOI /
cassegrain-design questions; land the fixed-mirror-figure scenarios; make
design exploration self-serve in the GUI.

Done (verified):
- **Chain drained clean, zero intervention**: prime_focus 12,096.3 MWh
  (0.7137), prime_focus_flat 8,333.5 (0.4917). Frozen UI chip ≠ dead
  process — log mtime + process list are ground truth.
- **AOI study**: annual DNI-weighted 33.3°; representative instant
  2026-02-20 09:27; axicon aim rays cross axis at tip+7.8 m mean —
  owner's "9 m" confirmed; F1=36,000 fixed for all pf/cassegrain work.
- **Cassegrain SETTLED** (blocking-comparability + 30 m cap + lowest
  feasible dish): rim 15 m @ z 30,000, F1 36,000, K −6.582109,
  |R_v| 31,548.867, vertex 27,151.783. Models script-built
  (build_cassegrain_model.py); verification pending seat.
- **Concentration adjudicated with traced data** (owner pushed back
  twice, correctly): cumulative field spot pf 474 mm / axicon 595 /
  cassegrain ≥636 ideal. Axicon wins on the ground; the relay
  magnification (≥1.7, forced by blocking floor + coverage) is why.
- **Axicon scan**: built 27k/20° is Pareto-efficient under the owner's
  sagittal cap; angle pinned ~18–21° by coverage. Estimator validated
  (+0.64% / +0.19% vs the two traced runs).
- **Fixed figures**: sweep --fixed-shapes + three pf36000 tables +
  4-day test family (not launched). mean_cos vs median differ 40% in c5.
- **GUI**: Design tab (design_eval.py single-source math, live
  cross-section, .optx exports) + plot_style.py (white/2pt/tight, 600 dpi
  PNG + PDF + CSV export on every tab). Both by the same Opus implementer,
  reviewed; it caught my numpy-2 repr bug in the CSV generator.

Owner decisions: sagittal correction ≤ built axicon's is a HARD rule;
30 m secondary diameter cap; comparable blocking beats +1.4% energy;
plain-language summaries; paper plot style (private memory updated).

Loose ends: verifications before cassegrain sweep (seat, minutes);
4-day run ordering owner's call; push blocked for agents — owner pushes;
cross-geometry report should carry a concentration column.
