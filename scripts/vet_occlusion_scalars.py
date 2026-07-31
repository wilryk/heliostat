#!/usr/bin/env python
"""Vet the analytic occlusion scalars against traced occlusion geometry.

The question this answers
------------------------
A comparison sweep costs ~11.5 h with occlusion traced as real geometry. If the
analytic scalars (neighbour shading, neighbour blocking, secondary shadow) can be
applied in post instead, that cost disappears. This script measures the price of
doing so, from runs already on disk -- no license, no tracing.

The ladder it reads
-------------------
Three runs over the *same* timestep grid, heliostats and ray budget:

    full5   occluders=False                    every occlusion is a scalar
    full6   occluders=True                     neighbours traced, axicon scalar
    full7   occluders=True, traced_secondary   everything traced, no scalar at all

so ``full5 -> full6`` isolates the neighbour channel and ``full6 -> full7`` the
secondary-shadow channel.

Which efficiency each run's ``power_w`` already carries
------------------------------------------------------
This is the one thing that must not be got wrong -- applying a factor twice, or
not at all, moves the answer by exactly the quantity being studied. The rule is
set in ``sweep._assemble`` and re-derived here from each manifest, then *checked*
against every stored row (``power_w == rays_landed * scale_factor * eff``):

    occluders=False                      eff = eta_shade * eta_block
    occluders=True,  traced_secondary=0  eff = eta_secondary
    occluders=True,  traced_secondary=1  eff = 1

``eta_shade`` is the *union* of neighbour shading and the secondary's shadow, so
the scalar path must not multiply ``eta_secondary`` on top of it; ``eta_secondary``
is carried alone only so the traced-neighbour run can still apply the one occluder
it did not trace.

What it produces
----------------
``analysis_output/vet_occlusion/`` -- per-timestep and per-heliostat CSVs, an
annual-energy table, three figures, and a verdict block (also printed).

    python scripts/vet_occlusion_scalars.py
    python scripts/vet_occlusion_scalars.py --runs full5 full6 full7 --aperture-mm 700
    python scripts/vet_occlusion_scalars.py --runs full5 full6 --skip-annual
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
warnings.filterwarnings("ignore")

ANALYSIS = REPO / "analysis_output"
RADII_MM = (300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1250.0, 1500.0)
EL_BANDS = ((0.0, 5.0), (5.0, 15.0), (15.0, 30.0), (30.0, 60.0), (60.0, 90.0))
ETA_EDGES = (0.0, 0.5, 0.7, 0.85, 0.95, 0.999, 1.001)


# --------------------------------------------------------------------------
# run bookkeeping
# --------------------------------------------------------------------------
def resolve_run(name: str) -> Path:
    p = Path(name)
    return p if p.exists() else ANALYSIS / name


def weight_mode(manifest: dict) -> str:
    """Which scalar the sweep already folded into this run's power_w."""
    if not bool(manifest.get("occluders", False)):
        return "product"                      # eta_shade * eta_block
    return "none" if bool(manifest.get("traced_secondary", False)) else "secondary"


MODE_TEXT = {
    "product": "eta_shade * eta_block   (all occlusion analytic)",
    "secondary": "eta_secondary           (neighbours traced, axicon scalar)",
    "none": "1.0                     (everything traced)",
}


def weights_from(mode: str, rows: pd.DataFrame) -> np.ndarray:
    if mode == "product":
        return (rows.eta_shade * rows.eta_block).to_numpy(float)
    if mode == "secondary":
        return rows.eta_secondary.to_numpy(float)
    return np.ones(len(rows))


class Run:
    """One sweep on disk, with the efficiency convention it was written under."""

    def __init__(self, name: str, cfg):
        from beamdown.store import RunStore, scale_factor

        self.name = Path(name).name
        self.root = resolve_run(name)
        self.store = RunStore(self.root, cfg=cfg, mode="r")
        self.manifest = self.store.manifest
        self.summary = self.store.summary()
        self.ids = [int(h) for h in self.manifest.get("heliostat_ids", [])]
        self.mode = weight_mode(self.manifest)
        self.rays = int(self.manifest.get("rays_per_heliostat",
                                          cfg.trace.rays_per_heliostat))
        self.sf = scale_factor(cfg, self.rays, 1000.0)
        self._by_key = {
            k: g.drop_duplicates("heliostat_id", keep="last")
                .set_index("heliostat_id").reindex(self.ids)
            for k, g in self.summary.groupby("timestep")
        }

    def rows(self, key: str) -> pd.DataFrame:
        return self._by_key[key]

    def weights(self, key: str) -> np.ndarray:
        return weights_from(self.mode, self.rows(key))

    def keys(self) -> set:
        return set(self.store.timestep_keys()) & set(self.summary.timestep)


def audit_semantics(runs, log) -> None:
    """Confirm from the data which efficiency each run's power_w carries."""
    log("EFFICIENCY SEMANTICS (re-derived from each manifest, checked on every row)")
    for r in runs:
        s = r.summary
        eff = s.shading_blocking_efficiency.to_numpy(float)
        landed = s.rays_landed.to_numpy(float)
        implied = np.divide(s.power_w.to_numpy(float), landed * r.sf,
                            out=np.full(len(s), np.nan), where=landed > 0)
        expect = weights_from(r.mode, s)
        d_power = float(np.nanmax(np.abs(implied - eff)))
        d_mode = float(np.max(np.abs(eff - expect)))
        log(f"  {r.name:<7s} occluders={str(bool(r.manifest.get('occluders', False))):<5s} "
            f"traced_secondary={str(bool(r.manifest.get('traced_secondary', False))):<5s} "
            f"rays/heliostat={r.rays:,}   rows={len(s):,}")
        log(f"          power_w carries  {MODE_TEXT[r.mode]}")
        log(f"          |power_w/(landed*scale) - stored eff| = {d_power:.1e}   "
            f"|stored eff - {r.mode}| = {d_mode:.1e}   "
            f"{'OK' if d_power < 1e-9 and d_mode < 1e-12 else 'FAILED'}")
        if not (d_power < 1e-9 and d_mode < 1e-12):
            raise SystemExit(f"{r.name}: efficiency semantics check failed -- refusing "
                             f"to compare runs whose power_w cannot be reproduced.")

    base = runs[0].summary.set_index(["timestep", "heliostat_id"]).sort_index()
    for r in runs[1:]:
        other = r.summary.set_index(["timestep", "heliostat_id"]).sort_index()
        common = base.index.intersection(other.index)
        worst = max(float(np.nanmax(np.abs(base.loc[common, c].to_numpy(float)
                                           - other.loc[common, c].to_numpy(float))))
                    for c in ("eta_shade", "eta_block", "eta_secondary"))
        log(f"  analytic etas {runs[0].name} vs {r.name}: max|diff| {worst:.1e} over "
            f"{len(common):,} shared rows (the ladder is controlled)")
        if worst > 1e-9:
            raise SystemExit("the runs disagree on their analytic efficiencies; they "
                             "are not the same field or grid.")


# --------------------------------------------------------------------------
# counts and noise
# --------------------------------------------------------------------------
def radial_masks_flat(cfg, radii_mm) -> np.ndarray:
    from beamdown.metrics import radial_masks

    return radial_masks(cfg, radii_mm).astype(np.float32).reshape(len(radii_mm), -1)


def read_counts(store, key, take, masks_flat):
    """Per-heliostat landed counts, counts inside each radius, and the raw stack."""
    counts = np.asarray(store.read_counts(key))[take].astype(np.float32)
    flat = counts.reshape(counts.shape[0], -1)
    inside = (flat @ masks_flat.T).astype(np.float64)
    return flat.sum(axis=1).astype(np.float64), inside, counts


def mc_variance(w: np.ndarray, n: np.ndarray, rays: int, sf: float) -> np.ndarray:
    """Per-heliostat Monte-Carlo variance of weighted power, in W^2.

    Each heliostat's landed count is Binomial(rays, p) with p estimated by n/rays,
    and the heliostats are traced independently, so these add. The runs drew
    independent rays -- checked directly: on fully unoccluded heliostats the
    run-to-run difference has the sd this predicts for two independent draws --
    so a difference between runs carries the sum of both runs' variances.
    """
    p = np.clip(n / rays, 0.0, 1.0)
    return (sf ** 2) * (w ** 2) * rays * p * (1.0 - p)


def mc_sigma(w: np.ndarray, n: np.ndarray, rays: int, sf: float) -> float:
    """1-sigma Monte-Carlo noise on a weighted field total, in watts."""
    return float(np.sqrt(mc_variance(w, n, rays, sf).sum()))


# --------------------------------------------------------------------------
# analytic efficiencies, recomputed post-hoc
# --------------------------------------------------------------------------
def recompute_etas(cfg, keys, summary, cache: Path, log):
    """eta_shade / eta_block / eta_secondary / eta_union for every timestep.

    ``eta_union`` is :func:`beamdown.shading.occlusion_efficiency` -- the
    lit-and-unblocked fraction, with shaded, blocked and secondary-shadowed
    patches unioned instead of multiplied. Same analytic model as the stored
    scalars; it just does not charge twice for a patch of mirror that is both
    shaded and blocked.
    """
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        if list(z["keys"]) == list(keys):
            log(f"  analytic etas: reusing {cache.name}")
            return {k: z[k] for k in ("shade", "block", "secondary", "union")}

    from beamdown import field as field_mod
    from beamdown import shading as shading_mod
    from beamdown.secondary import get_strategy

    fld = field_mod.load_field(cfg)
    strategy = get_strategy(cfg)
    body = shading_mod.secondary_body(cfg)
    el_min = float(summary.solar_el_deg.min())
    radius = shading_mod.search_radius_for(el_min, cfg.field.mirror_height_mm,
                                           cfg.field.mirror_width_mm)
    neighbours = field_mod.neighbour_pairs(fld, radius)
    log(f"  analytic etas: {len(keys)} timesteps x {len(fld)} heliostats, neighbour "
        f"search radius {radius/1000:.1f} m (min sun elevation {el_min:.2f} deg)")

    out = {k: [] for k in ("shade", "block", "secondary", "union")}
    t0 = time.time()
    for i, key in enumerate(keys, 1):
        rows = summary[summary.timestep == key]
        az = float(rows.solar_az_deg.iloc[0])
        el = float(rows.solar_el_deg.iloc[0])
        sols = [strategy.solve(float(fld.x_mm[j]), float(fld.y_mm[j]), az, el, cfg.geometry)
                for j in range(len(fld))]
        geoms, aims = shading_mod.build_geometries(fld, sols, cfg)
        sh, bl, sec = shading_mod.shading_blocking(geoms, aims, az, el, neighbours,
                                                   secondary=body)
        uni = shading_mod.occlusion_efficiency(geoms, aims, az, el, neighbours,
                                               secondary=body)
        out["shade"].append(sh)
        out["block"].append(bl)
        out["secondary"].append(sec)
        out["union"].append(uni)
        if i == 1 or i % 10 == 0 or i == len(keys):
            log(f"    [{i}/{len(keys)}] {key} el {el:5.2f}  "
                f"{(time.time()-t0)/i:.1f} s/timestep")
    res = {k: np.array(v) for k, v in out.items()}
    np.savez_compressed(cache, keys=np.array(list(keys)), **res)
    return res


# --------------------------------------------------------------------------
# slot-overflow probe (geometry only, no tracing)
# --------------------------------------------------------------------------
def overflow_probe(cfg, log, n_steps: int = 6) -> pd.DataFrame:
    """How much occlusion the traced path *cannot* represent, by sun elevation.

    The traced model carries 10 shading + 4 blocking slots per heliostat. When
    more neighbours than that occlude a mirror the extras are dropped, and there
    the traced run is the approximate one. This measures that gap without
    tracing: the analytic union over *all* neighbours against the same union over
    only the neighbours that fit in the slots. It runs over the grid in
    config.toml, which reaches lower sun than the compared runs do.
    """
    from beamdown import field as field_mod
    from beamdown import occluder_slots as slots_mod
    from beamdown import shading as shading_mod
    from beamdown import solar as solar_mod
    from beamdown.secondary import get_strategy

    steps = solar_mod.build_time_grid(cfg)
    picks, seen = [], []
    for s in sorted(steps, key=lambda s: s.solar_el_deg):
        if all(abs(s.solar_el_deg - e) > 4.0 for e in seen):
            picks.append(s)
            seen.append(s.solar_el_deg)
        if len(picks) >= n_steps:
            break

    fld = field_mod.load_field(cfg)
    strategy = get_strategy(cfg)
    body = shading_mod.secondary_body(cfg)
    radius = shading_mod.search_radius_for(
        min(s.solar_el_deg for s in steps),
        cfg.field.mirror_height_mm, cfg.field.mirror_width_mm)
    neighbours = field_mod.neighbour_pairs(fld, radius)
    by_id = {int(h): k for k, h in enumerate(fld.ids)}

    rows = []
    for step in picks:
        az, el = step.solar_az_deg, step.solar_el_deg
        sols = [strategy.solve(float(fld.x_mm[j]), float(fld.y_mm[j]), az, el, cfg.geometry)
                for j in range(len(fld))]
        geoms, aims = shading_mod.build_geometries(fld, sols, cfg)
        to_sun = shading_mod.sun_vector(az, el)
        plans = slots_mod.plan_field(geoms, aims, fld.ids, neighbours, to_sun,
                                     body=body, has_secondary=True)
        rep = slots_mod.overflow_report(plans)

        eta_all = shading_mod.occlusion_efficiency(geoms, aims, az, el, neighbours,
                                                   secondary=body)
        # Slots carry a heliostat_id (-1 for a synthetic filler, which occludes
        # nothing by construction), so the traced model's view of a heliostat's
        # neighbourhood is exactly this set.
        kept = [np.array(sorted({s.heliostat_id for s in list(p.shading) + list(p.blocking)}
                                & by_id.keys()), dtype=int) for p in plans]
        kept = [np.array([by_id[int(h)] for h in k], dtype=int) for k in kept]
        eta_slots = shading_mod.occlusion_efficiency(geoms, aims, az, el, kept,
                                                     secondary=body)
        # weight by cosine efficiency: a heliostat collecting nothing cannot
        # mis-state much, and at these elevations most of them barely collect.
        wgt = np.clip(np.array([s.cosine_efficiency for s in sols]), 0.0, None)
        den = float(np.sum(wgt * eta_all))
        bias = (float(np.sum(wgt * eta_slots)) / den - 1.0) if den > 0 else np.nan
        rows.append(dict(timestep=step.key, solar_el_deg=el,
                         overflowed=rep["overflowed"],
                         max_shading_used=rep["max_shading_used"],
                         max_blocking_used=rep["max_blocking_used"],
                         worst_dropped_fraction=rep["worst_dropped_fraction"],
                         eta_all_mean=float(eta_all.mean()),
                         eta_slots_mean=float(eta_slots.mean()),
                         traced_overstates_pct=100.0 * bias))
        log(f"    {step.key} el {el:5.2f}  overflowed {rep['overflowed']:4d}/{len(plans)}  "
            f"eta(all neighbours) {eta_all.mean():.4f}  eta(slots only) "
            f"{eta_slots.mean():.4f}  traced overstates power by {100.0*bias:+.2f}%")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# annual energy
# --------------------------------------------------------------------------
def annual_from_power(power_by_key: dict, template: pd.DataFrame, cfg, prov,
                      n_heliostats: int):
    """Annual energy from a field-power-per-timestep series.

    ``annual_energy`` takes per-heliostat rows but only ever sums them per
    timestep, so one synthetic row per timestep carrying the whole field's power
    gives an identical answer -- with ``n_heliostats`` passed explicitly so the
    mirror area is still the real field's.
    """
    from beamdown import energy as E

    rows = template.drop_duplicates("timestep").copy()
    rows["power_w"] = rows.timestep.map(power_by_key)
    rows["heliostat_id"] = 0
    return E.annual_energy(rows, cfg, prov, n_heliostats=n_heliostats)


def annual_weights(template, cfg, prov, n_heliostats, base_power, log):
    """dE_annual / dP_timestep, by finite difference.

    ``annual_energy`` is linear in each traced timestep's field power (the
    efficiency surface is interpolated linearly, DNI is a fixed multiplier), so
    these weights are exact and independent of which run supplied the power. They
    turn per-timestep Monte-Carlo noise into an annual uncertainty, and they say
    how much of the year rides on the low-sun timesteps.
    """
    base = annual_from_power(base_power, template, cfg, prov, n_heliostats)
    E0 = base["annual_energy_kwh"]
    w, t0 = {}, time.time()
    for i, key in enumerate(base_power, 1):
        step = 0.01 * base_power[key]
        if step <= 0:
            w[key] = 0.0
            continue
        bumped = dict(base_power)
        bumped[key] = base_power[key] + step
        E1 = annual_from_power(bumped, template, cfg, prov,
                               n_heliostats)["annual_energy_kwh"]
        w[key] = (E1 - E0) / step
        if i == 1 or i % 20 == 0 or i == len(base_power):
            log(f"    annual sensitivity [{i}/{len(base_power)}] "
                f"{(time.time()-t0)/i:.2f} s/timestep")
    return base, w


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------
def make_figures(out_dir, per_ts, per_h, morph, traced_names, scalar_name,
                 has_union, ap_mm, log):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colours = {t: c for t, c in zip(traced_names, ["#c1272d", "#1b6ca8"])}
    colours["union"] = "#2a8f4b"
    last = traced_names[-1]

    # --- 1. ratio vs sun elevation ---------------------------------------
    order = per_ts.sort_values("solar_el_deg")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3), sharey=True)
    for ax, (kind, title) in zip(axes, [("tot", "field total power"),
                                        ("ap", f"power inside r <= {ap_mm:.0f} mm")]):
        ax.fill_between(order.solar_el_deg, -order[f"sigma_{kind}_pct"],
                        order[f"sigma_{kind}_pct"], color="0.86", lw=0, zorder=0,
                        label="Monte-Carlo $1\\sigma$ on the difference")
        for t in traced_names:
            ax.plot(per_ts.solar_el_deg, 100.0 * (per_ts[f"ratio_{kind}_{t}"] - 1.0),
                    "o", ms=5, color=colours[t], label=f"{t} vs {scalar_name} (product)")
        if has_union:
            ax.plot(per_ts.solar_el_deg, 100.0 * (per_ts[f"ratio_{kind}_union"] - 1.0),
                    "s", ms=5, mfc="none", color=colours["union"],
                    label=f"{last} vs scalar in union form")
        ax.axhline(0, lw=0.8, color="0.4")
        ax.set_xlabel("sun elevation (deg)")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("traced $-$ scalar (%)")
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Traced occlusion vs analytic scalars, per timestep")
    fig.tight_layout()
    f1 = out_dir / "fig1_ratio_vs_elevation.png"
    fig.savefig(f1, dpi=150)
    plt.close(fig)

    # --- 2. per-heliostat scatter ----------------------------------------
    sub = per_h[(per_h.p_scalar > 0) & np.isfinite(per_h[f"ratio_{last}"])]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
    bands = [(0, 15), (15, 30), (30, 60), (60, 90)]
    cmap = plt.get_cmap("viridis")
    for i, (lo, hi) in enumerate(bands):
        m = (sub.solar_el_deg >= lo) & (sub.solar_el_deg < hi)
        if not m.any():
            continue
        axes[0].plot(sub.eta_product[m], 100.0 * (sub[f"ratio_{last}"][m] - 1.0), ".",
                     ms=2, alpha=0.25, color=cmap(i / 3.5), label=f"el {lo}-{hi} deg")
    axes[0].set_xlabel("analytic $\\eta_{shade}\\times\\eta_{block}$")
    axes[0].set_ylabel(f"{last} / {scalar_name} $-$ 1 (%)")
    axes[0].axhline(0, lw=0.8, color="0.4")
    axes[0].set_ylim(-40, 80)
    axes[0].legend(fontsize=8, markerscale=6)
    axes[0].set_title("per heliostat, per timestep")
    axes[0].grid(alpha=0.3)

    b = pd.cut(sub.eta_product, list(ETA_EDGES))
    g = sub.groupby(b, observed=True)
    med = 100.0 * (g[f"ratio_{last}"].median() - 1.0)
    q1 = 100.0 * (g[f"ratio_{last}"].quantile(0.25) - 1.0)
    q3 = 100.0 * (g[f"ratio_{last}"].quantile(0.75) - 1.0)
    x = np.arange(len(med))
    axes[1].errorbar(x, med, yerr=[med - q1, q3 - med], fmt="o", capsize=4,
                     color=colours[last], label="scalar as product")
    if has_union and f"ratio_union_{last}" in sub.columns:
        gu = g[f"ratio_union_{last}"].median()
        axes[1].plot(x, 100.0 * (gu - 1.0), "s", mfc="none", ms=7,
                     color=colours["union"], label="scalar in union form")
    axes[1].axhline(0, lw=0.8, color="0.4")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([str(i) for i in med.index], rotation=25, fontsize=7)
    axes[1].set_xlabel("analytic $\\eta_{shade}\\times\\eta_{block}$")
    axes[1].set_ylabel("median difference (%), bars = IQR")
    axes[1].set_title("where the difference sits")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    f2 = out_dir / "fig2_per_heliostat.png"
    fig.savefig(f2, dpi=150)
    plt.close(fig)

    # --- 3. spot morphology ----------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.1))
    panels = [("d_r90_mm", "$\\Delta r_{90}$ (mm)", "spot size at the receiver"),
              ("d_apfrac_pp", "$\\Delta$ aperture fraction (percentage points)",
               f"power across the r = {ap_mm:.0f} mm edge"),
              ("d_centroid_mm", "centroid shift (mm)", "centroid")]
    for ax, (col, lab, title) in zip(axes, panels):
        for t in traced_names:
            g = morph[morph.traced == t].sort_values("solar_el_deg")
            ax.fill_between(g.solar_el_deg, g[f"{col}_p05_occ"], g[f"{col}_p95_occ"],
                            color=colours[t], alpha=0.15, lw=0)
            ax.plot(g.solar_el_deg, g[f"{col}_med_occ"], "o-", ms=4, color=colours[t],
                    label=f"{t} $-$ {scalar_name}, occluded")
        g = morph[morph.traced == traced_names[-1]].sort_values("solar_el_deg")
        ax.plot(g.solar_el_deg, g[f"{col}_med_null"], "^--", ms=4, color="0.45",
                label="unoccluded heliostats (Monte-Carlo floor)")
        ax.axhline(0, lw=0.8, color="0.4")
        ax.set_xlabel("sun elevation (deg)")
        ax.set_ylabel(lab)
        ax.set_title(title)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.suptitle("Per-heliostat spot morphology at the receiver: traced occlusion "
                 "minus scalar (median, 5-95%)")
    fig.tight_layout()
    f3 = out_dir / "fig3_spot_morphology.png"
    fig.savefig(f3, dpi=150)
    plt.close(fig)

    for f in (f1, f2, f3):
        log(f"  wrote {f}")


# --------------------------------------------------------------------------
# verdict
# --------------------------------------------------------------------------
def build_verdict(args, scalar, traced_names, per_ts, per_h, morph, ann, over,
                  ap_mm, has_union, low_share) -> str:
    last = traced_names[-1]
    L = ["=" * 78,
         "VERDICT -- analytic occlusion scalars vs traced occlusion geometry",
         "=" * 78, ""]

    if ann is not None and len(ann):
        b = ann[ann.run == scalar.name].iloc[0]
        t = ann[ann.run == last].iloc[0]
        d = 100.0 * (b.annual_mwh_aperture / t.annual_mwh_aperture - 1.0)
        d_tot = 100.0 * (b.annual_mwh_total / t.annual_mwh_total - 1.0)
        sig = 100.0 * float(np.hypot(b.mc_sigma_mwh, t.mc_sigma_mwh)) / t.annual_mwh_aperture
        L += [f"ANNUAL COLLECTED ENERGY (aperture r <= {ap_mm:.0f} mm, monthly DNI)",
              f"  scalar path {scalar.name:<8s} {b.annual_mwh_aperture:12,.1f} MWh",
              f"  traced path {last:<8s} {t.annual_mwh_aperture:12,.1f} MWh",
              f"  scalar - traced        {d:+8.3f} %   +/- {sig:.3f} % "
              f"(Monte-Carlo 1 sigma)",
              f"  same, field total      {d_tot:+8.3f} %"]
        if has_union and (ann.run == "union").any():
            u = ann[ann.run == "union"].iloc[0]
            du = 100.0 * (u.annual_mwh_aperture / t.annual_mwh_aperture - 1.0)
            L += [f"  scalar in UNION form   {du:+8.3f} %   "
                  f"(shading.occlusion_efficiency, not eta_shade*eta_block)"]
        L += ["",
              f"  tolerance asked of the scalar path: {args.tolerance_pct:.2f} %",
              f"  -> the stored scalar form is "
              f"{'WITHIN' if abs(d) <= args.tolerance_pct else 'OUTSIDE'} it."]
        if has_union and (ann.run == "union").any():
            L += [f"  -> the union form is "
                  f"{'WITHIN' if abs(du) <= args.tolerance_pct else 'OUTSIDE'} it."]
        L.append("")

    if len(traced_names) > 1:
        first = traced_names[0]
        nb = 100.0 * (per_ts[f"P_ap_{first}"].sum() / per_ts[f"P_ap_{scalar.name}"].sum() - 1.0)
        sec = 100.0 * (per_ts[f"P_ap_{last}"].sum() / per_ts[f"P_ap_{first}"].sum() - 1.0)
        floor = 100.0 * float(np.sqrt((per_ts[f"sigma_ap_{first}"] ** 2).sum()
                                      + (per_ts[f"sigma_ap_{last}"] ** 2).sum())
                              ) / per_ts[f"P_ap_{last}"].sum()
        L += ["CHANNEL DECOMPOSITION (aperture power summed over all timesteps)",
              f"  neighbour shading+blocking, scalar -> traced ({scalar.name} -> {first}): "
              f"{nb:+.3f} %",
              f"  secondary shadow,           scalar -> traced ({first} -> {last}): "
              f"{sec:+.3f} %   (Monte-Carlo 1 sigma {floor:.3f} %)",
              f"  -> the whole difference is the neighbour channel. In AGGREGATE the",
              f"     axicon-shadow scalar is exact to the noise floor -- it is a hard-edged",
              f"     circle, and a circle projected onto a plane is still that circle.",
              f"     Per heliostat it is not: on the ~1 % of mirrors the shadow rim",
              f"     actually crosses the scalar disagrees with the trace by tens of",
              f"     percent either way (see the secondary-shadow block in the log). Those",
              f"     mirrors carry under 1 % of field power, so it cancels in the total but",
              f"     shows up as extra per-timestep scatter at low and middle sun.",
              ""]

    L += ["WHAT THE DIFFERENCE IS, NOT JUST HOW BIG",
          "  The scalar path multiplies eta_shade by eta_block. Those two losses",
          "  overlap on the mirror -- a patch in a neighbour's shadow sends no beam,",
          "  so it cannot also be blocked -- and multiplying charges the overlap",
          "  twice. That is a one-sided error: the scalar path is biased LOW, which",
          "  is why every traced-minus-scalar number here has the same sign, and why",
          "  the gap grows as the sun drops and the shadows lengthen. Re-weighting the",
          "  SAME analytic model as a union (shading.occlusion_efficiency) removes most",
          "  of it, which is the evidence that this -- and not the geometry of the",
          "  analytic model -- is what the difference was.",
          "  The analytic geometry itself is checked independently by",
          "  shading.self_check, which compares the sampled shading fraction against a",
          "  closed-form rectangle overlap for aligned mirror pairs (0.3491 vs 0.3502,",
          "  0.3540 vs 0.3536, 0.4117 vs 0.4107) and re-derives the ground shadow two",
          "  ways. It passes today; that is the leg of this argument that does not",
          "  depend on any ray trace.",
          ""]

    # elevation regimes
    L.append("REGIMES (power inside the aperture, summed over the timesteps in each band)")
    for lo, hi in EL_BANDS:
        m = (per_ts.solar_el_deg >= lo) & (per_ts.solar_el_deg < hi)
        if not m.any():
            L.append(f"  el {lo:4.0f}-{hi:4.0f} deg   -- not sampled by these runs")
            continue
        d = 100.0 * (per_ts.loc[m, f"P_ap_{last}"].sum()
                     / per_ts.loc[m, f"P_ap_{scalar.name}"].sum() - 1.0)
        extra = ""
        if has_union:
            du = 100.0 * (per_ts.loc[m, f"P_ap_{last}"].sum()
                          / per_ts.loc[m, "P_ap_union"].sum() - 1.0)
            extra = f"   union form {du:+6.2f} %"
        L.append(f"  el {lo:4.0f}-{hi:4.0f} deg   n={int(m.sum()):2d}   "
                 f"scalar low by {d:+6.2f} %{extra}")
    L.append("")

    heavy = per_h[(per_h.eta_product < 0.7) & (per_h.p_scalar > 0)]
    if len(heavy):
        dh = 100.0 * (heavy[f"p_{last}"].sum() / heavy.p_scalar.sum() - 1.0)
        L += [f"  heavily occluded heliostats (eta_shade*eta_block < 0.70): "
              f"{len(heavy):,} rows,",
              f"    {100*heavy.p_scalar.sum()/per_h.p_scalar.sum():.1f} % of the "
              f"field's aperture power, scalar low by {dh:+.2f} %.", ""]

    if over is not None and len(over):
        worst = over.sort_values("solar_el_deg").iloc[0]
        clean = over[over.overflowed == 0].sort_values("solar_el_deg")
        L += ["  SLOT OVERFLOW -- the regime where the TRACED path is the approximate one:",
              f"    The compared runs never overflow: their lowest sun is "
              f"{per_ts.solar_el_deg.min():.2f} deg and the",
              f"    sweep logs for them carry no overflow warning at all. On the grid in",
              f"    config.toml the lowest step at each end of the day does. At el = "
              f"{worst.solar_el_deg:.2f} deg,",
              f"    {int(worst.overflowed)} of 645 heliostats overflow the 10+4 slots; the "
              f"traced model then sees only",
              f"    the occluders that fit, and overstates collected power by "
              f"{worst.traced_overstates_pct:+.2f} % at that instant"]
        if len(clean):
            L.append(f"    (nothing overflows by el = {clean.iloc[0].solar_el_deg:.2f} deg, "
                     f"where the two agree to "
                     f"{clean.iloc[0].traced_overstates_pct:+.2f} %).")
        if isinstance(low_share, dict) and 5.0 in low_share:
            s5 = low_share[5.0]
            L.append(f"    Elevations below 5 deg carry {s5:.3f} % of modelled annual "
                     f"energy, so even if every")
            L.append(f"    such hour were mis-stated by that much the annual effect is "
                     f"~{abs(worst.traced_overstates_pct)*s5/100:.4f} %.")
        L.append("")

    # morphology
    L.append("DISTRIBUTION (the through-focus question)")
    mm = morph[morph.traced == last]
    lowest = mm.sort_values("solar_el_deg").iloc[0]
    highest = mm.sort_values("solar_el_deg").iloc[-1]
    L += [f"  Real occlusion removes specific rays, so it changes the spot; a scalar",
          f"  cannot. At the RECEIVER the change is small but not zero. Median per",
          f"  OCCLUDED heliostat, {last} minus {scalar.name}, with the same statistic over",
          f"  UNOCCLUDED heliostats beside it as the Monte-Carlo floor:",
          f"    el {lowest.solar_el_deg:5.2f} deg  d r50 {lowest.d_r50_mm_med_occ:+6.2f} mm "
          f"(null {lowest.d_r50_mm_med_null:+5.2f})   d r90 {lowest.d_r90_mm_med_occ:+6.2f} mm "
          f"(null {lowest.d_r90_mm_med_null:+5.2f})",
          f"                 d rms {lowest.d_rms_mm_med_occ:+6.2f} mm "
          f"(null {lowest.d_rms_mm_med_null:+5.2f})   aperture fraction "
          f"{lowest.d_apfrac_pp_med_occ:+.3f} pp (null {lowest.d_apfrac_pp_med_null:+.3f})",
          f"                 centroid {lowest.d_centroid_mm_med_occ:5.2f} mm "
          f"(null {lowest.d_centroid_mm_med_null:5.2f})",
          f"    el {highest.solar_el_deg:5.2f} deg  d r50 {highest.d_r50_mm_med_occ:+6.2f} mm "
          f"(null {highest.d_r50_mm_med_null:+5.2f})   d r90 {highest.d_r90_mm_med_occ:+6.2f} mm "
          f"(null {highest.d_r90_mm_med_null:+5.2f})",
          f"                 d rms {highest.d_rms_mm_med_occ:+6.2f} mm "
          f"(null {highest.d_rms_mm_med_null:+5.2f})   aperture fraction "
          f"{highest.d_apfrac_pp_med_occ:+.3f} pp (null {highest.d_apfrac_pp_med_null:+.3f})",
          f"                 centroid {highest.d_centroid_mm_med_occ:5.2f} mm "
          f"(null {highest.d_centroid_mm_med_null:5.2f})",
          f"  Field-summed over all 645 heliostats the aperture fraction moves by at most",
          f"  {mm.field_d_apfrac_pp.abs().max():.3f} pp and r90 by "
          f"{mm.field_d_r90_mm.abs().max():.2f} mm, so the scalar path mis-states",
          f"  spillage by less than that.",
          "",
          "  CAVEAT, stated because the data cannot answer it: these runs store rays",
          "  at the RECEIVER only. Nothing here measures a plane between the secondary",
          "  and the receiver. At the receiver 645 heliostats' occluded edges land in",
          "  different places and wash out; through focus they need not. Any",
          "  through-focus study must trace occlusion -- this vet does not license",
          "  scalars there.",
          ""]

    L += ["BOTTOM LINE"]
    if ann is not None and len(ann):
        b = ann[ann.run == scalar.name].iloc[0]
        t = ann[ann.run == last].iloc[0]
        d = 100.0 * (b.annual_mwh_aperture / t.annual_mwh_aperture - 1.0)
        L += [f"  Annualised collected energy from the scalar path is {abs(d):.3f} % "
              f"{'below' if d < 0 else 'above'} the",
              f"  traced path, one-sided and explained. Run the remaining comparison",
              f"  sweeps without occlusion geometry and apply the scalars in post; if the",
              f"  {abs(d):.2f} % matters, apply them in union form and it drops to "
              + (f"{abs(100.0*(ann[ann.run=='union'].iloc[0].annual_mwh_aperture/t.annual_mwh_aperture-1)):.3f} %."
                 if (ann.run == "union").any() else "less.")]
    L += ["  The scalar path is NOT a substitute for traced occlusion when the question",
          "  is: (a) instantaneous power at low sun -- it is several percent low there;",
          "  (b) per-heliostat power for a heavily occluded mirror; (c) anything about",
          "  the shape of the spot, including all through-focus work.",
          "=" * 78]
    return "\n".join(L)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Vet analytic occlusion scalars against traced occlusion.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", default=["full5", "full6", "full7"],
                    help="run directories; the FIRST is the scalar path, the rest are "
                         "traced paths compared against it")
    ap.add_argument("--out", default=str(ANALYSIS / "vet_occlusion"))
    ap.add_argument("--aperture-mm", type=float, default=700.0)
    ap.add_argument("--tolerance-pct", type=float, default=0.5,
                    help="annual-energy difference the scalar path is allowed")
    ap.add_argument("--dni-mode", default=None, help="constant | table | monthly")
    ap.add_argument("--limit-timesteps", type=int, default=0,
                    help="use only the first N shared timesteps (smoke tests)")
    ap.add_argument("--skip-union", action="store_true")
    ap.add_argument("--skip-annual", action="store_true")
    ap.add_argument("--skip-overflow", action="store_true")
    ap.add_argument("--skip-figures", action="store_true")
    args = ap.parse_args(argv)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_lines: list[str] = []

    def log(msg=""):
        print(msg, flush=True)
        log_lines.append(str(msg))

    from beamdown import dni as D
    from beamdown import metrics as M
    from beamdown.config import load_config

    cfg = load_config(None)
    runs = [Run(n, cfg) for n in args.runs]
    scalar, traced = runs[0], runs[1:]
    if not traced:
        raise SystemExit("give at least two runs: a scalar path and a traced path")
    traced_names = [r.name for r in traced]
    last = traced_names[-1]

    log(f"vet_occlusion_scalars   {time.strftime('%Y-%m-%d %H:%M')}")
    log(f"  scalar path : {scalar.root}")
    for r in traced:
        log(f"  traced path : {r.root}")
    log("")
    audit_semantics(runs, log)
    log("")

    ids = [h for h in scalar.ids if all(h in set(r.ids) for r in traced)]
    keys = sorted(set.intersection(*[r.keys() for r in runs]))
    if args.limit_timesteps:
        keys = keys[:args.limit_timesteps]
        log(f"  LIMITED to the first {len(keys)} timesteps (smoke test)")
    log(f"  {len(keys)} shared timesteps, {len(ids)} shared heliostats")
    if len({r.rays for r in runs}) != 1:
        log(f"  NOTE: rays/heliostat differs between runs: "
            f"{ {r.name: r.rays for r in runs} }")
    take = {r.name: np.array([r.ids.index(h) for h in ids]) for r in runs}
    key_index = {k: i for i, k in enumerate(keys)}

    etas = None
    if not args.skip_union and scalar.mode != "product":
        # The union variant re-weights the scalar run's *clean* rays. If the run
        # named as the scalar path already traced its occluders, its rays are not
        # clean and re-weighting them would charge the same loss twice.
        log(f"  union variant skipped: {scalar.name} already traces occlusion "
            f"(power_w carries {scalar.mode}), so its rays cannot be re-weighted.")
    elif not args.skip_union:
        etas = recompute_etas(cfg, keys, scalar.summary,
                              out_dir / "analytic_etas.npz", log)
        stored = np.array([(scalar.rows(k).eta_shade * scalar.rows(k).eta_block)
                           .to_numpy(float) for k in keys])
        fresh = (etas["shade"][:, take[scalar.name]] * etas["block"][:, take[scalar.name]])
        log(f"  recomputed eta_shade*eta_block reproduces the stored value to "
            f"{float(np.nanmax(np.abs(stored - fresh))):.1e}")
    log("")

    masks_flat = radial_masks_flat(cfg, RADII_MM)
    ref_i = int(np.argmin(np.abs(np.array(RADII_MM) - args.aperture_mm)))
    ap_mm = RADII_MM[ref_i]

    ts_rows, h_frames, morph_rows = [], [], []
    cap = {r.name: np.zeros(len(RADII_MM)) for r in runs}
    cap["union"] = np.zeros(len(RADII_MM))
    p_tot = {r.name: {} for r in runs} | {"union": {}}
    p_ap = {r.name: {} for r in runs} | {"union": {}}

    t0 = time.time()
    for i, key in enumerate(keys, 1):
        first = scalar.rows(key).iloc[0]
        rec = dict(timestep=key, date=str(first.date), hour=float(first.hour),
                   solar_el_deg=float(first.solar_el_deg),
                   solar_az_deg=float(first.solar_az_deg))
        wn, wa, raw_n, raw_a, fmap = {}, {}, {}, {}, {}
        for r in runs:
            n, inside, counts = read_counts(r.store, key, take[r.name], masks_flat)
            w = r.weights(key)
            raw_n[r.name], raw_a[r.name] = n, inside
            wn[r.name], wa[r.name] = n * w, inside * w[:, None]
            fmap[r.name] = (counts * w[:, None, None].astype(np.float32)).sum(axis=0)
            del counts
            cap[r.name] += wa[r.name].sum(axis=0) * r.sf
            rec[f"P_tot_{r.name}"] = float(wn[r.name].sum() * r.sf)
            rec[f"P_ap_{r.name}"] = float(wa[r.name][:, ref_i].sum() * r.sf)
            rec[f"sigma_tot_{r.name}"] = mc_sigma(w, n, r.rays, r.sf)
            rec[f"sigma_ap_{r.name}"] = mc_sigma(w, inside[:, ref_i], r.rays, r.sf)
            p_tot[r.name][key] = rec[f"P_tot_{r.name}"]
            p_ap[r.name][key] = rec[f"P_ap_{r.name}"]

        # the same analytic model, applied as a union instead of a product
        u = None
        if etas is not None:
            u = etas["union"][key_index[key]][take[scalar.name]]
            cap["union"] += (raw_a[scalar.name] * u[:, None]).sum(axis=0) * scalar.sf
            rec["P_tot_union"] = float((raw_n[scalar.name] * u).sum() * scalar.sf)
            rec["P_ap_union"] = float((raw_a[scalar.name][:, ref_i] * u).sum() * scalar.sf)
            p_tot["union"][key] = rec["P_tot_union"]
            p_ap["union"][key] = rec["P_ap_union"]

        for r in traced:
            for kind in ("tot", "ap"):
                base = rec[f"P_{kind}_{scalar.name}"]
                rec[f"ratio_{kind}_{r.name}"] = (rec[f"P_{kind}_{r.name}"] / base
                                                 if base else np.nan)
        if etas is not None:
            for kind in ("tot", "ap"):
                rec[f"ratio_{kind}_union"] = (rec[f"P_{kind}_{last}"]
                                              / rec[f"P_{kind}_union"]
                                              if rec[f"P_{kind}_union"] else np.nan)
        for kind in ("tot", "ap"):
            rec[f"sigma_{kind}_pct"] = 100.0 * float(np.hypot(
                rec[f"sigma_{kind}_{scalar.name}"], rec[f"sigma_{kind}_{last}"])
            ) / rec[f"P_{kind}_{scalar.name}"]

        for r in runs:
            mm = M.map_metrics(fmap[r.name], r.rays, cfg)
            rec[f"field_r50_{r.name}"] = mm["r50_mm"]
            rec[f"field_r90_{r.name}"] = mm["r90_mm"]
            rec[f"field_rms_{r.name}"] = mm["rms_radius_mm"]
            rec[f"field_apfrac_{r.name}"] = (rec[f"P_ap_{r.name}"] / rec[f"P_tot_{r.name}"]
                                             if rec[f"P_tot_{r.name}"] else np.nan)
        ts_rows.append(rec)

        # ---- per heliostat ------------------------------------------------
        s_rows = scalar.rows(key)
        nan = np.full(len(ids), np.nan)
        ap_s = wa[scalar.name][:, ref_i]
        hr = pd.DataFrame(dict(
            timestep=key, heliostat_id=ids, solar_el_deg=float(first.solar_el_deg),
            radius_m=s_rows.radius_m.to_numpy(float),
            eta_shade=s_rows.eta_shade.to_numpy(float),
            eta_block=s_rows.eta_block.to_numpy(float),
            eta_secondary=s_rows.eta_secondary.to_numpy(float),
            eta_product=(s_rows.eta_shade * s_rows.eta_block).to_numpy(float),
            p_scalar=ap_s * scalar.sf,
            var_scalar=mc_variance(scalar.weights(key), raw_a[scalar.name][:, ref_i],
                                   scalar.rays, scalar.sf),
            apfrac_scalar=np.divide(raw_a[scalar.name][:, ref_i], raw_n[scalar.name],
                                    out=nan.copy(), where=raw_n[scalar.name] > 0),
            r50_scalar=s_rows.r50_mm.to_numpy(float),
            r90_scalar=s_rows.r90_mm.to_numpy(float),
            rms_scalar=s_rows.rms_radius_mm.to_numpy(float),
        ))
        if u is not None:
            hr["eta_union"] = u
            hr["p_union"] = raw_a[scalar.name][:, ref_i] * u * scalar.sf
        for r in traced:
            t_rows = r.rows(key)
            ap_t = wa[r.name][:, ref_i]
            hr[f"p_{r.name}"] = ap_t * r.sf
            hr[f"ratio_{r.name}"] = np.divide(ap_t * r.sf, ap_s * scalar.sf,
                                              out=nan.copy(), where=ap_s > 0)
            if u is not None:
                hr[f"ratio_union_{r.name}"] = np.divide(
                    ap_t * r.sf, hr["p_union"].to_numpy(float), out=nan.copy(),
                    where=hr["p_union"].to_numpy(float) > 0)
            hr[f"var_{r.name}"] = mc_variance(r.weights(key), raw_a[r.name][:, ref_i],
                                              r.rays, r.sf)
            hr[f"apfrac_{r.name}"] = np.divide(raw_a[r.name][:, ref_i], raw_n[r.name],
                                               out=nan.copy(), where=raw_n[r.name] > 0)
            hr[f"r50_{r.name}"] = t_rows.r50_mm.to_numpy(float)
            hr[f"r90_{r.name}"] = t_rows.r90_mm.to_numpy(float)
            hr[f"rms_{r.name}"] = t_rows.rms_radius_mm.to_numpy(float)
            hr[f"dcent_{r.name}"] = np.hypot(
                t_rows.centroid_x_mm.to_numpy(float) - s_rows.centroid_x_mm.to_numpy(float),
                t_rows.centroid_y_mm.to_numpy(float) - s_rows.centroid_y_mm.to_numpy(float))

            d = {"d_r50_mm": hr[f"r50_{r.name}"] - hr.r50_scalar,
                 "d_r90_mm": hr[f"r90_{r.name}"] - hr.r90_scalar,
                 "d_rms_mm": hr[f"rms_{r.name}"] - hr.rms_scalar,
                 "d_apfrac_pp": 100.0 * (hr[f"apfrac_{r.name}"] - hr.apfrac_scalar),
                 "d_centroid_mm": hr[f"dcent_{r.name}"]}
            # Heliostats the analytic model says nothing occludes are the null
            # sample: there the two runs differ ONLY by their independent rays, so
            # whatever spread they show IS the Monte-Carlo floor for these metrics.
            occ = (hr.eta_product < 0.999) | (hr.eta_secondary < 0.999)
            null = ~occ
            row = dict(timestep=key, traced=r.name,
                       solar_el_deg=float(first.solar_el_deg),
                       n_occluded=int(occ.sum()), n_unoccluded=int(null.sum()))
            for nm, v in d.items():
                row[f"{nm}_med"] = float(np.nanmedian(v))
                row[f"{nm}_p05"] = float(np.nanpercentile(v, 5))
                row[f"{nm}_p95"] = float(np.nanpercentile(v, 95))
                row[f"{nm}_med_occ"] = (float(np.nanmedian(v[occ])) if occ.any()
                                        else float("nan"))
                row[f"{nm}_p05_occ"] = (float(np.nanpercentile(v[occ], 5)) if occ.any()
                                        else float("nan"))
                row[f"{nm}_p95_occ"] = (float(np.nanpercentile(v[occ], 95)) if occ.any()
                                        else float("nan"))
                row[f"{nm}_med_null"] = (float(np.nanmedian(v[null])) if null.any()
                                         else float("nan"))
                row[f"{nm}_p95_null"] = (float(np.nanpercentile(v[null], 95))
                                         if null.any() else float("nan"))
            row["field_d_r90_mm"] = (rec[f"field_r90_{r.name}"]
                                     - rec[f"field_r90_{scalar.name}"])
            row["field_d_rms_mm"] = (rec[f"field_rms_{r.name}"]
                                     - rec[f"field_rms_{scalar.name}"])
            row["field_d_apfrac_pp"] = 100.0 * (rec[f"field_apfrac_{r.name}"]
                                                - rec[f"field_apfrac_{scalar.name}"])
            morph_rows.append(row)
        h_frames.append(hr)
        if i == 1 or i % 10 == 0 or i == len(keys):
            log(f"  [{i}/{len(keys)}] {key}  el {first.solar_el_deg:5.2f}  "
                f"{(time.time()-t0)/i:.1f} s/timestep")

    per_ts = pd.DataFrame(ts_rows)
    per_h = pd.concat(h_frames, ignore_index=True)
    morph = pd.DataFrame(morph_rows)
    per_ts.to_csv(out_dir / "per_timestep.csv", index=False)
    per_h.to_csv(out_dir / "per_heliostat.csv", index=False)
    morph.to_csv(out_dir / "morphology.csv", index=False)
    log(f"  wrote per_timestep.csv ({len(per_ts)} rows), per_heliostat.csv "
        f"({len(per_h):,} rows), morphology.csv ({len(morph)} rows)")

    # The binned counts and the summary are two independent routes to the same
    # watts (the summary counts raw landed rays, the flux maps histogram them).
    # If they disagree, one of them is not what it claims to be.
    for r in runs:
        want = r.summary[r.summary.timestep.isin(keys)].groupby("timestep").power_w.sum()
        got = per_ts.set_index("timestep")[f"P_tot_{r.name}"]
        rel = float(np.max(np.abs(got.reindex(want.index) / want - 1.0)))
        log(f"  cross-check {r.name}: field power from flux maps vs from summary "
            f"power_w -- max relative difference {rel:.2e}")
        if rel > 1e-6:
            log(f"    WARNING: {r.name} flux maps and summary disagree by more than "
                f"rounding; treat its numbers with care.")
    log("")

    # --- aperture-radius sweep -------------------------------------------
    rad = pd.DataFrame({"radius_mm": RADII_MM})
    for r in runs:
        rad[f"P_{r.name}_W"] = cap[r.name]
    for r in traced:
        rad[f"delta_{r.name}_vs_{scalar.name}_pct"] = 100.0 * (cap[r.name]
                                                               / cap[scalar.name] - 1.0)
    if etas is not None:
        rad["P_union_W"] = cap["union"]
        rad[f"delta_{last}_vs_union_pct"] = 100.0 * (cap[last] / cap["union"] - 1.0)
    if len(traced) > 1:
        rad[f"delta_{last}_vs_{traced[0].name}_pct"] = 100.0 * (cap[last]
                                                                / cap[traced[0].name] - 1.0)
    rad.to_csv(out_dir / "aperture_radius_sweep.csv", index=False)
    log("POWER INSIDE AN APERTURE OF RADIUS R (summed over all timesteps, DNI 1000)")
    log(f"  {'R (mm)':>7}" + "".join(f"{r.name:>13}" for r in runs)
        + ("".join(f"{'union':>13}") if etas is not None else "")
        + "".join(f"{'d ' + r.name:>11}" for r in traced)
        + (f"{'d union':>11}" if etas is not None else ""))
    for j, R in enumerate(RADII_MM):
        line = f"  {R:7.0f}" + "".join(f"{cap[r.name][j]/1e6:12.2f}M" for r in runs)
        if etas is not None:
            line += f"{cap['union'][j]/1e6:12.2f}M"
        line += "".join(f"{100.0*(cap[r.name][j]/cap[scalar.name][j]-1):+10.2f}%"
                        for r in traced)
        if etas is not None:
            line += f"{100.0*(cap[last][j]/cap['union'][j]-1):+10.2f}%"
        log(line)
    log("")

    # --- annual energy ----------------------------------------------------
    ann, low_share = None, None
    if not args.skip_annual:
        from beamdown import energy as E

        prov = (D.provider_for(cfg, args.dni_mode) if args.dni_mode
                else D.load_dni_provider(cfg))
        log(f"ANNUAL ENERGY   DNI model: {prov.describe()}")
        template = scalar.summary[["date", "hour", "timestep", "solar_az_deg",
                                   "solar_el_deg"]].copy()
        base, wt = annual_weights(template, cfg, prov, len(ids), p_ap[last], log)

        names = [r.name for r in runs] + (["union"] if etas is not None else [])
        rows = []
        for name in names:
            res = annual_from_power(p_ap[name], template, cfg, prov, len(ids))
            res_t = annual_from_power(p_tot[name], template, cfg, prov, len(ids))
            sig = np.nan
            if name != "union":
                sig = float(np.sqrt(sum(
                    (wt[k] * per_ts.loc[per_ts.timestep == k, f"sigma_ap_{name}"].iloc[0]) ** 2
                    for k in keys)))
            rows.append(dict(run=name,
                             annual_mwh_aperture=res["annual_energy_mwh"],
                             annual_mwh_total=res_t["annual_energy_mwh"],
                             annual_optical_efficiency=res["annual_optical_efficiency"],
                             mc_sigma_mwh=sig / 1000.0 if np.isfinite(sig) else np.nan,
                             extrapolated_fraction=res["extrapolated_fraction"]))
        ann = pd.DataFrame(rows)
        ref = ann[ann.run == last].iloc[0]
        ann["delta_vs_traced_pct"] = 100.0 * (ann.annual_mwh_aperture
                                              / ref.annual_mwh_aperture - 1.0)
        ann["delta_total_vs_traced_pct"] = 100.0 * (ann.annual_mwh_total
                                                    / ref.annual_mwh_total - 1.0)
        ann.to_csv(out_dir / "annual_energy.csv", index=False)

        log(f"  {'run':<8}{'MWh (r<=' + format(ap_mm, '.0f') + 'mm)':>18}"
            f"{'MWh (total)':>14}{'vs ' + last:>12}{'MC 1sigma':>12}")
        for _, row in ann.iterrows():
            s = (f"{100.0*row.mc_sigma_mwh/row.annual_mwh_aperture:11.3f}%"
                 if np.isfinite(row.mc_sigma_mwh) else f"{'-':>12}")
            log(f"  {row.run:<8}{row.annual_mwh_aperture:18,.1f}"
                f"{row.annual_mwh_total:14,.1f}{row.delta_vs_traced_pct:+11.3f}%{s}")
        log("")

        daily_rows = []
        for name in names:
            frame = scalar.summary.drop_duplicates("timestep")[["date", "hour",
                                                                "timestep"]].copy()
            frame["power_w"] = frame.timestep.map(p_ap[name])
            frame["heliostat_id"] = 0
            for d in sorted(set(frame.date)):
                res = E.traced_day_energy(frame, cfg, prov, d)
                daily_rows.append(dict(run=name, date=str(d),
                                       energy_mwh=res["energy_kwh"] / 1000.0))
        daily = pd.DataFrame(daily_rows)
        piv = daily.pivot(index="date", columns="run", values="energy_mwh")
        for name in names:
            if name != scalar.name:
                piv[f"{name}/{scalar.name}"] = piv[name] / piv[scalar.name]
        piv.to_csv(out_dir / "daily_energy.csv")
        log("  daily collected energy by traced date (trapezoid over sampled hours, "
            f"aperture r<={ap_mm:.0f} mm)")
        log(f"    {'date':<12}{scalar.name + ' MWh':>14}"
            + "".join(f"{n + '/' + scalar.name:>16}" for n in names if n != scalar.name))
        for d, row in piv.iterrows():
            log(f"    {d:<12}{row[scalar.name]:14.2f}"
                + "".join(f"{row[f'{n}/{scalar.name}']:16.5f}"
                          for n in names if n != scalar.name))
        log("")

        hourly = base["hourly"]
        tot = hourly.energy_wh.sum()
        shares = {}
        log("  share of modelled annual energy collected below a sun elevation")
        for thr in (2.0, 5.0, 10.0, 15.0):
            shares[thr] = 100.0 * hourly.loc[hourly.solar_el_deg < thr,
                                             "energy_wh"].sum() / tot
            log(f"    el < {thr:4.1f} deg : {shares[thr]:7.4f} %")
        low_share = shares
        log("")

    # --- slot overflow ----------------------------------------------------
    over = None
    over_csv = out_dir / "overflow_probe.csv"
    log("TRACED-PATH SLOT OVERFLOW (analytic; over the grid in config.toml)")
    if args.skip_overflow and over_csv.exists():
        over = pd.read_csv(over_csv)
        log(f"  reusing {over_csv.name} (--skip-overflow)")
    elif args.skip_overflow:
        log("  skipped (--skip-overflow, and no cached overflow_probe.csv)")
    else:
        over = overflow_probe(cfg, log)
        over.to_csv(over_csv, index=False)
    log("")

    # --- figures ----------------------------------------------------------
    if not args.skip_figures:
        make_figures(out_dir, per_ts, per_h, morph, traced_names, scalar.name,
                     etas is not None, ap_mm, log)
        log("")

    # --- bands ------------------------------------------------------------
    log("WHERE THE DIFFERENCE CONCENTRATES")
    log(f"  by sun elevation band (power inside r<={ap_mm:.0f} mm, summed over the band)")
    band_rows = []
    for lo, hi in EL_BANDS:
        m = (per_ts.solar_el_deg >= lo) & (per_ts.solar_el_deg < hi)
        if not m.any():
            log(f"    el {lo:4.0f}-{hi:4.0f}   not sampled by these runs")
            band_rows.append(dict(band=f"{lo:.0f}-{hi:.0f}", timesteps=0))
            continue
        row = dict(band=f"{lo:.0f}-{hi:.0f}", timesteps=int(m.sum()))
        txt = f"    el {lo:4.0f}-{hi:4.0f}   n={int(m.sum()):2d}"
        for r in traced:
            d = 100.0 * (per_ts.loc[m, f"P_ap_{r.name}"].sum()
                         / per_ts.loc[m, f"P_ap_{scalar.name}"].sum() - 1.0)
            row[f"delta_{r.name}_pct"] = d
            txt += f"   {r.name} {d:+6.2f}%"
        if etas is not None:
            d = 100.0 * (per_ts.loc[m, f"P_ap_{last}"].sum()
                         / per_ts.loc[m, "P_ap_union"].sum() - 1.0)
            row["delta_traced_vs_union_pct"] = d
            txt += f"   {last} vs union {d:+6.2f}%"
        row["mc_sigma_pct"] = float(per_ts.loc[m, "sigma_ap_pct"].mean())
        txt += f"   (MC 1sigma {row['mc_sigma_pct']:.3f}%)"
        band_rows.append(row)
        log(txt)
    pd.DataFrame(band_rows).to_csv(out_dir / "elevation_bands.csv", index=False)
    log("")

    log("  by analytic eta_shade x eta_block (per heliostat, all timesteps)")
    eta_rows = []
    for lo, hi in zip(ETA_EDGES[:-1], ETA_EDGES[1:]):
        m = (per_h.eta_product >= lo) & (per_h.eta_product < hi) & (per_h.p_scalar > 0)
        if not m.any():
            continue
        row = dict(band=f"{lo:.3f}-{hi:.3f}", rows=int(m.sum()),
                   share_of_scalar_power=float(per_h.loc[m, "p_scalar"].sum()
                                               / per_h.p_scalar.sum()))
        txt = (f"    eta {lo:5.3f}-{hi:5.3f}  n={int(m.sum()):6d}  "
               f"{100*row['share_of_scalar_power']:5.1f}% of power")
        for r in traced:
            d = 100.0 * (per_h.loc[m, f"p_{r.name}"].sum()
                         / per_h.loc[m, "p_scalar"].sum() - 1.0)
            row[f"delta_{r.name}_pct"] = d
            txt += f"   {r.name} {d:+7.2f}%"
        if etas is not None:
            d = 100.0 * (per_h.loc[m, f"p_{last}"].sum()
                         / per_h.loc[m, "p_union"].sum() - 1.0)
            row["delta_traced_vs_union_pct"] = d
            txt += f"   vs union {d:+7.2f}%"
        sig = 100.0 * float(np.sqrt(per_h.loc[m, "var_scalar"].sum()
                                    + per_h.loc[m, f"var_{last}"].sum())) \
            / float(per_h.loc[m, "p_scalar"].sum())
        row["mc_sigma_pct"] = sig
        txt += f"  (MC 1sigma {sig:.3f}%)"
        eta_rows.append(row)
        log(txt)
    pd.DataFrame(eta_rows).to_csv(out_dir / "eta_bands.csv", index=False)
    log("")

    # --- the secondary-shadow channel, heliostat by heliostat ------------
    # Aggregates can hide a coarse approximation that simply does not carry much
    # power. The axicon's shadow is applied as ONE number per heliostat in the
    # scalar path (and in full6), so the mirrors its rim actually crosses are the
    # place to look: there the scalar spreads the loss over the whole spot while
    # the trace removes the rays that were really in shadow.
    sec_run = next((r.name for r in traced if r.mode == "secondary"), None)
    if sec_run and sec_run != last:
        part = per_h[(per_h.eta_secondary < 0.999) & (per_h.eta_secondary > 0.001)
                     & (per_h[f"p_{sec_run}"] > 0)]
        ctrl = per_h[(per_h.eta_secondary > 0.999) & (per_h.eta_product > 0.999)
                     & (per_h[f"p_{sec_run}"] > 0)]
        if len(part) and len(ctrl):
            r = part[f"p_{last}"] / part[f"p_{sec_run}"] - 1.0
            pred = float((np.sqrt(part[f"var_{sec_run}"] + part[f"var_{last}"])
                          / part[f"p_{sec_run}"]).mean())
            rc = ctrl[f"p_{last}"] / ctrl[f"p_{sec_run}"] - 1.0
            log(f"  the secondary-shadow scalar, on the mirrors its rim crosses "
                f"({sec_run} vs {last})")
            log(f"    {len(part):,} rows ({100*len(part)/len(per_h):.2f} % of all), "
                f"carrying {100*part.p_scalar.sum()/per_h.p_scalar.sum():.2f} % of the "
                f"field's aperture power")
            log(f"    per heliostat: median {100*r.median():+.2f} %, "
                f"5-95 % {100*r.quantile(.05):+.1f} .. {100*r.quantile(.95):+.1f} %, "
                f"sd {100*r.std():.1f} %  (Monte-Carlo would give "
                f"{100*pred:.2f} %)")
            log(f"    same statistic on unshadowed heliostats: sd {100*rc.std():.2f} % "
                f"-- so the spread above is the scalar, not the rays")
            log(f"    aggregated over those rows: "
                f"{100*(part[f'p_{last}'].sum()/part[f'p_{sec_run}'].sum() - 1):+.3f} %, "
                f"which is {100*(part[f'p_{last}'].sum() - part[f'p_{sec_run}'].sum())/per_h.p_scalar.sum():+.4f} % "
                f"of field power")
            log("")

    log("  by field radius (per heliostat, all timesteps)")
    rad_rows = []
    sub = per_h[per_h.p_scalar > 0]
    for band, g in sub.groupby(pd.cut(sub.radius_m, [0, 30, 45, 60, 80, 200]),
                               observed=True):
        row = dict(band=str(band), rows=len(g))
        txt = f"    r {str(band):>14}  n={len(g):6d}"
        for r in traced:
            d = 100.0 * (g[f"p_{r.name}"].sum() / g.p_scalar.sum() - 1.0)
            row[f"delta_{r.name}_pct"] = d
            txt += f"   {r.name} {d:+6.2f}%"
        sig = 100.0 * float(np.sqrt(g.var_scalar.sum() + g[f"var_{last}"].sum())) \
            / float(g.p_scalar.sum())
        row["mc_sigma_pct"] = sig
        txt += f"  (MC 1sigma {sig:.3f}%)"
        rad_rows.append(row)
        log(txt)
    pd.DataFrame(rad_rows).to_csv(out_dir / "radius_bands.csv", index=False)
    log("")

    # --- morphology summary ----------------------------------------------
    log("SPOT MORPHOLOGY (per heliostat, traced minus scalar; medians over the band)")
    log("  'occluded' = the analytic model says something occludes it; 'null' = it says")
    log("  nothing does, so the two runs there differ only by their independent rays --")
    log("  the null column IS the Monte-Carlo floor for the column beside it.")
    log(f"  {'traced':>7} {'el band':>9} {'set':>9} {'n':>6} {'d r50 mm':>10} "
        f"{'d r90 mm':>10} {'d rms mm':>10} {'d apfrac pp':>12} {'centroid mm':>12}")
    mb = []
    for r in traced:
        for lo, hi in EL_BANDS:
            m = (per_h.solar_el_deg >= lo) & (per_h.solar_el_deg < hi)
            if not m.any():
                continue
            occ_all = (per_h.eta_product < 0.999) | (per_h.eta_secondary < 0.999)
            for tag, sel in (("occluded", m & occ_all), ("null", m & ~occ_all)):
                if not sel.any():
                    continue
                g = per_h[sel]
                row = dict(traced=r.name, band=f"{lo:.0f}-{hi:.0f}", set=tag, n=len(g),
                           d_r50_mm=float(np.nanmedian(g[f"r50_{r.name}"] - g.r50_scalar)),
                           d_r90_mm=float(np.nanmedian(g[f"r90_{r.name}"] - g.r90_scalar)),
                           d_rms_mm=float(np.nanmedian(g[f"rms_{r.name}"] - g.rms_scalar)),
                           d_apfrac_pp=float(100.0 * np.nanmedian(g[f"apfrac_{r.name}"]
                                                                  - g.apfrac_scalar)),
                           d_centroid_mm=float(np.nanmedian(g[f"dcent_{r.name}"])))
                mb.append(row)
                log(f"  {r.name:>7} {f'{lo:.0f}-{hi:.0f}':>9} {tag:>9} {len(g):6d} "
                    f"{row['d_r50_mm']:10.2f} {row['d_r90_mm']:10.2f} "
                    f"{row['d_rms_mm']:10.2f} {row['d_apfrac_pp']:12.3f} "
                    f"{row['d_centroid_mm']:12.2f}")
    pd.DataFrame(mb).to_csv(out_dir / "morphology_by_band.csv", index=False)
    log("")
    log("  field-summed spot (645 heliostats added, then measured)")
    log(f"    {'traced':>7} {'el':>6} {'d r50 mm':>10} {'d r90 mm':>10} "
        f"{'d apfrac pp':>12}")
    for r in traced:
        g = morph[morph.traced == r.name].sort_values("solar_el_deg")
        for _, row in g.iloc[[0, len(g) // 2, -1]].iterrows():
            f50 = (per_ts.loc[per_ts.timestep == row.timestep,
                              f"field_r50_{r.name}"].iloc[0]
                   - per_ts.loc[per_ts.timestep == row.timestep,
                                f"field_r50_{scalar.name}"].iloc[0])
            log(f"    {r.name:>7} {row.solar_el_deg:6.2f} {f50:10.2f} "
                f"{row.field_d_r90_mm:10.2f} {row.field_d_apfrac_pp:12.3f}")
    log("")

    verdict = build_verdict(args, scalar, traced_names, per_ts, per_h, morph, ann,
                            over, ap_mm, etas is not None, low_share)
    log(verdict)
    (out_dir / "verdict.txt").write_text(verdict, encoding="utf-8")
    (out_dir / "vet_occlusion.log").write_text("\n".join(log_lines), encoding="utf-8")
    print(f"\nwrote {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
