"""Recompute shading and blocking on a finished sweep, without re-tracing.

Why this can exist at all
-------------------------
Shading and blocking never entered the ray trace. They are scalar efficiencies
computed analytically per (heliostat, timestep) and folded into watts-per-ray at
write time by :func:`beamdown.metrics.spot_metrics`. The stored ray counts are
therefore innocent of them, and a correction to the geometry is a multiply over
the summary table rather than hours of Quadoa.

That is the whole reason ``shading.py`` sits outside the trace, and this module
is what cashes it in.

What it does not touch
----------------------
``raw/`` and ``flux/`` are counts and stay exactly as traced. Only the summary's
weight-bearing columns move, and the original is copied aside first.

Safety
------
Pointing is re-solved from the strategy rather than read back, because the aim
points that blocking needs were never stored. The re-solved ``rot_az``/``rot_el``
are then checked against what the sweep recorded, and a mismatch aborts: it would
mean the config or the strategy has changed since the run, in which case the
whole summary -- not just its weights -- is stale, and quietly rescaling it would
paper over that.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

# Columns the sweep multiplied by shading x blocking before writing.
WEIGHTED = ("power_w", "peak_flux_w_m2")

# Largest pointing disagreement tolerated between the stored run and a fresh
# solve, in degrees. Anything above this is a changed model, not rounding.
POINTING_TOL_DEG = 1e-6


@dataclass
class RescaleReport:
    run: str
    timesteps: int = 0
    rows: int = 0
    power_before_w: float = 0.0
    power_after_w: float = 0.0
    eta_shade_before: float = 0.0
    eta_shade_after: float = 0.0
    eta_block_before: float = 0.0
    eta_block_after: float = 0.0
    eta_secondary: float = 0.0
    worst_pointing_deg: float = 0.0
    per_step: list = field(default_factory=list)
    applied: bool = False
    backup: str | None = None

    def describe(self) -> str:
        d = 100.0 * (self.power_after_w / self.power_before_w - 1.0) if self.power_before_w else 0.0
        lines = [
            f"{self.run}: {self.rows:,} rows over {self.timesteps} timesteps",
            "",
            f"  {'':22s} {'before':>10s} {'after':>10s}",
            f"  {'eta_shade (mean)':22s} {self.eta_shade_before:10.4f} {self.eta_shade_after:10.4f}",
            f"  {'eta_block (mean)':22s} {self.eta_block_before:10.4f} {self.eta_block_after:10.4f}",
            f"  {'eta_secondary (mean)':22s} {'-':>10s} {self.eta_secondary:10.4f}",
            f"  {'delivered power (MW)':22s} {self.power_before_w/1e6:10.3f} "
            f"{self.power_after_w/1e6:10.3f}   {d:+.2f}%",
            "",
            f"  worst pointing disagreement vs the stored run: {self.worst_pointing_deg:.2e} deg",
        ]
        if self.per_step:
            lines += ["", "  per timestep",
                      f"    {'timestep':>16s} {'sun el':>7s} {'eta_sec':>8s} {'hit':>4s} "
                      f"{'delta':>8s}"]
            for s in self.per_step:
                lines.append(f"    {s['timestep']:>16s} {s['el']:7.1f} {s['eta_sec']:8.4f} "
                             f"{s['n_hit']:4d} {s['delta_pct']:+7.2f}%")
        lines += ["", f"  {'APPLIED' if self.applied else 'DRY RUN -- nothing written'}"]
        if self.backup:
            lines.append(f"  original summary copied to {self.backup}")
        return "\n".join(lines)


def recompute_weights(run_dir, cfg, apply: bool = False,
                      progress=print) -> RescaleReport:
    """Recompute eta_shade / eta_block / eta_secondary and rescale the summary."""
    from . import field as field_mod
    from . import shading as shading_mod
    from .secondary import get_strategy
    from .store import RunStore

    store = RunStore(Path(run_dir), cfg=cfg, mode="r")
    summary = store.summary()
    report = RescaleReport(run=str(store.root))

    dup = int(summary.duplicated(subset=["timestep", "heliostat_id"]).sum())
    if dup:
        raise ValueError(
            f"{store.root} has {dup} duplicate (timestep, heliostat) row(s). "
            f"Rescaling would compound whatever produced them -- fix the run first."
        )

    full = field_mod.load_field(cfg)
    ids = [int(i) for i in store.manifest.get("heliostat_ids", [])]
    if not ids:
        raise ValueError(f"{store.root}/manifest.json has no heliostat_ids")
    index = {int(h): k for k, h in enumerate(full.ids)}
    fld = full.subset([index[h] for h in ids])

    strategy = get_strategy(cfg)
    cone = shading_mod.secondary_body(cfg)
    radius = shading_mod.search_radius_for(
        float(summary.solar_el_deg.min()),
        cfg.field.mirror_height_mm, cfg.field.mirror_width_mm,
    )
    neighbours = field_mod.neighbour_pairs(fld, radius)

    keys = sorted(summary.timestep.unique())
    updated = summary.set_index(["timestep", "heliostat_id"]).sort_index()
    worst_pointing = 0.0

    for n, key in enumerate(keys, 1):
        rows = summary[summary.timestep == key].set_index("heliostat_id").reindex(ids)
        az = float(rows.solar_az_deg.iloc[0])
        el = float(rows.solar_el_deg.iloc[0])

        solutions = [strategy.solve(float(fld.x_mm[i]), float(fld.y_mm[i]),
                                    az, el, cfg.geometry) for i in range(len(fld))]
        drift = max(
            max(abs(s.rot_az_deg - r), abs(s.rot_el_deg - e))
            for s, r, e in zip(solutions, rows.rot_az_deg, rows.rot_el_deg)
        )
        worst_pointing = max(worst_pointing, float(drift))
        if drift > POINTING_TOL_DEG:
            raise ValueError(
                f"{key}: re-solved pointing differs from the stored run by "
                f"{drift:.3e} deg. The geometry or strategy has changed since the "
                f"sweep, so the traced spots are stale too -- re-run rather than "
                f"rescale."
            )

        geoms, aims = shading_mod.build_geometries(fld, solutions, cfg)
        eta_shade, eta_block, eta_secondary = shading_mod.shading_blocking(
            geoms, aims, az, el, neighbours, secondary=cone
        )

        old_eff = (rows.eta_shade.to_numpy(float) * rows.eta_block.to_numpy(float))
        new_eff = eta_shade * eta_block
        # A timestep with the sun down has zero weight both before and after;
        # dividing would turn a legitimate zero into a nan.
        ratio = np.where(old_eff > 0, new_eff / np.where(old_eff > 0, old_eff, 1.0), 0.0)

        before = float((rows.power_w.to_numpy(float)).sum())
        after = float((rows.power_w.to_numpy(float) * ratio).sum())
        report.per_step.append({
            "timestep": key, "el": el,
            "eta_sec": float(eta_secondary.mean()),
            "n_hit": int((eta_secondary < 0.999).sum()),
            "delta_pct": 100.0 * (after / before - 1.0) if before else 0.0,
        })
        report.power_before_w += before
        report.power_after_w += after

        idx = pd.MultiIndex.from_product([[key], ids], names=["timestep", "heliostat_id"])
        for col in WEIGHTED:
            if col in updated.columns:
                updated.loc[idx, col] = updated.loc[idx, col].to_numpy(float) * ratio
        updated.loc[idx, "eta_shade"] = eta_shade
        updated.loc[idx, "eta_block"] = eta_block
        updated.loc[idx, "eta_secondary"] = eta_secondary
        updated.loc[idx, "shading_blocking_efficiency"] = new_eff

        if progress and (n == 1 or n % 10 == 0 or n == len(keys)):
            progress(f"  [{n}/{len(keys)}] {key}  el {el:5.1f}  "
                     f"eta_sec {eta_secondary.mean():.4f}")

    report.timesteps = len(keys)
    report.rows = len(summary)
    report.worst_pointing_deg = worst_pointing
    report.eta_shade_before = float(summary.eta_shade.mean())
    report.eta_block_before = float(summary.eta_block.mean())
    report.eta_shade_after = float(updated.eta_shade.mean())
    report.eta_block_after = float(updated.eta_block.mean())
    report.eta_secondary = float(updated.eta_secondary.mean())

    if apply:
        target = store.root / "summary.csv"
        backup = store.root / "summary_before_secondary_shading.csv"
        if not backup.exists():
            shutil.copy2(target, backup)
            report.backup = str(backup)
        out = updated.reset_index()
        # Restore the sweep's column order so the file stays diffable against
        # its own backup, with the new column appended rather than interleaved.
        cols = [c for c in summary.columns if c in out.columns]
        cols += [c for c in out.columns if c not in cols]
        out[cols].to_csv(target, index=False)
        report.applied = True

    return report
