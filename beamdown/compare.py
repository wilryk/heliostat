"""Compare two sweeps of the same field -- "is this one actually better?"

Written for the axicon shape-correction A/B, but nothing here knows about that:
it compares any two runs over the same heliostats and timesteps.

The headline metric is **power delivered inside an aperture of radius R**, swept
over R rather than evaluated at one radius. A single radius cannot settle the
question: a change that tightens the tail while shaving the peak wins at small R
and loses at large R, and quoting whichever radius suits is not an answer. If one
run dominates at *every* radius, that is a real result; if the curves cross, the
crossing radius is the useful number and it belongs in the report.

Both runs emit the same number of rays per heliostat and get the same DNI,
reflectivity and shading weights, so absolute watts are directly comparable and
no normalisation is needed -- and nothing is re-traced, since this reads stored
counts.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

DEFAULT_RADII_MM = (300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1250.0, 1500.0)


@dataclass
class ComparisonReport:
    label_a: str
    label_b: str
    keys: list
    radii_mm: np.ndarray
    captured_a: np.ndarray          # (R,) watts summed over timesteps
    captured_b: np.ndarray
    per_timestep: pd.DataFrame
    per_heliostat: pd.DataFrame
    crossover_mm: float | None
    reference_mm: float
    notes: list = field(default_factory=list)

    def describe(self) -> str:
        L = [f"{self.label_a}  vs  {self.label_b}",
             f"{len(self.keys)} shared timestep(s), "
             f"{self.per_heliostat.heliostat_id.nunique()} heliostats",
             ""]
        L.append("power inside an aperture of radius R (summed over timesteps)")
        L.append(f"  {'R (mm)':>8} {self.label_a[:12]:>14} {self.label_b[:12]:>14} "
                 f"{'delta':>12} {'':>8}")
        for r, a, b in zip(self.radii_mm, self.captured_a, self.captured_b):
            d = (b - a) / a * 100 if a else float("nan")
            L.append(f"  {r:8.0f} {a/1000:13.1f}k {b/1000:13.1f}k {d:+11.2f}% "
                     f"{'better' if d > 0 else 'worse' if d < 0 else '=':>8}")
        rel = np.max(np.abs(self.captured_b - self.captured_a)
                     / np.where(self.captured_a > 0, self.captured_a, 1.0))
        if rel < 1e-12:
            L.append("\n  the two runs are identical to machine precision")
        elif self.crossover_mm is None:
            better = self.captured_b[0] > self.captured_a[0]
            L.append(f"\n  no crossover -- {self.label_b if better else self.label_a} "
                     f"delivers more at every radius tested")
        else:
            L.append(f"\n  curves cross at R = {self.crossover_mm:.0f} mm")
        L += ["", f"per timestep, at R = {self.reference_mm:.0f} mm"]
        L.append(f"  {'timestep':>16} {'sun el':>7} {self.label_a[:10]:>12} "
                 f"{self.label_b[:10]:>12} {'delta':>9}")
        for _, row in self.per_timestep.iterrows():
            L.append(f"  {row.timestep:>16} {row.solar_el_deg:7.2f} "
                     f"{row.captured_a/1000:11.1f}k {row.captured_b/1000:11.1f}k "
                     f"{row.delta_pct:+8.2f}%")
        L += ["", *self.notes]
        return "\n".join(L)


def _radial_mask(cfg, radii_mm):
    """Boolean masks over the receiver grid, one per aperture radius.

    Kept as a thin alias so the GUI, the figures and this comparison all agree on
    what "inside the aperture" means -- a divergence here would show up as two
    tools reporting different captured power for the same run.
    """
    from .metrics import radial_masks

    return radial_masks(cfg, radii_mm)


def _weights(summary, key, ids, use_shading, notes=None, manifest=None):
    """Per-heliostat efficiency weights, aligned to the manifest's row order.

    Which columns carry the weight is :func:`beamdown.store.occlusion_weight_columns`'s
    call, from ``manifest``: a run that put its neighbours in the ray path has
    shading and blocking in the counts already and must not have them applied
    again, and a scalar run applies either the union column or the older
    eta_shade x eta_block. Getting this wrong charges the same loss twice and
    shows up as a uniform deficit at every aperture radius -- which reads
    exactly like a real optical result and is not one.

    The summary is appended to per timestep, so a resumed -- or accidentally
    duplicated -- run can hold more than one row for the same
    (timestep, heliostat). Keeping the last is right for a resume, where the
    later trace supersedes the earlier, but a duplicate is worth saying out loud
    rather than silently averaging away: it can also mean two processes wrote
    the same run, in which case the counts on disk are not trustworthy either.
    """
    rows = summary[summary.timestep == key]
    dup = int(rows.heliostat_id.duplicated().sum())
    if dup and notes is not None:
        notes.append(f"  WARNING: {key} has {dup} duplicate heliostat row(s) in the "
                     f"summary; kept the last of each. Check the run was not written "
                     f"by two processes.")
    rows = rows.drop_duplicates(subset="heliostat_id", keep="last").set_index("heliostat_id")
    w = np.ones(len(ids))
    if use_shading:
        from .store import occlusion_weight_columns

        for col in occlusion_weight_columns(manifest or {}, rows.columns):
            if col in rows.columns:
                w = w * rows[col].reindex(ids).fillna(1.0).to_numpy()
    return w


def compare_runs(store_a, store_b, cfg, label_a="A", label_b="B",
                 radii_mm=DEFAULT_RADII_MM, reference_mm=700.0,
                 dni_w_m2=1000.0, use_shading=True) -> ComparisonReport:
    """Compare two stores over every timestep and heliostat they share."""
    from .store import scale_factor

    radii_mm = np.asarray(radii_mm, float)
    ids_a = list(store_a.manifest.get("heliostat_ids", []))
    ids_b = list(store_b.manifest.get("heliostat_ids", []))
    ids = [h for h in ids_a if h in set(ids_b)]
    if not ids:
        raise ValueError("the two runs share no heliostats")
    row_a = {h: i for i, h in enumerate(ids_a)}
    row_b = {h: i for i, h in enumerate(ids_b)}
    take_a = np.array([row_a[h] for h in ids])
    take_b = np.array([row_b[h] for h in ids])

    sum_a, sum_b = store_a.summary(), store_b.summary()
    keys = sorted(set(store_a.timestep_keys()) & set(store_b.timestep_keys())
                  & set(sum_a.timestep) & set(sum_b.timestep))
    if not keys:
        raise ValueError("the two runs share no completed timesteps")

    masks = _radial_mask(cfg, radii_mm)
    ref_i = int(np.argmin(np.abs(radii_mm - reference_mm)))
    cap_a = np.zeros(len(radii_mm))
    cap_b = np.zeros(len(radii_mm))
    per_step, per_helio, notes = [], [], []

    for key in keys:
        step = {}
        for tag, store, take, summ in (("a", store_a, take_a, sum_a),
                                       ("b", store_b, take_b, sum_b)):
            counts = np.asarray(store.read_counts(key))[take].astype(np.float64)
            w = _weights(summ, key, ids, use_shading, notes,
                         manifest=store.manifest)
            rays = int(store.manifest.get("rays_per_heliostat",
                                          cfg.trace.rays_per_heliostat))
            watts = scale_factor(cfg, rays, dni_w_m2)
            wc = counts * w[:, None, None]
            step[tag] = dict(
                field=wc.sum(axis=0) * watts,
                per_h=np.einsum("hij,rij->hr", wc, masks.astype(np.float64)) * watts,
            )
        fa = np.array([step["a"]["field"][m].sum() for m in masks])
        fb = np.array([step["b"]["field"][m].sum() for m in masks])
        cap_a += fa
        cap_b += fb

        first = sum_a[sum_a.timestep == key].iloc[0]
        per_step.append(dict(timestep=key,
                             solar_az_deg=float(first.solar_az_deg),
                             solar_el_deg=float(first.solar_el_deg),
                             captured_a=fa[ref_i], captured_b=fb[ref_i],
                             delta_pct=(fb[ref_i] - fa[ref_i]) / fa[ref_i] * 100))

        ha, hb = step["a"]["per_h"][:, ref_i], step["b"]["per_h"][:, ref_i]
        per_helio.append(pd.DataFrame(dict(
            timestep=key, heliostat_id=ids, captured_a=ha, captured_b=hb,
            delta_w=hb - ha,
            delta_pct=np.where(ha > 0, (hb - ha) / np.where(ha > 0, ha, 1) * 100, np.nan),
        )))

    # Crossover: the radius where the sign of the difference flips.
    diff = cap_b - cap_a
    cross = None
    sign = np.sign(diff)
    flip = np.nonzero(np.diff(sign) != 0)[0]
    if len(flip):
        i = int(flip[0])
        t = abs(diff[i]) / (abs(diff[i]) + abs(diff[i + 1]))
        cross = float(radii_mm[i] + t * (radii_mm[i + 1] - radii_mm[i]))

    return ComparisonReport(
        label_a=label_a, label_b=label_b, keys=keys, radii_mm=radii_mm,
        captured_a=cap_a, captured_b=cap_b,
        per_timestep=pd.DataFrame(per_step),
        per_heliostat=pd.concat(per_helio, ignore_index=True),
        crossover_mm=cross, reference_mm=float(radii_mm[ref_i]),
        notes=sorted(set(notes)),
    )


def _require_axicon(cfg, what: str) -> None:
    """Refuse to compute an axicon-only quantity for another layout.

    ``foreshortening`` and everything built on it measure the size of the
    *axicon's* sagittal correction. For prime focus and Cassegrain that correction
    is identically zero, so these functions would not merely be inapplicable --
    they would happily return 0.0 and a flat attribution table, which reads as a
    real (null) result rather than as the wrong question.
    """
    layout = cfg.optics.secondary
    if layout != "axicon":
        raise ValueError(
            f"{what} is axicon-only: it measures the size of the axicon's "
            f"sagittal shape correction, which layout {layout!r} does not have "
            f"(its heliostats carry no secondary correction at all). Point cfg at "
            f"an axicon run, or drop this analysis for {layout!r}."
        )


def foreshortening(cfg, x_mm, y_mm, solar_az_deg, solar_el_deg) -> float:
    """``L**2`` for one heliostat -- the axicon shape correction's size.

    Exposed here so a comparison can check *which* heliostats moved, not just
    that the total moved. If the run-to-run improvement really comes from the
    shape correction, it has to concentrate where ``1/L**2`` is largest; if it
    shows up somewhere else, the change is doing something other than advertised.
    """
    from .secondary import axicon as AX

    _require_axicon(cfg, "compare.foreshortening")
    g = cfg.geometry
    sec_h = g.secondary_height_mm
    drop = sec_h - g.receiver_height_mm
    radius = float(np.hypot(x_mm, y_mm))
    offset, height, _, _ = AX.receiver_correction(radius, sec_h, drop, g.axicon_angle_deg)
    aim = np.array([x_mm / radius * offset, y_mm / radius * offset, sec_h + height])
    mirror = np.array([float(x_mm), float(y_mm), 0.0])

    *_, u, v = AX.heliostat_orientation(aim, mirror, solar_az_deg, solar_el_deg)
    to_aim = aim - mirror
    direction = to_aim / np.linalg.norm(to_aim)
    horizontal = np.array([direction[0], direction[1], 0.0])
    horizontal /= np.linalg.norm(horizontal)
    sag = np.cross(horizontal, np.array([0.0, 0.0, 1.0]))
    sag /= np.linalg.norm(sag)
    return float(np.dot(sag, u) ** 2 + np.dot(sag, v) ** 2)


def predicted_shape_change(cfg, x_mm, y_mm, solar_az_deg, solar_el_deg):
    """How far the foreshortening moves this heliostat's coefficients.

    ``1/L**2`` on its own is *not* the right predictor, which is worth being
    explicit about: the coefficient change is ``(1/L**2 - 1)`` times the axicon
    term, and that term's own weight varies across the field independently of
    ``L``. A heliostat can have a large ``1/L**2`` and still barely move because
    its axicon correction was small to begin with. Measured over one timestep,
    banding by this quantity correlates with delivered power at r = +0.83 versus
    +0.52 for ``1/L**2`` -- and only this one comes out monotonic.
    """
    import functools

    from .secondary import axicon as AX
    from .secondary import get_strategy

    _require_axicon(cfg, "compare.predicted_shape_change")
    strategy = get_strategy(cfg)
    real = AX.axicon_shape_correction
    corrected = strategy.solve(x_mm, y_mm, solar_az_deg, solar_el_deg, cfg.geometry)
    AX.axicon_shape_correction = functools.partial(real, foreshorten=1.0)
    try:
        plain = strategy.solve(x_mm, y_mm, solar_az_deg, solar_el_deg, cfg.geometry)
    finally:
        AX.axicon_shape_correction = real
    return (float(np.hypot(corrected.c3 - plain.c3, corrected.c5 - plain.c5)),
            float(abs(corrected.c4 - plain.c4)))


def attribute_by_foreshortening(report, cfg, field, n_bins: int = 5) -> pd.DataFrame:
    """Group the per-heliostat change by how far the correction moved it.

    The consistency check: power gained must concentrate where the coefficients
    actually changed, and must be ~0 where they did not. A flat trend, or a gain
    in the untouched band, means something other than the shape correction moved.
    """
    _require_axicon(cfg, "compare.attribute_by_foreshortening")
    rows_by_id = {int(h): i for i, h in enumerate(field.ids)}
    per = report.per_heliostat
    rows = []
    for key in report.keys:
        az, el = _sun_for_key(report, key)
        sub = per[per.timestep == key]
        for hid, dw in zip(sub.heliostat_id, sub.delta_w):
            i = rows_by_id[int(hid)]
            x, y = float(field.x_mm[i]), float(field.y_mm[i])
            d_astig, _ = predicted_shape_change(cfg, x, y, az, el)
            rows.append(dict(heliostat_id=int(hid), timestep=key,
                             inv_L2=1.0 / foreshortening(cfg, x, y, az, el),
                             d_astig=d_astig, delta_w=dw))
    df = pd.DataFrame(rows)
    df["band"] = pd.qcut(df.d_astig, n_bins, duplicates="drop")
    out = (df.groupby("band", observed=True)
             .agg(heliostats=("heliostat_id", "size"),
                  d_astig_mean=("d_astig", "mean"),
                  inv_L2_mean=("inv_L2", "mean"),
                  delta_w_mean=("delta_w", "mean"),
                  delta_w_total=("delta_w", "sum"))
             .reset_index())
    out.attrs["r_d_astig"] = float(np.corrcoef(df.d_astig, df.delta_w)[0, 1])
    out.attrs["r_inv_L2"] = float(np.corrcoef(df.inv_L2, df.delta_w)[0, 1])
    out.attrs["monotonic"] = bool(out.delta_w_mean.is_monotonic_increasing)
    return out


def _sun_for_key(report, key):
    row = report.per_timestep[report.per_timestep.timestep == key].iloc[0]
    return float(row.solar_az_deg), float(row.solar_el_deg)
