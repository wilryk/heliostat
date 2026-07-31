"""Export a Quadoa model set up for specific heliostats at a specific instant.

Closes the loop between analysis and inspection. When the summary table or the
explorer flags a heliostat -- worst power, biggest spot, most blocked -- this
writes an ``.optx`` with that heliostat's exact pointing, shape and sun position
already loaded, so it can be opened in the Quadoa GUI and examined directly:
3D view, spot diagram, whatever the question needs.

Several heliostats can be exported at once, one per configuration, so they can
be flipped through in the GUI's multiconfig selector and compared.

All four sequences remain available in the exported model; the GUI's sequence
selector picks between them. Sequence 2 (index 1) has the smallest ray count and
redraws the 3D view fastest, which is usually what you want for looking around.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class ExportReport:
    path: Path
    timestep: str
    solar_az_deg: float
    solar_el_deg: float
    heliostats: list[dict]
    columns: int

    def describe(self) -> str:
        lines = [
            f"wrote {self.path}",
            f"  timestep {self.timestep}   sun az {self.solar_az_deg:.2f}  "
            f"el {self.solar_el_deg:.2f}",
            f"  {len(self.heliostats)} heliostat(s) in {self.columns} configuration(s)",
            "",
            f"  {'cfg':>3} {'id':>5} {'x (m)':>9} {'y (m)':>9} {'rot_az':>9} "
            f"{'rot_el':>8} {'aoi':>7} {'peak kW/m2':>11} {'r90 mm':>8}",
        ]
        for h in self.heliostats:
            lines.append(
                f"  {h['config']:>3} {h['heliostat_id']:>5} {h['x_m']:>9.2f} "
                f"{h['y_m']:>9.2f} {h['rot_az_deg']:>9.3f} {h['rot_el_deg']:>8.3f} "
                f"{h['aoi_deg']:>7.2f} "
                + (f"{h['peak_flux_kw_m2']:>11.1f}" if h.get('peak_flux_kw_m2') is not None else f"{'-':>11}")
                + (f" {h['r90_mm']:>8.1f}" if h.get('r90_mm') is not None else f" {'-':>8}")
            )
        lines += [
            "",
            "  Open in Quadoa and use the multiconfig selector to switch heliostats.",
            "  Sequence 2 (index 1) redraws the 3D view fastest.",
        ]
        return "\n".join(lines)


def _timestep_sun(summary, timestep: str):
    rows = summary[summary.timestep == timestep]
    if not len(rows):
        raise KeyError(f"timestep {timestep!r} not in the store")
    first = rows.iloc[0]
    return float(first["solar_az_deg"]), float(first["solar_el_deg"]), rows


def export_for_inspection(
    cfg,
    heliostat_ids,
    timestep: str = None,
    summary=None,
    solar_az_deg: float = None,
    solar_el_deg: float = None,
    out_path=None,
    strategy=None,
) -> ExportReport:
    """Write an .optx holding the given heliostats at the given instant.

    Either pass ``timestep`` together with a store ``summary`` (the usual route,
    so the sun position and expected metrics come from the actual run), or pass
    ``solar_az_deg``/``solar_el_deg`` directly for a hypothetical instant.
    """
    from . import field as F
    from .model_edit import expand_multiconfig
    from .secondary import get_strategy
    from .session import QuadoaSession

    heliostat_ids = [int(h) for h in np.atleast_1d(heliostat_ids)]
    strategy = strategy or get_strategy(cfg)

    metrics_by_id: dict[int, dict] = {}
    if timestep is not None:
        if summary is None:
            raise ValueError("timestep requires the store summary")
        solar_az_deg, solar_el_deg, rows = _timestep_sun(summary, timestep)
        for hid in heliostat_ids:
            match = rows[rows.heliostat_id == hid]
            if len(match):
                r = match.iloc[0]
                metrics_by_id[hid] = {
                    "peak_flux_kw_m2": float(r.get("peak_flux_w_m2", np.nan)) / 1000.0,
                    "r90_mm": float(r.get("r90_mm", np.nan)),
                    "eta_shade": float(r.get("eta_shade", np.nan)),
                    "eta_block": float(r.get("eta_block", np.nan)),
                }
    if solar_az_deg is None or solar_el_deg is None:
        raise ValueError("need either a timestep or explicit sun angles")

    label = timestep or f"az{solar_az_deg:.0f}_el{solar_el_deg:.0f}"
    if out_path is None:
        stem = "_".join(str(h) for h in heliostat_ids[:4])
        suffix = "" if len(heliostat_ids) <= 4 else f"_plus{len(heliostat_ids)-4}"
        out_path = cfg.path(f"models/inspect_h{stem}{suffix}_{label}.optx")
    out_path = Path(out_path)

    # The exported model needs at least one column per heliostat. Never fewer
    # than the source has -- shrinking would discard its stored values, and the
    # spare columns are harmless.
    from .model_edit import current_columns

    n = len(heliostat_ids)
    have = current_columns(cfg.model_path.read_text(encoding="utf-8"))
    report = expand_multiconfig(cfg.model_path, out_path, max(n, have))

    full = F.load_field(cfg)
    id_to_row = {int(i): k for k, i in enumerate(full.ids)}

    session = QuadoaSession(cfg)
    entries = []
    try:
        # Work on the exported copy, not the source model.
        session.core.loadModelFile(str(out_path))
        session.core.applyChangesAndInitModel()
        session.set_global_geometry()
        session.set_sun(solar_az_deg, solar_el_deg)

        for config, hid in enumerate(heliostat_ids):
            if hid not in id_to_row:
                raise KeyError(f"heliostat id {hid} is not in the field")
            row = id_to_row[hid]
            x, y = float(full.x_mm[row]), float(full.y_mm[row])
            sol = strategy.solve(x, y, solar_az_deg, solar_el_deg, cfg.geometry)
            session.set_heliostat(x, y, sol, config)

            entry = {
                "config": config,
                "heliostat_id": hid,
                "x_m": x / 1000.0,
                "y_m": y / 1000.0,
                "rot_az_deg": sol.rot_az_deg,
                "rot_el_deg": sol.rot_el_deg,
                "aoi_deg": sol.aoi_deg,
                "c3": sol.c3, "c4": sol.c4, "c5": sol.c5,
            }
            entry.update(metrics_by_id.get(hid, {}))
            entries.append(entry)

        session.core.setConfig(0)
        session.core.applyChangesAndInitModel()
        session.core.saveModelFile(str(out_path))
    finally:
        session.close()

    return ExportReport(
        path=out_path,
        timestep=label,
        solar_az_deg=solar_az_deg,
        solar_el_deg=solar_el_deg,
        heliostats=entries,
        columns=report["columns_after"],
    )


def pick_heliostats(summary, timestep: str = None, by: str = "power_w",
                    n: int = 1, worst: bool = True) -> list[int]:
    """Select heliostats by a summary metric, for exporting.

    With a ``timestep`` the ranking is at that instant; without one it is
    averaged over the whole run, which is usually the fairer question -- a
    heliostat can look fine at noon and be badly blocked morning and evening.
    """
    data = summary[summary.timestep == timestep] if timestep else summary
    if timestep:
        ranked = data.sort_values(by, ascending=worst)
    else:
        agg = "mean" if by in {
            "transmission", "eta_shade", "eta_block", "peak_flux_w_m2",
            "r90_mm", "r50_mm", "cosine_efficiency", "spillage",
        } else "sum"
        ranked = (data.groupby("heliostat_id", as_index=False)
                  .agg(**{by: (by, agg)})
                  .sort_values(by, ascending=worst))
    return [int(h) for h in ranked.heliostat_id.head(n)]
