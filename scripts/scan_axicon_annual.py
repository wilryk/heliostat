"""Axicon geometry scan: does (tip height, half-angle) beat the built design?

The cassegrain got a design scan (scan_cassegrain_annual.py); same fairness
for the axicon.  Same validated estimator -- cos(AOI) x exact UNION occlusion
(neighbours + the cone, per candidate) x rho^2 on the sweeps' own 94-step
grid, annualised by beamdown.energy.annual_energy.  The estimator reproduced
traced full8 to +0.64%, so trust DIFFERENCES between rows, not tenths.

Free knobs: cone tip height and half-angle.  Held fixed: receiver at
z = 7,000 mm (receiver_offset follows the tip so the receiver never moves),
aperture radius 15,000 mm (the 30 m manufacturability cap -- rim height
follows as tip + 15,000*tan(angle)).  Feasibility: every heliostat's beam
must meet the cone inside the aperture, above the tip radius.

Concentration is reported two ways, because the axicon's dominant real
aberration is NOT in the sun-disk term:

  ideal_r90    cosine-and-DNI-weighted ideal sun image r90 (0.5 deg over the
               true mirror->cone->receiver path, arrival stretch included).
  crowding     beam half-footprint on the cone / cone hit radius.  The cone's
               sagittal curvature goes like 1/radius, the heliostat's
               pre-correction is exact only at the footprint centre, so the
               residual blows up as hits crowd the axis (owner's finding).
               Calibration from traced full8: the built design's INNER
               mirrors hit at only 1,943 mm (max crowding ~0.34) and trace
               at r90 774 mm = 3.1x their sun limit, while its outer mirrors
               (crowding ~0.08) trace at 1.3x.  Treat any row whose
               max_crowding exceeds the built design's as a concentration
               REGRESSION whatever ideal_r90 says, and rows that push hits
               outward as the real opportunity.

Usage::

    python scripts/scan_axicon_annual.py
    python scripts/scan_axicon_annual.py --tips 26000 27000 28000 --angles 18 20 22
"""

from __future__ import annotations

import argparse
import datetime as _dt
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

TRACED_DATES = ["2026-12-21", "2026-01-21", "2026-02-20", "2026-03-20",
                "2026-04-21", "2026-05-21", "2026-06-21"]
F2_MM, AP_R_MM = 7000.0, 15000.0
BUILT = (27000.0, 20.0)          # the design full8 traced
FULL8_MWH = 10152.2


def lower_priority() -> None:
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--tips", type=float, nargs="+",
                   default=[24000, 27000, 30000, 33000])
    p.add_argument("--angles", type=float, nargs="+",
                   default=[15.0, 20.0, 25.0])
    p.add_argument("--nu", type=int, default=15)
    p.add_argument("--nv", type=int, default=9)
    a = p.parse_args(argv)
    lower_priority()

    from beamdown import energy as energy_mod
    from beamdown import field as field_mod
    from beamdown import shading as shading_mod
    from beamdown.config import load_config
    from beamdown.dni import provider_for
    from beamdown.secondary import get_strategy
    from beamdown.secondary.axicon import receiver_correction
    from beamdown.solar import build_time_grid

    cfg = load_config(None)
    full = field_mod.load_field(cfg)
    x = np.asarray(full.x_mm, float)
    y = np.asarray(full.y_mm, float)
    n = len(x)
    R = np.hypot(x, y)
    area, rho = cfg.field.mirror_area_m2, cfg.optics.mirror_reflectivity
    provider = provider_for(cfg)
    theta_half = np.deg2rad(0.25)                     # 0.5 deg full sun

    dates = [_dt.date.fromisoformat(d) for d in TRACED_DATES]
    steps = build_time_grid(cfg, dates)
    print(f"grid: {len(steps)} timesteps; sampling {a.nu}x{a.nv}; "
          f"DNI {provider.describe()}; estimator bias on the built design "
          f"was +0.64% vs traced {FULL8_MWH:,.1f} MWh")

    nb_cache: dict[float, list] = {}
    step_nb = []
    for s in steps:
        r = shading_mod.search_radius_for(
            s.solar_el_deg, cfg.field.mirror_height_mm, cfg.field.mirror_width_mm)
        key = round(r, 0)
        if key not in nb_cache:
            nb_cache[key] = field_mod.neighbour_pairs(full, r)
        step_nb.append(nb_cache[key])

    def annualise(power_per_step):
        rows = [{"date": str(s.date), "hour": s.hour, "heliostat_id": 0,
                 "power_w": float(pw), "solar_az_deg": s.solar_az_deg,
                 "solar_el_deg": s.solar_el_deg}
                for s, pw in zip(steps, power_per_step)]
        return energy_mod.annual_energy(pd.DataFrame(rows), cfg, provider,
                                        n_heliostats=n)

    rows = []
    for tip in a.tips:
        for ang in a.angles:
            drop = tip - F2_MM
            x_r, y_r, x_a, y_a = receiver_correction(R, tip, drop, ang)
            rim_z = tip + AP_R_MM * np.tan(np.deg2rad(ang))
            hit_ok = (x_a > 500.0) & (x_a < AP_R_MM)   # on the cone, above tip radius
            if not hit_ok.all():
                print(f"  tip {tip:,.0f} angle {ang:g}: INFEASIBLE -- "
                      f"{int((~hit_ok).sum())} beams miss the cone "
                      f"(hit radius {x_a.min():,.0f}..{x_a.max():,.0f})")
                continue

            # Sagittal-crowding metric, sun-independent like the geometry.
            aim = np.stack([x / R * x_r, y / R * x_r,
                            np.full(n, tip) + y_r], axis=1)
            mirror = np.stack([x, y, np.zeros(n)], axis=1)
            path_total = np.linalg.norm(aim - mirror, axis=1)
            path_after = np.hypot(x_a, (tip + y_a) - F2_MM)
            w_half = 0.5 * max(cfg.field.mirror_width_mm,
                               cfg.field.mirror_height_mm)
            crowding = (w_half * path_after / path_total) / x_a

            object.__setattr__(cfg.geometry, "secondary_height_mm", float(tip))
            object.__setattr__(cfg.geometry, "receiver_offset_mm", float(F2_MM - tip))
            object.__setattr__(cfg.geometry, "axicon_angle_deg", float(ang))
            strategy = get_strategy(cfg)
            body = shading_mod.secondary_body(cfg)

            t0 = time.monotonic()
            power, wsum, r90sum = [], 0.0, 0.0
            for k, s in enumerate(steps):
                sols = [strategy.solve(float(x[i]), float(y[i]), s.solar_az_deg,
                                       s.solar_el_deg, cfg.geometry)
                        for i in range(n)]
                geoms, aims = shading_mod.build_geometries(full, sols, cfg)
                occ = shading_mod.occlusion_efficiency(
                    geoms, aims, s.solar_az_deg, s.solar_el_deg, step_nb[k],
                    nu=a.nu, nv=a.nv, secondary=body)
                cosv = np.array([so.cosine_efficiency for so in sols])
                power.append(1000.0 * area * rho ** 2 * float((cosv * occ).sum()))
                # concentration proxy accumulates with the same weights
                path = np.array([so.focal_dist_mm for so in sols])
                stretch = 1.0 / np.cos(np.arctan2(x_a, (tip + y_a) - F2_MM))
                r90 = np.sqrt(0.9) * theta_half * path * np.sqrt(stretch)
                w = cosv * provider.dni(s.date, s.hour)
                r90sum += float((r90 * w).sum())
                wsum += float(w.sum())

            res = annualise(power)
            rows.append({
                "tip_m": tip / 1000, "angle_deg": ang, "rim_z_m": rim_z / 1000,
                "annual_MWh": res["annual_energy_mwh"],
                "eta": res["annual_optical_efficiency"],
                "ideal_r90_mm": r90sum / wsum,
                "min_hit_m": x_a.min() / 1000,
                "max_crowding": float(crowding.max()),
                "cone_depth_m": (rim_z - tip) / 1000,
            })
            print(f"  tip {tip:,.0f} angle {ang:g}: {res['annual_energy_mwh']:,.1f} MWh, "
                  f"ideal r90 {r90sum / wsum:,.0f} mm, inner hit "
                  f"{x_a.min():,.0f} mm, max crowding {crowding.max():.2f}  "
                  f"({time.monotonic() - t0:.0f}s)")

    tab = pd.DataFrame(rows).sort_values("annual_MWh", ascending=False)
    built = tab[(tab.tip_m == BUILT[0] / 1000) & (tab.angle_deg == BUILT[1])]
    if len(built):
        b = float(built.annual_MWh.iloc[0])
        tab["vs_built_pct"] = 100 * (tab.annual_MWh / b - 1)
        tab["r90_vs_built_pct"] = 100 * (tab.ideal_r90_mm
                                         / float(built.ideal_r90_mm.iloc[0]) - 1)
    with pd.option_context("display.width", 200, "display.float_format",
                           lambda v: f"{v:,.2f}"):
        print("\n" + tab.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
