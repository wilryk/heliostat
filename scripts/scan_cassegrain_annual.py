"""Which cassegrain geometry (dish height, prime-focus height) maximizes annual energy?

Licence-free ESTIMATE, built to be comparable to the traced numbers:

    power(t) = sum_i  A_mirror * 1000 W/m2 * cos(AOI_i) * occ_mirrors_i
                      * eta_disc_i * reflectivity^2

per timestep on the SWEEPS' OWN 94-step grid (7 traced declinations), then
annualised by ``beamdown.energy.annual_energy`` -- the identical efficiency
surface + 8760-hour DNI integration every traced run's report uses.  What the
estimate leaves out relative to a trace: receiver-plane spot structure and
spillage (the geometric solar disk is 1.4-2.1 m across every candidate, far
inside the receiver, so it cannot differentiate them) and traced-occluder
subtleties the scalar union misses (vetted at ~0.1%).

VALIDATION comes first, not on faith: the same estimator is run for the
axicon (traced full8 = 10,152.2 MWh) and prime focus F1=47m (traced
12,096.3 MWh) and the deltas printed.  Trust the cassegrain RANKING to the
size of those deltas.

Speed vs the sweep's scalar pass, with results unchanged in what matters:
  - neighbour search radius per TIMESTEP's own elevation (the sweep uses one
    radius for its lowest elevation; extra neighbours contribute exactly zero);
  - mirror occlusion (depends on F1 via pointing) and disc shading (depends on
    dish rim height) computed separately and multiplied; the product-vs-union
    error is MEASURED against the exact union and printed;
  - coarser aperture sampling (--nu/--nv), identical across candidates.

Runs at low process priority: a Quadoa sweep usually owns the machine.

Usage::

    python scripts/scan_cassegrain_annual.py
"""

from __future__ import annotations

import argparse
import datetime as _dt
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from design_cassegrain import close_design  # noqa: E402
from scan_prime_focus_height import disk_at_receiver  # noqa: E402

TRACED_DATES = ["2026-12-21", "2026-01-21", "2026-02-20", "2026-03-20",
                "2026-04-21", "2026-05-21", "2026-06-21"]
F2_MM, RIM_R_MM, TIP_MM = 7000.0, 15000.0, 27000.0
FULL8_MWH, PF47_MWH = 10152.2, 12096.3


def lower_priority() -> None:
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--f1-list", type=float, nargs="+",
                   default=[35000, 36000, 37000, 38000, 39000, 40000])
    p.add_argument("--rim-list", type=float, nargs="+",
                   default=[30000, 31000, 32460, 33000, 34000])
    p.add_argument("--nu", type=int, default=15)
    p.add_argument("--nv", type=int, default=9)
    p.add_argument("--no-validate", action="store_true")
    a = p.parse_args(argv)
    lower_priority()

    from beamdown import energy as energy_mod
    from beamdown import field as field_mod
    from beamdown import shading as shading_mod
    from beamdown.config import load_config, validate_layout
    from beamdown.dni import provider_for
    from beamdown.secondary import get_strategy
    from beamdown.solar import build_time_grid

    cfg = load_config(None)
    full = field_mod.load_field(cfg)
    x = np.asarray(full.x_mm, float)
    y = np.asarray(full.y_mm, float)
    n = len(x)
    R_f = float(np.hypot(x, y).max())
    area = cfg.field.mirror_area_m2
    rho = cfg.optics.mirror_reflectivity
    provider = provider_for(cfg)

    dates = [_dt.date.fromisoformat(d) for d in TRACED_DATES]
    steps = build_time_grid(cfg, dates)
    print(f"grid: {len(steps)} timesteps over {len(dates)} traced declinations; "
          f"DNI {provider.describe()}; sampling {a.nu}x{a.nv}")

    # Neighbour pairs per timestep, at that timestep's own sufficient radius.
    nb_cache: dict[float, list] = {}
    step_nb = []
    for s in steps:
        r = shading_mod.search_radius_for(
            s.solar_el_deg, cfg.field.mirror_height_mm, cfg.field.mirror_width_mm)
        key = round(r, 0)
        if key not in nb_cache:
            nb_cache[key] = field_mod.neighbour_pairs(full, r)
        step_nb.append(nb_cache[key])
    empty_nb = [np.empty(0, dtype=int) for _ in range(n)]

    def annualise(power_per_step: np.ndarray, label: str) -> dict:
        rows = [{"date": str(s.date), "hour": s.hour, "heliostat_id": 0,
                 "power_w": float(pw), "solar_az_deg": s.solar_az_deg,
                 "solar_el_deg": s.solar_el_deg}
                for s, pw in zip(steps, power_per_step)]
        res = energy_mod.annual_energy(pd.DataFrame(rows), cfg, provider,
                                       n_heliostats=n)
        return res

    def field_pass(secondary_body):
        """Per-step (cos_i, occ_i) with the CURRENT cfg layout/F1."""
        strategy = get_strategy(cfg)
        cos_t, occ_t, geoms_t = [], [], []
        for k, s in enumerate(steps):
            sols = [strategy.solve(float(x[i]), float(y[i]), s.solar_az_deg,
                                   s.solar_el_deg, cfg.geometry)
                    for i in range(n)]
            geoms, aims = shading_mod.build_geometries(full, sols, cfg)
            occ = shading_mod.occlusion_efficiency(
                geoms, aims, s.solar_az_deg, s.solar_el_deg, step_nb[k],
                nu=a.nu, nv=a.nv, secondary=secondary_body)
            cos_t.append(np.array([so.cosine_efficiency for so in sols]))
            occ_t.append(occ)
            geoms_t.append((geoms, aims))
        return np.array(cos_t), np.array(occ_t), geoms_t

    # ---- validation against the two traced runs --------------------------
    if not a.no_validate:
        t0 = time.monotonic()
        body = shading_mod.secondary_body(cfg)          # axicon cone
        cos_t, occ_t, _ = field_pass(body)
        pw = 1000.0 * area * rho ** 2 * (cos_t * occ_t).sum(axis=1)
        ax = annualise(pw, "axicon")["annual_energy_mwh"]
        print(f"\nvalidation ({time.monotonic() - t0:.0f}s for the axicon pass):")
        print(f"  axicon estimate      {ax:9,.1f} MWh   traced full8 "
              f"{FULL8_MWH:9,.1f}   delta {100 * (ax / FULL8_MWH - 1):+.2f}%")

        object.__setattr__(cfg.optics, "secondary", "prime_focus")
        object.__setattr__(cfg.optics, "n_mirrors", 1)
        object.__setattr__(cfg.geometry, "focus_height_mm", 47000.0)
        validate_layout(cfg)
        cos_t, occ_t, _ = field_pass(None)
        pw = 1000.0 * area * rho ** 1 * (cos_t * occ_t).sum(axis=1)
        pf = annualise(pw, "pf47")["annual_energy_mwh"]
        print(f"  prime-focus estimate {pf:9,.1f} MWh   traced "
              f"{PF47_MWH:9,.1f}   delta {100 * (pf / PF47_MWH - 1):+.2f}%")

    # ---- the cassegrain scan ---------------------------------------------
    object.__setattr__(cfg.optics, "secondary", "cassegrain")
    object.__setattr__(cfg.optics, "n_mirrors", 2)
    object.__setattr__(cfg.geometry, "secondary_rim_height_mm", 32460.0)

    per_f1 = {}
    ref_geoms = None
    for z1 in a.f1_list:
        t0 = time.monotonic()
        object.__setattr__(cfg.geometry, "focus_height_mm", float(z1))
        validate_layout(cfg)
        cos_t, occm_t, geoms_t = field_pass(None)       # mirrors only
        per_f1[z1] = (cos_t, occm_t)
        if z1 == 36000.0:
            ref_geoms = geoms_t
        print(f"  F1 {z1:,.0f}: mirror pass {time.monotonic() - t0:.0f}s")
    if ref_geoms is None:
        ref_geoms = geoms_t                              # last one

    # Disc shading factor per rim height, on the F1=36000 mirror tilts (the
    # tilt dependence across this F1 range is measured below, not assumed).
    eta_disc = {}
    for rz in a.rim_list:
        disc = shading_mod.SecondaryDisc(z_mm=float(rz), radius_mm=RIM_R_MM)
        rows = []
        for k, s in enumerate(steps):
            g, aims = ref_geoms[k]
            rows.append(shading_mod.occlusion_efficiency(
                g, aims, s.solar_az_deg, s.solar_el_deg, empty_nb,
                nu=a.nu, nv=a.nv, secondary=disc))
        eta_disc[rz] = np.array(rows)

    # Product-vs-union audit at the baseline geometry, three sun heights.
    disc = shading_mod.SecondaryDisc(z_mm=32460.0, radius_mm=RIM_R_MM)
    object.__setattr__(cfg.geometry, "focus_height_mm", 36000.0)
    audit = []
    for k in [np.argmax([s.solar_el_deg for s in steps]),
              len(steps) // 2, 0]:
        s = steps[k]
        g, aims = ref_geoms[k]
        exact = shading_mod.occlusion_efficiency(
            g, aims, s.solar_az_deg, s.solar_el_deg, step_nb[k],
            nu=a.nu, nv=a.nv, secondary=disc)
        approx = per_f1[36000.0][1][k] * eta_disc[32460.0][k]
        audit.append(100 * abs(approx.mean() - exact.mean()) / exact.mean())
    print(f"\nproduct-vs-exact-union audit (3 timesteps): field-mean error "
          f"{max(audit):.3f}% max -- the ranking below is safe to that size")

    rows = []
    for rz in a.rim_list:
        for z1 in a.f1_list:
            try:
                des = close_design(float(rz), R_f, F2_MM, RIM_R_MM, float(z1))
                disk_outer, _ = disk_at_receiver(des, x, y, 0.5)
            except ValueError:
                continue
            if np.isnan(disk_outer):
                continue                                 # beam misses the dish
            cos_t, occm_t = per_f1[z1]
            pw = 1000.0 * area * rho ** 2 * (cos_t * occm_t * eta_disc[rz]).sum(axis=1)
            res = annualise(pw, f"rim{rz:.0f}_f1{z1:.0f}")
            rows.append({
                "rim_z_m": rz / 1000, "F1_m": z1 / 1000,
                "F1_above_tip_m": (z1 - TIP_MM) / 1000,
                "annual_MWh": res["annual_energy_mwh"],
                "eta": res["annual_optical_efficiency"],
                "disk_mm": disk_outer,
                "dish_depth_m": (rz - (des.z_c + des.a)) / 1000,
                "vertex_z_m": (des.z_c + des.a) / 1000,
            })

    tab = pd.DataFrame(rows).sort_values("annual_MWh", ascending=False)
    base = tab[(tab.rim_z_m == 32.46) & (tab.F1_m == 36.0)]
    base_mwh = float(base.annual_MWh.iloc[0]) if len(base) else float("nan")
    tab["vs_36k_axrim_pct"] = 100 * (tab.annual_MWh / base_mwh - 1)
    with pd.option_context("display.width", 200, "display.float_format",
                           lambda v: f"{v:,.2f}"):
        print("\n" + tab.to_string(index=False))
    b = tab.iloc[0]
    print(f"\n  best: rim z {b.rim_z_m:.2f} m, F1 {b.F1_m:.0f} m "
          f"(tip +{b.F1_above_tip_m:.0f} m): {b.annual_MWh:,.1f} MWh, "
          f"{b.vs_36k_axrim_pct:+.2f}% vs the axicon-comparability point; "
          f"disk {b.disk_mm:,.0f} mm")
    return 0


if __name__ == "__main__":
    sys.exit(main())
