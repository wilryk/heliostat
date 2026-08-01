"""Fast axicon iterator: fine (tip, angle) grid, owner's sagittal-correction cap.

The owner's design rule, encoded: a candidate is admissible only if the
sagittal curvature correction it demands of the WORST heliostat is no larger
than the built design's (tip 27,000, 20 deg).  That correction is pure
geometry -- the axicon solve's ``1/focal_dist_s - 1/focal_dist`` -- so the
whole evaluation is analytic and the grid runs in seconds:

  optics   from the field extremes (inner ring bears the correction and the
           crowding; outer ring bears coverage and the largest sun image);
  energy   cos(AOI) summed over the sweeps' 94-step grid x DNI, times an
           occlusion-vs-mean-aim-height curve interpolated from the exact
           union computations at the representative instant
           (scan_prime_focus_height.py); the axicon's own exact point sits
           on that curve to ~0.1%, which is what licenses the shortcut.
           RELATIVE index, built design = 1.0000 -- do not read as MWh.
           scan_axicon_annual.py is the exact (slow) cross-check.

Usage::

    python scripts/quick_axicon_iterate.py
    python scripts/quick_axicon_iterate.py --top 20
"""

from __future__ import annotations

import argparse
import datetime as _dt
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import numpy as np  # noqa: E402

TRACED_DATES = ["2026-12-21", "2026-01-21", "2026-02-20", "2026-03-20",
                "2026-04-21", "2026-05-21", "2026-06-21"]
F2_MM, AP_R_MM, TIP_RAD_MM = 7000.0, 15000.0, 500.0
BUILT = (27000.0, 20.0)

# Exact union occlusion at the representative instant vs shared-focus aim
# height (scan_prime_focus_height.py, disc at 32,460).  Monotone; used as a
# RELATIVE transfer curve over mean aim height.
OCC_H = np.array([30000., 32000., 33000., 34000., 35000., 36000., 38000.,
                  40000., 43000., 47000.])
OCC_V = np.array([0.9226, 0.9307, 0.9345, 0.9382, 0.9412, 0.9445, 0.9494,
                  0.9534, 0.9582, 0.9626])


def geometry_terms(R, tip, ang_deg):
    """Vectorised copy of the axicon solve's sun-independent geometry."""
    from beamdown.secondary.axicon import receiver_correction

    drop = tip - F2_MM
    x_r, y_r, x_a, y_a = receiver_correction(R, tip, drop, ang_deg)
    alpha = np.deg2rad(ang_deg)
    cone_dist = np.hypot(x_a, y_a)
    s_prime = -np.hypot(drop + y_a, x_a)
    radius_axicon = cone_dist / np.tan(alpha)
    axicon_aoi = np.arctan2(x_a, drop + y_a) + alpha
    s = 1.0 / (2.0 * np.cos(axicon_aoi) / radius_axicon - 1.0 / s_prime)
    return x_r, y_r, x_a, y_a, s, s_prime


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--tips", type=float, nargs=3, default=[24000, 33000, 500],
                   metavar=("LO", "HI", "STEP"))
    p.add_argument("--angles", type=float, nargs=3, default=[13.0, 22.0, 0.5],
                   metavar=("LO", "HI", "STEP"))
    p.add_argument("--top", type=int, default=12)
    a = p.parse_args(argv)

    from beamdown import field as F
    from beamdown.config import load_config
    from beamdown.dni import provider_for
    from beamdown.solar import build_time_grid

    cfg = load_config(None)
    full = F.load_field(cfg)
    x = np.asarray(full.x_mm, float)
    y = np.asarray(full.y_mm, float)
    R = np.hypot(x, y)
    n = len(x)
    theta_half = np.deg2rad(0.25)

    dates = [_dt.date.fromisoformat(d) for d in TRACED_DATES]
    steps = build_time_grid(cfg, dates)
    prov = provider_for(cfg)
    azr = np.deg2rad([s.solar_az_deg for s in steps])
    elr = np.deg2rad([s.solar_el_deg for s in steps])
    sun = np.stack([np.cos(elr) * np.cos(np.pi / 2 - azr),
                    np.cos(elr) * np.sin(np.pi / 2 - azr), np.sin(elr)], axis=1)
    dni = np.array([prov.dni(s.date, s.hour) for s in steps])

    def evaluate(tip, ang):
        x_r, y_r, x_a, y_a, s, s_prime = geometry_terms(R, tip, ang)
        aim = np.stack([x / R * x_r, y / R * x_r, np.full(n, tip + y_r)], axis=1)
        mirror = np.stack([x, y, np.zeros(n)], axis=1)
        to_aim = aim - mirror
        f = np.linalg.norm(to_aim, axis=1)
        f_s = f + s + s_prime
        dpow = np.abs(1.0 / f_s - 1.0 / f)          # extra sagittal power, 1/mm

        cov_ok = bool(((x_a > TIP_RAD_MM) & (x_a < AP_R_MM)).all()
                      and (f_s > 0).all())

        t = to_aim / f[:, None]
        cosaoi = np.cos(0.5 * np.arccos(np.clip(t @ sun.T, -1.0, 1.0)))
        cos_annual = float((cosaoi * dni[None, :]).sum())

        z_axis = aim[:, 2] * R / (R - x_r)          # x_r < 0: crossing height
        occ = float(np.interp(np.clip(z_axis.mean(), OCC_H[0], OCC_H[-1]),
                              OCC_H, OCC_V))

        path_after = np.hypot(x_a, (tip + y_a) - F2_MM)
        stretch = 1.0 / np.cos(np.arctan2(x_a, (tip + y_a) - F2_MM))
        r90 = np.sqrt(0.9) * theta_half * f * np.sqrt(stretch)
        w = (cosaoi * dni[None, :]).sum(axis=1)
        return {
            "cov_ok": cov_ok,
            "max_dpow": float(dpow.max()),
            "inner_hit_mm": float(x_a.min()),
            "energy_raw": cos_annual * occ,
            "r90_pw": float(np.average(r90, weights=w)),
            "rim_z": tip + AP_R_MM * np.tan(np.deg2rad(ang)),
        }

    built = evaluate(*BUILT)
    print(f"built design (tip 27,000, 20 deg): max sagittal correction "
          f"{built['max_dpow']:.3e} /mm, inner hit {built['inner_hit_mm']:,.0f} mm, "
          f"pw ideal r90 {built['r90_pw']:,.0f} mm")
    print(f"admissibility: max correction <= built's, full coverage, "
          f"F1-floor n/a (axicon)\n")

    rows = []
    for tip in np.arange(a.tips[0], a.tips[1] + 1, a.tips[2]):
        for ang in np.arange(a.angles[0], a.angles[1] + 1e-9, a.angles[2]):
            e = evaluate(float(tip), float(ang))
            if not e["cov_ok"] or e["max_dpow"] > built["max_dpow"] * (1 + 1e-9):
                continue
            rows.append({
                "tip_m": tip / 1000, "angle": ang, "rim_z_m": e["rim_z"] / 1000,
                "energy_idx": e["energy_raw"] / built["energy_raw"],
                "r90_pw_mm": e["r90_pw"],
                "r90_vs_built": e["r90_pw"] / built["r90_pw"] - 1,
                "inner_hit_m": e["inner_hit_mm"] / 1000,
                "corr_vs_built": e["max_dpow"] / built["max_dpow"] - 1,
            })

    import pandas as pd
    tab = pd.DataFrame(rows)
    print(f"{len(tab)} admissible candidates of the grid")
    for name, srt in (("BEST ENERGY", tab.sort_values("energy_idx", ascending=False)),
                      ("BEST CONCENTRATION (smallest pw ideal r90)",
                       tab.sort_values("r90_pw_mm"))):
        print(f"\n-- {name} --")
        with pd.option_context("display.width", 200, "display.float_format",
                               lambda v: f"{v:,.4f}"):
            print(srt.head(a.top).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
