"""Assert beamdown.secondary.axicon reproduces heliostat/heliostat_shape_solve.py
apart from one deliberate correction.

The axicon math was ported, not rewritten, so any *unexplained* difference here
is a porting bug. The one explained difference is the foreshortening factor in
``axicon_shape_correction``: the axicon's cylinder axis is not the mirror's own
sagittal axis, so the plain ``2 f cos(aoi)`` relation over-strengthens the
correction by 1/L**2. Forcing that factor to 1 must reproduce the legacy values
exactly, which is what keeps this a parity test rather than a loosened one.

Run with `python tests/test_axicon_parity.py` (no pytest required).
"""

import functools
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "legacy"))  # original scripts, for parity comparison

from beamdown.config import load_config          # noqa: E402
from beamdown import field as F                  # noqa: E402
from beamdown.secondary import axicon as AX      # noqa: E402
from beamdown.secondary import get_strategy      # noqa: E402
from heliostat.heliostat_shape_solve import get_heliostat_axicon_shape as legacy  # noqa: E402

KEYS = ["rot_az", "rot_el", "c3", "c4", "c5"]


def main() -> int:
    cfg = load_config()
    fld = F.load_field(cfg)
    # Explicitly "axicon", NOT cfg.optics.secondary: this test is about the axicon
    # port being faithful to the legacy solver, and it must keep testing that even
    # when config.toml has been switched to prime_focus or cassegrain for a run.
    strategy = get_strategy("axicon")
    geom = cfg.geometry

    # Sun positions spanning the whole traced envelope.
    sun_cases = [
        (77.9, 50.1), (41.9, 76.7), (294.3, 66.7), (271.7, 8.8),
        (65.6, 3.1), (12.2, 55.8), (112.1, 11.8), (155.9, 75.2),
        (247.7, 10.0), (330.0, 45.0),
    ]

    def sweep():
        """Largest |new - legacy| per coefficient, over every heliostat."""
        worst = dict.fromkeys(KEYS, 0.0)
        n = 0
        for az, el in sun_cases:
            for i in range(0, len(fld), 7):
                x, y = float(fld.x_mm[i]), float(fld.y_mm[i])
                old = legacy(
                    x, y, az, el,
                    geom.secondary_height_mm, geom.receiver_height_mm, geom.axicon_angle_deg,
                )
                new = strategy.solve(x, y, az, el, geom)
                got = (new.rot_az_deg, new.rot_el_deg, new.c3, new.c4, new.c5)
                for key, o, g in zip(KEYS, old, got):
                    worst[key] = max(worst[key], abs(float(o) - float(g)))
                n += 1
        return worst, n

    # With the foreshortening disabled the port must be bit-for-bit faithful.
    real = AX.axicon_shape_correction
    AX.axicon_shape_correction = functools.partial(real, foreshorten=1.0)
    try:
        flat, n = sweep()
    finally:
        AX.axicon_shape_correction = real
    live, _ = sweep()

    print(f"axicon parity: {n} cases ({len(sun_cases)} sun positions x "
          f"{len(range(0, len(fld), 7))} heliostats)")
    print(f"  {'':8s} {'L^2 forced to 1':>18s} {'as corrected':>18s}")
    for key in KEYS:
        print(f"  {key:8s} {flat[key]:18.3e} {live[key]:18.3e}")

    tol = 1e-12
    ported = all(v <= tol for v in flat.values())          # the port itself
    pointing = all(live[k] <= tol for k in ("rot_az", "rot_el"))
    corrected = all(live[k] > 0.0 for k in ("c3", "c4", "c5"))
    ok = ported and pointing and corrected
    print(f"\n  port faithful with L^2=1 : {ported}")
    print(f"  pointing still identical : {pointing}")
    print(f"  shape deliberately moved : {corrected}")
    print("PASS" if ok else f"FAIL (tolerance {tol:g})")

    # Diagnostics the legacy function threw away.
    sol = strategy.solve(float(fld.x_mm[0]), float(fld.y_mm[0]), 77.9, 50.1, geom)
    print(f"\nsample diagnostics: aoi {sol.aoi_deg:.2f} deg, "
          f"focal {sol.focal_dist_mm/1000:.1f} m, cos_eff {sol.cosine_efficiency:.4f}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
