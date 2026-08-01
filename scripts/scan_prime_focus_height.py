"""Scan cassegrain prime-focus height (and dish height): disk size + occlusion.

Two questions, both licence-free, both answered on the sweeps' own math:

1. SOLAR DISK SIZE at the receiver, as a function of prime-focus height z1
   and dish rim height z_r (receiver F2 and 15 m rim radius fixed).  The
   heliostat focuses the sun (default 0.5 deg full angle) into a disk of
   diameter theta * slant at F1; the hyperboloid relays that image to F2
   magnified by |P->F2| / |P->F1|, with P the exact ray/quadric intercept
   from ``design_cassegrain.Cassegrain.intersect``.  No paraxial guessing,
   no aberrations: this is the pure geometric solar image.  Configurations
   where any heliostat's beam misses the dish are marked infeasible.

2. SHADING / BLOCKING at one instant (default 2026-02-20 09:27, the
   average-AOI instant) as a function of z1, using exactly the machinery the
   sweep's scalar branch uses: ``shading.shading_blocking`` and the UNION
   ``shading.occlusion_efficiency`` over the full 645-mirror field, with the
   cassegrain disc (15 m at z_r) as the shading body.  The axicon at the
   same instant is printed as the reference row.  Pointing for a cassegrain
   field is the shared-focus solve, so blocking depends on z1 through the
   beam's upward angle -- lower F1, flatter beams, more mirror-on-mirror
   blocking.  This quantifies where "too low" starts.

Usage::

    python scripts/scan_prime_focus_height.py
    python scripts/scan_prime_focus_height.py --f1-min 33000 --sun-deg 0.53
"""

from __future__ import annotations

import argparse
import datetime as _dt
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import numpy as np  # noqa: E402

from design_cassegrain import close_design, rays_to_f1  # noqa: E402

F2_MM = 7000.0            # receiver, fixed (27000 - 20000; config echoes this)
RIM_R_MM = 15000.0        # aperture radius, fixed (the axicon's)
TIP_MM = 27000.0          # axicon tip, for the "+x m" labels


def disk_at_receiver(des, x_mm, y_mm, sun_full_deg: float):
    """(outer-heliostat disk, field-max disk) diameter at F2, mm; NaN if any miss.

    diameter = theta * slant(heliostat->F1) * |P->F2| / |P->F1|.
    """
    theta = np.deg2rad(sun_full_deg)
    o, d = rays_to_f1(x_mm, y_mm, des)
    t, hit, ok, _ = des.intersect(o, d)
    r_hit = np.hypot(hit[:, 0], hit[:, 1])
    if not (ok.all() and (r_hit <= des.rim_r * (1.0 + 1e-9)).all()):
        return float("nan"), float("nan")
    slant = np.linalg.norm(des.F1[None, :] - o, axis=1)
    u = np.linalg.norm(des.F1[None, :] - hit, axis=1)
    v = np.linalg.norm(des.F2[None, :] - hit, axis=1)
    disk = theta * slant * v / u
    R = np.hypot(x_mm, y_mm)
    return float(disk[np.argmax(R)]), float(disk.max())


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--sun-deg", type=float, default=0.5,
                   help="solar full field of view, deg (default 0.5)")
    p.add_argument("--f1-min", type=float, default=35000.0,
                   help="lowest prime focus to scan, mm (default 35000 = tip+8m,"
                        " the owner's blocking floor)")
    p.add_argument("--f1-max", type=float, default=47000.0)
    p.add_argument("--date", default="2026-02-20")
    p.add_argument("--hour", type=float, default=9.45401542341586,
                   help="local decimal hour (default: the average-AOI instant)")
    p.add_argument("--rim-z", type=float, default=32460.0,
                   help="dish rim height used in part 2, mm")
    a = p.parse_args(argv)

    from beamdown import field as field_mod
    from beamdown import shading as shading_mod
    from beamdown.config import load_config, validate_layout
    from beamdown.secondary import get_strategy
    from beamdown.solar import sun_position

    cfg = load_config(None)
    full = field_mod.load_field(cfg)
    x = np.asarray(full.x_mm, float)
    y = np.asarray(full.y_mm, float)

    # ---- part 1: disk size over (rim height, F1) -------------------------
    rims = [30000.0, 31000.0, 32460.0, 33000.0, 34000.0]
    f1s = np.arange(a.f1_min, a.f1_max + 1.0, 1000.0)
    print(f"solar disk diameter at the receiver, OUTER heliostat, mm "
          f"(sun {a.sun_deg:g} deg; '--' = beam misses the dish, i.e. "
          f"F1 above the aperture-fill limit for that rim height)")
    print(f"receiver F2 z={F2_MM:,.0f}, rim radius {RIM_R_MM:,.0f}, both fixed\n")
    hdr = "  rim z (mm) |" + "".join(f" {int(z1 / 1000):>5} km"[:-2] + "m"
                                     for z1 in f1s)
    hdr = "  rim z (mm) |" + "".join(f" F1={z1 / 1000:>4.0f}m" for z1 in f1s)
    print(hdr)
    best = (float("inf"), None, None)
    for rz in rims:
        cells = []
        for z1 in f1s:
            try:
                des = close_design(rz, float(np.hypot(x, y).max()), F2_MM,
                                   RIM_R_MM, float(z1))
                outer, worst = disk_at_receiver(des, x, y, a.sun_deg)
            except ValueError:
                outer = float("nan")
            if np.isnan(outer):
                cells.append("      --")
            else:
                cells.append(f" {outer:>7,.0f}")
                if outer < best[0]:
                    best = (outer, rz, float(z1))
        tag = " (axicon rim)" if rz == 32460.0 else ""
        print(f"  {rz:>10,.0f} |" + "".join(cells) + tag)
    o, rz, z1 = best
    des = close_design(rz, float(np.hypot(x, y).max()), F2_MM, RIM_R_MM, z1)
    _, worst = disk_at_receiver(des, x, y, a.sun_deg)
    print(f"\n  smallest outer-heliostat disk in the feasible wedge: "
          f"{o:,.0f} mm at rim z {rz:,.0f}, F1 {z1:,.0f} "
          f"(tip {z1 - TIP_MM:+,.0f}); field-max disk there {worst:,.0f} mm")
    print(f"  K={-(des.c / des.a) ** 2:.4f}  |R_v|={des.b2 / des.a:,.1f}  "
          f"vertex z={des.z_c + des.a:,.1f}")

    # ---- part 2: occlusion vs F1 at the representative instant ----------
    date = _dt.date.fromisoformat(a.date)
    az, el = sun_position(cfg.site.latitude, cfg.site.longitude,
                          cfg.site.timezone, date.year, date.month, date.day,
                          a.hour)
    print(f"\n\nshading/blocking at {date} {a.hour:.3f} h "
          f"(sun az {az:.2f}, el {el:.2f}), full 645-mirror field,")
    print(f"cassegrain disc 15 m at z {a.rim_z:,.0f}; UNION occlusion "
          f"(lit AND unblocked), field means:\n")

    radius = shading_mod.search_radius_for(
        el, cfg.field.mirror_height_mm, cfg.field.mirror_width_mm)
    neighbours = field_mod.neighbour_pairs(full, radius)

    def occlusion_row(cfg):
        strategy = get_strategy(cfg)
        sols = [strategy.solve(float(x[i]), float(y[i]), az, el, cfg.geometry)
                for i in range(len(full))]
        geoms, aims = shading_mod.build_geometries(full, sols, cfg)
        body = shading_mod.secondary_body(cfg)
        shade, block, sec = shading_mod.shading_blocking(
            geoms, aims, az, el, neighbours, secondary=body)
        occ = shading_mod.occlusion_efficiency(
            geoms, aims, az, el, neighbours, secondary=body)
        return shade, block, sec, occ

    print(f"  {'F1 (mm)':>10} {'tip+':>6} | {'shade':>6} {'disc':>6} "
          f"{'block':>6} {'UNION':>6} | {'worst block':>11} {'n<0.9':>5}")

    object.__setattr__(cfg.optics, "secondary", "cassegrain")
    object.__setattr__(cfg.optics, "n_mirrors", 2)
    object.__setattr__(cfg.geometry, "secondary_rim_height_mm", float(a.rim_z))
    for z1 in [30000, 32000, 33000, 34000, 35000, 36000, 38000, 40000, 43000,
               47000]:
        object.__setattr__(cfg.geometry, "focus_height_mm", float(z1))
        validate_layout(cfg)
        shade, block, sec, occ = occlusion_row(cfg)
        print(f"  {z1:>10,} {(z1 - TIP_MM) / 1000:>5.1f}m | {shade.mean():.4f} "
              f"{sec.mean():.4f} {block.mean():.4f} {occ.mean():.4f} | "
              f"{block.min():>11.4f} {(block < 0.9).sum():>5}")

    # The axicon, same instant, as the comparability anchor.
    axicon = load_config(None)
    shade, block, sec, occ = occlusion_row(axicon)
    print(f"  {'axicon ref':>10} {'':>6} | {shade.mean():.4f} {sec.mean():.4f} "
          f"{block.mean():.4f} {occ.mean():.4f} | {block.min():>11.4f} "
          f"{(block < 0.9).sum():>5}")
    print("\n  shade = neighbours+disc unioned; disc = disc alone; "
          "UNION = lit AND unblocked (the sweep's scalar weight)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
