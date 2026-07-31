"""List the neighbours that occlude one heliostat, ready to build in Quadoa.

Why this exists
---------------
Shading and blocking are computed analytically, outside the ray trace. That is
what makes them revisable without re-tracing, but it also means they have never
been checked against Quadoa on the same geometry. This emits the occluding
rectangles for one (heliostat, timestep) in exactly the parameterisation the
model already uses for a heliostat, so they can be added by hand and traced.

The parameterisation
--------------------
A heliostat in ``heliostat_field_model_mcfg.optx`` is::

    pos  heli_pos          x = posx, y = posy, z = 0          (no rotation)
      pos  heliostat       all zero
        pos  heli_coord_shift   ry = -rot_el, rz = rot_az     (order xyz)
          surf helio_surf       rx = -90, rz = 90             (order xyz)
            aperture rect  s_x = 2500, s_y = 1500             (half-sizes)

so an occluder is that same assembly with four numbers changed: ``posx``,
``posy``, ``rot_az``, ``rot_el``. Nothing else differs -- an occluding neighbour
*is* a heliostat, pointed at the same aim point by the same solver.

Two things to change on the copy: drop the Zernike form (an occluder only needs
to be opaque, and its curvature is irrelevant at these separations -- 1.5 mm of
sag across the aperture), and make it absorbing rather than ``ideal_mirror``, or
its reflections will land somewhere on the receiver and flatter the result.

The rectangle axes were checked against the model's own transform chain rather
than assumed: composing ``Rz(rot_az) Ry(-rot_el)`` with ``Rz(90) Rx(-90)`` sends
local +x to ``normalize(z x n)`` and local +y to ``n x u``, which is the basis
:mod:`beamdown.shading` uses, to twelve decimal places.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Occluder:
    heliostat_id: int
    posx_mm: float
    posy_mm: float
    rot_az_deg: float
    rot_el_deg: float
    distance_m: float
    shades_pct: float          # of the target's aperture, this neighbour alone
    blocks_pct: float
    corners_mm: np.ndarray     # (4, 3) world coordinates, as a cross-check


def occluders_for(cfg, summary, heliostat_id: int, timestep: str,
                  min_effect_pct: float = 0.05):
    """Every neighbour that shades or blocks the given heliostat, with its effect.

    Effects are reported *alone* -- what that one neighbour costs if it were the
    only other mirror in the field. They do not add to the totals, because two
    neighbours can shadow the same patch; the totals come from the union and are
    returned separately.
    """
    from . import field as field_mod
    from . import shading as shading_mod
    from .secondary import get_strategy

    rows = summary[summary.timestep == timestep]
    if not len(rows):
        raise KeyError(f"timestep {timestep!r} is not in this run")
    rows = rows.set_index("heliostat_id")
    if heliostat_id not in rows.index:
        raise KeyError(f"heliostat {heliostat_id} is not in this run")

    az = float(rows.solar_az_deg.iloc[0])
    el = float(rows.solar_el_deg.iloc[0])

    fld = field_mod.load_field(cfg)
    strategy = get_strategy(cfg)
    solutions = [strategy.solve(float(fld.x_mm[i]), float(fld.y_mm[i]),
                                az, el, cfg.geometry) for i in range(len(fld))]
    geoms, aims = shading_mod.build_geometries(fld, solutions, cfg)

    radius = shading_mod.search_radius_for(el, cfg.field.mirror_height_mm,
                                           cfg.field.mirror_width_mm)
    neighbours = field_mod.neighbour_pairs(fld, radius)

    i = int(np.flatnonzero(fld.ids == heliostat_id)[0])
    target = geoms[i]
    pts = target.sample_points(101, 101)
    to_sun = shading_mod.sun_vector(az, el)
    beam = aims[i] - pts

    found = []
    for j in neighbours[i]:
        nb = geoms[int(j)]
        s = float(shading_mod._blocked_mask(pts, to_sun, [nb]).mean()) * 100.0
        b = float(shading_mod._blocked_mask(pts, beam, [nb]).mean()) * 100.0
        if max(s, b) < min_effect_pct:
            continue
        sol = solutions[int(j)]
        corners = np.array([
            nb.centre + su * nb.half_width * nb.u + sv * nb.half_height * nb.v
            for su, sv in ((-1, -1), (1, -1), (1, 1), (-1, 1))
        ])
        found.append(Occluder(
            heliostat_id=int(fld.ids[int(j)]),
            posx_mm=float(fld.x_mm[int(j)]), posy_mm=float(fld.y_mm[int(j)]),
            rot_az_deg=float(sol.rot_az_deg), rot_el_deg=float(sol.rot_el_deg),
            distance_m=float(np.linalg.norm(nb.centre - target.centre) / 1000.0),
            shades_pct=s, blocks_pct=b, corners_mm=corners,
        ))

    found.sort(key=lambda o: -(o.shades_pct + o.blocks_pct))

    # Prime focus has no body over the field, so nothing here is shaded by one.
    cone = shading_mod.secondary_body(cfg)
    by_cone = (cone.occludes(pts, to_sun) if cone is not None
               else np.zeros(len(pts), dtype=bool))
    nbrs = [geoms[int(j)] for j in neighbours[i]]
    shaded = shading_mod._blocked_mask(pts, to_sun, nbrs) | by_cone

    totals = {
        "eta_shade": float(1.0 - shaded.mean()),
        "eta_secondary": float(1.0 - by_cone.mean()),
        "eta_block": float(1.0 - shading_mod._blocked_mask(pts, beam, nbrs).mean()),
        "solar_az_deg": az, "solar_el_deg": el,
        "rot_az_deg": float(solutions[i].rot_az_deg),
        "rot_el_deg": float(solutions[i].rot_el_deg),
        "posx_mm": float(fld.x_mm[i]), "posy_mm": float(fld.y_mm[i]),
        "aim_mm": aims[i],
        "n_neighbours_searched": int(len(neighbours[i])),
        "search_radius_m": radius / 1000.0,
    }
    return found, totals


def describe(found, totals, cfg, heliostat_id: int, timestep: str) -> str:
    L = [
        f"heliostat {heliostat_id} at {timestep}",
        f"  sun            az {totals['solar_az_deg']:8.3f}   el {totals['solar_el_deg']:7.3f} deg",
        f"  target         posx {totals['posx_mm']:10.1f}  posy {totals['posy_mm']:10.1f} mm",
        f"                 rot_az {totals['rot_az_deg']:8.3f}  rot_el {totals['rot_el_deg']:7.3f} deg",
        f"  aim point      {totals['aim_mm'][0]:8.1f} {totals['aim_mm'][1]:9.1f} "
        f"{totals['aim_mm'][2]:9.1f} mm",
        "",
        "  predicted by beamdown (101x101 sampling of the aperture)",
        f"    eta_shade      {totals['eta_shade']:.4f}   "
        f"({100*(1-totals['eta_shade']):.2f}% of the aperture in shadow)",
        f"    eta_secondary  {totals['eta_secondary']:.4f}   (the axicon's share of that)",
        f"    eta_block      {totals['eta_block']:.4f}   "
        f"({100*(1-totals['eta_block']):.2f}% of the beam intercepted)",
        "",
        f"  {len(found)} neighbour(s) contribute, of "
        f"{totals['n_neighbours_searched']} within {totals['search_radius_m']:.1f} m",
        "",
        "  Build each as a copy of the heliostat assembly with these four values.",
        "  Keep aperture rect s_x = 2500, s_y = 1500 (half-sizes); drop the Zernike",
        "  form and make the surface absorbing, not ideal_mirror.",
        "",
        f"    {'id':>5s} {'posx (mm)':>12s} {'posy (mm)':>12s} {'rot_az':>9s} {'rot_el':>8s} "
        f"{'dist m':>7s} {'shades':>7s} {'blocks':>7s}",
    ]
    for o in found:
        L.append(f"    {o.heliostat_id:5d} {o.posx_mm:12.1f} {o.posy_mm:12.1f} "
                 f"{o.rot_az_deg:9.3f} {o.rot_el_deg:8.3f} {o.distance_m:7.2f} "
                 f"{o.shades_pct:6.2f}% {o.blocks_pct:6.2f}%")
    L += [
        "",
        "  'shades'/'blocks' are that neighbour acting alone. They do not sum to the",
        "  totals above -- two neighbours can shadow the same patch, so the totals are",
        "  the union, which is what a trace with all of them present will show.",
        "",
        "  Corner coordinates (world mm), if you would rather place vertices directly:",
    ]
    for o in found:
        L.append(f"    heliostat {o.heliostat_id}")
        for k, c in enumerate(o.corners_mm):
            L.append(f"      corner {k}  {c[0]:11.1f} {c[1]:11.1f} {c[2]:9.1f}")
    return "\n".join(L)
