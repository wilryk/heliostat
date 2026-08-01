#!/usr/bin/env python
"""Design calculator for a CASSEGRAIN (hyperboloid) beam-down secondary.

Licence-free: numpy + pandas, plus :mod:`beamdown.design_eval` for the geometry
itself. It imports nothing from ``quadoa`` and touches no package state beyond
that one pure-math module, reads only ``config.toml`` and the field position
file, and writes no files -- so it is safe to run while a trace is in progress.
Nothing here is written back into ``config.toml``; the lines the package would
need are printed as a suggestion for the user to paste.

The conic closure itself (:class:`~beamdown.design_eval.Cassegrain`,
:func:`~beamdown.design_eval.close_design`,
:func:`~beamdown.design_eval.rays_to_f1`) lived in this file first and now lives
in ``beamdown/design_eval.py``, so the GUI's Design tab and this calculator
cannot drift apart. This script's own behaviour and printed output are
unchanged; what follows describes the math wherever it is written down.


The geometry
------------
Everything is a surface of revolution about the tower axis ``z``, with the
heliostat field on the ground plane ``z = 0`` and the final receiver aperture
(F2) on the axis at ``z = receiver_height_mm`` (7,000 mm for the current
config: ``secondary_height_mm + receiver_offset_mm = 27000 - 20000``).

The secondary is a **hyperboloid of revolution** with its two foci on the axis:

    F1 = (0, 0, z1)     the upper / "prime" focus, above the reflector
    F2 = (0, 0, z2)     the lower focus, z2 = 7,000 mm, the real receiver

Every heliostat aims *and focuses* at F1, exactly as it would in a plain
prime-focus tower. The hyperboloid is inserted **below** F1, intercepting the
converging bundle before it reaches F1, and relays it to F2. This is a textbook
Cassegrain, turned upside down: the field is the primary, F1 is the prime focus
that is never actually formed, and the hyperboloid's conjugate-focus property
(``|d1 - d2| = 2a`` and the reflection law that a ray aimed at one focus leaves
aimed at the other) does the relay with zero added aberration *for the ideal
point-source case*. Because the relay is exact for every ray aimed at F1
regardless of where on the sheet it lands, the secondary contributes no
spherical aberration of its own -- unlike the axicon it replaces, which has
power in the sagittal direction only and needs a per-heliostat astigmatism
correction.

Two-sheet / branch choice
~~~~~~~~~~~~~~~~~~~~~~~~~
A hyperboloid of revolution has two sheets, one wrapped around each focus, and
picking the wrong one silently produces a design that "closes" arithmetically
but reflects light nowhere near F2. With the axis quantities

    c      = half the focus separation = (z1 - z2) / 2
    z_c    = midpoint of the foci      = (z1 + z2) / 2      (the centre)
    a      = semi-transverse axis, 0 < a < c
    b^2    = c^2 - a^2

the two sheets have vertices at ``z_c + a`` (the sheet enclosing F1, opening
upward) and ``z_c - a`` (the sheet enclosing F2, opening downward). The
reflector we want is **the sheet surrounding F1**, i.e. vertex at ``z_c + a``,
which sits *below* F1 and *above* the centre. It is a bowl opening upward toward
F1, and it is silvered on its **underside** -- the side the field sees, on which
it is convex, exactly like a Cassegrain secondary. This module does not assume
that: :func:`reflect_check` intersects the real quadric, reflects about the real
surface normal, and asserts the reflected ray hits F2, which is the only proof
that the branch and every sign are right.

Conic / sag convention
~~~~~~~~~~~~~~~~~~~~~~
Written as a rotationally symmetric conic sag measured from the vertex, with the
local axis pointing *out of the dish* (up, toward F1), so the sag is positive:

    z(r) = r^2 / ( R_v * (1 + sqrt(1 - (1 + K) r^2 / R_v^2)) )

    R_v = b^2 / a        vertex radius of curvature
    K   = -e^2           conic constant,  e = c / a  > 1  (hyperbola)

This is checked in code against the algebraic form ``a*(sqrt(1 + r^2/b^2) - 1)``.


The design closure -- the point of this script
----------------------------------------------
The receiver F2 (7,000 mm) and the rim radius (15,000 mm, matching the axicon
aperture) are fixed by the existing plant. That leaves exactly two free
parameters: the prime-focus height ``z1`` and the hyperboloid's ``2a``. They are
closed in three steps:

1. The user picks the **rim height** ``z_r`` -- the height of the 15 m-radius rim
   circle, i.e. how high the dish's edge hangs. ``--rim-height-mm``. This is the
   one genuine engineering degree of freedom (tower/truss height, wind load,
   how much of the sky the dish blocks, how deep the dish is).

   **The chosen design is z_r = 32,460 mm**, which is the axicon's own rim
   height (``27000 + 15000*tan(20 deg)``). Picking it means the hyperboloid's
   rim circle sits exactly where the cone's rim circle sat: same aperture, same
   height, therefore the same silhouette against the sky and the same shading
   footprint on the field, so the shading model needs no re-derivation and the
   two secondaries are compared on optics alone. That is the default here; the
   trade table still spans 24-34 m so the cost of the choice is visible.

2. **Aperture fill** pins ``z1``. The ray from the *outermost* heliostat, at the
   true maximum field radius ``R_f``, must pass exactly through the rim point
   ``(15000, z_r)``. That ray runs from ``(R_f, 0, 0)`` to ``(0, 0, z1)``; on it,
   ``z = t*z1`` and ``r = R_f*(1 - t)``, so at ``z = z_r`` the radius is
   ``r = R_f * (1 - z_r/z1)``. Setting that equal to the rim radius ``r_rim``:

       R_f (1 - z_r/z1) = r_rim
       1 - z_r/z1       = r_rim / R_f
       z_r / z1         = (R_f - r_rim) / R_f
       z1               = z_r * R_f / (R_f - r_rim)

   *Rationale*: the spec is "same aperture as the axicon", 15 m radius. That is
   only physically true if the field's bundle actually fills 15 m at the rim. If
   ``z1`` were lower, the bundle would be narrower than the dish at the rim
   height and the outer annulus of the secondary would be dead metal; if ``z1``
   were higher, the outermost rings of heliostats would throw their beams past
   the rim and be wasted mirror -- 645 heliostats' worth of capital, spilling.
   Aperture fill is the unique choice where the last heliostat *grazes* the rim,
   so no mirror is wasted and no reflector is wasted. Everything inside the
   outermost ring then lands strictly inside the rim automatically, because
   ``r`` at the rim plane scales linearly with the heliostat's own radius.

3. **The rim point lying on the surface** pins ``a``. For any point P on a
   hyperboloid, ``| |P-F1| - |P-F2| | = 2a``. Applying it at the rim point
   ``P = (15000, 0, z_r)`` gives ``2a`` directly, and with ``c`` already known
   from ``z1`` and ``z2`` the whole conic follows.

``--f1-height-mm`` overrides step 2: the user supplies ``z1`` and aperture fill
is skipped, but ``a`` still comes from step 3, so the rim stays at
``(15000, z_r)`` by construction and the design still closes.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

# The geometry, single-sourced. ``RIM_RADIUS_MM`` is the axicon's own aperture;
# ``DEFAULT_RIM_HEIGHT_MM`` (32,460 = 27000 + 15000*tan(20 deg)) is the axicon's
# own rim height, so the replacement occupies the same rim circle in space as the
# cone it removes -- same aperture AND same height, hence the same shading
# footprint on the field. ``TRADE_RIM_HEIGHTS_MM`` is the 24-34 m trade sweep
# with that chosen height spliced in so it appears in the table.
from beamdown.design_eval import (  # noqa: E402
    DEFAULT_RIM_HEIGHT_MM,
    RIM_RADIUS_MM,
    TRADE_RIM_HEIGHTS_MM,
    Cassegrain,
    close_design,
    rays_to_f1,
)


# ----------------------------------------------------------------------------
# inputs: read config.toml and the field file, read-only
# ----------------------------------------------------------------------------

def read_config(path: Path) -> dict:
    """The handful of values needed, straight out of ``config.toml``.

    Read with ``tomllib`` rather than through ``beamdown.config`` so this script
    stays importable with a trace running and cannot touch package state.
    """
    try:
        import tomllib
    except ModuleNotFoundError:                          # py<3.11
        import tomli as tomllib                          # type: ignore

    with open(path, "rb") as fh:
        raw = tomllib.load(fh)

    geom = raw["geometry"]
    secondary_height = float(geom["secondary_height_mm"])
    receiver_offset = float(geom["receiver_offset_mm"])
    return {
        "secondary_height_mm": secondary_height,
        "receiver_offset_mm": receiver_offset,
        "receiver_height_mm": secondary_height + receiver_offset,
        "axicon_angle_deg": float(geom["axicon_angle_deg"]),
        "axicon_aperture_radius_mm": float(geom.get("axicon_aperture_radius_mm", 15000.0)),
        "positions_file": str(raw["field"]["positions_file"]),
    }


def read_field_mm(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Heliostat centres in **mm**. Source file is in metres (columns 'X (m)').

    Mirrors the unit handling of ``beamdown/field.py`` (metres in, mm out)
    without importing it.
    """
    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)

    def pick(aliases: tuple[str, ...], which: str) -> str:
        low = {str(c).strip().lower(): c for c in df.columns}
        for alias in aliases:
            if alias in low:
                return low[alias]
        raise KeyError(f"no {which} column in {path.name}; has {list(df.columns)}")

    xcol = pick(("x (m)", "x_s (m)", "x", "x (mm)"), "x")
    ycol = pick(("y (m)", "y_s (m)", "y", "y (mm)"), "y")
    x = df[xcol].to_numpy(float)
    y = df[ycol].to_numpy(float)
    scale = 1.0 if "mm" in str(xcol).lower() else 1000.0
    return x * scale, y * scale


# ----------------------------------------------------------------------------
# checks and diagnostics
# ----------------------------------------------------------------------------

def coverage(x_mm, y_mm, des: Cassegrain) -> dict:
    """Does every real heliostat's beam land on the dish, inside the rim?

    Two independent tests. The cheap one crosses the rim *plane*; the exact one
    solves the quadric. They must agree on who misses.
    """
    o, d = rays_to_f1(x_mm, y_mm, des)
    R = np.hypot(o[:, 0], o[:, 1])

    # rim-plane crossing: on the ray z = t*z1 and r = R*(1-t), so at z = rim_z
    r_plane = R * (1.0 - des.rim_z / des.z1)

    t, hit, ok, roots = des.intersect(o, d)
    r_hit = np.hypot(hit[:, 0], hit[:, 1])

    tol = 1e-6 * des.rim_r                       # 0.015 mm: the grazing ray
    inside = ok & (r_hit <= des.rim_r + tol)
    miss_plane = r_plane > des.rim_r + tol
    graze = ok & (np.abs(r_hit - des.rim_r) <= tol)

    # reflect and measure the miss distance at F2
    dr = des.reflect(d[ok], hit[ok])
    v = des.F2[None, :] - hit[ok]
    t2 = np.sum(v * dr, axis=1)
    miss_f2 = np.linalg.norm(v - t2[:, None] * dr, axis=1)

    cos_i = np.abs(np.sum(d[ok] * des.normal(hit[ok]), axis=1))
    aoi = np.degrees(np.arccos(np.clip(cos_i, -1.0, 1.0)))

    path = t[ok] + np.linalg.norm(des.F2[None, :] - hit[ok], axis=1)

    return {
        "n": int(R.size),
        "R_mm": R,
        "r_plane_mm": r_plane,
        "r_hit_mm": r_hit,
        "z_hit_mm": hit[:, 2],
        "ok": ok,
        "inside": inside,
        "n_miss_plane": int(miss_plane.sum()),
        "n_no_intersect": int((~ok).sum()),
        "n_outside_rim": int((ok & ~inside).sum()),
        "n_graze": int(graze.sum()),
        "aoi_deg": aoi,
        "path_mm": path,
        "miss_f2_mm": miss_f2,
        "roots": roots,
        "t_hit": t,
        "hit": hit,
        "agree": bool(np.array_equal(miss_plane, ~inside)),
    }


def reflect_check(x_mm, y_mm, des: Cassegrain, tol_mm: float = 1e-6):
    """Exact end-to-end optical proof for a handful of heliostats.

    Intersect heliostat->F1 with the quadric exactly, reflect about the exact
    analytic normal, and measure the perpendicular distance from F2 to the
    reflected ray. Anything but ~0 means the branch, the sign of ``a``, or the
    normal is wrong. Returns ``(DataFrame, all_passed)``.
    """
    o, d = rays_to_f1(x_mm, y_mm, des)
    t, hit, ok, roots = des.intersect(o, d)
    if not ok.all():
        raise RuntimeError("a sampled heliostat's ray does not meet the used sheet")

    dr = des.reflect(d, hit)
    v = des.F2[None, :] - hit
    t2 = np.sum(v * dr, axis=1)
    miss = np.linalg.norm(v - t2[:, None] * dr, axis=1)

    n = des.normal(hit)
    aoi = np.degrees(np.arccos(np.clip(np.abs(np.sum(d * n, axis=1)), -1.0, 1.0)))
    r_hit = np.hypot(hit[:, 0], hit[:, 1])

    # the defining hyperbola relation, at the actual intercept
    two_a = np.abs(np.linalg.norm(hit - des.F1[None, :], axis=1)
                   - np.linalg.norm(hit - des.F2[None, :], axis=1))

    df = pd.DataFrame({
        "R_helio_mm": np.hypot(o[:, 0], o[:, 1]),
        "r_hit_mm": r_hit,
        "z_hit_mm": hit[:, 2],
        "|f(P)|": np.abs(des.implicit(hit)),
        "2a_at_hit": two_a,
        "aoi_deg": aoi,
        "refl_dz": dr[:, 2],
        "t_to_F2_mm": t2,
        "miss_F2_mm": miss,
        "pass": (miss < tol_mm) & (dr[:, 2] < 0.0) & (t2 > 0.0),
    })
    df["2a_err_mm"] = np.abs(two_a - 2.0 * des.a)
    return df, bool(df["pass"].all())


def sample_indices(R: np.ndarray) -> np.ndarray:
    """Origin-nearest, three quartiles, outermost -- by radius."""
    order = np.argsort(R, kind="stable")
    n = order.size
    picks = [0, n // 4, n // 2, (3 * n) // 4, n - 1]
    return order[np.array(sorted(set(picks)))]


# ----------------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------------

def trade_table(rim_heights, field_radius_mm, z2, x_mm, y_mm,
                rim_r=RIM_RADIUS_MM) -> pd.DataFrame:
    rows = []
    for rz in np.asarray(rim_heights, float):
        des = close_design(float(rz), field_radius_mm, z2, rim_r)
        cov = coverage(x_mm, y_mm, des)
        rows.append({
            "rim_z_m": rz / 1000.0,
            "F1_z1_m": des.z1 / 1000.0,
            "a_mm": des.a,
            "c_mm": des.c,
            "b2_mm2": des.b2,
            "e": des.e,
            "K": des.K,
            "|R_v|_mm": abs(des.R_v),
            "vertex_z_m": des.vertex_z / 1000.0,
            "sag_mm": des.sag_mm,
            "f/D_prime": des.z1 / (2.0 * field_radius_mm),
            "aoi_min_deg": cov["aoi_deg"].min(),
            "aoi_max_deg": cov["aoi_deg"].max(),
            "path_ratio": cov["path_mm"].max() / cov["path_mm"].min(),
            "n_miss": cov["n_miss_plane"],
            "max_miss_F2_mm": cov["miss_f2_mm"].max(),
        })
    return pd.DataFrame(rows)


def print_design(des: Cassegrain, x_mm, y_mm, cfg: dict) -> bool:
    R = np.hypot(x_mm, y_mm)
    cov = coverage(x_mm, y_mm, des)
    idx = sample_indices(R)
    chk, passed = reflect_check(x_mm[idx], y_mm[idx], des)

    def mm_m(v):
        return f"{v:>14,.3f} mm  ({v/1000.0:>9.4f} m)"

    print()
    print("=" * 78)
    print(f"CASSEGRAIN SECONDARY DESIGN   rim height z_r = {des.rim_z:,.1f} mm "
          f"({des.rim_z/1000:.3f} m)")
    print("=" * 78)
    print(f"  closure           : {'APERTURE FILL (z1 derived)' if des.aperture_filled else 'F1 HEIGHT OVERRIDE (z1 given)'}")
    print(f"  field radius R_f  : {mm_m(des.field_radius_mm)}   (max over "
          f"{R.size} real heliostats)")
    print(f"  rim radius        : {mm_m(des.rim_r)}   (same aperture as the axicon)")
    print(f"  receiver F2       : {mm_m(des.z2)}   on axis")

    print()
    print("-- FOCI ---------------------------------------------------------------")
    print(f"  F1 (prime focus)  : (0, 0, {des.z1:,.3f}) mm   = {des.z1/1000:.4f} m")
    print("      F1 is ALSO the receiver location for the separate 'prime focus'")
    print("      layout: remove the secondary entirely and put the receiver here,")
    print("      and the same heliostat aiming/focusing solution is unchanged.")
    print(f"  F2 (real receiver): (0, 0, {des.z2:,.3f}) mm   = {des.z2/1000:.4f} m")
    print(f"  focus separation  : {mm_m(des.z1 - des.z2)}")

    print()
    print("-- CONIC --------------------------------------------------------------")
    print(f"  c  (half focus sep)     : {mm_m(des.c)}")
    print(f"  a  (semi-transverse)    : {mm_m(des.a)}")
    print(f"  2a                      : {mm_m(2.0 * des.a)}")
    print(f"  b^2 = c^2 - a^2         : {des.b2:>18,.3f} mm^2   (b = {des.b:,.3f} mm)")
    print(f"  eccentricity e = c/a    : {des.e:>18.9f}")
    print(f"  conic constant K = -e^2 : {des.K:>18.9f}")
    print(f"  |R_v| = b^2/a           : {abs(des.R_v):>18,.6f} mm   "
          f"({abs(des.R_v)/1000:.6f} m)")
    print("  sag convention: z(r) = r^2 / ( R_v (1 + sqrt(1 - (1+K) r^2 / R_v^2)) ),")
    print("    measured from the vertex with the local axis pointing OUT of the dish")
    print("    (upward, toward F1); with that orientation R_v is POSITIVE.")
    print("  QUADOA: enter |R_v| and K above. The SIGN of the radius depends on the")
    print("    orientation you give that surface's local z axis, which you set by")
    print("    hand -- CHECK that the dish is concave toward F1 (bowl opening up,")
    print("    reflective underside facing the field); if the sag comes out with the")
    print("    vertex above the rim instead of below, flip the radius sign.")

    print()
    print("-- SHAPE --------------------------------------------------------------")
    print(f"  centre (focus midpoint) : {mm_m(des.z_c)}")
    print(f"  vertex z = z_c + a      : {mm_m(des.vertex_z)}   (used sheet, encloses F1)")
    print(f"  other sheet vertex      : {mm_m(des.z_c - des.a)}   (encloses F2 -- NOT used)")
    print(f"  rim z (echo)            : {mm_m(des.rim_z)}")
    print(f"  sag vertex->rim (depth) : {mm_m(des.sag_mm)}")
    print(f"  vertex is below F1 by   : {mm_m(des.z1 - des.vertex_z)}")
    sag_alg = float(des.sag_at(des.rim_r))
    sag_con = float(des.sag_at_conic(des.rim_r))
    print(f"  sag cross-check         : algebraic {sag_alg:,.9f} vs conic formula "
          f"{sag_con:,.9f} mm  (|d| = {abs(sag_alg-sag_con):.3e})")

    print()
    print("-- DIAGNOSTICS ACROSS THE FIELD ---------------------------------------")
    ok = cov["ok"]
    print(f"  intercept radius on dish : {cov['r_hit_mm'][ok].min():>12,.3f} .. "
          f"{cov['r_hit_mm'][ok].max():>12,.3f} mm")
    print(f"  intercept height on dish : {cov['z_hit_mm'][ok].min():>12,.3f} .. "
          f"{cov['z_hit_mm'][ok].max():>12,.3f} mm")
    print(f"  incidence angle on dish  : {cov['aoi_deg'].min():>12.6f} .. "
          f"{cov['aoi_deg'].max():>12.6f} deg")
    print(f"  path helio->dish->F2     : {cov['path_mm'].min():>12,.1f} .. "
          f"{cov['path_mm'].max():>12,.1f} mm")
    print(f"  path length ratio max/min: {cov['path_mm'].max()/cov['path_mm'].min():>12.6f}")
    d1 = np.linalg.norm(cov["hit"] - des.F1[None, :], axis=1)
    d2 = np.linalg.norm(cov["hit"] - des.F2[None, :], axis=1)
    ratio = d1[ok] / d2[ok]
    print(f"  |P-F1|/|P-F2| at hit     : {ratio.min():>12.6f} .. {ratio.max():>12.6f}")

    print()
    tol_outer = 1e-6 * des.field_radius_mm
    print("-- COVERAGE CHECK over all real heliostats ----------------------------")
    print(f"  heliostats                                    : {cov['n']}")
    print(f"  rays crossing rim plane at r > {des.rim_r:,.0f} mm (MISS) : "
          f"{cov['n_miss_plane']}")
    print(f"  rays with no intersection on the used sheet    : {cov['n_no_intersect']}")
    print(f"  rays intersecting outside the rim radius       : {cov['n_outside_rim']}")
    n_outer = int(np.sum(np.abs(R - des.field_radius_mm) <= tol_outer))
    print(f"  rays grazing the rim (r_hit == rim, 1e-6 rel)  : {cov['n_graze']}"
          f"   (the outermost ring holds {n_outer} heliostats at exactly R_f, so by "
          "construction every one of them grazes)")
    print(f"  plane test and exact quadric test agree        : {cov['agree']}")
    print(f"  rim-plane crossing radius r_i                  : "
          f"{cov['r_plane_mm'].min():,.3f} .. {cov['r_plane_mm'].max():,.3f} mm")
    print(f"  worst |F2 miss| over ALL {cov['n']} heliostats       : "
          f"{cov['miss_f2_mm'].max():.3e} mm")
    if cov["n_miss_plane"] or cov["n_outside_rim"] or cov["n_no_intersect"]:
        worst = np.argsort(cov["r_plane_mm"])[::-1][:5]
        print("  worst-offending heliostats (index, R_mm, r_plane_mm, r_hit_mm):")
        for i in worst:
            print(f"    {i:>4d}  R={cov['R_mm'][i]:>10,.1f}  "
                  f"r_plane={cov['r_plane_mm'][i]:>10,.1f}  "
                  f"r_hit={cov['r_hit_mm'][i]:>10,.1f}")

    print()
    print("-- EXACT END-TO-END REFLECTION CHECK (5 heliostats) -------------------")
    print("   ray helio->F1  x  quadric (exact)  ->  reflect about exact normal")
    print("   ->  perpendicular distance from F2 to the reflected ray, tol 1e-6 mm")
    with pd.option_context("display.width", 200, "display.max_columns", 40,
                           "display.float_format", lambda v: f"{v:,.6g}"):
        print(chk.to_string(index=False))
    print(f"   ALL PASS: {passed}   (2a reproduced at every intercept to "
          f"{chk['2a_err_mm'].max():.3e} mm)")
    # Make the branch discrimination explicit rather than implied: show where the
    # rejected root of the same quadratic lands, and that it is on the other sheet.
    o_s, d_s = rays_to_f1(x_mm[idx], y_mm[idx], des)
    t_s, hit_s, _, roots_s = des.intersect(o_s, d_s)
    other = np.where(np.isclose(roots_s[:, 0], t_s), roots_s[:, 1], roots_s[:, 0])
    p_other = o_s + other[:, None] * d_s
    r_other = np.hypot(p_other[:, 0], p_other[:, 1])
    why = ["F2 sheet (z<z_c)" if z < des.z_c else "F1 sheet, far outside rim"
           for z in p_other[:, 2]]
    print("   branch discrimination -- the rejected root of the same quadratic:")
    for z, r, w in zip(p_other[:, 2], r_other, why):
        print(f"     z = {z:>13,.0f} mm  r = {r:>13,.0f} mm  ->  {w}")
    print(f"   accepted roots z = [" + ", ".join(f"{z:,.0f}" for z in hit_s[:, 2])
          + f"] mm, every one inside [vertex {des.vertex_z:,.0f}, rim "
          f"{des.rim_z:,.0f}] and r <= {des.rim_r:,.0f} mm: the used sheet is the "
          "one enclosing F1, first surface met going up.")

    print()
    print("-- config.toml SUGGESTION (not written -- paste by hand) ---------------")
    print("  these are the only keys beamdown/config.py's Geometry/OpticsSpec")
    print("  dataclasses actually accept for this layout -- Config.load_config()")
    print("  builds Geometry(**raw['geometry']), so any OTHER key under [geometry]")
    print("  (e.g. a cassegrain_* key) is a stray dataclass kwarg and load_config()")
    print("  raises TypeError immediately, it does not warn-and-ignore.")
    print("  [geometry]")
    print(f"  # secondary_height_mm and receiver_offset_mm are UNCHANGED "
          f"({cfg['secondary_height_mm']:.1f} / {cfg['receiver_offset_mm']:.1f}); "
          f"receiver stays at {cfg['receiver_height_mm']:.1f} mm (F2).")
    print(f"  focus_height_mm          = {des.z1:.6f}   # F1, prime focus: what "
          "every heliostat aims and focuses at (z1 of this design)")
    print(f"  secondary_rim_height_mm  = {des.rim_z:.1f}   # rim height of the "
          "hyperboloid, for the shading silhouette")
    print(f"  # axicon_aperture_radius_mm ({cfg['axicon_aperture_radius_mm']:.1f}) "
          "already doubles as the cassegrain disc radius -- no change needed.")
    print("  [optics]")
    print('  secondary                = "cassegrain"   # strategy in '
          "beamdown/secondary/cassegrain.py")
    print("  n_mirrors                = 2   # already correct: heliostat + "
          "secondary, same as axicon")
    print(f"  # existing axicon keys (axicon_angle_deg = {cfg['axicon_angle_deg']:.1f}) "
          "become unused but are harmless to leave in place.")
    print()
    print("-- HYPERBOLOID CONSTANTS for Quadoa (enter these when building the ---")
    print("-- surface there -- the Python side above never consumes them) --------")
    print(f"  # vertex z              = {des.vertex_z:.9f} mm")
    print(f"  # |R_v| (vertex radius) = {des.R_v:.9f} mm   (sign per Quadoa local axis)")
    print(f"  # K (conic constant)    = {des.K:.9f}")
    print(f"  # rim radius            = {des.rim_r:.1f} mm")
    print(f"  # dish sag (vertex->rim): {des.sag_mm:.6f} mm")
    print(f"  # a (semi-transverse)   = {des.a:.9f} mm")
    print(f"  # c (half focus sep)    = {des.c:.9f} mm")
    print(f"  # centre z (foci mid)   = {des.z_c:.9f} mm")
    print("=" * 78)
    return passed


# ----------------------------------------------------------------------------

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Design calculator for a Cassegrain (hyperboloid) beam-down secondary.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--rim-height-mm", type=float, default=DEFAULT_RIM_HEIGHT_MM,
                   help="height of the 15 m-radius rim circle (default 32460, the "
                        "axicon's own rim height -- the chosen design)")
    p.add_argument("--f1-height-mm", type=float, default=None,
                   help="override z1 directly; skips the aperture-fill condition "
                        "(a is still set by the rim point)")
    p.add_argument("--rim-radius-mm", type=float, default=RIM_RADIUS_MM,
                   help="aperture rim radius (default 15000, matching the axicon)")
    p.add_argument("--config", type=Path, default=REPO / "config.toml")
    p.add_argument("--no-trade-table", action="store_true")
    args = p.parse_args(argv)

    cfg = read_config(args.config)
    field_path = REPO / cfg["positions_file"]
    x_mm, y_mm = read_field_mm(field_path)
    R = np.hypot(x_mm, y_mm)
    R_f = float(R.max())
    z2 = cfg["receiver_height_mm"]

    axicon_rim = (cfg["secondary_height_mm"]
                  + cfg["axicon_aperture_radius_mm"]
                  * np.tan(np.deg2rad(cfg["axicon_angle_deg"])))

    print("=" * 78)
    print("CASSEGRAIN SECONDARY DESIGN CALCULATOR")
    print("=" * 78)
    print(f"  config              : {args.config}")
    print(f"  field file          : {field_path.name}")
    print(f"  heliostats          : {R.size}")
    print(f"  field radius R_f    : {R_f:,.3f} mm = {R_f/1000:.6f} m   "
          f"(min {R.min():,.1f} mm = {R.min()/1000:.3f} m)")
    print(f"  receiver F2 (fixed) : {z2:,.1f} mm = {z2/1000:.3f} m "
          f"(= secondary_height_mm {cfg['secondary_height_mm']:,.0f} "
          f"+ receiver_offset_mm {cfg['receiver_offset_mm']:,.0f})")
    print(f"  rim radius (fixed)  : {args.rim_radius_mm:,.1f} mm = "
          f"{args.rim_radius_mm/1000:.1f} m")
    print(f"  axicon it replaces  : tip z = {cfg['secondary_height_mm']:,.0f} mm, "
          f"half-angle {cfg['axicon_angle_deg']:.1f} deg, "
          f"rim z = {axicon_rim:,.1f} mm")
    print(f"  aperture fill       : z1 = z_r * R_f / (R_f - r_rim) = z_r * "
          f"{R_f/(R_f-args.rim_radius_mm):.9f}")

    if not args.no_trade_table:
        print()
        print("-- TRADE TABLE over rim height (aperture-fill closure) ----------------")
        tt = trade_table(TRADE_RIM_HEIGHTS_MM, R_f, z2, x_mm, y_mm, args.rim_radius_mm)
        with pd.option_context("display.width", 250, "display.max_columns", 40):
            print(tt.to_string(index=False, float_format=lambda v: f"{v:,.4f}"))
        print("  rim_z_m   rim height (the choice).  F1_z1_m  prime focus / prime-focus receiver.")
        print("  sag_mm    dish depth.  n_miss  heliostats whose beam clears the rim (want 0).")
        print("  path_ratio  max/min of heliostat->dish->F2 optical path across the field.")

    des = close_design(args.rim_height_mm, R_f, z2, args.rim_radius_mm, args.f1_height_mm)
    passed = print_design(des, x_mm, y_mm, cfg)
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
