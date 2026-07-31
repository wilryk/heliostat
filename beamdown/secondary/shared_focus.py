"""Layouts where the whole field aims at one point on the axis.

Both :mod:`beamdown.secondary.prime_focus` and
:mod:`beamdown.secondary.cassegrain` share this solver, and the sharing is the
point: they are optically identical *upstream of F1*. Every heliostat in the
field aims at and focuses on the single point

    F1 = (0, 0, focus_height_mm)

on the tower axis. That is the substantive contrast with the axicon, which has no
single aim point at all -- it computes one per heliostat as a function of radial
position, because a cone has no focus and each heliostat has to hit the part of
the cone that will send it down to the receiver.

What the two layouts then do with the converging bundle differs, but not here:

**prime_focus** puts a downward-facing receiver aperture physically at F1. One
reflection, no secondary mirror, nothing above the field to cast a shadow.

**cassegrain** puts a hyperboloid secondary below F1, using F1 as its far
(virtual) focus, and relays the bundle down to the existing receiver. Two
reflections, and a circular horizontal shadow silhouette over the field.

Neither needs a shape correction on the heliostat. Prime focus has no second
optic to correct for. The Cassegrain's relay is the hyperboloid's job -- its
conic constants are chosen in Quadoa to image F1 onto the receiver, and the
Python side never needs them for pointing. So the astigmatism correction that
``axicon.axicon_shape_correction`` supplies is simply zero for both, and
:func:`beamdown.secondary.mirror.to_quadoa_zernike` is called with the base shape
alone.
"""

from __future__ import annotations

import numpy as np

from .base import HeliostatSolution, SecondaryStrategy
from .mirror import heliostat_orientation, heliostat_shape, to_quadoa_zernike


def focus_point(geometry) -> np.ndarray:
    """``F1``, the single aim point shared by the whole field."""
    height = getattr(geometry, "focus_height_mm", None)
    if height is None:
        raise ValueError(
            "geometry.focus_height_mm is not set. Layouts 'prime_focus' and "
            "'cassegrain' aim the whole field at (0, 0, focus_height_mm), so "
            "that height has no default -- add focus_height_mm to the "
            "[geometry] section of config.toml."
        )
    return np.array([0.0, 0.0, float(height)], dtype=float)


class SharedFocusStrategy(SecondaryStrategy):
    """Aim every heliostat at one on-axis point, with no shape correction.

    Deliberately has no ``name`` and is not registered: it is the shared body,
    and the two concrete layouts below register themselves against it. Their
    :meth:`solve` outputs are identical for identical inputs -- ``tests/verify.py``
    asserts that -- because everything that distinguishes them lives in the
    *other* seams: which body (if any) :func:`beamdown.shading.secondary_body`
    puts over the field, how many reflections ``optics.n_mirrors`` counts, and
    which ``.optx`` gets loaded.
    """

    def solve(self, x_mm, y_mm, solar_az_deg, solar_el_deg, geometry) -> HeliostatSolution:
        aim = focus_point(geometry)
        mirror_pos = np.array([float(x_mm), float(y_mm), 0.0], dtype=float)

        # Unlike the axicon there is no radial direction to be undefined, so a
        # heliostat at the field origin is perfectly well posed: it simply looks
        # straight up the axis. Nothing here divides by the field radius.
        focal_dist = float(np.linalg.norm(aim - mirror_pos))

        (
            rot_az_deg,
            rot_el_deg,
            rot_astig_deg,
            rad_s,
            rad_t,
            aoi_deg,
            _u_mirror,
            _v_mirror,
        ) = heliostat_orientation(aim, mirror_pos, solar_az_deg, solar_el_deg)

        _c0, c3, c4, c5 = heliostat_shape(rot_astig_deg, rad_s, rad_t)
        qc3, qc4, qc5 = to_quadoa_zernike(c3, c4, c5)

        return HeliostatSolution(
            rot_az_deg=float(rot_az_deg),
            rot_el_deg=float(rot_el_deg),
            c3=float(qc3),
            c4=float(qc4),
            c5=float(qc5),
            aoi_deg=float(aoi_deg),
            focal_dist_mm=focal_dist,
            cosine_efficiency=float(np.cos(np.deg2rad(aoi_deg))),
            extras={
                "rot_astig_deg": float(rot_astig_deg),
                "rad_s_mm": float(rad_s),
                "rad_t_mm": float(rad_t),
                # No secondary-induced sagittal shift, so the sagittal focal
                # distance is the plain one. Reported anyway so the summary
                # table has the same columns whichever layout produced it.
                "focal_dist_s_mm": focal_dist,
                # The aim-point contract -- see HeliostatSolution. Constant
                # across the field for these layouts, unlike the axicon.
                "aim_x_mm": float(aim[0]),
                "aim_y_mm": float(aim[1]),
                "aim_z_mm": float(aim[2]),
            },
        )
