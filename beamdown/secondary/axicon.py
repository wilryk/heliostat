"""Axicon secondary: heliostat pointing and astigmatic shape correction.

Ported from ``heliostat/heliostat_shape_solve.py``, then corrected in one place.
Pointing is unchanged and still matches the original to machine precision;
``tests/test_axicon_parity.py`` asserts that, and also that the shape matches
once the foreshortening factor in :func:`axicon_shape_correction` is forced back
to 1, so the departure is exactly one deliberate factor and nothing else. Also
changed: naming, structure, and the fact that intermediate quantities (angle of
incidence, focal distance, cosine efficiency) are now returned instead of
discarded, so they can land in the summary table.

Physics summary
---------------
A heliostat is aimed by the bisector of the sun vector and the vector to its
aim point, which for the beam-down layout is not the receiver itself but the
point where the ray meets the axicon cone (:func:`receiver_correction`). An
off-axis mirror used at incidence angle ``aoi`` is astigmatic: the sagittal and
tangential focal radii differ by ``cos(aoi)`` and ``1/cos(aoi)``, which is what
``c3``/``c4``/``c5`` correct.

The axicon adds a second astigmatism: it has optical power in the sagittal
direction only (its tangential focal length is effectively infinite), so each
heliostat gets an extra sagittal-only correction on top of the standard one.
The two corrections are summed in the Quadoa normalised Zernike convention.

Sagittal for the *axicon* is not sagittal for the *mirror*, though, and the
``cos(aoi)`` relation only holds for the mirror's own sagittal axis. The extra
correction therefore carries a foreshortening factor -- see
:func:`axicon_shape_correction`. It matters most where the two sagittal axes are
furthest apart, which correlates with, but is not caused by, being near the
tower.
"""

from __future__ import annotations

import numpy as np

from .base import HeliostatSolution, SecondaryStrategy, register

# Layout-agnostic mirror math, shared with the prime-focus and Cassegrain
# strategies. Re-exported here because ``beamdown.compare`` and the parity tests
# reach for them as ``axicon.heliostat_orientation`` / ``axicon.heliostat_shape``.
from .mirror import (  # noqa: F401
    heliostat_orientation,
    heliostat_shape,
    to_quadoa_zernike,
)


def axicon_shape_correction(u, v, sagittal_vector, focal_dist, focal_dist_s, aoi_rad,
                            foreshorten=None):
    """Extra sagittal-only correction contributed by the axicon.

    The axicon has no tangential power, represented here by a very large
    tangential focal length.

    The correction is a cylinder, but its axis is the *axicon's* sagittal
    direction, which is not the mirror's own -- so the plain sagittal relation
    ``rad = 2 f cos(aoi)`` does not apply to it unmodified. Projecting that
    direction onto the mirror yields two things: the angle the cylinder must run
    (used below, and correct as it stands), and the projection's length, which
    says how much of the direction lies in the mirror surface at all. A step
    across the mirror buys only that fraction of travel along the direction
    being corrected, and wavefront error grows as the square of distance, so the
    required curvature carries the squared length. It is 1 exactly when the
    axicon's sagittal direction lies in the mirror plane, where the plain
    relation is already right; ignoring it makes the correction 1/L**2 too
    strong, which is a factor of ~2 for the inner field.

    ``foreshorten`` overrides that factor. Passing ``1.0`` restores the earlier
    behaviour, which is how the parity tests confirm this is the only
    intentional departure from ``legacy/heliostat/heliostat_shape_solve.py``.
    """
    d_power = 1.0 / focal_dist_s - 1.0 / focal_dist
    f_tangential = 1e20
    f_sagittal = 1.0 / d_power

    sag_u = np.dot(sagittal_vector, u)
    sag_v = np.dot(sagittal_vector, v)
    if foreshorten is None:
        foreshorten = (sag_u ** 2 + sag_v ** 2) / np.dot(sagittal_vector, sagittal_vector)

    rad_s = f_sagittal * 2.0 * np.cos(aoi_rad) / foreshorten
    rad_t = f_tangential * 2.0

    rot_astig = np.arctan2(sag_v, sag_u)
    return heliostat_shape(np.rad2deg(rot_astig), rad_s, rad_t)


def receiver_correction(mirror_radial_position, axicon_height, receiver_offset, axicon_angle_deg):
    """Where a ray from a heliostat meets the axicon cone, and the aim offset.

    Solves the intersection of the mirror-to-receiver line with the cone
    surface, in the radial/height plane containing the heliostat.
    """
    alpha = np.deg2rad(axicon_angle_deg)

    x_r = -receiver_offset * np.sin(2 * alpha)
    y_r = receiver_offset * np.cos(2 * alpha)

    x_m = mirror_radial_position
    y_m = -axicon_height

    slope = (y_r - y_m) / (x_r - x_m)
    x_a = (slope * x_m - y_m) / (slope - np.tan(alpha))
    y_a = np.tan(alpha) * x_a

    return x_r, y_r, x_a, y_a


@register
class AxiconStrategy(SecondaryStrategy):
    """Conical (axicon) secondary reflector."""

    name = "axicon"

    def global_params(self, geometry) -> dict[str, float]:
        """The axicon model also carries the cone's half-angle."""
        params = super().global_params(geometry)
        params["axi_angle"] = float(geometry.axicon_angle_deg)
        return params

    def solve(self, x_mm, y_mm, solar_az_deg, solar_el_deg, geometry) -> HeliostatSolution:
        secondary_height = geometry.secondary_height_mm
        receiver_height = geometry.receiver_height_mm
        axicon_angle_deg = geometry.axicon_angle_deg

        drop = secondary_height - receiver_height
        field_radius = np.hypot(x_mm, y_mm)
        if field_radius == 0.0:
            raise ValueError("Heliostat at the field origin has no defined radial direction")

        (
            receiver_radial_offset,
            receiver_height_offset,
            axicon_radial_intersection,
            axicon_height_intersection,
        ) = receiver_correction(
            mirror_radial_position=field_radius,
            axicon_height=secondary_height,
            receiver_offset=drop,
            axicon_angle_deg=axicon_angle_deg,
        )

        # Push the aim point out along this heliostat's own radial direction.
        aim = np.array([
            x_mm / field_radius * receiver_radial_offset,
            y_mm / field_radius * receiver_radial_offset,
            secondary_height + receiver_height_offset,
        ], dtype=float)
        mirror_pos = np.array([x_mm, y_mm, 0.0], dtype=float)

        alpha = np.deg2rad(axicon_angle_deg)
        cone_dist = np.hypot(axicon_radial_intersection, axicon_height_intersection)
        s_prime = -np.hypot(drop + axicon_height_intersection, axicon_radial_intersection)
        radius_axicon = cone_dist / np.tan(alpha)

        axicon_aoi = np.rad2deg(
            np.arctan2(axicon_radial_intersection, drop + axicon_height_intersection) + alpha
        )
        s = 1.0 / (2.0 * np.cos(np.deg2rad(axicon_aoi)) / radius_axicon - 1.0 / s_prime)

        to_aim = aim - mirror_pos
        focal_dist = float(np.linalg.norm(to_aim))
        focal_dist_s = focal_dist + (s + s_prime)

        (
            rot_az_deg,
            rot_el_deg,
            rot_astig_deg,
            rad_s,
            rad_t,
            aoi_deg,
            u_mirror,
            v_mirror,
        ) = heliostat_orientation(aim, mirror_pos, solar_az_deg, solar_el_deg)

        c0, c3, c4, c5 = heliostat_shape(rot_astig_deg, rad_s, rad_t)

        # Sagittal direction of the beam, in the horizontal plane.
        focus_dir = to_aim / focal_dist
        focus_xy = np.array([focus_dir[0], focus_dir[1], 0.0])
        focus_xy /= np.linalg.norm(focus_xy)
        sagittal_vector = np.cross(focus_xy, np.array([0.0, 0.0, 1.0]))

        c0_c, c3_c, c4_c, c5_c = axicon_shape_correction(
            u_mirror, v_mirror, sagittal_vector, focal_dist, focal_dist_s, np.deg2rad(aoi_deg)
        )

        # Convert to Quadoa's normalised Zernike convention and sum the two
        # contributions.
        qc3, qc4, qc5 = to_quadoa_zernike(c3, c4, c5, c3_c, c4_c, c5_c)

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
                "axicon_aoi_deg": float(axicon_aoi),
                "focal_dist_s_mm": float(focal_dist_s),
                "aim_x_mm": float(aim[0]),
                "aim_y_mm": float(aim[1]),
                "aim_z_mm": float(aim[2]),
            },
        )
