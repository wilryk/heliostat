"""Mirror math that every secondary layout shares.

These three functions know nothing about what sits at the top of the tower. They
answer "given an aim point, how is this heliostat pointed and what shape does it
need", which is the same question for the axicon, for prime focus, and for a
Cassegrain relay. Only the *choice of aim point* and the *extra correction on
top* are layout-specific, and those live in the individual strategy modules.

:func:`to_quadoa_zernike` is the one place the sign and normalisation conversion
into Quadoa's normalised Zernike convention happens. It used to be inline in
``axicon.py``; it is here so a layout with no correction term can reuse it
instead of copying three lines of sign flips that are easy to get subtly wrong.
The arithmetic is written out exactly as it was, term by term, so the axicon's
``c3``/``c4``/``c5`` are bit-identical to before the factoring --
``tests/verify.py`` stage "axicon solve() regression" pins that against a literal
table captured from the pre-refactor code.
"""

from __future__ import annotations

import numpy as np

# Normalisation constants for Quadoa's Zernike convention. Computed with
# ``np.sqrt`` rather than ``math.sqrt`` purely so the factored helper reproduces
# the previous inline arithmetic bit for bit.
_SQRT3 = np.sqrt(3.0)
_SQRT6 = np.sqrt(6.0)


def heliostat_orientation(receiver_pos, mirror_pos, solar_az_deg, solar_el_deg):
    """Aim a flat mirror to send the sun to ``receiver_pos``.

    Returns pointing angles, the astigmatic rotation, the sagittal/tangential
    focal radii, the angle of incidence, and the mirror's local basis vectors.

    The mirror normal is the *bisector* of the direction to the sun and the
    direction to the aim point, both taken as unit vectors from the mirror
    centre. Anything reconstructing the normal from ``rot_az``/``rot_el`` (the
    shading code, the GUI, the new-strategy tests) has to use exactly this
    convention or it will disagree with the ray trace.
    """
    receiver_pos = np.asarray(receiver_pos, dtype=float)
    mirror_pos = np.asarray(mirror_pos, dtype=float)

    solar_az = np.deg2rad(solar_az_deg)
    solar_el = np.deg2rad(solar_el_deg)

    # Unit vector from mirror toward the sun. Azimuth is compass bearing, hence
    # the pi/2 - az conversion into standard math convention.
    to_sun = np.array([
        np.cos(solar_el) * np.cos(np.pi / 2 - solar_az),
        np.cos(solar_el) * np.sin(np.pi / 2 - solar_az),
        np.sin(solar_el),
    ])
    to_sun /= np.linalg.norm(to_sun)

    to_target = receiver_pos - mirror_pos
    focal_length = np.linalg.norm(to_target)
    to_target = to_target / focal_length

    normal = to_sun + to_target
    normal /= np.linalg.norm(normal)

    rot_el = np.arcsin(normal[2])
    rot_az = np.arctan2(normal[1], normal[0])

    aoi = 0.5 * np.arccos(np.clip(np.dot(to_sun, to_target), -1.0, 1.0))

    up = np.array([0.0, 0.0, 1.0])
    u = np.cross(up, normal)
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    v /= np.linalg.norm(v)

    rot_astig = np.arctan2(np.dot(to_sun, v), np.dot(to_sun, u))

    radius = focal_length * 2.0
    rad_s = radius * np.cos(aoi)
    rad_t = radius / np.cos(aoi)

    return (
        np.rad2deg(rot_az),
        np.rad2deg(rot_el),
        np.rad2deg(rot_astig),
        rad_s,
        rad_t,
        np.rad2deg(aoi),
        u,
        v,
    )


def heliostat_shape(rot_astig_deg, rad_s, rad_t):
    """Curvature radii -> Zernike-like shape coefficients."""
    rot_astig = np.deg2rad(rot_astig_deg)
    curv_t = 1.0 / rad_t
    curv_s = 1.0 / rad_s

    c0 = 0.125 * (curv_s + curv_t)
    c3 = 0.25 * (curv_t - curv_s) * np.sin(2 * rot_astig)
    c4 = 0.125 * (curv_t + curv_s)
    c5 = 0.25 * (curv_t - curv_s) * np.cos(2 * rot_astig)
    return c0, c3, c4, c5


def to_quadoa_zernike(c3, c4, c5, c3_corr=0.0, c4_corr=0.0, c5_corr=0.0):
    """Sum base shape and secondary correction in Quadoa's Zernike convention.

    The sign flips are Quadoa's orientation conventions, not physics. A layout
    whose secondary contributes no extra astigmatism (prime focus; Cassegrain,
    whose relay is handled by the hyperboloid itself rather than by bending the
    heliostat) passes the correction terms as zero and gets the base shape
    converted alone.
    """
    qc3 = c3 / _SQRT6 + (-c3_corr / _SQRT6)
    qc4 = -c4 / _SQRT3 + (-c4_corr / _SQRT3)
    qc5 = -c5 / _SQRT6 + (c5_corr / _SQRT6)
    return qc3, qc4, qc5
