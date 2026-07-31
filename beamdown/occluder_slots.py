"""Choose which neighbours occupy the model's occluder slots, per heliostat.

Why slots
---------
Shading and blocking have always been scalars applied after the trace. That is
right for total power and wrong for spot *shape*: a shadow removes a particular
patch of the mirror, and that patch images to a particular part of the spot. To
see the real spot -- especially defocused -- the occluders have to be in the ray
path.

Quadoa geometry cannot be created through the API, so the model carries a fixed
number of occluder surfaces and each trace moves them into place.

There is no "off" position -- see :func:`_rank`. A sequential trace visits every
listed surface, so an unused slot cannot be removed, only moved, and moving it
somewhere silly costs rays. Every slot therefore holds a real neighbour, and the
ones past the genuine occluders hold neighbours that occlude nothing.

How many slots
--------------
Measured over the field at five representative timesteps: shading needs at most
6 neighbours (only at the lowest sun; 0-2 is typical) and blocking at most 2,
never more. :data:`N_SHADE` and :data:`N_BLOCK` leave headroom on both, and
:func:`plan_field` reports any heliostat that would have overflowed rather than
silently dropping an occluder. Every heliostat in this field has at least 24
neighbours within the search radius, so the slots always fill.

Ranking
-------
Slots go first to the neighbours that block the most, measured the same way
:mod:`beamdown.shading` measures them, then to non-occluding neighbours ordered
by how far ahead of the mirror they sit.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

N_SHADE = 10
N_BLOCK = 4

# A neighbour that hides less than this of the aperture is not worth a slot.
# Well below the 1/375 quantisation of the analytic sampling grid.
MIN_FRACTION = 5e-4

PARK_MM = 200000.0

# Where a synthetic filler goes when no real neighbour lies ahead: far enough
# forward that the plane is well conditioned, far enough sideways that its
# rectangle cannot touch a bundle a few metres across.
FILLER_AHEAD_MM = 30000.0
FILLER_LATERAL_MM = 25000.0


@dataclass
class Slot:
    heliostat_id: int
    x_mm: float
    y_mm: float
    rot_az_deg: float
    rot_el_deg: float
    fraction: float


# Height of the axicon-shadow plane in the model. Must match the literal z on
# the ax0 <pos> element; it is not a parameter there.
AXICON_PLANE_Z_MM = 13500.0


@dataclass
class HeliostatPlan:
    row: int
    heliostat_id: int
    shading: list = field(default_factory=list)
    blocking: list = field(default_factory=list)
    dropped_shading: int = 0
    dropped_blocking: int = 0
    dropped_fraction: float = 0.0
    # Where the secondary's shadow falls on the ax0 plane. One value per
    # timestep -- it depends on the sun alone -- carried here so a worker gets
    # it with the rest of the heliostat's setup.
    #
    # ``None`` means this layout has no secondary body and so no shadow plane to
    # position: prime focus. ``session.set_occluders`` skips the write, which is
    # also what happens anyway because a prime-focus model has no ``ax0_x``
    # parameter -- the ``None`` just makes the intent explicit rather than
    # relying on the probe.
    axicon_xy: tuple | None = (0.0, 0.0)


def axicon_shadow_centre(body, to_sun, plane_z_mm: float = AXICON_PLANE_Z_MM):
    """Centre of the secondary's shadow circle on a horizontal plane.

    Takes any secondary body exposing ``rim_height_mm``, which is
    :class:`~beamdown.shading.SecondaryCone` and
    :class:`~beamdown.shading.SecondaryDisc`.

    The shadow *is* a circle of the rim's radius, exactly, whenever the sun is
    high enough -- because a horizontal circle projected along a fixed direction
    onto a horizontal plane is a congruent circle, merely translated. For the
    Cassegrain's disc that is exact for every sun direction, full stop, since the
    rim circle is the whole silhouette.

    For the axicon's cone it is an approximation, because the cone's shadow is the
    convex hull of that circle and the *vertex's* shadow point, and the point
    falls inside the circle only while

        (rim height - vertex height) / tan(elevation)  <=  rim radius

    which holds above about 20 degrees. Below that the shadow is already being
    thrown clear of the field -- 158 m at 9.7 degrees against a 90 m field
    radius -- so the circle is exact everywhere it matters. Checked against the
    exact ray-cone test over all 44 timesteps of a full sweep: zero disagreeing
    samples except at 16 degrees, where 3 of 40,635 differ.
    """
    ts = np.asarray(to_sun, float)
    if ts[2] <= 1e-6:
        return (float(PARK_MM), float(PARK_MM))
    drop = (body.rim_height_mm - plane_z_mm) / ts[2]
    return (float(-ts[0] * drop), float(-ts[1] * drop))


def _rank(geoms, ids, candidates, points, direction, limit):
    """Fill every slot with a real neighbour, the strongest occluders first.

    Every slot is filled even when fewer neighbours actually occlude, because
    there is no safe way to switch a slot off. The model's sequence is fixed, so
    a listed surface must be intersected on every trace; "parking" it far away
    does not remove it, it just moves where the intersection happens, and moving
    an infinite plane a thousand kilometres out cost 34% of the rays to
    precision alone. Distances of 1e5 to 1e7 mm were no better -- erratic
    between 0.72 and 1.01 of the unobstructed count.

    A real neighbour that happens to occlude nothing is the honest filler: its
    plane sits a few metres away where the intersection is well conditioned, and
    its rectangle is somewhere the rays are not, so it blocks exactly nothing.
    Duplicates are harmless too -- two coincident obscurations stop the same
    rays as one -- which covers a heliostat with fewer neighbours than slots.
    """
    from .shading import _blocked_mask

    # A filler must be AHEAD of the heliostat along the ray, not merely nearby.
    # Quadoa marches the sequence in order, so every listed plane is crossed
    # before helio_surf; a plane sitting behind the mirror makes the ray double
    # back and it is lost. A heliostat on the up-sun rim has all its neighbours
    # down-sun, and filling its slots with the nearest ones cost 28% of its rays
    # while the analytic answer was a clean 1.000.
    origin = np.mean(points, axis=0)
    unit = np.asarray(direction, float)
    unit = unit.mean(axis=0) if unit.ndim == 2 else unit
    unit = unit / np.linalg.norm(unit)

    occluding, ahead = [], []
    for j in candidates:
        g = geoms[int(j)]
        f = float(_blocked_mask(points, direction, [g]).mean())
        reach = float((g.centre - origin) @ unit)
        if f > MIN_FRACTION:
            occluding.append((f, reach, int(j)))
        elif reach > 0.0:
            ahead.append((f, reach, int(j)))

    occluding.sort(key=lambda t: -t[0])
    ahead.sort(key=lambda t: -t[1])        # furthest ahead first: best separated
    order = occluding + ahead

    kept = order[:limit]
    dropped = [t for t in occluding[limit:]]
    slots = [Slot(heliostat_id=int(ids[j]),
                  x_mm=float(geoms[j].centre[0]), y_mm=float(geoms[j].centre[1]),
                  rot_az_deg=_az_of(geoms[j]), rot_el_deg=_el_of(geoms[j]),
                  fraction=f)
             for f, _r, j in kept]

    # A heliostat on the up-sun rim of the field can have no neighbour ahead of
    # it at all, and then there is nothing real to put in the spare slots. Place
    # a synthetic one: ahead along the ray so the plane is crossed in the right
    # order, and displaced far enough sideways that its 5 x 3 m rectangle cannot
    # reach a bundle only metres wide. Its orientation is copied from the
    # heliostat, which guarantees the plane is not parallel to the ray.
    while len(slots) < limit:
        centre = origin + unit * FILLER_AHEAD_MM
        side = np.cross(unit, np.array([0.0, 0.0, 1.0]))
        norm = np.linalg.norm(side)
        side = (side / norm) if norm > 1e-9 else np.array([1.0, 0.0, 0.0])
        centre = centre + side * FILLER_LATERAL_MM
        slots.append(Slot(heliostat_id=-1,
                          x_mm=float(centre[0]), y_mm=float(centre[1]),
                          rot_az_deg=slots[0].rot_az_deg if slots else 0.0,
                          rot_el_deg=slots[0].rot_el_deg if slots else 45.0,
                          fraction=0.0))
    return slots, len(dropped), sum(t[0] for t in dropped)


def _az_of(geom) -> float:
    n = geom.normal
    return float(np.degrees(np.arctan2(n[1], n[0])))


def _el_of(geom) -> float:
    n = geom.normal
    return float(np.degrees(np.arcsin(np.clip(n[2], -1.0, 1.0))))


def plan_field(geoms, aims, ids, neighbours, to_sun, nu: int = 25, nv: int = 15,
               n_shade: int = N_SHADE, n_block: int = N_BLOCK, body=None,
               *, has_secondary: bool = True):
    """One :class:`HeliostatPlan` per heliostat, for one instant.

    ``body`` is whatever :func:`beamdown.shading.secondary_body` returned, and
    ``has_secondary`` says whether this layout has one *at all*. The two are
    distinct: ``body=None, has_secondary=True`` is "there is a shadow plane in the
    model but the sun is not usefully placed, so park it", which is the existing
    behaviour for a model whose ax0 slot is not being traced;
    ``has_secondary=False`` is prime focus, where there is no shadow plane and the
    plan carries no ``axicon_xy`` for anyone to write.
    """
    if not has_secondary:
        ax_xy = None
    elif body is not None:
        ax_xy = axicon_shadow_centre(body, to_sun)
    else:
        ax_xy = (float(PARK_MM), float(PARK_MM))

    plans = []
    for i, geom in enumerate(geoms):
        pts = geom.sample_points(nu, nv)
        cand = neighbours[i]
        s, ds, fs = _rank(geoms, ids, cand, pts, to_sun, n_shade)
        b, db, fb = _rank(geoms, ids, cand, pts, aims[i] - pts, n_block)
        plans.append(HeliostatPlan(row=i, heliostat_id=int(ids[i]),
                                   shading=s, blocking=b,
                                   dropped_shading=ds, dropped_blocking=db,
                                   dropped_fraction=float(fs + fb),
                                   axicon_xy=ax_xy))
    return plans


def overflow_report(plans) -> dict:
    """How often the slot count was not enough, and by how much."""
    over = [p for p in plans if p.dropped_shading or p.dropped_blocking]
    return {
        "heliostats": len(plans),
        "overflowed": len(over),
        "max_shading_used": max((len(p.shading) for p in plans), default=0),
        "max_blocking_used": max((len(p.blocking) for p in plans), default=0),
        "worst_dropped_fraction": max((p.dropped_fraction for p in plans), default=0.0),
    }
