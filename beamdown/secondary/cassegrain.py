"""Cassegrain: a hyperboloid secondary relays the F1 bundle to the receiver.

Every heliostat aims at ``F1 = (0, 0, focus_height_mm)``, exactly as for prime
focus. A hyperboloid mirror sits below F1, intercepts the converging bundle before
it gets there, and re-images it onto the receiver -- F1 is the hyperboloid's far
(virtual) focus and the receiver is its near one. Two reflections, like the
axicon.

Why the Python side never sees the conic constants
--------------------------------------------------
Pointing a heliostat only requires knowing where the bundle is *supposed* to
converge, and that is F1 by construction: a hyperboloid's defining property is
that rays headed for one focus leave headed for the other, so the heliostat's job
is unchanged by the relay sitting in the way. The surface itself is built by hand
in Quadoa (conic constant, vertex radius, position) and the trace does the rest.
Nothing here needs those numbers, which is why this file is as short as prime
focus's.

For the same reason there is no shape correction. The axicon needs one because a
cone has no focus and contributes sagittal-only power that the heliostat has to
pre-compensate; a hyperboloid is stigmatic between its two foci and contributes
none.

Consequences handled in the seams outside this file:

* **A circular shadow.** Unlike the axicon's cone, the silhouette is a horizontal
  circle of radius ``axicon_aperture_radius_mm`` at ``secondary_rim_height_mm``,
  so :func:`beamdown.shading.secondary_body` returns a
  :class:`beamdown.shading.SecondaryDisc`. Shading only, never blocking -- the
  same argument as for the cone, and it still holds because F1 sits *above* the
  disc: a beam reaching the disc is the beam arriving.
* **Two reflections.** ``optics.n_mirrors`` should be 2, as for the axicon.

Pointing is the shared single-focus solve -- see
:mod:`beamdown.secondary.shared_focus`.
"""

from __future__ import annotations

from .base import register
from .shared_focus import SharedFocusStrategy


@register
class CassegrainStrategy(SharedFocusStrategy):
    """Hyperboloid secondary relaying the common focus onto the receiver."""

    name = "cassegrain"
