"""Prime focus: no secondary at all, receiver sitting at F1.

The simplest layout and the natural baseline for the other two, because it is the
only one that costs a single reflection. The receiver is a horizontal aperture
facing *down* placed physically at ``F1 = (0, 0, focus_height_mm)``, and every
heliostat aims straight at it.

Consequences a reader should hold onto, all of which are handled in seams outside
this file:

* **Nothing above the field.** There is no secondary body, so
  :func:`beamdown.shading.secondary_body` returns ``None`` for this layout and
  ``eta_secondary`` stays 1.0 everywhere. ``shading_blocking`` already tolerates
  ``secondary=None``.
* **No ax0 slot to plan.** With no body there is no shadow circle to place, so
  ``plan_field(has_secondary=False)`` leaves ``axicon_xy`` as ``None`` and
  ``session.set_occluders`` skips the write. (It would skip it anyway, because a
  prime-focus model has no ``ax0_x`` parameter for ``has_axicon_slot`` to find --
  the ``None`` makes the intent explicit rather than relying on that probe.)
* **Its own ``.optx``.** ``models/heliostat_field_prime_focus.optx``, built from
  the axicon model by ``scripts/build_prime_focus_model.py``: sequence 3 becomes
  ``sun -> helio_surf -> prime_focus`` and the detector's height becomes the
  ``pf_height`` parameter :meth:`PrimeFocusStrategy.global_params` writes. Point
  ``[trace] model_file`` at it; ``analysis_seq = 3`` is unchanged.
* **One reflection.** ``optics.n_mirrors`` should be 1, not 2; the receiver
  aperture is not a mirror. :func:`beamdown.config.load_config` warns if it
  disagrees rather than correcting it, because throughput is applied when stored
  runs are *read*.

Pointing is the shared single-focus solve -- see
:mod:`beamdown.secondary.shared_focus`.
"""

from __future__ import annotations

from .base import register
from .shared_focus import SharedFocusStrategy


@register
class PrimeFocusStrategy(SharedFocusStrategy):
    """Receiver at the field's common focus; no secondary mirror."""

    name = "prime_focus"

    def global_params(self, geometry) -> dict[str, float]:
        """Also positions the detector, so it cannot drift off the aim point.

        ``pf_height`` is the ``single_param`` that
        ``models/heliostat_field_prime_focus.optx`` puts on the ``prime_focus``
        surface's ``z`` -- the same mechanism ``sec_height`` uses on
        ``secondary`` and ``rec_offset`` on ``receiver``. Writing it here is what
        makes the model's detector and this layout's aim point the *same number*
        rather than two numbers that happen to agree: :func:`shared_focus.focus_point`
        aims every heliostat at ``(0, 0, focus_height_mm)``, and this puts the
        plane the rays are counted on at exactly that height.

        Get this wrong and nothing fails. The trace runs, the spot is round, the
        centroid is on axis, and the spot size is simply the wrong one -- it is
        the defocus at whatever height the file happened to store. That failure
        mode is why the height is a parameter at all; see
        ``scripts/build_prime_focus_model.py``.

        ``sec_height`` and ``rec_offset`` are inherited and still written. The
        prime-focus model is a copy of the axicon one, so it still carries both
        parameters and both surfaces -- ``secondary`` and ``receiver`` are simply
        no longer in sequence 3's path. Writing them is therefore harmless and
        keeps the two models' parameter sets identical, which is what lets
        ``analysis_seq = 3`` mean the same thing in both.

        Against a model that has no ``pf_height`` -- the stock
        ``heliostat_field_model_mcfg.optx``, whose ``prime_focus`` surface still
        has a literal ``z = 27000`` -- Quadoa silently ignores the write and the
        detector stays where it was. ``scripts/verify_prime_focus_model.py``
        reads the parameter back for exactly that reason; run it once when a
        licence seat is free, before trusting a prime-focus sweep.
        """
        params = super().global_params(geometry)
        # focus_point() raises with a full explanation when this is None, and it
        # is reached on the very first solve() -- so a None here would fail the
        # run either way. Deferring to it keeps one error message rather than
        # two that could drift apart.
        height = getattr(geometry, "focus_height_mm", None)
        if height is not None:
            params["pf_height"] = float(height)
        return params
