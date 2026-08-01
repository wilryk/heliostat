"""Secondary-reflector strategies.

A strategy turns (heliostat position, sun direction, geometry) into the pointing
and shape coefficients that get written to the Quadoa model. Everything
downstream -- tracing, storage, metrics, figures, annual energy -- is written
against the :class:`~beamdown.secondary.base.SecondaryStrategy` interface and
does not know or care which secondary is in play.

Three layouts are registered:

``axicon``       conical secondary; the aim point is computed per heliostat from
                 its radial position, because a cone has no focus.
``prime_focus``  no secondary; the receiver sits at the field's common focus F1.
``cassegrain``   hyperboloid secondary relaying F1 down to the receiver.

The last two share one aim point for the whole field and one solver -- see
:mod:`beamdown.secondary.shared_focus`. The strategy name is recorded in each
run's manifest so results from different layouts stay directly comparable.

Orthogonal to the layout there is a second axis, the heliostat FIGURE, with
three points on it. Each is a property of the *run*, not of a layout, so each is
applied once by a wrapper around whichever strategy
:func:`~beamdown.secondary.base.get_strategy` was asked for, and each is recorded
in the manifest next to the layout name:

``focused``      the default: the mirror is re-figured every timestep for that
                 instant's angle of incidence and slant range.
``flat_mirrors`` :class:`~beamdown.secondary.base.FlatHeliostats` -- pointing
                 kept, optical power dropped entirely.
``fixed_shapes`` :class:`~beamdown.secondary.base.FixedShapeHeliostats` --
                 pointing kept, figure frozen per heliostat to a table, which is
                 what a mirror ground once actually does.

The last two both write ``c3``/``c4``/``c5``, so ``get_strategy`` refuses them
together rather than ordering them.

Everything downstream -- tracing, storage, metrics, figures, annual energy -- is
written against the :class:`~beamdown.secondary.base.SecondaryStrategy`
interface. Note that the interface is wider than the ABC's single abstract
method: the ``extras["aim_*_mm"]`` keys documented on
:class:`~beamdown.secondary.base.HeliostatSolution` are part of the contract, and
:meth:`~beamdown.secondary.base.SecondaryStrategy.global_params` decides which
model-wide parameters get written into the ``.optx``.
"""

from .base import (
    FixedShapeError,
    FixedShapeHeliostats,
    FlatHeliostats,
    HeliostatSolution,
    SecondaryStrategy,
    available,
    get_strategy,
    load_fixed_shapes,
)

__all__ = [
    "FixedShapeError",
    "FixedShapeHeliostats",
    "FlatHeliostats",
    "HeliostatSolution",
    "SecondaryStrategy",
    "available",
    "get_strategy",
    "load_fixed_shapes",
]
