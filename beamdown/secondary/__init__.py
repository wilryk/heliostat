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

Orthogonal to the layout there is a second axis: ``[optics] flat_mirrors``.
Flat heliostats keep their pointing and lose their optical power, and that is a
property of the *run*, not of a layout -- so it is applied once, by
:class:`~beamdown.secondary.base.FlatHeliostats`, wrapping whichever strategy
:func:`~beamdown.secondary.base.get_strategy` was asked for. It is recorded in
the manifest next to the layout name.

Everything downstream -- tracing, storage, metrics, figures, annual energy -- is
written against the :class:`~beamdown.secondary.base.SecondaryStrategy`
interface. Note that the interface is wider than the ABC's single abstract
method: the ``extras["aim_*_mm"]`` keys documented on
:class:`~beamdown.secondary.base.HeliostatSolution` are part of the contract, and
:meth:`~beamdown.secondary.base.SecondaryStrategy.global_params` decides which
model-wide parameters get written into the ``.optx``.
"""

from .base import (
    FlatHeliostats,
    HeliostatSolution,
    SecondaryStrategy,
    available,
    get_strategy,
)

__all__ = [
    "FlatHeliostats",
    "HeliostatSolution",
    "SecondaryStrategy",
    "available",
    "get_strategy",
]
