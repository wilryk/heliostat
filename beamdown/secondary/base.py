"""Interface every secondary-reflector strategy implements."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field, replace


@dataclass(frozen=True)
class HeliostatSolution:
    """Pointing and shape for one heliostat at one instant.

    ``rot_az_deg``/``rot_el_deg`` and ``c3``/``c4``/``c5`` are exactly the five
    per-config multiconfig parameters the Quadoa model expects. The remaining
    fields are diagnostics carried into the summary table -- they cost nothing to
    compute and answer most "why is this heliostat bad?" questions without a
    re-trace.

    The ``extras`` aim-point contract
    --------------------------------
    ``extras`` is free-form *except* for three keys that every strategy MUST
    set::

        extras["aim_x_mm"], extras["aim_y_mm"], extras["aim_z_mm"]

    the world-coordinate point this heliostat is aimed at and focused on.
    :func:`beamdown.shading.build_geometries` reads them to build the outgoing
    beam direction that the *blocking* test is measured along, and it falls back
    to ``(0, 0, secondary_height_mm)`` when a key is missing. That fallback is
    silent and it is wrong for every layout whose aim point is not the secondary
    vertex on the axis -- the axicon pushes the aim point out radially, prime
    focus and Cassegrain put it at ``focus_height_mm`` rather than
    ``secondary_height_mm``. A strategy that forgets these keys therefore does
    not fail, it just reports plausible and incorrect blocking. Set them.

    Beyond those three, ``rot_astig_deg``, ``rad_s_mm``, ``rad_t_mm`` and
    ``focal_dist_s_mm`` are set by every layout that has a meaningful value for
    them; layout-specific diagnostics (the axicon's ``axicon_aoi_deg``) are the
    strategy's own business. Nothing downstream requires the optional ones -- the
    summary table takes whatever it is given.
    """

    rot_az_deg: float
    rot_el_deg: float
    c3: float
    c4: float
    c5: float
    aoi_deg: float = float("nan")
    focal_dist_mm: float = float("nan")
    cosine_efficiency: float = float("nan")
    extras: dict = field(default_factory=dict)


class SecondaryStrategy(ABC):
    """Maps heliostat position + sun direction to pointing and shape."""

    name: str = "base"

    @abstractmethod
    def solve(
        self,
        x_mm: float,
        y_mm: float,
        solar_az_deg: float,
        solar_el_deg: float,
        geometry,
    ) -> HeliostatSolution:
        """Solve one heliostat. ``geometry`` is a :class:`beamdown.config.Geometry`.

        Implementations must populate the ``aim_*_mm`` extras documented on
        :class:`HeliostatSolution`.
        """

    def global_params(self, geometry) -> dict[str, float]:
        """Model-wide ``single_param`` values this layout's ``.optx`` expects.

        Written once per session by
        :meth:`beamdown.session.QuadoaSession.set_global_geometry`. The base set
        is what every beam-down model has; a layout whose model carries extra
        shared geometry (the axicon's ``axi_angle``) overrides and adds to it.

        Kept on the strategy rather than in ``session`` because writing a
        parameter a model does not have is silently ignored by Quadoa -- so a
        prime-focus model would not error on a stray ``axi_angle``, it would just
        leave an unexplained write in the log. Naming only what the layout
        actually has keeps that honest.
        """
        return {
            "sec_height": float(geometry.secondary_height_mm),
            "rec_offset": float(geometry.receiver_offset_mm),
        }

    def describe(self) -> str:
        return self.name


class FlatHeliostats(SecondaryStrategy):
    """Wraps any layout and removes the heliostats' optical power.

    "Flat mirrors" is a property of the RUN, not of the secondary: the paper
    compares focused against flat heliostats *for each* of the three layouts, so
    this composes with all of them rather than being a fourth registry entry.

    What it changes and what it must not
    ------------------------------------
    Exactly three numbers: ``c3``, ``c4``, ``c5`` are forced to ``0.0``. Those are
    the ``z3``/``z4``/``z5`` coefficients of the ONE active ``<form
    type="zernike">`` on ``helio_surf``, whose base ``radius`` is ``inf``. All
    other coefficients on that form (``z0``-``z2``) are literal zeros in the
    ``.optx``, and every other form on the surface (a second zernike, two cosine
    ripples, a biconic) is ``active="false"``. So zeroing these three leaves a
    Zernike form that contributes zero sag and zero slope everywhere -- a plane,
    identical to deactivating the form, without having to edit the model file per
    run. See ``models/heliostat_field_model_mcfg.optx``.

    Pointing is untouched. ``rot_az_deg``/``rot_el_deg`` come straight through, as
    do the diagnostics (``aoi_deg``, ``focal_dist_mm``, ``cosine_efficiency``) and
    the whole ``extras`` dict -- including the ``aim_*_mm`` keys the blocking test
    reads, which would silently change the *shading* answer if they moved.
    ``rad_s_mm``/``rad_t_mm``/``rot_astig_deg`` stay in ``extras`` on purpose: they
    are what the mirror *would* have been bent to, which is the useful diagnostic
    for a flat run, and nothing downstream turns them back into a shape.

    Expect the physics: a flat 5 m x 3 m mirror at 30-120 m throws a spot far
    larger than a focused one, so transmission into the aperture and collected
    energy drop a long way. That is the measurement, not a bug.
    """

    def __init__(self, inner: SecondaryStrategy):
        self.inner = inner
        # The layout's own name, unchanged, so anything comparing against
        # ``cfg.optics.secondary`` still matches. ``describe()`` is where the
        # flatness shows up.
        self.name = inner.name

    def solve(self, x_mm, y_mm, solar_az_deg, solar_el_deg, geometry) -> HeliostatSolution:
        solution = self.inner.solve(x_mm, y_mm, solar_az_deg, solar_el_deg, geometry)
        # ``replace`` rather than rebuilding the dataclass field by field: every
        # field this does not name is carried over by identity, so the non-shape
        # outputs are bit-identical to the focused solve by construction rather
        # than by a copy that could drift.
        return replace(solution, c3=0.0, c4=0.0, c5=0.0)

    def global_params(self, geometry) -> dict[str, float]:
        """Delegated: flatness is per-heliostat shape, not model-wide geometry."""
        return self.inner.global_params(geometry)

    def describe(self) -> str:
        return f"{self.inner.describe()} (flat heliostats, no focusing)"


_REGISTRY: dict[str, type[SecondaryStrategy]] = {}


def register(cls: type[SecondaryStrategy]) -> type[SecondaryStrategy]:
    _REGISTRY[cls.name] = cls
    return cls


def _load_strategies() -> None:
    """Import every strategy module so the registry is complete."""
    from . import axicon, cassegrain, prime_focus  # noqa: F401


def available() -> list[str]:
    """Every registered layout name, for error messages and validation."""
    _load_strategies()
    return sorted(_REGISTRY)


def _layout_and_flat(source) -> tuple[str, bool]:
    """``(layout name, flat?)`` from a config, an ``[optics]`` section, or a name.

    A bare string is a layout name and carries no opinion about flatness, which
    is what keeps ``get_strategy("axicon")`` meaning exactly what it always did
    in the parity and regression tests.
    """
    if isinstance(source, str):
        return source, False
    optics = getattr(source, "optics", source)
    name = getattr(optics, "secondary", None)
    if not isinstance(name, str):
        raise TypeError(
            f"get_strategy() wants a layout name, a Config, or an OpticsSpec; "
            f"got {type(source).__name__}"
        )
    return name, bool(getattr(optics, "flat_mirrors", False))


def get_strategy(source, *, flat: bool | None = None) -> SecondaryStrategy:
    """Build the strategy this run needs, flat mirrors included.

    ``source`` is normally the whole :class:`beamdown.config.Config`, and that is
    deliberate. ``[optics] flat_mirrors`` is a property of the run that must
    compose with all three layouts, so it has to be applied somewhere no call site
    can forget it -- and the only way for this function to apply it is to be given
    the object that carries it. Every caller in the package therefore passes
    ``cfg`` (or ``cfg.optics``) rather than ``cfg.optics.secondary``, and the
    zeroing itself exists in exactly one place, :class:`FlatHeliostats`.

    Passing a bare layout name still works and means "this layout, focused"; the
    parity tests use it to pin the axicon whatever config.toml selects.
    ``flat=`` overrides whatever the config said, in either direction.
    """
    _load_strategies()

    name, flat_from_config = _layout_and_flat(source)
    if flat is None:
        flat = flat_from_config

    try:
        strategy = _REGISTRY[name]()
    except KeyError:
        raise ValueError(
            f"Unknown secondary {name!r}. Available: {sorted(_REGISTRY)}"
        ) from None

    return FlatHeliostats(strategy) if flat else strategy
