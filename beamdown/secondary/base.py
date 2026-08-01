"""Interface every secondary-reflector strategy implements."""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, replace
from pathlib import Path


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


class FixedShapeError(RuntimeError):
    """A heliostat has no row in the fixed-figure table.

    Hard rather than a fallback to the solved shape: a run that silently mixed
    ground mirrors with re-figured ones would report an annual energy between the
    two and look entirely plausible. Absent means the table is wrong for this
    field, which is a build-the-table problem, not a run-time one.
    """


# Position -> figure lookups are keyed to 1e-3 mm (one nanometre of field
# position), which is far finer than any real placement difference and far
# coarser than float noise between the CSV's repr round-trip and the field
# file's own parse. Both the loader and the lookup go through this function so
# the two roundings cannot drift apart.
_SHAPE_KEY_DP = 3


def shape_key(x_mm: float, y_mm: float) -> tuple[float, float]:
    """The ``(x, y)`` key a fixed-figure table is indexed by."""
    return (round(float(x_mm), _SHAPE_KEY_DP), round(float(y_mm), _SHAPE_KEY_DP))


def load_fixed_shapes(path) -> dict[tuple[float, float], tuple[float, float, float]]:
    """Read a ``heliostat,x_mm,y_mm,c3,c4,c5`` CSV into a position -> figure map.

    ``#`` lines anywhere in the file are metadata (which layout, which weighting,
    when it was generated) and are skipped, so the table stays self-describing
    without a sidecar. Written by ``scripts/build_fixed_shapes.py``::

        # fixed mirror figure, mode=mean_cos, layout=pf36000
        heliostat,x_mm,y_mm,c3,c4,c5
        0,22300.0,-60000.0,1.234e-05,-8.5e-06,3.1e-07

    Repeated positions are accepted only when they carry the SAME figure. The
    field file genuinely contains two coincident pairs (heliostats 144 = 192 and
    241 = 289 are byte-identical positions -- see CLAUDE.md), and a figure that is
    a function of position must come out equal for both, so that case is normal.
    Two rows for one position that DISAGREE are an error rather than last-wins:
    the generator contradicted itself and picking either one would hide it.
    """
    import pandas as pd

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"fixed-shape table {path} does not exist. Build it with "
            f"scripts/build_fixed_shapes.py before the sweep -- a missing table "
            f"is not the same as no fixed figure, so it is not treated as one."
        )

    # float_precision="round_trip": pandas' default C parser is fast, not
    # correctly rounded, and returns a double up to an ULP away from the one the
    # generator wrote. That is harmless for the positions (the key rounds to
    # 1e-3 mm) but it means a table does NOT round-trip, so a figure read back
    # would not be the figure computed. The exact parse costs nothing on a
    # 645-row file.
    table_df = pd.read_csv(path, comment="#", float_precision="round_trip")
    wanted = ["x_mm", "y_mm", "c3", "c4", "c5"]
    missing = [c for c in wanted if c not in table_df.columns]
    if missing:
        raise ValueError(
            f"fixed-shape table {path} is missing column(s) {missing}. Expected "
            f"header 'heliostat,x_mm,y_mm,c3,c4,c5' (leading '#' lines are "
            f"comments); got {list(table_df.columns)}"
        )

    # Coerced explicitly so a column that is not numeric is named here, with the
    # value that broke it. The failure this catches is real and silent-looking:
    # under numpy 2 a scalar's repr is "np.float64(30000.0)", so a generator that
    # formats positions with !r writes a whole column of strings and pandas types
    # it as object -- and the first sign of it otherwise would be a float()
    # traceback with no filename in it.
    numeric = table_df[wanted].apply(pd.to_numeric, errors="coerce")
    bad = numeric.isna().any(axis=1)
    if bad.any():
        row = table_df[wanted][bad].iloc[0].to_dict()
        raise ValueError(
            f"fixed-shape table {path} has non-numeric values; first bad row "
            f"(data row {int(bad.idxmax())}): {row}. Columns must be plain "
            f"numbers -- numpy 2 reprs like 'np.float64(30000.0)' are the usual "
            f"cause, so format with float(v) rather than repr(v) when writing."
        )

    table: dict[tuple[float, float], tuple[float, float, float]] = {}
    for x, y, c3, c4, c5 in numeric.itertuples(index=False, name=None):
        key = shape_key(x, y)
        figure = (float(c3), float(c4), float(c5))
        # Agreement to 1e-9 relative, not bit-equality: a generator that
        # accumulates a year-weighted mean over 4,876 timesteps reaches the two
        # coincident heliostats in different orders, so their coefficients differ
        # in the last few ULPs. That is float noise, not two different mirrors --
        # the measured spread on the real tables is ~1e-12 relative. Anything
        # larger means the two rows really do disagree. First row wins; at this
        # tolerance the choice cannot matter.
        if key in table and not all(
                math.isclose(a, b, rel_tol=1e-9, abs_tol=0.0)
                for a, b in zip(table[key], figure)):
            raise ValueError(
                f"fixed-shape table {path} gives position ({key[0]}, {key[1]}) mm "
                f"two DIFFERENT figures, {table[key]} and {figure}. Coincident "
                f"heliostats are expected in this field (144 = 192, 241 = 289) "
                f"but they must be ground the same, since the figure is a "
                f"function of position."
            )
        table.setdefault(key, figure)

    if not table:
        raise ValueError(f"fixed-shape table {path} has no data rows")
    return table


class FixedShapeHeliostats(SecondaryStrategy):
    """Wraps any layout and freezes each heliostat's FIGURE to a ground shape.

    The focused solve re-figures every mirror at every timestep -- the astigmatic
    Zernike that exactly corrects that instant's angle of incidence and slant
    range. Real glass is ground once. This models that: pointing still tracks the
    sun exactly as the focused run does, but ``c3``/``c4``/``c5`` come from a
    table computed once per heliostat and never move again. The result sits
    between the focused idealisation and :class:`FlatHeliostats`, which is the
    point of measuring it.

    What it changes and what it must not
    ------------------------------------
    Exactly the same three numbers ``FlatHeliostats`` touches, and for the same
    reason: they are the ``z3``/``z4``/``z5`` coefficients of the one active
    ``<form type="zernike">`` on ``helio_surf``. Everything else -- pointing, the
    diagnostics, and the whole ``extras`` dict including the ``aim_*_mm`` keys the
    blocking test is measured along -- comes through by ``replace``, so it is
    bit-identical to the focused solve by construction.

    ``extras`` therefore still carries ``rad_s_mm``/``rad_t_mm``/``rot_astig_deg``:
    what the mirror *would* have been bent to at this instant, which is the
    diagnostic that says how far this timestep is from the frozen figure. Nothing
    downstream turns those back into a shape.

    A heliostat with no row raises :class:`FixedShapeError`. See its docstring for
    why that is not a fallback.
    """

    def __init__(self, inner: SecondaryStrategy, table: dict, source: str = ""):
        self.inner = inner
        self.table = table
        # Only for describe() and the error message -- the run's identity in the
        # manifest is written by store.write_manifest from the config, not read
        # back off the strategy.
        self.source = str(source)
        # The layout's own name, unchanged, so anything comparing against
        # ``cfg.optics.secondary`` still matches; describe() carries the figure.
        self.name = inner.name

    def solve(self, x_mm, y_mm, solar_az_deg, solar_el_deg, geometry) -> HeliostatSolution:
        # Looked up BEFORE the inner solve so a table that does not cover this
        # field fails on the first heliostat of the first timestep, not after a
        # trace has been paid for.
        key = shape_key(x_mm, y_mm)
        figure = self.table.get(key)
        if figure is None:
            raise FixedShapeError(
                f"heliostat at ({key[0]}, {key[1]}) mm has no row in the "
                f"fixed-shape table {self.source or '<in memory>'} "
                f"({len(self.table)} rows). The table must cover the WHOLE field, "
                f"not just this run's traced subset -- the sweep solves every "
                f"heliostat to build the occlusion geometry. Falling back to the "
                f"solved shape would mix ground and re-figured mirrors in one "
                f"annual number."
            )
        solution = self.inner.solve(x_mm, y_mm, solar_az_deg, solar_el_deg, geometry)
        c3, c4, c5 = figure
        return replace(solution, c3=c3, c4=c4, c5=c5)

    def global_params(self, geometry) -> dict[str, float]:
        """Delegated: the figure is per-heliostat shape, not model-wide geometry."""
        return self.inner.global_params(geometry)

    def describe(self) -> str:
        return (f"{self.inner.describe()} (fixed figure, {len(self.table)} "
                f"heliostats from {self.source or '<in memory>'})")


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


def _layout_flat_shapes(source) -> tuple[str, bool, str]:
    """``(layout, flat?, fixed-shape path)`` from a config, an ``[optics]``, or a name.

    A bare string is a layout name and carries no opinion about flatness or about
    a fixed figure, which is what keeps ``get_strategy("axicon")`` meaning exactly
    what it always did in the parity and regression tests.
    """
    if isinstance(source, str):
        return source, False, ""
    optics = getattr(source, "optics", source)
    name = getattr(optics, "secondary", None)
    if not isinstance(name, str):
        raise TypeError(
            f"get_strategy() wants a layout name, a Config, or an OpticsSpec; "
            f"got {type(source).__name__}"
        )
    return (name,
            bool(getattr(optics, "flat_mirrors", False)),
            str(getattr(optics, "fixed_shapes", "") or ""))


def _fixed_shapes_path(source, name: str) -> Path:
    """Resolve a fixed-shape path the way the config resolves every other path.

    A whole ``Config`` knows its ``repo_root``, so a relative path in config.toml
    or on the command line means the same thing whatever directory the sweep was
    launched from -- and each worker resolves it identically, because it rebuilds
    the same Config from the same file. A bare ``OpticsSpec`` has no root and
    falls back to the process's cwd.
    """
    resolve = getattr(source, "path", None)
    return Path(resolve(name)) if callable(resolve) else Path(name)


def get_strategy(source, *, flat: bool | None = None,
                 fixed_shapes: str | None = None) -> SecondaryStrategy:
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

    ``[optics] fixed_shapes`` rides the same seam for the same reason:
    :class:`FixedShapeHeliostats` is the only place a frozen figure is applied,
    and this is the only place a strategy is built. ``fixed_shapes=""``
    explicitly means "no frozen figure" and beats a config that names a table,
    the way ``flat=False`` beats ``flat_mirrors = true``.

    Flat and fixed are refused together rather than ordered. Both write
    ``c3``/``c4``/``c5``, so a precedence rule would mean one of the two flags a
    user passed did nothing while the run reported both.
    """
    _load_strategies()

    name, flat_from_config, shapes_from_config = _layout_flat_shapes(source)
    if flat is None:
        flat = flat_from_config
    if fixed_shapes is None:
        fixed_shapes = shapes_from_config

    if flat and fixed_shapes:
        raise ValueError(
            f"flat mirrors and a fixed figure ({fixed_shapes}) both set c3/c4/c5 "
            f"and cannot both apply. Flat means no figure at all; fixed means one "
            f"ground figure. Pick one -- on the command line, drop --flat-mirrors, "
            f"or add --focused-mirrors if config.toml sets flat_mirrors = true."
        )

    try:
        strategy = _REGISTRY[name]()
    except KeyError:
        raise ValueError(
            f"Unknown secondary {name!r}. Available: {sorted(_REGISTRY)}"
        ) from None

    if flat:
        return FlatHeliostats(strategy)
    if fixed_shapes:
        path = _fixed_shapes_path(source, fixed_shapes)
        return FixedShapeHeliostats(strategy, load_fixed_shapes(path), source=str(path))
    return strategy
