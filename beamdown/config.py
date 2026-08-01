"""Typed configuration loaded from config.toml.

Replaces the constants that were hardcoded at the top of main.py and
heliostat_geometry.py -- including the ones that disagreed between the two.
"""

from __future__ import annotations

import datetime as _dt
import math
import os
import warnings
from dataclasses import dataclass, field as _dc_field, fields as _dc_fields
from pathlib import Path
from typing import Sequence

import tomli

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = REPO_ROOT / "config.toml"


@dataclass(frozen=True)
class Site:
    latitude: float
    longitude: float
    timezone: int


@dataclass(frozen=True)
class Geometry:
    secondary_height_mm: float
    receiver_offset_mm: float
    axicon_angle_deg: float
    # Rim of the conical secondary, matching the model's circular aperture on the
    # ``secondary`` surface. It is a 30 m wide opaque body directly over a field
    # whose innermost heliostats are at 30 m radius, so it shades them.
    #
    # Shared with the Cassegrain layout, whose hyperboloid is the same 30 m
    # across -- hence the name staying ``axicon_`` even though it is no longer
    # axicon-only. Renaming it would invalidate every stored run's config copy.
    axicon_aperture_radius_mm: float = 15000.0

    # -- shared-focus layouts (prime_focus, cassegrain) -------------------------
    #
    # F1 = (0, 0, focus_height_mm) is the ONE point the whole field aims at and
    # focuses on. Unused by the axicon, which has no single aim point: it derives
    # one per heliostat from radial position. No default, because guessing it
    # would silently mis-aim 645 heliostats; ``load_config`` requires it when the
    # layout needs it.
    focus_height_mm: float | None = None

    # Height of the Cassegrain hyperboloid's rim, i.e. of the horizontal circle
    # its silhouette projects from. Only the *shadow* geometry needs it, so it is
    # a plain height rather than anything derived from the conic constants.
    #
    # There is deliberately no axicon equivalent here: the cone's rim height is
    # derived from vertex height + radius * tan(angle) by
    # ``shading.SecondaryCone.rim_height_mm``, and duplicating it as a config
    # value would let the two disagree.
    secondary_rim_height_mm: float | None = None

    @property
    def receiver_height_mm(self) -> float:
        """Absolute receiver height.

        The model stores ``rec_offset`` relative to the secondary, while
        ``get_heliostat_axicon_shape`` wants an absolute receiver height.
        """
        return self.secondary_height_mm + self.receiver_offset_mm


@dataclass(frozen=True)
class FieldSpec:
    positions_file: str
    downselect_file: str
    n_configs: int
    mirror_width_mm: float
    mirror_height_mm: float
    # Height of the mirror's centre of rotation above grade.
    #
    # Defaults to 0 because that is what the optical model uses: the multiconfig
    # block carries ``posx``/``posy`` and no ``posz``, so Quadoa traces every
    # heliostat from z = 0. Raising it here alone would tilt the beam-to-aim
    # vector that blocking is measured along without the ray trace agreeing --
    # a silent disagreement between the two halves of the pipeline.
    #
    # Shading between mirrors that all share a height is invariant to that
    # height, and ``shading.self_check`` asserts it, so this knob is inert for
    # mutual shading. It exists for the day the model gains a real pedestal, and
    # for anything that needs the mirror to actually sit above ground -- at zero
    # the lower half of a steeply tilted mirror is underground.
    pedestal_height_mm: float = 0.0
    # Height used only when *drawing* the field from above, so cast shadows land
    # beside the mirrors instead of straddling them. Safe to set independently
    # of the line above precisely because mutual shading is invariant to a height
    # every mirror shares -- the picture shows the same occlusion relationships
    # the numbers were computed from, just translated to where a real pedestal
    # would put them.
    draw_pedestal_height_mm: float = 5000.0

    @property
    def mirror_area_m2(self) -> float:
        return (self.mirror_width_mm * self.mirror_height_mm) / 1e6

    @property
    def mirror_half_diagonal_mm(self) -> float:
        return 0.5 * math.hypot(self.mirror_width_mm, self.mirror_height_mm)


@dataclass(frozen=True)
class OpticsSpec:
    secondary: str
    mirror_reflectivity: float
    n_mirrors: int

    # Flat heliostats: keep the pointing, drop the focusing. A second comparison
    # axis orthogonal to ``secondary`` -- each of the three layouts can be run
    # with focused or flat mirrors, which is why this is a run-level flag rather
    # than a fourth layout.
    #
    # Applied in exactly one place, ``secondary.get_strategy``, which wraps the
    # chosen strategy in ``FlatHeliostats`` and forces c3 = c4 = c5 = 0. Nothing
    # downstream branches on it, so no code path can keep curvature by accident.
    #
    # Expect the physical consequence: a flat 5 m x 3 m mirror at 30-120 m throws
    # a far larger spot, so transmission and collected energy drop a long way.
    # That is the measurement.
    flat_mirrors: bool = False

    # Fixed mirror figure: path to a CSV of per-heliostat c3/c4/c5, empty for
    # none. A THIRD point on the same comparison axis as ``flat_mirrors`` --
    # focused re-figures the mirror every timestep (an idealisation: no ground
    # glass changes shape hourly), flat has no figure at all, and this freezes
    # each mirror to one figure it was ground to once. Pointing still tracks.
    #
    # Mutually exclusive with ``flat_mirrors``: both write c3/c4/c5, so letting
    # them compose would mean one silently winning. ``secondary.get_strategy``
    # refuses the pair rather than picking.
    #
    # A path, not a table, because the value has to survive the trip through
    # ``apply_overrides`` into every sweep worker, which re-reads config.toml
    # from disk and replays the override dict.
    fixed_shapes: str = ""

    @property
    def throughput(self) -> float:
        """Combined reflectivity, 0.9^2 = 0.81 by default.

        The model's coatings are ``ideal_mirror`` (100%), so this factor is
        applied in post-processing rather than by the ray trace.
        """
        return self.mirror_reflectivity ** self.n_mirrors


@dataclass(frozen=True)
class SourceSpec:
    power_w: float
    aperture_radius_mm: float
    angular_spread_rad: float

    @property
    def aperture_area_m2(self) -> float:
        return math.pi * (self.aperture_radius_mm / 1000.0) ** 2

    @property
    def dni_w_m2(self) -> float:
        """Implied DNI. Works out to exactly 1000 W/m^2 for the shipped model."""
        return self.power_w / self.aperture_area_m2

    def watts_per_ray(self, rays_emitted: int) -> float:
        """Power carried by one emitted ray.

        The source is never modified by this package, so this is a constant for
        a given ray budget and one ray always means the same thing.
        """
        if rays_emitted <= 0:
            raise ValueError("rays_emitted must be positive")
        return self.power_w / rays_emitted


@dataclass(frozen=True)
class DNISpec:
    mode: str
    constant_w_m2: float
    table_file: str

    def __post_init__(self) -> None:
        allowed = {"constant", "table", "monthly"}
        if self.mode not in allowed:
            raise ValueError(f"dni.mode must be one of {sorted(allowed)}, got {self.mode!r}")


def chunk_plan(rays_per_heliostat: int, rays_per_trace: int) -> list[int]:
    """Ray counts per ``traceRays`` call, summing **exactly** to the budget.

    One heliostat is not one ``traceRays`` call: the budget is split into chunks
    of ``rays_per_trace`` because a single enormous call allocates the whole ray
    set at once. So ``len(chunk_plan(...))`` is the number of
    ``setRayDistributionCount1`` + ``traceRays`` + ``getRayPos`` round trips each
    heliostat costs -- the count ``--rays-per-trace`` exists to control.

    The exact sum is the point of putting this in one place: it is used by
    :attr:`TraceSpec.chunk_sizes`, by :meth:`beamdown.session.QuadoaSession.trace`
    and by the GUI's derived "N traceRays calls per heliostat" label, and a
    remainder dropped in any one of them would silently emit fewer rays than the
    run says it emitted, scaling every reported watt with it. A budget that does
    not divide evenly gets a short final chunk:
    ``chunk_plan(100000, 30000) == [30000, 30000, 30000, 10000]``.

    A chunk larger than the whole budget is clamped to it rather than rejected,
    so this stays a pure arithmetic helper; refusing that combination is
    :func:`validate_trace`'s job, and only the command line can hit it (see
    :meth:`TraceSpec.__post_init__`).
    """
    total = int(rays_per_heliostat)
    if total <= 0:
        return []
    per = int(rays_per_trace)
    if per <= 0 or per > total:
        per = total
    full, rem = divmod(total, per)
    return [per] * full + ([rem] if rem else [])


@dataclass(frozen=True)
class TraceSpec:
    model_file: str
    quadoa_folder: str
    analysis_seq: int
    analysis_surface: int
    figure_seq: int
    bulk_config: int
    rays_per_heliostat: int
    rays_per_trace: int
    n_workers: int
    max_retries: int

    def __post_init__(self) -> None:
        if not 1 <= self.n_workers <= 4:
            raise ValueError(
                f"n_workers must be 1-4 (USB HASP license limit), got {self.n_workers}"
            )
        # A chunk bigger than the budget is meaningless, and clamping it silently
        # is safe for a value read from the FILE: config.toml's rays_per_trace is
        # a ceiling ("split into chunks of at most this many"), not a request.
        #
        # Overrides are different, and deliberately do not come through here --
        # apply_overrides writes the field directly. A command line that says
        # "--rays 6000 --rays-per-trace 30000" has asked for two contradictory
        # things explicitly, so :func:`validate_trace` refuses it instead of
        # quietly running something else.
        if self.rays_per_trace > self.rays_per_heliostat:
            object.__setattr__(self, "rays_per_trace", self.rays_per_heliostat)

    @property
    def n_chunks(self) -> int:
        """``traceRays`` calls per heliostat -- what ``--rays-per-trace`` buys."""
        return len(self.chunk_sizes)

    @property
    def chunk_sizes(self) -> list[int]:
        """Ray counts per traceRays call, summing exactly to rays_per_heliostat."""
        return chunk_plan(self.rays_per_heliostat, self.rays_per_trace)


@dataclass(frozen=True)
class SweepSpec:
    dates: tuple[_dt.date, ...]
    sunrise_margin_min: float
    hour_step: float


@dataclass(frozen=True)
class ReceiverSpec:
    window_mm: float
    grid_size: int

    @property
    def bin_size_mm(self) -> float:
        return 2.0 * self.window_mm / self.grid_size

    @property
    def bin_area_m2(self) -> float:
        return (self.bin_size_mm / 1000.0) ** 2

    @property
    def edges(self):
        import numpy as np

        return np.linspace(-self.window_mm, self.window_mm, self.grid_size + 1)


@dataclass(frozen=True)
class StorageSpec:
    root: str
    raw_rays: str

    def __post_init__(self) -> None:
        allowed = {"all", "downselected", "none"}
        if self.raw_rays not in allowed:
            raise ValueError(f"raw_rays must be one of {sorted(allowed)}, got {self.raw_rays!r}")


@dataclass(frozen=True)
class Config:
    site: Site
    geometry: Geometry
    field: FieldSpec
    optics: OpticsSpec
    source: SourceSpec
    dni: DNISpec
    trace: TraceSpec
    sweep: SweepSpec
    receiver: ReceiverSpec
    storage: StorageSpec
    repo_root: Path = _dc_field(default=REPO_ROOT)

    # -- path helpers -----------------------------------------------------
    def path(self, relative: str) -> Path:
        p = Path(relative)
        return p if p.is_absolute() else self.repo_root / p

    @property
    def positions_path(self) -> Path:
        return self.path(self.field.positions_file)

    @property
    def downselect_path(self) -> Path:
        return self.path(self.field.downselect_file)

    @property
    def model_path(self) -> Path:
        return self.path(self.trace.model_file)

    @property
    def fixed_shapes_path(self) -> Path | None:
        """The fixed-figure table, or ``None`` when the run has no fixed figure."""
        return self.path(self.optics.fixed_shapes) if self.optics.fixed_shapes else None

    @property
    def output_root(self) -> Path:
        return self.path(self.storage.root)


# How many reflections each layout puts in the path, and therefore what
# ``optics.n_mirrors`` should be. Prime focus has no secondary mirror; the other
# two bounce off the heliostat and then off the secondary.
_REFLECTIONS = {"axicon": 2, "cassegrain": 2, "prime_focus": 1}

# The layouts ``[optics] secondary`` and ``beamdown sweep --secondary`` accept.
# One source of truth so the CLI's choices cannot drift from the rules below;
# ``beamdown.secondary.available()`` is the registry's own answer but importing it
# just to build an argument parser would pull in the whole strategy stack.
SECONDARY_LAYOUTS = tuple(_REFLECTIONS)

# Reflections in the path per layout, exposed so a CLI or GUI can offer the
# consistent ``--n-mirrors`` for a chosen layout instead of making the user
# remember which layout has a secondary mirror in it.
def reflections_for(layout: str) -> int | None:
    return _REFLECTIONS.get(layout)

# Geometry keys each layout cannot do without, in the order they are reported.
_REQUIRED_GEOMETRY = {
    "prime_focus": ("focus_height_mm",),
    "cassegrain": ("focus_height_mm", "secondary_rim_height_mm"),
}


def _validate_layout(cfg: "Config") -> None:
    """Cross-check ``[optics] secondary`` against the ``[geometry]`` keys it needs.

    Separate from ``Geometry.__post_init__`` because which keys are required
    depends on the *optics* section, and a frozen dataclass validating a sibling
    section would have to reach outside itself.
    """
    layout = cfg.optics.secondary

    missing = [
        key for key in _REQUIRED_GEOMETRY.get(layout, ())
        if getattr(cfg.geometry, key, None) is None
    ]
    if missing:
        raise ValueError(
            f"[optics] secondary = {layout!r} requires "
            + ", ".join(f"[geometry] {k}" for k in missing)
            + f" to be set in config.toml, but {'they are' if len(missing) > 1 else 'it is'} "
            f"missing. {layout!r} aims the whole field at the single on-axis point "
            f"F1 = (0, 0, focus_height_mm), so that height has no sensible default"
            + (
                "; secondary_rim_height_mm is the height of the hyperboloid rim "
                "whose circular silhouette shades the field."
                if "secondary_rim_height_mm" in missing else "."
            )
        )

    expected = _REFLECTIONS.get(layout)
    if expected is not None and cfg.optics.n_mirrors != expected:
        warnings.warn(
            f"[optics] secondary = {layout!r} has {expected} reflection"
            f"{'s' if expected != 1 else ''} in the path "
            f"({'heliostat only' if expected == 1 else 'heliostat + secondary'}), "
            f"but n_mirrors = {cfg.optics.n_mirrors}, giving throughput "
            f"{cfg.optics.throughput:.4f} instead of "
            f"{cfg.optics.mirror_reflectivity ** expected:.4f}. "
            f"n_mirrors is NOT changed automatically: optics.throughput is applied "
            f"when a stored run is READ, not when it is written, so flipping it "
            f"here would silently rescale the numbers reported for every existing "
            f"run in analysis_output/. Edit n_mirrors yourself once you have "
            f"decided which runs that should affect.",
            stacklevel=3,
        )


def validate_layout(cfg: "Config") -> None:
    """Public entry point for :func:`_validate_layout`.

    ``load_config`` already runs it, but the CLI's ``--secondary`` /
    ``--n-mirrors`` / ``--focus-height-mm`` overrides land *after* the file has
    been read, so the combination that will actually run has to be checked once
    more -- with the same message, from the same place, rather than a second
    half-copy of the rules in ``beamdown.cli``.
    """
    _validate_layout(cfg)


def validate_trace(cfg: "Config") -> None:
    """Check the ray budget and the per-call chunk against each other.

    Runs *after* overrides land, for the same reason :func:`validate_layout`
    does: ``TraceSpec.__post_init__`` only sees the file, and
    :func:`apply_overrides` writes fields directly, so a combination that only a
    command line can produce is only checkable here.

    The rule the CLI implements on top of this:

    * ``--rays`` alone keeps the historical clamp -- the chunk is config.toml's
      value capped at the new budget, so a small test run does not emit a full
      60,000-ray chunk.
    * ``--rays`` **and** ``--rays-per-trace`` are both honoured literally, and a
      chunk larger than the budget is an error rather than a silent clamp,
      because the two flags then disagree about what was asked for.
    """
    t = cfg.trace
    if t.rays_per_heliostat <= 0:
        raise ValueError(
            f"[trace] rays_per_heliostat must be positive, got {t.rays_per_heliostat}"
        )
    if t.rays_per_trace <= 0:
        raise ValueError(
            f"[trace] rays_per_trace must be positive, got {t.rays_per_trace}. "
            f"It is the ray count of one traceRays call, not a divisor."
        )
    if t.rays_per_trace > t.rays_per_heliostat:
        raise ValueError(
            f"[trace] rays_per_trace = {t.rays_per_trace:,} is larger than "
            f"rays_per_heliostat = {t.rays_per_heliostat:,}. A chunk cannot be "
            f"bigger than the whole budget: the chunks sum to the budget, so the "
            f"trace would emit {t.rays_per_heliostat:,} rays in one call and the "
            f"chunk size would be a fiction. Ask for --rays-per-trace <= --rays "
            f"(one call per heliostat is --rays-per-trace equal to --rays), or "
            f"drop --rays-per-trace and let it be clamped to the budget."
        )


def apply_overrides(cfg: "Config", overrides: dict | None) -> None:
    """Force ``{section: {field: value}}`` onto an already-loaded config.

    The sections are frozen dataclasses because a config read from disk should
    not drift while a sweep runs, and ``object.__setattr__`` is the escape hatch
    ``beamdown.cli`` and ``scripts/report_energy.py`` already use to let a
    command-line flag beat the file.

    It lives here, next to the definitions, because the override set has to be
    replayed in **two** places: on the driver's config, and again on the copy
    each sweep worker loads from disk for itself
    (:func:`beamdown.sweep._init_worker`). A value set only on the driver's copy
    does not reach the trace -- which is exactly how ``--rays`` came to report a
    ray budget the workers never used.

    Deliberately does not validate: which checks matter depends on what moved,
    so the caller decides (the CLI calls :func:`validate_layout` once, after
    applying everything).
    """
    for section, values in (overrides or {}).items():
        target = getattr(cfg, section, None)
        if target is None:
            raise ValueError(
                f"unknown config section {section!r}; expected one of "
                f"{sorted(f.name for f in _dc_fields(cfg))}"
            )
        for key, value in (values or {}).items():
            if not hasattr(target, key):
                raise ValueError(f"[{section}] has no field {key!r}")
            object.__setattr__(target, key, value)


def _as_date(value) -> _dt.date:
    if isinstance(value, _dt.datetime):
        return value.date()
    if isinstance(value, _dt.date):
        return value
    return _dt.date.fromisoformat(str(value))


def load_config(path: str | os.PathLike | None = None) -> Config:
    """Load and validate config.toml."""
    cfg_path = Path(path) if path is not None else DEFAULT_CONFIG
    with open(cfg_path, "rb") as fh:
        raw = tomli.load(fh)

    sweep_raw = raw["sweep"]
    cfg = Config(
        site=Site(**raw["site"]),
        geometry=Geometry(**raw["geometry"]),
        field=FieldSpec(**raw["field"]),
        optics=OpticsSpec(**raw["optics"]),
        source=SourceSpec(**raw["source"]),
        dni=DNISpec(**raw["dni"]),
        trace=TraceSpec(**raw["trace"]),
        sweep=SweepSpec(
            dates=tuple(_as_date(d) for d in sweep_raw["dates"]),
            sunrise_margin_min=float(sweep_raw["sunrise_margin_min"]),
            hour_step=float(sweep_raw["hour_step"]),
        ),
        receiver=ReceiverSpec(**raw["receiver"]),
        storage=StorageSpec(**raw["storage"]),
        repo_root=cfg_path.resolve().parent,
    )
    _validate_layout(cfg)
    return cfg
