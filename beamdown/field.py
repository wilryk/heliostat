"""Heliostat field: positions, downselection, and neighbour queries.

Positions are held in **millimetres** internally, matching Quadoa and
``get_heliostat_axicon_shape``. The source files are in metres; conversion
happens once, here, rather than being repeated as ``xpos_s *= 1000`` at every
call site.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

_X_ALIASES = ("x (m)", "x_s (m)", "x", "x (mm)")
_Y_ALIASES = ("y (m)", "y_s (m)", "y", "y (mm)")


def _pick_column(df: pd.DataFrame, aliases: tuple[str, ...], which: str) -> str:
    lowered = {str(c).strip().lower(): c for c in df.columns}
    for alias in aliases:
        if alias in lowered:
            return lowered[alias]
    raise KeyError(
        f"No {which} column found. Looked for {aliases}, file has {list(df.columns)}"
    )


def _read_xy(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read an x/y position table from .xlsx or .csv, in metres."""
    path = Path(path)
    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)
    xcol = _pick_column(df, _X_ALIASES, "x")
    ycol = _pick_column(df, _Y_ALIASES, "y")
    x = df[xcol].to_numpy(dtype=float)
    y = df[ycol].to_numpy(dtype=float)
    if "mm" in str(xcol).lower():
        x, y = x / 1000.0, y / 1000.0
    return x, y


@dataclass(frozen=True)
class HeliostatField:
    """Heliostat centre positions, in millimetres."""

    x_mm: np.ndarray
    y_mm: np.ndarray
    ids: np.ndarray
    source: str = ""

    def __len__(self) -> int:
        return int(self.x_mm.size)

    @property
    def x_m(self) -> np.ndarray:
        return self.x_mm / 1000.0

    @property
    def y_m(self) -> np.ndarray:
        return self.y_mm / 1000.0

    @property
    def radius_mm(self) -> np.ndarray:
        return np.hypot(self.x_mm, self.y_mm)

    @property
    def azimuth_deg(self) -> np.ndarray:
        """Compass bearing of each heliostat from the tower, degrees CW from +y."""
        return np.mod(np.degrees(np.arctan2(self.x_mm, self.y_mm)), 360.0)

    @property
    def xy_mm(self) -> np.ndarray:
        return np.column_stack((self.x_mm, self.y_mm))

    def subset(self, indices) -> "HeliostatField":
        idx = np.asarray(indices, dtype=int)
        return HeliostatField(
            x_mm=self.x_mm[idx],
            y_mm=self.y_mm[idx],
            ids=self.ids[idx],
            source=f"{self.source}[subset n={idx.size}]",
        )

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "heliostat_id": self.ids,
                "x_m": self.x_m,
                "y_m": self.y_m,
                "radius_m": self.radius_mm / 1000.0,
                "azimuth_deg": self.azimuth_deg,
            }
        )

    def describe(self) -> str:
        r = self.radius_mm / 1000.0
        return (
            f"{len(self)} heliostats from {self.source}\n"
            f"  x {self.x_m.min():.1f}..{self.x_m.max():.1f} m, "
            f"y {self.y_m.min():.1f}..{self.y_m.max():.1f} m\n"
            f"  radius {r.min():.1f}..{r.max():.1f} m"
        )


def load_field(cfg) -> HeliostatField:
    """Load the full heliostat field (645 positions)."""
    path = cfg.positions_path
    x, y = _read_xy(path)
    fld = HeliostatField(
        x_mm=x * 1000.0,
        y_mm=y * 1000.0,
        ids=np.arange(x.size, dtype=int),
        source=path.name,
    )
    warn_coincident(fld)
    return fld


def coincident_pairs(field: HeliostatField, tol_mm: float = 1.0):
    """Positions shared by two or more heliostats.

    Two heliostats at the same point cannot both exist. They also fail to shade
    each other -- :func:`beamdown.shading._fraction_unoccluded` requires the
    occluder to be strictly ahead (``t > 1e-6``), and a coincident mirror is at
    ``t = 0`` -- so each is traced and summed at full power, double-counting that
    position. Cheap to detect, and invisible in every downstream total, so it is
    checked at load rather than left to be noticed.
    """
    tree = cKDTree(field.xy_mm)
    return sorted({(int(min(i, j)), int(max(i, j)))
                   for i, j in tree.query_pairs(r=float(tol_mm))})


def warn_coincident(field: HeliostatField, tol_mm: float = 1.0) -> list:
    """Report coincident positions without removing them.

    Dropping a heliostat would change every published total, so which one goes
    -- if either -- is a decision for whoever owns the position file.
    """
    import warnings

    pairs = coincident_pairs(field, tol_mm)
    if pairs:
        listed = ", ".join(f"{i}={j}" for i, j in pairs)
        warnings.warn(
            f"{field.source}: {len(pairs)} coincident heliostat position(s) "
            f"({listed}). They do not shade each other and each is traced at full "
            f"power, so those positions are double-counted.",
            stacklevel=2,
        )
    return pairs


def match_to_field(field: HeliostatField, x_m, y_m, tol_m: float = 1e-3) -> np.ndarray:
    """Map external x/y positions onto indices in ``field``.

    Used to turn a hand-made downselect file into field indices, so a
    downselected run and the full run share one heliostat numbering.
    """
    tree = cKDTree(field.xy_mm)
    query = np.column_stack((np.asarray(x_m, float) * 1000.0, np.asarray(y_m, float) * 1000.0))
    dist, idx = tree.query(query)
    bad = dist > tol_m * 1000.0
    if np.any(bad):
        raise ValueError(
            f"{int(bad.sum())} of {len(query)} positions are not in the field "
            f"(worst mismatch {dist.max()/1000.0:.4f} m)"
        )
    return idx.astype(int)


def downselect(
    field: HeliostatField,
    n: int,
    method: str = "farthest_point",
    seed: int = 0,
) -> np.ndarray:
    """Choose ``n`` representative heliostats. Returns indices into ``field``.

    ``farthest_point``
        Greedy max-dispersion, the same approach as
        ``heliostat/read_helio_pos.py`` and what produced
        ``downselected_x,y centers.xlsx``. Reproducible here: it starts from the
        heliostat nearest the field centroid rather than a random one. Heavily
        weights the perimeter, which is good for showing the extremes.

    ``uniform``
        Stratified over radius rings and azimuth sectors, giving coverage that
        tracks heliostat density instead of over-representing the boundary.
        Better when the downselect is meant to stand in for the whole field.
    """
    if n >= len(field):
        return np.arange(len(field))
    if n <= 0:
        raise ValueError("n must be positive")

    if method == "farthest_point":
        return _farthest_point(field.xy_mm, n)
    if method == "uniform":
        return _stratified_uniform(field, n, seed)
    raise ValueError(f"unknown downselect method {method!r}")


def _farthest_point(points: np.ndarray, n: int) -> np.ndarray:
    centroid = points.mean(axis=0)
    start = int(np.argmin(np.linalg.norm(points - centroid, axis=1)))

    selected = [start]
    min_dist = np.linalg.norm(points - points[start], axis=1)
    for _ in range(n - 1):
        nxt = int(np.argmax(min_dist))
        selected.append(nxt)
        min_dist = np.minimum(min_dist, np.linalg.norm(points - points[nxt], axis=1))
    return np.array(sorted(selected), dtype=int)


def _stratified_uniform(field: HeliostatField, n: int, seed: int) -> np.ndarray:
    """Split into radius rings of equal population, then spread over azimuth."""
    rng = np.random.default_rng(seed)
    r = field.radius_mm
    az = field.azimuth_deg
    order = np.argsort(r)

    n_rings = max(1, int(round(np.sqrt(n))))
    rings = np.array_split(order, n_rings)

    base, extra = divmod(n, n_rings)
    chosen: list[int] = []
    for i, ring in enumerate(rings):
        take = base + (1 if i < extra else 0)
        if take <= 0 or ring.size == 0:
            continue
        take = min(take, ring.size)
        # Spread the picks evenly in azimuth within this ring.
        ring_sorted = ring[np.argsort(az[ring])]
        offset = rng.integers(0, max(1, ring_sorted.size // take))
        picks = (np.linspace(0, ring_sorted.size, take, endpoint=False).astype(int) + offset)
        chosen.extend(ring_sorted[np.clip(picks, 0, ring_sorted.size - 1)].tolist())

    chosen = sorted(set(chosen))
    # Backfill if de-duplication lost any.
    if len(chosen) < n:
        remaining = [i for i in np.argsort(r) if i not in set(chosen)]
        chosen.extend(remaining[: n - len(chosen)])
    return np.array(sorted(chosen[:n]), dtype=int)


def load_or_build_downselect(cfg, field: HeliostatField, method: str = "farthest_point"):
    """Use the configured downselect file if present, else compute one.

    Returns ``(indices, provenance)``.
    """
    path = cfg.downselect_path
    if path.exists():
        x, y = _read_xy(path)
        idx = match_to_field(field, x, y)
        if idx.size != cfg.field.n_configs:
            raise ValueError(
                f"{path.name} has {idx.size} positions but field.n_configs = "
                f"{cfg.field.n_configs}. Update config.toml or the file."
            )
        return idx, f"file:{path.name}"
    idx = downselect(field, cfg.field.n_configs, method=method)
    return idx, f"computed:{method}"


def neighbour_pairs(field: HeliostatField, search_radius_mm: float):
    """Neighbours within ``search_radius_mm`` of each heliostat.

    Returns a list of index arrays, one per heliostat, excluding itself. Used by
    :mod:`beamdown.shading` to limit shading/blocking tests to plausible
    occluders instead of all 645.
    """
    tree = cKDTree(field.xy_mm)
    groups = tree.query_ball_point(field.xy_mm, r=search_radius_mm)
    return [np.array([j for j in g if j != i], dtype=int) for i, g in enumerate(groups)]
