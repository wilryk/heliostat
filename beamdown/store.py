"""On-disk store for sweep results.

Layout under ``storage.root``::

    manifest.json                run metadata + quantisation scale
    summary.csv                  one row per (timestep, heliostat)
    raw/<key>_rays.npy           int16 (N, 2)  receiver x/y, all heliostats concatenated
    raw/<key>_index.npy          int64 (H, 3)  [heliostat_id, start, count]
    flux/<key>.npy               uint32 (H, G, G)  per-heliostat bin counts

Design notes
------------
**Counts are stored, never scaled flux.** Watts-per-ray, mirror reflectivity,
shading/blocking, and DNI are all applied at read time by
:func:`scale_factor`. Every one of those can therefore be revised without
re-tracing -- which is the whole point, given a sweep costs hours.

**Raw rays are the source of truth; flux maps are a cache.** They are written
during the sweep because binning during the trace is nearly free, but they are
fully reconstructible from the raw rays via :meth:`RunStore.rebin`.

**int16 quantisation** over the +/-window_mm receiver window gives 0.03 mm
resolution at half the size of float32 -- far finer than the 15.6 mm flux bins,
and irrelevant next to Monte-Carlo noise.

CSV (not Parquet) for the summary: 29k rows is trivially small, it needs no
extra dependency, and it opens in Excel.
"""

from __future__ import annotations

import datetime as _dt
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

INT16_MAX = 32767
SUMMARY_NAME = "summary.csv"
MANIFEST_NAME = "manifest.json"


@dataclass
class TimestepResult:
    """Everything one traced timestep contributes to the store."""

    key: str
    date: str
    hour: float
    solar_az_deg: float
    solar_el_deg: float
    heliostat_ids: np.ndarray
    rays_emitted: int
    counts: np.ndarray          # uint32 (H, G, G)
    rays: np.ndarray | None     # int16 (N, 2), or None when raw retention is off
    index: np.ndarray | None    # int64 (H, 3)
    rows: pd.DataFrame          # per-heliostat summary rows


def scale_factor(cfg, rays_emitted: int, dni_w_m2: float = 1000.0) -> float:
    """Watts per landed ray, including reflectivity and DNI.

    Multiply a bin count by this and divide by bin area to get W/m^2.
    """
    return (
        cfg.source.watts_per_ray(rays_emitted)
        * cfg.optics.throughput
        * (dni_w_m2 / 1000.0)
    )


def occlusion_weight_columns(manifest: dict, columns=()) -> tuple[str, ...]:
    """Summary columns a reader must multiply into a run's stored ray counts.

    One place decides, because getting it wrong is invisible: every wrong
    answer is a plausible number a few percent from the right one. Four
    readers ask (the GUI, ``compare``, ``cli figures`` and ``rescale``), and
    the answer depends on how the run was written:

    ``occluders`` true
        Neighbour shading and blocking are already in the ray counts.
        Only the secondary's shadow is left as a scalar -- and not even that
        once ``traced_secondary`` says its shadow plane is in the model.

    ``occlusion_form == "union"`` (runs from 2026-07-31 on)
        The sweep applied :func:`beamdown.shading.occlusion_efficiency` and
        stored it as ``eta_occlusion``. Reproduce it from that column; it is
        NOT recoverable from eta_shade and eta_block, whose product deletes
        the overlap twice.

    absent / ``"product"`` (every earlier run)
        ``eta_shade x eta_block``, which is what those runs' ``power_w``
        already carries. Left exactly as it was so old runs keep reading the
        way they always did.

    ``columns`` is the summary's column index when the caller has it; a
    union-form run whose summary somehow lacks ``eta_occlusion`` falls back to
    the product columns rather than silently weighting by nothing.
    """
    if manifest.get("occluders", False):
        return () if manifest.get("traced_secondary", False) else ("eta_secondary",)
    if manifest.get("occlusion_form", "product") == "union":
        if not len(columns) or "eta_occlusion" in columns:
            return ("eta_occlusion",)
    return ("eta_shade", "eta_block")


class RunStore:
    """Reader/writer for a sweep output directory."""

    def __init__(self, root, cfg=None, mode: str = "r"):
        self.root = Path(root)
        self.cfg = cfg
        self.mode = mode
        self.raw_dir = self.root / "raw"
        self.flux_dir = self.root / "flux"
        if mode == "w":
            for d in (self.root, self.raw_dir, self.flux_dir):
                d.mkdir(parents=True, exist_ok=True)
        elif not self.root.exists():
            raise FileNotFoundError(f"No store at {self.root}")
        self._manifest: dict | None = None

    # -- manifest ---------------------------------------------------------
    @property
    def manifest(self) -> dict:
        if self._manifest is None:
            path = self.root / MANIFEST_NAME
            self._manifest = json.loads(path.read_text()) if path.exists() else {}
        return self._manifest

    def write_manifest(self, extra: dict | None = None) -> None:
        cfg = self.cfg
        payload = {
            "created": _dt.datetime.now().isoformat(timespec="seconds"),
            "quantisation_scale_mm": cfg.receiver.window_mm / INT16_MAX,
            "receiver_window_mm": cfg.receiver.window_mm,
            "grid_size": cfg.receiver.grid_size,
            "rays_per_heliostat": cfg.trace.rays_per_heliostat,
            "source_power_w": cfg.source.power_w,
            "throughput": cfg.optics.throughput,
            "secondary": cfg.optics.secondary,
            # The other half of the run's optical identity. Written here rather
            # than by run_sweep so EVERY writer of a store records it, and read
            # back by the GUI's run readout: a finished run has to be able to say
            # whether its 645 heliostats were focused or flat, because the two
            # differ by a factor in collected energy and nothing else in the
            # store distinguishes them.
            "flat_mirrors": bool(getattr(cfg.optics, "flat_mirrors", False)),
            "site": asdict(cfg.site),
            "geometry": asdict(cfg.geometry),
            "raw_rays": cfg.storage.raw_rays,
        }
        # The third point on the figure axis, written ONLY when it applies. An
        # absent key means the historical behaviour -- the mirror was re-figured
        # every timestep -- which is what every run before this option was true
        # of, so a reader must not have to distinguish "no key" from "not fixed".
        # The value is the table's path as the run was given it: the CSV is not
        # copied into the store, and the field positions in it are what make a
        # run reproducible.
        fixed_shapes = str(getattr(cfg.optics, "fixed_shapes", "") or "")
        if fixed_shapes:
            payload["fixed_shapes"] = fixed_shapes
        payload.update(extra or {})
        self._manifest = payload
        (self.root / MANIFEST_NAME).write_text(json.dumps(payload, indent=2))

    @property
    def quant_scale(self) -> float:
        return float(self.manifest.get("quantisation_scale_mm", 1.0))

    # -- quantisation -----------------------------------------------------
    @staticmethod
    def inside_window(xy_mm: np.ndarray, window_mm: float) -> np.ndarray:
        """Mask of rays within the storable receiver window."""
        return (np.abs(xy_mm[:, 0]) <= window_mm) & (np.abs(xy_mm[:, 1]) <= window_mm)

    @staticmethod
    def quantise(xy_mm: np.ndarray, window_mm: float) -> np.ndarray:
        """Float mm -> int16.

        Rays must already be inside the window -- use :meth:`inside_window` to
        filter first. Clipping here instead would pile out-of-window rays onto
        the boundary, inventing a hot ring at the receiver edge and making the
        raw store disagree with the binned counts.
        """
        scaled = np.clip(xy_mm / window_mm, -1.0, 1.0) * INT16_MAX
        return np.rint(scaled).astype(np.int16)

    def dequantise(self, raw: np.ndarray) -> np.ndarray:
        return raw.astype(np.float32) * np.float32(self.quant_scale)

    # -- writing ----------------------------------------------------------
    def write_timestep(self, result: TimestepResult) -> None:
        np.save(self.flux_dir / f"{result.key}.npy", result.counts)
        if result.rays is not None and result.index is not None:
            np.save(self.raw_dir / f"{result.key}_rays.npy", result.rays)
            np.save(self.raw_dir / f"{result.key}_index.npy", result.index)
        self.append_summary(result.rows)

    def append_summary(self, rows: pd.DataFrame) -> None:
        path = self.root / SUMMARY_NAME
        rows.to_csv(path, mode="a" if path.exists() else "w",
                    header=not path.exists(), index=False)

    # -- reading ----------------------------------------------------------
    def timestep_keys(self) -> list[str]:
        return sorted(p.stem for p in self.flux_dir.glob("*.npy"))

    def has_timestep(self, key: str) -> bool:
        """Used to make a sweep resumable after an interruption."""
        return (self.flux_dir / f"{key}.npy").exists()

    def read_counts(self, key: str, mmap: bool = True) -> np.ndarray:
        """Per-heliostat bin counts, uint32 (H, G, G)."""
        return np.load(self.flux_dir / f"{key}.npy", mmap_mode="r" if mmap else None)

    def read_index(self, key: str) -> np.ndarray:
        return np.load(self.raw_dir / f"{key}_index.npy")

    def read_rays(self, key: str, heliostat_id: int | None = None) -> np.ndarray:
        """Receiver x/y in mm. One heliostat, or all of them concatenated.

        Memory-maps the file and slices, so reading one heliostat out of a
        200 MB timestep does not read the whole thing.
        """
        rays_path = self.raw_dir / f"{key}_rays.npy"
        if not rays_path.exists():
            raise FileNotFoundError(
                f"No raw rays for {key} (storage.raw_rays was "
                f"{self.manifest.get('raw_rays')!r} for this run)"
            )
        raw = np.load(rays_path, mmap_mode="r")
        if heliostat_id is None:
            return self.dequantise(np.asarray(raw))

        index = self.read_index(key)
        match = index[index[:, 0] == heliostat_id]
        if match.size == 0:
            raise KeyError(f"heliostat {heliostat_id} not in timestep {key}")
        _, start, count = match[0]
        return self.dequantise(np.asarray(raw[start:start + count]))

    def rebin(self, key: str, grid_size: int, window_mm: float,
              heliostat_id: int | None = None) -> np.ndarray:
        """Re-histogram raw rays at a different resolution or window."""
        xy = self.read_rays(key, heliostat_id)
        edges = np.linspace(-window_mm, window_mm, grid_size + 1)
        counts, _, _ = np.histogram2d(xy[:, 1], xy[:, 0], bins=[edges, edges])
        return counts

    # -- summary ----------------------------------------------------------
    def summary(self) -> pd.DataFrame:
        path = self.root / SUMMARY_NAME
        if not path.exists():
            raise FileNotFoundError(f"No summary at {path}")
        df = pd.read_csv(path)
        df["date"] = pd.to_datetime(df["date"]).dt.date
        return df

    # -- aggregation ------------------------------------------------------
    def field_flux(self, key: str, cfg=None, dni_w_m2: float = 1000.0,
                   efficiency: np.ndarray | None = None) -> np.ndarray:
        """Whole-field receiver flux in W/m^2 for one timestep.

        Flux maps add linearly, so the combined field map is a weighted sum over
        per-heliostat maps -- no re-trace, and per-heliostat efficiency factors
        (shading, blocking) drop straight in as weights.
        """
        cfg = cfg or self.cfg
        counts = np.asarray(self.read_counts(key)).astype(np.float64)
        if efficiency is not None:
            counts = counts * np.asarray(efficiency, float)[:, None, None]
        total = counts.sum(axis=0)
        rays_emitted = int(self.manifest.get("rays_per_heliostat",
                                             cfg.trace.rays_per_heliostat))
        return total * scale_factor(cfg, rays_emitted, dni_w_m2) / cfg.receiver.bin_area_m2

    def heliostat_flux(self, key: str, heliostat_row: int, cfg=None,
                       dni_w_m2: float = 1000.0, efficiency: float = 1.0) -> np.ndarray:
        """Single-heliostat receiver flux in W/m^2."""
        cfg = cfg or self.cfg
        counts = np.asarray(self.read_counts(key)[heliostat_row]).astype(np.float64) * efficiency
        rays_emitted = int(self.manifest.get("rays_per_heliostat",
                                             cfg.trace.rays_per_heliostat))
        return counts * scale_factor(cfg, rays_emitted, dni_w_m2) / cfg.receiver.bin_area_m2
