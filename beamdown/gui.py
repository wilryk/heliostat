"""Desktop GUI for examining a sweep, heliostat by heliostat.

    python -m beamdown.gui --output analysis_output/demo25

The summary table answers *which* heliostat is worst; this answers *why*, and
hands off to Quadoa when the answer needs geometry. Six views over the same
selection -- field map, receiver spot, through-the-day curves, the distribution
across the field, annual energy, and the raw table -- with one heliostat and one
timestep selected at a time, so switching view never loses your place.

No *analysis* here re-traces: shading/blocking and DNI are weights applied to
stored counts, so toggling them is instant, and the "Open in Quadoa" button
exports an ``.optx`` with the selected heliostat's pointing, shape and sun
already loaded.

The **Design** tab is the odd one out: it is not a view of the loaded run at
all. It evaluates candidate secondaries -- cone, dish, or no secondary --
straight from :mod:`beamdown.design_eval`, draws the cross-section that says
where the light goes, and can hand the current sliders to the ``.optx``
builders. All of that is geometry over config plus the field file, so it costs
no licence seat and works with no run loaded.

The **Trace** tab is the exception, and the only place a licence seat is spent:
it composes a ``beamdown sweep`` command line from a form, shows exactly what it
will run, and either copies it or launches it as a detached subprocess with a log
tail and a Stop button. It edits no configuration -- every option is a flag,
because config.toml is shared state that a running sweep re-reads.

Built on Tkinter (standard library) with embedded matplotlib, so there is
nothing to install. A sweep that is still running can be opened read-only --
press Reload to pick up timesteps as they land.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import queue
import sys
import threading
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import (  # noqa: E402
    FigureCanvasTkAgg,
    NavigationToolbar2Tk,
)
from matplotlib.figure import Figure  # noqa: E402

# The layouts the Trace tab offers and the reflection count each implies, taken
# from the config module so the tab cannot offer a layout the CLI would reject.
# chunk_plan comes from the same place for the same reason: the tab's
# "N traceRays calls per heliostat" label must be the count the trace will
# actually make, not a second opinion about it.
from .config import SECONDARY_LAYOUTS, chunk_plan, reflections_for  # noqa: E402

# The one shared paper style. Applied in BeamdownGUI.__init__ before any Figure
# is built, because constrained_layout is read at figure-construction time --
# see plot_style.apply's docstring. Every "Save figure…" button goes back out
# through plot_style.save_figure, so what the GUI shows and what a script writes
# are the same figure.
from . import plot_style  # noqa: E402

# The DNI models the Energy tab's selector offers, alongside "whatever
# config.toml says" -- kept in sync with beamdown.dni's supported modes by
# hand, since that module has no public list of them to import.
DNI_MODE_DEFAULT = "(config default)"
DNI_MODES = ("constant", "table", "monthly")

# Columns offered for colouring the field map and for the through-day plot.
# Anything numeric in the summary works; these are just the useful ones first.
PREFERRED_COLUMNS = [
    "power_w", "power_in_aperture_w", "aperture_spillage", "transmission",
    "eta_shade", "eta_block", "eta_occlusion", "cosine_efficiency", "peak_flux_w_m2",
    "r90_mm", "r50_mm", "rms_radius_mm", "spillage", "aoi_deg", "rays_landed",
]

# Not in the summary -- recomputed from the stored counts against the aperture
# currently set, which is the whole point of the aperture control.
DERIVED_COLUMNS = ("power_in_aperture_w", "aperture_spillage")

# Summary columns that already carry shading x blocking and DNI, folded into
# watts-per-ray at sweep time by metrics.spot_metrics. Turning the weights off,
# or changing DNI, has to divide these back out -- otherwise the table quietly
# disagrees with the spot image drawn beside it, which is the bug that made the
# weights look like they were not being applied at all.
WEIGHTED_COLUMNS = ("power_w", "peak_flux_w_m2", "power_in_aperture_w")

BG = "#f7f7f7"

# Overlay colours, chosen to sit on the pale metric ramp below rather than on
# viridis: a translucent wash needs a light, low-chroma background to read.
SHADOW_COLOUR = (0.13, 0.35, 0.72)   # shading loss, on the mirror -- blue
BLOCK_COLOUR = (0.78, 0.13, 0.10)    # blocking loss, on the mirror -- red
GROUND_COLOUR = (0.32, 0.45, 0.62)   # the shadow pattern on the ground
CONE_EDGE = "#1b4f9c"                # outline of the secondary's own shadow
BODY_EDGE = "#b8541f"                # the secondary's physical footprint on the
                                     # axis -- deliberately a different hue from
                                     # CONE_EDGE, because the two circles are the
                                     # same size and only one of them moves

_METRIC_CMAP = None


def _metric_colormap():
    """Greys, truncated well short of black.

    The field map has to carry a metric *and* two translucent overlays. A full
    ramp to black leaves the overlays invisible at one end, so this stops at a
    mid grey and starts above white, keeping the whole range legible under a
    wash.
    """
    global _METRIC_CMAP
    if _METRIC_CMAP is None:
        import matplotlib as mpl
        from matplotlib.colors import ListedColormap

        base = mpl.colormaps["Greys"]
        _METRIC_CMAP = ListedColormap(base(np.linspace(0.10, 0.62, 128)),
                                      name="beamdown_metric")
    return _METRIC_CMAP


def _fmt(value) -> str:
    """Readable fixed-width-ish formatting for the readout panel."""
    if isinstance(value, (int, np.integer)):
        return f"{int(value):,}"
    if isinstance(value, (float, np.floating)):
        v = float(value)
        if not np.isfinite(v):
            return "-"
        if v == 0:
            return "0"
        if abs(v) >= 10000 or abs(v) < 0.01:
            return f"{v:.4g}"
        return f"{v:.4g}"
    return str(value)


# --------------------------------------------------------------------------
# Trace tab: setting up a sweep
# --------------------------------------------------------------------------
#
# Every option the Trace tab offers is a `beamdown sweep` flag, and that is a
# constraint rather than a convenience. Selecting a layout or a ray budget by
# editing config.toml would be shared mutable state: a sweep already running
# re-reads config.toml when it writes its end-of-run report, so an edit made
# mid-run corrupts that report, and two people cannot set up two different runs
# at once. So the tab builds a command line and nothing else -- which is also
# why the command is shown rather than hidden, and why it can simply be copied
# into a shell instead of launched from here.

# One heliostat trace, measured over full7/full8 (~800-900 ms/heliostat at
# 120,000 rays, 645 heliostats). Used only for the "this will take about..."
# estimate, the same 0.55 s/trace figure `beamdown info` quotes.
TRACE_SECONDS_PER_TRACE = 0.55

# Why the worker spinbox is pinned to 1. Measured, see scripts/run_full8.sh's
# header: asking for 2 got "asked for 2, 1 got a license seat", and the failed
# second request then leaked the first seat, so the NEXT single-worker run found
# zero seats and sat waiting ten minutes. Each failed request also pops a modal
# H0038 dialog that a human has to dismiss.
WORKERS_WARNING = (
    "1 worker. There is ONE HASP licence seat: a failed request for a second\n"
    "seat has been measured to leak the first, leaving the next run waiting\n"
    "with zero seats and a modal dialog to dismiss. Quadoa's tracer already\n"
    "runs at ~4x parallelism inside one session, so a second worker buys\n"
    "almost nothing here."
)


def build_sweep_argv(options: dict, python: str | None = None) -> list[str]:
    """Turn the Trace tab's form into the exact argv that will be run.

    Pure -- dict in, list out, no Tk and no filesystem -- because the command
    preview and the Run button must agree to the character. If they could
    disagree, the preview would be decoration rather than a record of what
    happened.

    ``--config`` is on the top-level parser, not on the sweep subcommand, so it
    has to come before the word ``sweep`` (see ``beamdown.cli.build_parser``).
    Flags whose value is ``None`` or blank are omitted, so the command shows only
    what actually departs from config.toml -- except ``--output``, ``--rays`` and
    ``--workers``, which are always written out because they decide where the run
    lands, what it costs, and whether it competes for the licence seat.
    """
    argv = [python or sys.executable, "-u", "-m", "beamdown"]
    if options.get("config"):
        argv += ["--config", str(options["config"])]
    argv += ["sweep"]

    if options.get("output"):
        argv += ["--output", str(options["output"])]

    # --dates wins over --suggest-dates in beamdown.cli._resolve_dates, so only
    # ever emit one of them; the form makes them mutually exclusive too.
    if options.get("dates"):
        argv += ["--dates", *[str(d) for d in options["dates"]]]
    elif options.get("suggest_dates"):
        argv += ["--suggest-dates", str(int(options["suggest_dates"]))]

    if options.get("all_heliostats"):
        argv += ["--all-heliostats"]
    if options.get("rays"):
        argv += ["--rays", str(int(options["rays"]))]
    # NOT one of the always-written flags: the chunk size is a departure from
    # config.toml like any other, and where it agrees with the file, saying so
    # would be a restatement rather than a decision. _trace_options only puts it
    # in the dict when it differs.
    if options.get("rays_per_trace"):
        argv += ["--rays-per-trace", str(int(options["rays_per_trace"]))]
    argv += ["--workers", str(int(options.get("workers") or 1))]

    if options.get("hour_step") is not None:
        argv += ["--hour-step", _trace_num(options["hour_step"])]
    if options.get("sunrise_margin_min") is not None:
        argv += ["--sunrise-margin-min", _trace_num(options["sunrise_margin_min"])]

    if options.get("secondary"):
        argv += ["--secondary", str(options["secondary"])]
    if options.get("focus_height_mm") is not None:
        argv += ["--focus-height-mm", _trace_num(options["focus_height_mm"])]
    if options.get("rim_height_mm") is not None:
        argv += ["--rim-height-mm", _trace_num(options["rim_height_mm"])]
    if options.get("n_mirrors"):
        argv += ["--n-mirrors", str(int(options["n_mirrors"]))]
    # Tri-state, like the CLI flag pair it maps to: absent from the dict means
    # "config.toml decides", and only a genuine disagreement with the file emits
    # anything -- in whichever direction the disagreement runs.
    if options.get("flat_mirrors") is not None:
        argv += ["--flat-mirrors" if options["flat_mirrors"] else "--focused-mirrors"]

    if options.get("occluders"):
        argv += ["--occluders"]
    elif options.get("model_file"):
        # --occluders selects its own model and the CLI refuses both together.
        argv += ["--model-file", str(options["model_file"])]

    if not options.get("resume", False):
        argv += ["--no-resume"]
    return argv


def describe_call_plan(rays, rays_per_trace) -> str:
    """The derived "N traceRays calls per heliostat" line, as one string.

    Pure, and built on ``config.chunk_plan`` -- the same function the trace
    itself splits on -- so the label cannot describe a split that will not
    happen. It is the read-only half of the rays/chunk pair: the two entries are
    inputs, the call count is what they mean.
    """
    sizes = chunk_plan(rays, rays_per_trace)
    if not sizes:
        return "→ no rays to trace"
    n = len(sizes)
    if n == 1:
        detail = f"{sizes[0]:,}"
    elif len(set(sizes)) == 1:
        detail = f"{n} x {sizes[0]:,}"
    else:                                    # an uneven budget: short last chunk
        detail = f"{n - 1} x {sizes[0]:,} + {sizes[-1]:,}"
    return (f"→ {n} traceRays call{'s' if n != 1 else ''} per heliostat"
            f"   ({detail} = {sum(sizes):,} rays)")


def _trace_num(value) -> str:
    """Numbers as a human would type them: 1.0 -> "1", 0.5 -> "0.5"."""
    f = float(value)
    return str(int(f)) if f == int(f) else repr(f)


def format_command(argv) -> str:
    """The argv as a single copy-pasteable line, quoting only what needs it."""
    return " ".join(f'"{a}"' if (" " in a or not a) else a for a in argv)


def trace_log_path(output_dir) -> Path:
    """``analysis_output/<name>.log`` beside the run directory, as run_full8.sh."""
    output_dir = Path(output_dir)
    return output_dir.parent / f"{output_dir.name}.log"


def trace_lock_dir(output_dir) -> Path:
    """``analysis_output/.<name>.lock``, the same lock the run scripts take."""
    output_dir = Path(output_dir)
    return output_dir.parent / f".{output_dir.name}.lock"


def scan_locks(roots) -> list[dict]:
    """Every ``.<name>.lock`` directory under ``roots``, with its pid file.

    The lock is a *directory*, because ``mkdir`` is an atomic claim, and it holds
    a ``pid`` file -- exactly the convention in ``scripts/run_full*.sh``, whose
    comment explains what it is for: two sweeps writing one run directory
    produced duplicated summary rows and a raw-ray index that disagreed with the
    ray file. A lock anywhere also means something may be holding the single
    licence seat, so all of them are reported, not just this run's.
    """
    found: list[dict] = []
    seen: set = set()
    for root in roots:
        root = Path(root)
        if not root.is_dir():
            continue
        for path in sorted(root.glob(".*.lock")):
            resolved = path.resolve()
            if not path.is_dir() or resolved in seen:
                continue
            seen.add(resolved)
            pid_file = path / "pid"
            try:
                pid = pid_file.read_text(encoding="utf-8", errors="replace").strip()
            except OSError:
                pid = ""
            found.append({"path": path, "pid": pid or "(no pid file)"})
    return found


def launch_refusal(output_dir, resume: bool, locks) -> str | None:
    """Why this run must not start, or ``None`` if it may.

    Pure so it can be tested without launching anything, and so the Run button
    and the status line cannot disagree about the reason.
    """
    output_dir = Path(output_dir)

    if locks:
        held = "\n".join(f"    {d['path']}   pid {d['pid']}" for d in locks)
        return (
            "A sweep lock is being held, so a sweep is probably running:\n"
            f"{held}\n\n"
            "There is ONE Quadoa HASP licence seat. Starting a second sweep "
            "would fail to get a seat and, measured, would leak the seat the "
            "running sweep holds -- killing a run that may be 24 hours in. Wait "
            "for it to finish, or delete the lock directory by hand if you are "
            "certain that pid is gone."
        )

    if not output_dir.parent.is_dir():
        return (f"{output_dir.parent} does not exist, so the run directory and "
                f"its log file cannot be written there.")

    if output_dir.exists() and not resume:
        return (
            f"{output_dir} already exists and 'resume' is off. Writing a second "
            f"sweep into one run directory duplicates summary rows and leaves the "
            f"raw-ray index disagreeing with the ray file. Tick 'resume' to keep "
            f"its finished timesteps and continue, or choose another name."
        )
    return None


def parse_progress(text: str) -> dict:
    """Pull ``[i/N] ... eta M min`` out of a sweep log tail.

    Returns ``{}`` when the log has not reached its first timestep yet. The
    ``[i/N]`` counter and the ETA are matched separately on purpose: a resumed
    timestep logs ``[3/44] ... already done, skipping`` with no ETA at all, and
    the last real ETA is still the best estimate available.
    """
    import re

    steps = re.findall(r"\[(\d+)/(\d+)\]", text)
    etas = re.findall(r"eta\s+([\d.]+)\s+min", text)
    out: dict = {}
    if steps:
        out["done"], out["total"] = int(steps[-1][0]), int(steps[-1][1])
    if etas:
        out["eta_min"] = float(etas[-1])
    if "Sweep complete" in text:
        out["complete"] = True
    return out


def tail_file(path, max_bytes: int = 16384) -> str:
    """The last ``max_bytes`` of a file, decoded loosely.

    Seeks rather than reading the whole thing: a 24-hour sweep's log is polled
    every couple of seconds and grows all day.
    """
    path = Path(path)
    try:
        size = path.stat().st_size
        with open(path, "rb") as fh:
            if size > max_bytes:
                fh.seek(size - max_bytes)
            raw = fh.read()
    except OSError as exc:
        return f"(cannot read {path}: {exc})"
    text = raw.decode("utf-8", errors="replace")
    return text.split("\n", 1)[-1] if len(raw) >= max_bytes else text


def detached_popen_kwargs() -> dict:
    """Popen flags that keep a launched sweep alive after this window closes.

    A multi-hour trace must not die because the GUI was closed, and it has to be
    killable as a group -- ``workers > 1`` starts a real process pool. On Windows
    a child inherits its parent's console and dies with it, so
    ``DETACHED_PROCESS`` gives it none and ``CREATE_NEW_PROCESS_GROUP`` makes it
    the root of its own group. ``start_new_session`` does both jobs elsewhere.
    """
    import subprocess

    if sys.platform == "win32":
        return {"creationflags": (subprocess.CREATE_NEW_PROCESS_GROUP
                                 | subprocess.DETACHED_PROCESS)}
    return {"start_new_session": True}


class _Tooltip:
    """Minimal hover tooltip. Tkinter has none, and the licence warning needs
    to be readable without permanently spending a paragraph of panel height on
    it (the short version stays visible beside the spinbox regardless)."""

    def __init__(self, widget, text: str):
        self.widget, self.text, self.tip = widget, text, None
        widget.bind("<Enter>", self._show, add="+")
        widget.bind("<Leave>", self._hide, add="+")

    def _show(self, _event=None) -> None:
        if self.tip is not None:
            return
        x = self.widget.winfo_rootx() + 20
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 6
        self.tip = tk.Toplevel(self.widget)
        self.tip.wm_overrideredirect(True)
        self.tip.wm_geometry(f"+{x}+{y}")
        tk.Label(self.tip, text=self.text, justify="left", background="#ffffe0",
                 relief="solid", borderwidth=1, font=("", 8)).pack()

    def _hide(self, _event=None) -> None:
        if self.tip is not None:
            self.tip.destroy()
            self.tip = None


class BeamdownGUI:
    def __init__(self, root: tk.Tk, cfg, store, config_path: str | None = None):
        self.root = root
        self.cfg = cfg
        self.store = store
        # Kept so the Trace tab's command line can carry the same --config this
        # window was opened with; None means "beamdown's default", which is what
        # omitting the flag already does.
        self.config_path = config_path

        self._jobs: "queue.Queue[tuple]" = queue.Queue()
        self._field_cache: dict = {}
        self._rows_cache: dict = {}
        self._aperture_cache: dict = {}
        self._poly_cache: dict = {}
        # The secondary's physical footprint depends only on cfg, so unlike
        # _poly_cache this is not keyed on the timestep. A 1-tuple rather than the
        # array itself, so "computed, and it is None" (prime focus) is
        # distinguishable from "not computed yet".
        self._footprint_cache: tuple | None = None
        self._dirty: set = set()
        self._syncing = False  # true while a widget is being updated in code

        # Annual energy is an 8760-hour walk (energy.annual_energy) -- too slow
        # to redo on every redraw, so it is computed once in a background
        # thread on first view of the Energy tab and cached here. Keyed on
        # whether shading x blocking is applied, since that changes which
        # power feeds the integral; the GUI's per-view DNI override does not
        # apply (see _energy_input_summary), so it is not part of the key.
        self._energy_cache: dict = {}
        self._energy_job_key = None
        # Keyed on the selected DNI mode (None = config.toml's default) so a
        # live switch never rereads a table file it already parsed, but each
        # mode is only built the first time it is actually selected.
        self._dni_provider_cache: dict = {}

        # -- Trace tab state ------------------------------------------------
        # The Popen handle of the sweep this window launched, if any. Only one at
        # a time, because there is only one licence seat -- and it is deliberately
        # detached, so it outlives this process (see detached_popen_kwargs).
        self._trace_proc = None
        self._trace_log_path: Path | None = None
        self._trace_lock_dir: Path | None = None
        self._trace_argv: list[str] = []
        self._trace_estimate_cache: dict = {}

        # -- Design tab state -----------------------------------------------
        # The Design tab evaluates geometry, not the loaded run: it reads the
        # field file and config through beamdown.design_eval and never touches
        # the store, so it works with no run loaded and spends no licence seat.
        self._design_after = None      # pending debounce callback id
        self._design_busy = False      # guards scale <-> spinbox write-back
        self._design_drawn = False     # has the picture ever been drawn
        self._design_result: dict = {}
        self._design_section = None    # the drawn cross-section, for CSV export
        self._design_proc_busy = False  # an export subprocess is running

        # Before _build_layout, which constructs every Figure: constrained
        # layout is a figure-construction-time rcParam, so applying the style
        # afterwards would leave every tab on matplotlib's default margins.
        plot_style.apply()

        self.root.title("beamdown explorer")
        self.root.geometry("1500x900")

        self._load_run()
        self._build_layout()
        self._poll_jobs()
        self._refresh_all()

    # ------------------------------------------------------------------
    # data
    # ------------------------------------------------------------------
    def _load_run(self) -> None:
        """Read the store. Tolerates a sweep still in progress."""
        self.summary = self.store.summary()

        # A running sweep writes the flux file before the summary row, so only
        # trust timesteps that appear in both.
        on_disk = set(self.store.timestep_keys())
        in_summary = set(self.summary.timestep.unique())
        self.keys = sorted(on_disk & in_summary)
        if not self.keys:
            raise SystemExit(f"No completed timesteps in {self.store.root}")

        # Time is two-dimensional -- a date and an hour within it -- and a single
        # flat list of 44 keys hides that. Kept as an explicit index so the date
        # and hour selectors, and stepping, all agree on the ordering.
        self.keys_by_date = {}
        for key in self.keys:
            self.keys_by_date.setdefault(key.split("_")[0], []).append(key)
        self.dates = sorted(self.keys_by_date)
        self.dates_as_dates = [_dt.datetime.strptime(d, "%Y%m%d").date() for d in self.dates]

        # Row order in the counts arrays comes from the manifest, not from
        # sorting ids -- they need not agree, and guessing would silently show
        # the wrong heliostat's spot.
        manifest_ids = self.store.manifest.get("heliostat_ids")
        if manifest_ids:
            self.ids = [int(i) for i in manifest_ids]
        else:
            self.ids = sorted(int(i) for i in self.summary.heliostat_id.unique())
        self.row_of = {hid: i for i, hid in enumerate(self.ids)}

        numeric = set(self.summary.select_dtypes("number").columns)
        ordered = [c for c in PREFERRED_COLUMNS
                   if c in numeric or c in DERIVED_COLUMNS]
        self.columns = ordered + sorted(c for c in numeric if c not in ordered
                                        and c not in ("heliostat_id", "hour"))
        self._field_cache.clear()
        self._rows_cache.clear()
        self._aperture_cache.clear()
        self._poly_cache.clear()
        # Depends only on cfg, but a newly opened run may have been produced by a
        # different layout, so it is not carried across a reload either.
        self._footprint_cache = None
        # A new run means a new summary, so any cached annual energy is stale.
        # The DNI provider only depends on cfg, not on the run, but is cheap to
        # reload here too in case the table file changed underfoot.
        self._energy_cache.clear()
        self._energy_job_key = None
        self._dni_provider_cache = {}

    def _step_rows(self, key=None):
        """Raw summary rows for one timestep, exactly as the sweep wrote them."""
        key = key or self.key
        return self.summary[self.summary.timestep == key].set_index("heliostat_id")

    def _rows_for_display(self, key=None):
        """Summary rows adjusted to the weights and aperture currently selected.

        Every view goes through here. The sweep folded shading x blocking and a
        1000 W/m^2 DNI into ``power_w`` and ``peak_flux_w_m2`` before writing
        them, so a view that reads the summary directly is showing weighted
        numbers whatever the controls say -- which is why the weights used to
        look as though they only reached the spot images.
        """
        key = key or self.key
        ck = (key, bool(self.var_shading.get()), self._dni(), self._aperture())
        if ck not in self._rows_cache:
            self._rows_cache[ck] = self._build_rows(key)
            if len(self._rows_cache) > 60:
                self._rows_cache.pop(next(iter(self._rows_cache)))
        return self._rows_cache[ck]

    def _landed_and_inside(self, key: str):
        """Rays landed, and rays inside the aperture, per heliostat row.

        Cached on (timestep, radius) and nothing else: these are ray *counts*, so
        weights and DNI do not enter, and one 42 MB counts file need only be read
        once per timestep however often the weights are toggled. The through-day
        plot walks every timestep, which without this would re-read the whole run
        on each redraw.
        """
        from .metrics import radial_mask

        radius = self._aperture()
        ck = (key, radius)
        if ck not in self._aperture_cache:
            counts = np.asarray(self.store.read_counts(key, mmap=False),
                                dtype=np.float64).reshape(len(self.ids), -1)
            landed = counts.sum(axis=1)
            inside = (counts[:, radial_mask(self.cfg, radius).ravel()].sum(axis=1)
                      if radius else landed)
            self._aperture_cache[ck] = (landed, inside, counts.max(axis=1))
        return self._aperture_cache[ck]

    def _build_rows(self, key: str):
        from .store import scale_factor

        rows = self._step_rows(key).copy()
        landed, inside, peak = self._landed_and_inside(key)
        emitted = int(self.store.manifest.get("rays_per_heliostat",
                                              self.cfg.trace.rays_per_heliostat))
        watts = scale_factor(self.cfg, emitted, self._dni())
        weight = self._eta_by_row(rows) if self.var_shading.get() else np.ones(len(self.ids))

        if not self.var_shading.get():
            # Rebuild from the stored counts rather than dividing the sweep's
            # weights back out. The axicon shades some heliostats completely, so
            # their recorded power is exactly zero with a weight of exactly zero
            # -- and 0/0 is not the unshaded power, it is a nan. Counts never
            # carried the weights, so they always have the answer.
            by_row = pd.Series(landed * watts, index=self.ids)
            rows["power_w"] = by_row.reindex(rows.index)
            if "peak_flux_w_m2" in rows.columns:
                rows["peak_flux_w_m2"] = pd.Series(
                    peak * watts / self.cfg.receiver.bin_area_m2,
                    index=self.ids).reindex(rows.index)
            if "shading_blocking_efficiency" in rows.columns:
                rows["shading_blocking_efficiency"] = 1.0
        else:
            dni = self._dni() / 1000.0
            for col in ("power_w", "peak_flux_w_m2"):
                if col in rows.columns:
                    rows[col] = rows[col].to_numpy(float) * dni

        # Aperture columns come from the stored counts, not the summary: the
        # sweep had no aperture to report against.

        power = pd.Series(inside * watts * weight, index=self.ids)
        spill = pd.Series(np.where(landed > 0, 1.0 - inside / np.maximum(landed, 1), np.nan),
                          index=self.ids)
        rows["power_in_aperture_w"] = power.reindex(rows.index)
        rows["aperture_spillage"] = spill.reindex(rows.index)
        return rows

    def _field_polygons(self, key=None):
        """Mirror outlines seen from above, and the shadows they cast, in metres.

        Both come straight from the summary's ``rot_az``/``rot_el`` -- the sweep
        already recorded the pointing, so nothing has to be re-solved and the
        drawing cannot disagree with the shading that was actually applied.

        Returns ``(ids, outlines, shadows, secondary_shadow)``. The mirror
        polygon arrays are ``(N, 4, 2)``: the outline is the rectangle projected
        straight down, so a steeply tilted mirror correctly draws as a narrow
        sliver; the shadow is the same rectangle projected along the sun onto the
        ground. ``secondary_shadow`` is the axicon's silhouette.

        Drawn from ``draw_pedestal_height_mm``, not the calculation's height. The
        model traces from z = 0, where a tilted mirror is half underground and
        its shadow straddles it rather than falling beside it. Lifting the field
        for the picture is exact, not cosmetic: mutual shading is invariant to a
        height every mirror shares, so the geometry drawn is the geometry the
        numbers came from.
        """
        from .shading import MirrorGeometry, corner_shadow, sun_vector

        key = key or self.key
        if key in self._poly_cache:
            return self._poly_cache[key]

        rows = self._step_rows(key)
        half_w = self.cfg.field.mirror_width_mm / 2.0
        half_h = self.cfg.field.mirror_height_mm / 2.0
        z = getattr(self.cfg.field, "draw_pedestal_height_mm", 5000.0)
        to_sun = sun_vector(float(rows.solar_az_deg.iloc[0]),
                            float(rows.solar_el_deg.iloc[0]))

        outlines, shadows = [], []
        for hid, r in rows.iterrows():
            g = MirrorGeometry.build(r.x_m * 1000.0, r.y_m * 1000.0,
                                     r.rot_az_deg, r.rot_el_deg, half_w, half_h, z)
            corners = np.array([
                g.centre + su * g.half_width * g.u + sv * g.half_height * g.v
                for su, sv in ((-1, -1), (1, -1), (1, 1), (-1, 1))
            ])
            outlines.append(corners[:, :2] / 1000.0)
            # Below the horizon the sun casts nothing; corner_shadow would divide
            # by a vanishing z-component and fling the polygon to infinity.
            shadows.append(corner_shadow(g, to_sun)[:, :2] / 1000.0
                           if to_sun[2] > 1e-3 else corners[:, :2] / 1000.0)

        # The secondary is projected onto the mirror plane, not the ground: the
        # whole body is ~27 m up, so its shadow there is exactly the set of
        # heliostats it shades, with none of the straddling a mirror's own
        # shadow suffers from.
        #
        # Whatever the layout put up there -- the axicon's cone, the Cassegrain's
        # disc, or nothing at all for prime focus, which yields an empty polygon
        # that every consumer below already handles via ``if len(cone)``.
        from .shading import secondary_body, secondary_shadow

        cone = secondary_shadow(secondary_body(self.cfg), to_sun, ground_z=z)

        out = (rows.index.to_numpy(), np.array(outlines), np.array(shadows),
               cone / 1000.0)
        self._poly_cache[key] = out
        if len(self._poly_cache) > 12:
            self._poly_cache.pop(next(iter(self._poly_cache)))
        return out

    def secondary_footprint(self, n: int = 128):
        """Where the secondary physically IS, in metres, or ``None``.

        A circle of radius ``geometry.axicon_aperture_radius_mm`` on the tower
        axis -- the body's own rim, taken from
        :func:`beamdown.shading.secondary_body` so the cone and the disc agree
        with the shading code instead of with a second copy of the number. Prime
        focus has no body and gets ``None``; the caller draws nothing.

        This is deliberately NOT the shadow. The shadow moves with the sun and
        lands on whichever heliostats are currently occluded; the footprint is
        fixed at (0, 0) and says what is hanging over the field. Both can be on
        screen at once, so :meth:`_draw_field` gives them different linestyles and
        separate legend entries.

        Independent of the timestep, so it is computed once per run rather than
        per redraw -- the Field tab redraws on every step of the time slider.

        Returned as a CLOSED polyline: the last point repeats the first, so
        ``ax.plot`` draws a closed circle without the caller having to wrap it.
        """
        if self._footprint_cache is None:
            from .shading import secondary_body

            body = secondary_body(self.cfg)
            if body is None:
                self._footprint_cache = (None,)
            else:
                theta = np.linspace(0.0, 2.0 * np.pi, n)
                r = float(body.aperture_radius_mm) / 1000.0
                self._footprint_cache = (
                    np.column_stack([r * np.cos(theta), r * np.sin(theta)]),
                )
        return self._footprint_cache[0]

    def _eta_series(self, rows) -> np.ndarray:
        """The weight the sweep actually applied, in row order."""
        cols = self._weight_columns
        eta = np.ones(len(rows))
        for col in cols:
            if col in rows.columns:
                eta = eta * rows[col].to_numpy(float)
        return eta

    @property
    def traced_occluders(self) -> bool:
        """Did this run put its neighbours in the ray path?

        If so the stored counts already carry shading and blocking, and applying
        those columns again would charge the same loss twice.
        """
        return bool(self.store.manifest.get("occluders", False))

    @property
    def _weight_columns(self) -> tuple:
        """Which summary columns the sweep actually multiplied into the counts.

        Delegated to :func:`beamdown.store.occlusion_weight_columns` so the
        GUI, ``compare`` and ``cli figures`` cannot drift apart on it -- a
        union-form run weighted by the old product columns reads as a small
        uniform deficit, which is indistinguishable from a real result.
        """
        from .store import occlusion_weight_columns

        return occlusion_weight_columns(self.store.manifest,
                                        self.summary.columns)

    def _eta_by_row(self, rows) -> np.ndarray:
        """Weights in manifest row order, which is the order the counts are in."""
        eta = np.ones(len(self.ids))
        for col in self._weight_columns:
            if col in rows.columns:
                eta = eta * rows[col].reindex(self.ids).fillna(1.0).to_numpy(float)
        return eta

    def _selected_row(self):
        rows = self._rows_for_display()
        if self.selected in rows.index:
            return rows.loc[self.selected]
        return None

    def _efficiency(self, key=None):
        """Per-heliostat shading x blocking weights, in manifest row order."""
        if not self.var_shading.get():
            return None
        return self._eta_by_row(self._step_rows(key))

    def _dni(self) -> float:
        try:
            return float(self.var_dni.get())
        except (ValueError, tk.TclError):
            return 1000.0

    # ------------------------------------------------------------------
    # annual energy (Energy tab)
    # ------------------------------------------------------------------
    def _current_date(self) -> _dt.date:
        return _dt.datetime.strptime(self.date_key, "%Y%m%d").date()

    def _dni_mode(self) -> str | None:
        """The Energy tab's selected DNI model, or ``None`` for whatever
        ``[dni] mode`` in config.toml currently says -- the same sentinel
        :func:`beamdown.dni.provider_for` uses to mean "use the config"."""
        value = self.var_dni_mode.get()
        return None if value == DNI_MODE_DEFAULT else value

    def _dni_provider(self):
        """The DNI model selected for annual energy -- constant, table,
        monthly, or (default) whatever ``[dni] mode`` in config.toml says.

        Independent of the "DNI W/m²" entry in the display panel: that entry
        is a what-if override for a single instant's spot/table view, whereas
        the annual integral needs a real hour-by-hour series. Conflating the
        two would rescale all 8760 hours by whatever the spot happens to be
        showing.

        Cached per mode rather than as a single slot, so flipping the
        selector back and forth never rereads a table file it has already
        parsed, but a mode that has not been selected yet is never built.
        """
        mode = self._dni_mode()
        if mode not in self._dni_provider_cache:
            from . import dni as D

            self._dni_provider_cache[mode] = D.provider_for(self.cfg, mode)
        return self._dni_provider_cache[mode]

    def _energy_input_summary(self):
        """Summary rows for the annual energy calc, honouring the shading toggle.

        With shading applied (the default) the stored ``power_w`` is exactly
        right -- it already carries the real shading x blocking, which is why
        ``energy.annual_energy`` reads it directly elsewhere in the codebase.
        With the toggle off this rebuilds the unweighted power from the stored
        ray counts, the same way ``_build_rows`` does for the other views, but
        always at the trace's native 1000 W/m² -- the annual integral applies
        its own hour-by-hour DNI, so the display DNI override must not also
        leak in here.
        """
        if self.var_shading.get():
            return self.summary

        from .store import scale_factor

        emitted = int(self.store.manifest.get("rays_per_heliostat",
                                              self.cfg.trace.rays_per_heliostat))
        watts = scale_factor(self.cfg, emitted, 1000.0)
        frames = []
        for key in self.keys:
            rows = self._step_rows(key).copy()
            landed, _inside, _peak = self._landed_and_inside(key)
            rows["power_w"] = pd.Series(landed * watts, index=self.ids).reindex(rows.index)
            frames.append(rows.reset_index())
        return pd.concat(frames, ignore_index=True)

    def _energy_key(self):
        # DNI mode belongs in the key alongside shading: it changes the
        # annual total exactly as much as shading does, and without it a
        # mode switch would silently keep serving the previous model's
        # cached MWh.
        return (bool(self.var_shading.get()), self._dni_mode())

    def _ensure_energy(self) -> None:
        """Kick off the background compute if the current key is not cached.

        annual_energy walks all 8760 hours of the year -- fast enough on
        today's 44-timestep sweep, but the whole point of this path is not to
        assume that stays true once the 12-date, finer-sampled sweep lands, so
        it always runs off the UI thread and reports through ``self._jobs``
        exactly like the Quadoa export button does.
        """
        key = self._energy_key()
        if key in self._energy_cache or self._energy_job_key == key:
            return
        self._energy_job_key = key
        self._status("computing annual energy (background, walks 8760 hours)…")

        try:
            provider = self._dni_provider()
        except Exception as exc:
            # Building the provider can now fail on the UI thread -- "table"
            # or "monthly" needs a downloaded CSV that may not exist -- so this
            # is routed through the exact same reporting path a failure inside
            # the background thread uses (see _poll_jobs' "energy_fail"
            # branch): status line gets the message, the job key is freed so
            # the next attempt (e.g. after fetching the file) is not silently
            # ignored as "already running", and the GUI keeps working.
            self._jobs.put(("energy_fail", {"key": key, "exc": exc}))
            return

        summary = self._energy_input_summary()
        cfg = self.cfg

        def work():
            try:
                from . import energy as E

                annual = E.annual_energy(summary, cfg, provider)
                checks = E.cross_check_daily_energy(summary, cfg, provider, annual=annual)
                sine = E.fit_annual_sine(annual["daily"])
                payload = {"key": key, "annual": annual, "checks": checks, "sine": sine}
                self._jobs.put(("energy_ok", payload))
            except Exception as exc:
                self._jobs.put(("energy_fail", {"key": key, "exc": exc}))

        threading.Thread(target=work, daemon=True).start()

    def _bins(self) -> int:
        try:
            return max(4, int(self.var_bins.get()))
        except (ValueError, tk.TclError):
            return int(self.cfg.receiver.grid_size)

    @property
    def bin_area_m2(self) -> float:
        """Area of one displayed bin. Follows the bin count, not the config."""
        g = self._bins()
        return (2.0 * self.cfg.receiver.window_mm / g) ** 2 / 1e6

    def _weighted_counts(self, key: str, efficiency) -> np.ndarray:
        counts = np.asarray(self.store.read_counts(key)).astype(np.float64)
        if efficiency is not None:
            counts = counts * np.asarray(efficiency, float)[:, None, None]
        return counts

    def _regrid(self, counts: np.ndarray, key: str, heliostat_id=None) -> np.ndarray:
        """Re-histogram a stored map to the displayed bin count.

        Coarsening is a block sum over the stored bins -- exact, instant, and it
        needs no raw rays, so it works on a run that did not keep them. Refining
        genuinely needs the rays back, which is cheap for one heliostat and a
        whole-file read for the field, so the cost is reported rather than hidden.
        """
        stored = int(self.cfg.receiver.grid_size)
        want = self._bins()
        if want == stored:
            return counts
        if stored % want == 0:
            f = stored // want
            return counts.reshape(want, f, want, f).sum(axis=(1, 3))
        try:
            self._status(f"rebinning to {want}x{want} from raw rays…")
            return self.store.rebin(key, want, self.cfg.receiver.window_mm,
                                    heliostat_id=heliostat_id)
        except FileNotFoundError:
            self._status(f"{want}x{want} needs raw rays, which this run did not "
                         f"keep -- showing {stored}x{stored}")
            self.var_bins.set(str(stored))
            return counts

    def _field_counts_from_rays(self, key: str, weights) -> np.ndarray:
        """Field map re-histogrammed finer than the stored grid, weights intact.

        The rays of all heliostats are stored concatenated, so a single
        histogram over the lot would lose the per-heliostat shading weights --
        they differ from mirror to mirror, which is the whole reason the field
        spot changes shape when the weights are toggled. Each heliostat's slice
        is therefore binned separately and accumulated. One file read, 645 small
        histograms.
        """
        want = self._bins()
        edges = np.linspace(-self.cfg.receiver.window_mm,
                            self.cfg.receiver.window_mm, want + 1)
        rays = self.store.read_rays(key)
        index = self.store.read_index(key)
        total = np.zeros((want, want))
        for hid, start, count in index:
            if count <= 0:
                continue
            xy = rays[start:start + count]
            hist, _, _ = np.histogram2d(xy[:, 1], xy[:, 0], bins=[edges, edges])
            total += hist * (weights[self.row_of[int(hid)]] if weights is not None else 1.0)
        return total

    def _field_flux(self):
        from .store import scale_factor

        ck = (self.key, self._dni(), bool(self.var_shading.get()), self._bins())
        if ck not in self._field_cache:
            eff = self._efficiency()
            stored = int(self.cfg.receiver.grid_size)
            if self._bins() > stored and stored % self._bins() != 0:
                try:
                    self._status(f"rebinning the field to {self._bins()}² from raw rays…")
                    total = self._field_counts_from_rays(self.key, eff)
                except FileNotFoundError:
                    self._status(f"{self._bins()}² needs raw rays, which this run did "
                                 f"not keep -- showing {stored}²")
                    self.var_bins.set(str(stored))
                    total = self._weighted_counts(self.key, eff).sum(axis=0)
                else:
                    self._status("ready")
            else:
                total = self._regrid(self._weighted_counts(self.key, eff).sum(axis=0),
                                     self.key)
            emitted = int(self.store.manifest.get("rays_per_heliostat",
                                                  self.cfg.trace.rays_per_heliostat))
            self._field_cache[ck] = (total * scale_factor(self.cfg, emitted, self._dni())
                                     / self.bin_area_m2)
            # Summing 645 x 128 x 128 is cheap but not free; a handful of
            # cached timesteps keeps stepping through time snappy.
            if len(self._field_cache) > 12:
                self._field_cache.pop(next(iter(self._field_cache)))
        return self._field_cache[ck]

    def _heliostat_flux(self):
        from .store import scale_factor

        row = self.row_of.get(self.selected)
        if row is None:
            return None
        eff = 1.0
        if self.var_shading.get():
            r = self._selected_row()
            if r is not None:
                # Whatever the sweep applied -- product, union, or the
                # secondary alone -- rather than a second opinion about it.
                eff = float(np.prod([float(r[c]) for c in self._weight_columns
                                     if c in r.index] or [1.0]))
        counts = np.asarray(self.store.read_counts(self.key)[row]).astype(np.float64)
        counts = self._regrid(counts, self.key, heliostat_id=self.selected) * eff
        emitted = int(self.store.manifest.get("rays_per_heliostat",
                                              self.cfg.trace.rays_per_heliostat))
        return counts * scale_factor(self.cfg, emitted, self._dni()) / self.bin_area_m2

    @property
    def key(self) -> str:
        return self.keys[self.step_i]

    @property
    def date_key(self) -> str:
        return self.key.split("_")[0]

    def hour_labels(self, date_key: str) -> list:
        """Short 'HH:MM  el NN' labels for one date's timesteps."""
        out = []
        for key in self.keys_by_date[date_key]:
            rows = self.summary[self.summary.timestep == key]
            hhmm = key.split("_")[1]
            if len(rows):
                out.append(f"{hhmm[:2]}:{hhmm[2:]}   el {rows.iloc[0]['solar_el_deg']:5.1f}°")
            else:
                out.append(f"{hhmm[:2]}:{hhmm[2:]}")
        return out

    def date_labels(self) -> list:
        out = []
        for d in self.dates:
            rows = self.summary[self.summary.timestep == self.keys_by_date[d][0]]
            label = str(rows.iloc[0]["date"]) if len(rows) else d
            out.append(f"{label}   ({len(self.keys_by_date[d])} steps)")
        return out

    def step_label(self, key: str) -> str:
        rows = self.summary[self.summary.timestep == key]
        if not len(rows):
            return key
        r = rows.iloc[0]
        return (f"{r['date']}  {r['hour']:04.1f}h   "
                f"el {r['solar_el_deg']:5.1f}°  az {r['solar_az_deg']:6.1f}°")

    # ------------------------------------------------------------------
    # layout
    # ------------------------------------------------------------------
    def _build_layout(self) -> None:
        self.step_i = len(self.keys) // 2
        self.selected = self.ids[0]

        self.var_colour = tk.StringVar(value=self.columns[0])
        self.var_metric = tk.StringVar(value="power_w")
        self.var_shading = tk.BooleanVar(value=True)
        self.var_dni = tk.StringVar(value="1000")
        # Annual-energy DNI model (Energy tab only) -- deliberately a separate
        # variable from var_dni above, which is the single-instant display
        # override; see _dni_provider's docstring.
        self.var_dni_mode = tk.StringVar(value=DNI_MODE_DEFAULT)
        self.var_aperture = tk.StringVar(value="700")
        self.var_spotview = tk.StringVar(value="image")
        self.var_bins = tk.StringVar(value=str(self.cfg.receiver.grid_size))
        self.var_show_shadow = tk.BooleanVar(value=True)
        self.var_show_block = tk.BooleanVar(value=True)
        self.var_log = tk.BooleanVar(value=False)
        self.var_helio = tk.StringVar(value=str(self.selected))
        self.var_nexport = tk.StringVar(value="1")
        self.var_extreme = tk.StringVar(value="max")
        self.var_status = tk.StringVar(value="ready")

        outer = ttk.Frame(self.root, padding=6)
        outer.pack(fill="both", expand=True)

        self.panel = ttk.Frame(outer, width=330)
        self.panel.pack(side="left", fill="y", padx=(0, 8))
        self.panel.pack_propagate(False)
        self._build_panel()

        self.book = ttk.Notebook(outer)
        self.book.pack(side="right", fill="both", expand=True)
        self.book.bind("<<NotebookTabChanged>>", lambda e: self._redraw_current())

        self.figures, self.canvases = {}, {}
        for name in ("Field", "Spot", "Through day", "Distribution"):
            self.figures[name], self.canvases[name] = self._add_figure_tab(name)
        self._build_energy_tab()
        self._build_design_tab()
        self._build_table_tab()
        self._build_trace_tab()

        self.root.bind("<Left>", lambda e: self._arrow(self._step, -1))
        self.root.bind("<Right>", lambda e: self._arrow(self._step, +1))
        self.root.bind("<Up>", lambda e: self._arrow(self._step_date, -1))
        self.root.bind("<Down>", lambda e: self._arrow(self._step_date, +1))

        status = ttk.Label(self.root, textvariable=self.var_status, anchor="w",
                           relief="sunken", padding=(6, 2))
        status.pack(side="bottom", fill="x")

        # A sweep launched from the Trace tab outlives this window on purpose, so
        # closing it says so rather than silently orphaning a run and its lock.
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

    def _add_figure_tab(self, name: str):
        frame = ttk.Frame(self.book)
        self.book.add(frame, text=name)
        fig = Figure(figsize=(9, 6.4), dpi=100, facecolor="white")
        canvas = FigureCanvasTkAgg(fig, master=frame)
        self._export_bar(frame, name)          # packed to the bottom, below the toolbar
        canvas.get_tk_widget().pack(fill="both", expand=True)
        NavigationToolbar2Tk(canvas, frame).update()
        if name == "Field":
            canvas.mpl_connect("button_press_event", self._on_field_click)
        return fig, canvas

    def _export_bar(self, parent, name: str, figure: bool = True):
        """The same two actions, in the same place, on every tab that has data.

        Bottom-right of the tab, under the plot: "Save figure…" and "Save data
        (CSV)…". Identical wording and order everywhere so there is nothing to
        learn per tab -- the Table tab, which has no figure, simply omits the
        first rather than showing a dead button.
        """
        bar = ttk.Frame(parent, padding=(6, 2, 6, 4))
        bar.pack(side="bottom", fill="x")
        ttk.Label(bar, text=plot_style.describe(), foreground="#888888"
                  ).pack(side="left")
        ttk.Button(bar, text="Save data (CSV)…", width=18,
                   command=lambda n=name: self._save_data_dialog(n)
                   ).pack(side="right")
        if figure:
            ttk.Button(bar, text="Save figure…", width=14,
                       command=lambda n=name: self._save_figure_dialog(n)
                       ).pack(side="right", padx=(0, 4))
        return bar

    def _build_energy_tab(self) -> None:
        """Headline numbers above the annual-energy figure, same figure/canvas
        machinery as ``_add_figure_tab`` but with a label strip on top for the
        numbers the user actually asked for -- MWh today is the headline, not
        buried in a plot title."""
        frame = ttk.Frame(self.book)
        self.book.add(frame, text="Energy")

        top = ttk.Frame(frame, padding=(8, 8, 8, 2))
        top.pack(fill="x")

        # Annual-energy DNI model -- a separate control from the "DNI W/m²"
        # display override in the side panel (see _dni_provider), so this is
        # labelled explicitly rather than reusing that entry's wording.
        row = ttk.Frame(top); row.pack(fill="x")
        ttk.Label(row, text="annual-energy DNI model", width=22
                 ).pack(side="left")
        self.cmb_dni_mode = ttk.Combobox(row, state="readonly", width=18,
                                         textvariable=self.var_dni_mode,
                                         values=[DNI_MODE_DEFAULT, *DNI_MODES])
        self.cmb_dni_mode.pack(side="left")
        self.cmb_dni_mode.bind("<<ComboboxSelected>>", self._on_dni_mode)

        self.lbl_energy_headline = ttk.Label(top, text="", font=("", 14, "bold"))
        self.lbl_energy_headline.pack(anchor="w", pady=(6, 0))
        self.lbl_energy_annual = ttk.Label(top, text="", foreground="#444444")
        self.lbl_energy_annual.pack(anchor="w", pady=(2, 0))
        self.lbl_energy_check = ttk.Label(top, text="", foreground="#666666")
        self.lbl_energy_check.pack(anchor="w", pady=(2, 0))

        fig = Figure(figsize=(9, 5.6), dpi=100, facecolor="white")
        canvas = FigureCanvasTkAgg(fig, master=frame)
        self._export_bar(frame, "Energy")
        canvas.get_tk_widget().pack(fill="both", expand=True)
        NavigationToolbar2Tk(canvas, frame).update()
        self.figures["Energy"] = fig
        self.canvases["Energy"] = canvas

    # ------------------------------------------------------------------
    # Design tab
    # ------------------------------------------------------------------
    #
    # Everything below is licence-free and run-independent. It answers "what
    # would this secondary do?" from geometry alone (beamdown.design_eval), so
    # it works while a sweep holds the seat and works with no run loaded at all.
    # The two things it can produce are a picture and an .optx to validate in
    # Quadoa; it never writes into analysis_output/.

    # Slider ranges, in the units the labels show. Chosen to bracket the built
    # design and the settled candidates with room either side, not to encode a
    # feasibility limit -- the readout says when a setting is infeasible, the
    # slider does not silently prevent it.
    DESIGN_PARAMS = {
        "axicon": [
            ("tip", "Cone tip height", 24.0, 34.0, 0.1, "m", 1),
            ("angle", "Cone half-angle", 13.0, 25.0, 0.1, "°", 1),
        ],
        "cassegrain": [
            ("rim", "Dish rim height", 28.0, 35.0, 0.1, "m", 1),
            ("f1", "Prime focus height", 34.0, 42.0, 0.1, "m", 1),
        ],
        "prime_focus": [
            ("pf", "Focus height", 34.0, 50.0, 0.1, "m", 1),
        ],
    }

    DESIGN_LAYOUTS = [
        ("axicon", "Cone (axicon) — the built design"),
        ("cassegrain", "Dish (cassegrain) — hyperboloid relay"),
        ("prime_focus", "Prime focus — no secondary at all"),
    ]

    def _build_design_tab(self) -> None:
        frame = ttk.Frame(self.book)
        self.book.add(frame, text="Design")

        # Defaults are the built axicon and the settled cassegrain/prime-focus
        # candidates, so the tab opens on numbers the owner already knows.
        self._design_vars = {
            "tip": tk.DoubleVar(value=27.0),
            "angle": tk.DoubleVar(value=20.0),
            "rim": tk.DoubleVar(value=30.0),
            "f1": tk.DoubleVar(value=36.0),
            "pf": tk.DoubleVar(value=36.0),
        }
        self.var_design_layout = tk.StringVar(value="axicon")
        self.var_design_date = tk.StringVar(value="2026-02-20")
        self.var_design_hour = tk.StringVar(value="9.454")

        left = ttk.Frame(frame, width=430, padding=(8, 8, 4, 8))
        left.pack(side="left", fill="y")
        left.pack_propagate(False)

        # -- which secondary ------------------------------------------------
        box = ttk.LabelFrame(left, text="What sits above the field", padding=6)
        box.pack(fill="x")
        for value, label in self.DESIGN_LAYOUTS:
            ttk.Radiobutton(box, text=label, value=value,
                            variable=self.var_design_layout,
                            command=self._design_layout_changed).pack(anchor="w")

        # -- the knobs -------------------------------------------------------
        self.frm_design_params = ttk.LabelFrame(left, text="Shape it", padding=6)
        self.frm_design_params.pack(fill="x", pady=(6, 0))

        # -- what it does ----------------------------------------------------
        box = ttk.LabelFrame(left, text="What that gives you", padding=6)
        box.pack(fill="x", pady=(6, 0))
        self._design_lines = []
        for _ in range(9):
            lbl = ttk.Label(box, text="", wraplength=390, justify="left",
                            anchor="w")
            lbl.pack(fill="x", anchor="w")
            self._design_lines.append(lbl)

        from .design_eval import HONESTY_FOOTER

        ttk.Label(left, text=HONESTY_FOOTER, wraplength=410, justify="left",
                  foreground="#777777").pack(fill="x", pady=(6, 0))

        # -- export ----------------------------------------------------------
        box = ttk.LabelFrame(left, text="Export a model to check in Quadoa",
                             padding=6)
        box.pack(fill="x", pady=(6, 0))
        row = ttk.Frame(box); row.pack(fill="x")
        ttk.Label(row, text="date").pack(side="left")
        ttk.Entry(row, textvariable=self.var_design_date, width=11
                  ).pack(side="left", padx=(3, 8))
        ttk.Label(row, text="hour").pack(side="left")
        ttk.Entry(row, textvariable=self.var_design_hour, width=8
                  ).pack(side="left", padx=3)
        self.btn_design_figure = ttk.Button(
            box, text="Write 25-heliostat model (.optx) at date/hour…",
            command=self._design_export_figure)
        self.btn_design_figure.pack(fill="x", pady=(5, 0))
        self.btn_design_full = ttk.Button(
            box, text="Write full-field cassegrain model (.optx)",
            command=self._design_export_cassegrain)
        self.btn_design_full.pack(fill="x", pady=(3, 0))
        # No --force is ever passed: if the target exists the builder refuses,
        # and its refusal is what shows up in the box below, word for word.
        self.txt_design_log = tk.Text(box, height=9, wrap="none", font=("Consolas", 8))
        self.txt_design_log.pack(fill="both", expand=True, pady=(5, 0))

        # The picture and its export bar share a right-hand frame, so the bar
        # lands under the plot exactly as it does on every other tab rather than
        # under the controls column.
        right = ttk.Frame(frame)
        right.pack(side="right", fill="both", expand=True)
        fig = Figure(figsize=(9, 6.4), dpi=100, facecolor="white")
        canvas = FigureCanvasTkAgg(fig, master=right)
        self._export_bar(right, "Design")
        canvas.get_tk_widget().pack(fill="both", expand=True)
        NavigationToolbar2Tk(canvas, right).update()
        self.figures["Design"] = fig
        self.canvases["Design"] = canvas

        self._design_build_params()

    def _design_build_params(self) -> None:
        """Rebuild the slider rows for the selected layout."""
        for child in self.frm_design_params.winfo_children():
            child.destroy()
        for key, label, lo, hi, step, unit, dec in \
                self.DESIGN_PARAMS[self.var_design_layout.get()]:
            var = self._design_vars[key]
            row = ttk.Frame(self.frm_design_params)
            row.pack(fill="x", pady=(4, 0))
            ttk.Label(row, text=label).pack(side="left")
            ttk.Label(row, text=unit, foreground="#666666", width=2
                      ).pack(side="right")
            spin = ttk.Spinbox(row, textvariable=var, from_=lo, to=hi,
                               increment=step, width=8, format=f"%.{dec}f",
                               command=self._design_changed)
            spin.pack(side="right")
            spin.bind("<Return>", lambda e: self._design_changed())
            scale = ttk.Scale(self.frm_design_params, from_=lo, to=hi,
                              variable=var,
                              command=lambda v, k=key, s=step:
                                  self._design_scale(k, s, v))
            scale.pack(fill="x")

    def _design_scale(self, key: str, step: float, value) -> None:
        """Snap a dragged scale onto the spinbox's own increment.

        A ttk.Scale writes continuous floats into the shared variable, which
        would show as ``29.372814...`` in the Spinbox beside it. Snapping writes
        back into the same variable, which fires this callback again -- hence
        the guard rather than a second variable, which would let the two widgets
        disagree about what the design currently is.
        """
        if self._design_busy:
            return
        try:
            v = float(value)
        except (TypeError, ValueError):
            return
        snapped = round(round(v / step) * step, 6)
        self._design_busy = True
        try:
            if abs(snapped - self._design_vars[key].get()) > 1e-12:
                self._design_vars[key].set(snapped)
        except tk.TclError:
            pass
        finally:
            self._design_busy = False
        self._design_changed()

    def _design_changed(self, *_args) -> None:
        """Coalesce a burst of slider motion into one evaluation."""
        if self._design_after is not None:
            try:
                self.root.after_cancel(self._design_after)
            except Exception:
                pass
        self._design_after = self.root.after(100, self._design_refresh)

    def _design_layout_changed(self) -> None:
        self._design_build_params()
        self._design_changed()

    def _design_value(self, key: str) -> float:
        """The current value of one knob, tolerating a half-typed Spinbox."""
        try:
            return float(self._design_vars[key].get())
        except (tk.TclError, ValueError):
            return float("nan")

    def _design_evaluate(self) -> dict:
        """Evaluate the current design. Pure geometry -- no run, no licence."""
        from . import design_eval as DE

        # Explicitly the config this window was opened with, not design_eval's
        # default: a GUI started with --config other.toml must not evaluate one
        # field while drawing another.
        field = DE.load_field_data(self.config_path)
        layout = self.var_design_layout.get()
        if layout == "axicon":
            return DE.eval_axicon(self._design_value("tip") * 1000.0,
                                  self._design_value("angle"), field)
        if layout == "cassegrain":
            return DE.eval_cassegrain(self._design_value("rim") * 1000.0,
                                      self._design_value("f1") * 1000.0, field)
        return DE.eval_prime_focus(self._design_value("pf") * 1000.0, field)

    def _design_refresh(self) -> None:
        self._design_after = None
        try:
            res = self._design_evaluate()
        except Exception as exc:
            self._status(f"Design: {type(exc).__name__}: {exc}")
            traceback.print_exc()
            return
        self._design_result = res
        self._design_readout(res)
        try:
            self._draw_design(res)
        except Exception as exc:
            self._status(f"Design picture: {type(exc).__name__}: {exc}")
            traceback.print_exc()
        self._design_drawn = True
        # The full-field cassegrain builder only has a dish to build for the
        # cassegrain layout; disabling it is honest, a refusal dialog is not.
        self.btn_design_full.config(
            state=("normal" if res.get("layout") == "cassegrain"
                   and res.get("feasible") else "disabled"))

    # -- the words ------------------------------------------------------
    def _design_readout(self, res: dict) -> None:
        """Plain language only: no repo jargon, no symbol names without a gloss."""
        lines: list[tuple[str, str]] = []          # (text, colour)
        grey, ink, good, warn, bad = "#777777", "#222222", "#1a7f37", "#b26a00", "#b3261e"

        if not res.get("feasible"):
            lines.append(("This shape does not work:", bad))
            for note in res.get("notes", []):
                lines.append(("  " + note, bad))
            if res.get("layout") == "cassegrain" and "K" in res:
                lines.append((f"Dish shape would be: K {res['K']:.3f}, vertex "
                              f"radius {res['R_v_mm']/1000:.2f} m", grey))
        else:
            lines.append((f"Yearly energy (relative to the built axicon): "
                          f"{res['energy_index']:.3f}×", ink))
            lines.append((f"Estimated field blocking+shading kept: "
                          f"{res['occlusion']*100:.1f}%", ink))
            lines.append((f"Sun image at receiver (ideal): {res['r90_mm']:.0f} mm "
                          f"r90; inner {res['r90_inner_mm']:.0f} / outer "
                          f"{res['r90_outer_mm']:.0f}", ink))
            lines.append(("", ink))

            layout = res["layout"]
            if layout == "axicon":
                ratio = res["correction_ratio"]
                colour = good if ratio <= 0.95 else (warn if ratio <= 1.0 else bad)
                lines.append((f"Inner beams hit the cone "
                              f"{res['inner_hit_mm']/1000:.2f} m from the axis — "
                              f"sagittal correction at {ratio*100:.0f}% of the "
                              f"built design's limit", colour))
                lines.append((f"Cone: tip {res['tip_mm']/1000:.1f} m, rim "
                              f"{res['rim_z_mm']/1000:.2f} m, "
                              f"{res['cone_depth_mm']/1000:.2f} m tall", ink))
            elif layout == "cassegrain":
                lines.append((f"Hyperboloid: K {res['K']:.3f}, vertex radius "
                              f"{res['R_v_mm']/1000:.2f} m, vertex at "
                              f"{res['vertex_z_mm']/1000:.2f} m, dish "
                              f"{res['sag_mm']/1000:.2f} m deep", ink))
                fill = res["fill_fraction"]
                lines.append((f"Beam fills {fill*100:.1f}% of the dish",
                              good if fill > 0.9 else warn))
                lines.append((f"Relay magnification {res['magnification_min']:.2f}"
                              f"–{res['magnification_max']:.2f}×; solar disk at "
                              f"the receiver {res['disk_outer_mm']:.0f} mm from "
                              f"the outer ring", ink))
                lines.append((res["blocking_note"], grey))
            else:
                gain = (res["bounce_gain"] - 1.0) * 100.0
                lines.append((f"Receiver hangs at the focus, "
                              f"{res['f1_mm']/1000:.1f} m up — no secondary "
                              f"mirror, so one reflection instead of two "
                              f"(+{gain:.0f}% energy for free)", good))
                lines.append((f"Solar disk: {res['disk_inner_mm']:.0f} mm from "
                              f"the inner ring, {res['disk_outer_mm']:.0f} mm "
                              f"from the outer", ink))
                lines.append((f"Beams arrive {res['tilt_inner_deg']:.0f}°–"
                              f"{res['tilt_outer_deg']:.0f}° off vertical, "
                              f"stretching the spot {res['stretch_inner']:.2f}×–"
                              f"{res['stretch_outer']:.2f}×", ink))
                lines.append((res["blocking_note"], grey))

        for lbl, (text, colour) in zip(self._design_lines,
                                       lines + [("", ink)] * len(self._design_lines)):
            lbl.config(text=text, foreground=colour)

    # -- the picture ----------------------------------------------------
    def _draw_design(self, res: dict) -> None:
        """A radial cut through the tower axis: where the light actually goes.

        Half a section (radius >= 0 on the right, with a little of the far side
        shown where the axicon's aim image needs it) rather than a full one:
        the field is 90 m wide and the optics are 40 m tall, so a full section
        halves the scale for a mirror image that adds nothing. Drawn at equal
        aspect, because the whole point is that the angles are real.
        """
        import numpy as _np

        from . import design_eval as DE

        fig = self.figures["Design"]
        fig.clear()
        ax = fig.add_subplot(111)
        self._style(ax)

        field = DE.load_field_data(self.config_path)
        r_in, r_out = field.R_min / 1000.0, field.R_max / 1000.0
        f2 = DE.F2_MM / 1000.0
        rim_r = DE.RIM_RADIUS_MM / 1000.0
        try:
            ap = abs(float(self.var_aperture.get())) / 1000.0
        except (ValueError, tk.TclError):
            ap = 0.7
        ap = max(ap, 0.4)                       # or it is invisible at this scale

        C_RAY, C_REFL, C_SEC = "#c2701c", "#1b4f9c", "#333333"
        C_GHOST = "#aaaaaa"
        x_lo = -1.0
        top = 40.0

        # -- ground and the two heliostats we follow -------------------------
        ax.axhline(0.0, color="#999999", lw=1.0, zorder=1)
        for r, name in ((r_in, f"inner heliostat\n{r_in:,.0f} m"),
                        (r_out, f"outer heliostat\n{r_out:,.1f} m")):
            ax.plot([r, r], [0.0, 2.2], color="#444444", lw=3, solid_capstyle="butt",
                    zorder=4)
            ax.annotate(name, (r, 0.0), xytext=(0, -8), textcoords="offset points",
                        ha="center", va="top", fontsize=8, color="#444444")

        layout = res.get("layout")
        rays: list[tuple[float, float, float, float]] = []   # heliostat -> hit
        refl: list[tuple[float, float, float, float]] = []   # hit -> receiver
        ghost: list[tuple[float, float, float, float]] = []  # dashed continuation
        # The secondary's own profile, sampled as drawn -- kept so "Save data
        # (CSV)" exports the curve that is on screen rather than re-deriving a
        # second version of it that could differ.
        profile: list[tuple[float, float]] = []

        if layout == "axicon" and res.get("feasible"):
            tip = res["tip_mm"] / 1000.0
            ang = res["angle_deg"]
            rim_z = res["rim_z_mm"] / 1000.0
            top = max(top, rim_z + 6.0)

            R = _np.array([field.R_min, field.R_max])
            x_r, y_r, x_a, y_a, _s, _sp = DE.geometry_terms(R, res["tip_mm"], ang)
            aim = (float(x_r) / 1000.0, (res["tip_mm"] + float(y_r)) / 1000.0)
            x_lo = min(x_lo, aim[0] - 4.0)
            top = max(top, aim[1] + 4.0)

            # The cone flank, tip on the axis, rim at 15 m.
            profile = [(0.0, tip), (rim_r, rim_z)]
            ax.plot([0.0, rim_r], [tip, rim_z], color=C_SEC, lw=3, zorder=5)
            ax.annotate(f"cone flank, {ang:.1f}° half-angle",
                        (rim_r * 0.55, tip + rim_r * 0.55 * _np.tan(_np.deg2rad(ang))),
                        xytext=(6, 10), textcoords="offset points", fontsize=8,
                        color=C_SEC)
            ax.annotate(f"tip {tip:.1f} m", (0.0, tip), xytext=(-6, 2),
                        textcoords="offset points", fontsize=8, color=C_SEC,
                        ha="right")

            for i, r0 in enumerate((r_in, r_out)):
                hx, hz = float(x_a[i]) / 1000.0, (res["tip_mm"] + float(y_a[i])) / 1000.0
                rays.append((r0, 0.0, hx, hz))
                # Proven in design_eval: the aim point is the receiver mirrored
                # in the flank, so the reflected ray goes exactly to F2.
                refl.append((hx, hz, 0.0, f2))
                ghost.append((hx, hz, aim[0], aim[1]))

            ax.plot([aim[0]], [aim[1]], marker="*", ms=13, color=C_GHOST, zorder=6)
            ax.annotate("aim image\n(the receiver reflected in the cone —\n"
                        "every heliostat aims here)", (aim[0], aim[1]),
                        xytext=(0, -12), textcoords="offset points", fontsize=8,
                        color="#666666", ha="center", va="top")

        elif layout == "cassegrain" and res.get("feasible"):
            des = DE.close_design(res["rim_z_mm"], field.R_max, DE.F2_MM,
                                  DE.RIM_RADIUS_MM, res["f1_mm"])
            f1 = res["f1_mm"] / 1000.0
            top = max(top, f1 + 5.0)

            rr = _np.linspace(0.0, rim_r, 200)
            zz = des.surface_z(rr * 1000.0) / 1000.0
            profile = list(zip(rr.tolist(), _np.asarray(zz, float).tolist()))
            ax.plot(rr, zz, color=C_SEC, lw=3, zorder=5)
            ax.annotate(f"hyperboloid dish, K {res['K']:.2f}\n"
                        f"{res['sag_mm']/1000:.2f} m deep, 15 m rim",
                        (rim_r * 0.6, des.surface_z(rim_r * 600.0) / 1000.0),
                        xytext=(8, -26), textcoords="offset points", fontsize=8,
                        color=C_SEC)

            o, d = DE.rays_to_f1(_np.array([field.R_min, field.R_max]),
                                 _np.zeros(2), des)
            _t, hit, _ok, _ = des.intersect(o, d)
            for i, r0 in enumerate((r_in, r_out)):
                hx = float(_np.hypot(hit[i, 0], hit[i, 1])) / 1000.0
                hz = float(hit[i, 2]) / 1000.0
                rays.append((r0, 0.0, hx, hz))
                refl.append((hx, hz, 0.0, f2))
                ghost.append((hx, hz, 0.0, f1))

            ax.plot([0.0], [f1], marker="x", ms=10, mew=2, color=C_GHOST, zorder=6)
            ax.annotate(f"F1, the prime focus at {f1:.1f} m\n"
                        f"(the image the dish relays down)", (0.0, f1),
                        xytext=(8, 6), textcoords="offset points", fontsize=8,
                        color="#666666")

        elif layout == "prime_focus":
            f1 = res["f1_mm"] / 1000.0
            top = max(top, f1 + 5.0)
            for r0 in (r_in, r_out):
                rays.append((r0, 0.0, 0.0, f1))
            ax.plot([-ap, ap], [f1, f1], color="#b3261e", lw=5,
                    solid_capstyle="butt", zorder=6)
            ax.annotate(f"receiver, up at the focus, {f1:.1f} m\n"
                        f"(nothing above the field but this)", (0.0, f1),
                        xytext=(12, -14), textcoords="offset points", fontsize=8,
                        color="#b3261e")
            ax.annotate("no secondary mirror in this layout —\n"
                        "the beams simply keep going to F1",
                        (rim_r * 0.6, f1 * 0.55), fontsize=8, color=C_GHOST,
                        style="italic")

        for x0, z0, x1, z1 in rays:
            ax.plot([x0, x1], [z0, z1], color=C_RAY, lw=1.4, zorder=3)
        for x0, z0, x1, z1 in refl:
            ax.plot([x0, x1], [z0, z1], color=C_REFL, lw=1.4, zorder=3)
        for x0, z0, x1, z1 in ghost:
            ax.plot([x0, x1], [z0, z1], color=C_GHOST, lw=0.9, ls=":", zorder=2)
        if rays:
            x0, z0, x1, z1 = rays[-1]
            ax.annotate("beam up from the mirror",
                        (0.55 * x0 + 0.45 * x1, 0.55 * z0 + 0.45 * z1),
                        xytext=(0, 8), textcoords="offset points", fontsize=8,
                        color=C_RAY, ha="center")
        if refl:
            x0, z0, x1, z1 = refl[-1]
            ax.annotate("reflected down to the receiver",
                        (0.5 * (x0 + x1), 0.5 * (z0 + z1)),
                        xytext=(8, 0), textcoords="offset points", fontsize=8,
                        color=C_REFL, va="center")

        if not res.get("feasible"):
            # Say why the picture is empty, in the picture -- an unexplained
            # blank plot beside a red readout reads as a broken tab.
            import textwrap as _tw

            why = "\n".join(_tw.fill(n, 64) for n in res.get("notes", []))
            ax.text(0.5, 0.62, "this shape does not work\n\n" + why,
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color="#b3261e")

        # -- the receiver aperture at z = 7 m --------------------------------
        if layout == "prime_focus":
            ax.plot([-ap, ap], [f2, f2], color=C_GHOST, lw=4,
                    solid_capstyle="butt", zorder=5)
            ax.annotate("existing receiver level, 7 m — empty in this layout",
                        (ap, f2), xytext=(8, -3), textcoords="offset points",
                        fontsize=8, color=C_GHOST)
        else:
            ax.plot([-ap, ap], [f2, f2], color="#b3261e", lw=5,
                    solid_capstyle="butt", zorder=6)
            ax.annotate(f"receiver aperture, {f2:.0f} m\n(±{ap*1000:,.0f} mm)",
                        (ap, f2), xytext=(10, -2), textcoords="offset points",
                        fontsize=8, color="#b3261e")

        top += 4.0                      # headroom, so no label runs into the title
        ax.axvline(0.0, color="#dddddd", lw=1.0, ls="--", zorder=0)
        ax.annotate("tower axis", (0.0, top - 0.5), xytext=(4, -2),
                    textcoords="offset points", fontsize=8, color="#bbbbbb",
                    va="top")

        ax.set_xlim(x_lo - 2.0, r_out * 1.06)
        ax.set_ylim(-5.0, top)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("distance from the tower axis (m)", fontsize=9)
        ax.set_ylabel("height (m)", fontsize=9)
        title = dict(self.DESIGN_LAYOUTS)[layout] if layout else "design"
        ax.set_title(title, fontsize=10, pad=8)

        # Everything drawn, as points, for the CSV export: one row per vertex,
        # tagged by which part of the picture it belongs to. Two-point entries
        # are straight segments; "secondary" is the sampled profile.
        section = [("secondary", r, z) for r, z in profile]
        for tag, segs in (("beam_up", rays), ("beam_down", refl),
                          ("construction", ghost)):
            for i, (x0, z0, x1, z1) in enumerate(segs):
                section.append((f"{tag}_{i}", x0, z0))
                section.append((f"{tag}_{i}", x1, z1))
        self._design_section = pd.DataFrame(
            section, columns=["part", "radius_m", "height_m"])

        plot_style.finish(fig)
        self.canvases["Design"].draw_idle()

    # -- exporting a model ----------------------------------------------
    def _design_log(self, text: str) -> None:
        self.txt_design_log.delete("1.0", "end")
        self.txt_design_log.insert("1.0", text)
        self.txt_design_log.see("end")

    def _design_instant(self):
        """``(date, hour)`` from the two entry boxes, or None after complaining."""
        date = self.var_design_date.get().strip()
        hour = self.var_design_hour.get().strip()
        try:
            _dt.date.fromisoformat(date)
            float(hour)
        except ValueError as exc:
            self._design_log(f"date/hour: {exc}\n"
                             f"expected YYYY-MM-DD and a decimal hour like 9.454")
            return None
        return date, hour

    def _design_export_figure(self) -> None:
        """25-heliostat model for the current design, at the given instant.

        Runs the existing builder scripts as subprocesses rather than importing
        them: they are argparse programs with their own verification and their
        own refusals, and the point of this button is to get exactly what the
        command line gets, including the refusal to overwrite.
        """
        res = self._design_result or {}
        if not res.get("feasible"):
            self._design_log("this design is not feasible — nothing to export.\n"
                             "Move the sliders until the readout goes green.")
            return
        instant = self._design_instant()
        if instant is None:
            return
        date, hour = instant

        layout = res["layout"]
        argv = [sys.executable, "scripts/build_figure_model.py",
                "--date", date, "--hour", hour]
        if layout == "axicon":
            argv += ["--tip-height-mm", f"{res['tip_mm']:.1f}",
                     "--axicon-angle-deg", f"{res['angle_deg']:g}"]
        else:
            argv += ["--secondary",
                     "prime_focus" if layout == "prime_focus" else "cassegrain",
                     "--focus-height-mm", f"{res['f1_mm']:.1f}"]

        # For the cassegrain the figure model only carries the POINTING -- which
        # is the shared-focus solve, identical to prime focus at the same F1 --
        # so the dish itself is a second pass over the file the first one wrote.
        # Hence a chain rather than one command.
        chain = None
        if layout == "cassegrain":
            argv += ["--rim-height-mm", f"{res['rim_z_mm']:.1f}"]
            chain = [sys.executable, "scripts/build_cassegrain_model.py",
                     "--rim-z-mm", f"{res['rim_z_mm']:.1f}",
                     "--f1-mm", f"{res['f1_mm']:.1f}"]
        self._design_run([argv], chain_on_output=chain,
                         chain_suffix=f"_dish{res.get('rim_z_mm', 0)/1000:g}")


    def _design_export_cassegrain(self) -> None:
        """Full-field cassegrain sweep model at the current sliders' geometry."""
        res = self._design_result or {}
        if res.get("layout") != "cassegrain" or not res.get("feasible"):
            self._design_log("select a feasible Dish (cassegrain) design first.")
            return
        self._design_run([[sys.executable, "scripts/build_cassegrain_model.py",
                           "--rim-z-mm", f"{res['rim_z_mm']:.1f}",
                           "--f1-mm", f"{res['f1_mm']:.1f}"]])

    def _design_run(self, commands, chain_on_output=None,
                    chain_suffix: str = "_dish") -> None:
        """Run builder scripts off the UI thread and report their own words back.

        ``chain_on_output`` is a partial command line completed with
        ``--base <first command's output> --out <that, tagged>`` -- the output
        path is taken from the first script's own "wrote ..." line rather than
        recomputed here, so the two can never disagree about the file name.
        """
        import subprocess

        if self._design_proc_busy:
            self._design_log("a model build from this tab is still running.")
            return
        self._design_proc_busy = True
        for btn in (self.btn_design_figure, self.btn_design_full):
            btn.config(state="disabled")
        self._design_log("building…\n" + "\n".join(
            format_command(c) for c in commands))
        self._status("building a model (background, no licence seat)…")
        cwd = str(self.cfg.repo_root)

        def work():
            chunks, rc, chained = [], 0, False
            verdict = ""
            try:
                pending = list(commands)
                while pending:
                    argv = pending.pop(0)
                    proc = subprocess.run(argv, cwd=cwd, capture_output=True,
                                          text=True, errors="replace")
                    out = proc.stdout or ""
                    err = (proc.stderr or "").strip()
                    chunks.append(f"$ {format_command(argv)}\n{out}"
                                  + (f"\n[stderr]\n{err}\n" if err else ""))
                    rc = proc.returncode
                    # The verdict comes from STDOUT only. The builders load the
                    # field file, which warns on stderr about the two coincident
                    # heliostats every single time -- letting that be the last
                    # line would report a known field-file quirk as the result of
                    # the build.
                    tail = [ln.strip() for ln in out.splitlines() if ln.strip()]
                    said = [ln for ln in tail
                            if ln.startswith(("wrote ", "refusing ", "PASS",
                                              "FAIL", "cannot "))]
                    if said or tail:
                        verdict = (said or tail)[-1]
                    if rc != 0:
                        break
                    if chain_on_output is not None and not chained and not pending:
                        chained = True
                        wrote = [ln for ln in out.splitlines()
                                 if ln.startswith("wrote ")]
                        if not wrote:
                            chunks.append("cannot chain the dish onto it: the "
                                          "builder printed no 'wrote <path>' line")
                            rc = 1
                            break
                        base = Path(wrote[-1].split("wrote ", 1)[1]
                                    .split("  ")[0].strip())
                        pending.append(list(chain_on_output) + [
                            "--base", str(base),
                            "--out", str(base.with_name(base.stem + chain_suffix
                                                        + base.suffix))])
                self._jobs.put(("design_build", {"text": "\n".join(chunks),
                                                 "rc": rc, "verdict": verdict}))
            except Exception as exc:
                self._jobs.put(("design_build",
                                {"text": f"{type(exc).__name__}: {exc}",
                                 "rc": 1, "verdict": f"{type(exc).__name__}: {exc}"}))

        threading.Thread(target=work, daemon=True).start()

    def _build_panel(self) -> None:
        p = self.panel

        # -- run ---------------------------------------------------------
        box = ttk.LabelFrame(p, text="run", padding=6)
        box.pack(fill="x", pady=(0, 6))
        self.lbl_run = ttk.Label(box, text=self.store.root.name, font=("", 9, "bold"))
        self.lbl_run.pack(anchor="w")
        self.lbl_run_info = ttk.Label(box, text="", foreground="#666666")
        self.lbl_run_info.pack(anchor="w")
        row = ttk.Frame(box); row.pack(fill="x", pady=(4, 0))
        ttk.Button(row, text="Open…", command=self._open_run, width=9).pack(side="left")
        ttk.Button(row, text="Reload", command=self._reload, width=9).pack(side="left", padx=4)

        # -- time --------------------------------------------------------
        # Date and hour are separate selectors rather than one flat list of every
        # timestep: "the same hour on each date" and "through this day" are the
        # two questions actually asked of a sweep, and a single 44-entry combobox
        # answers neither without scrolling.
        box = ttk.LabelFrame(p, text="timestep", padding=6)
        box.pack(fill="x", pady=(0, 6))

        row = ttk.Frame(box); row.pack(fill="x")
        ttk.Label(row, text="date", width=5).pack(side="left")
        self.cmb_date = ttk.Combobox(row, state="readonly", width=27,
                                     values=self.date_labels())
        self.cmb_date.pack(side="left", fill="x", expand=True)
        self.cmb_date.bind("<<ComboboxSelected>>",
                           lambda e: self._set_date(self.cmb_date.current()))

        row = ttk.Frame(box); row.pack(fill="x", pady=(3, 0))
        ttk.Label(row, text="hour", width=5).pack(side="left")
        self.cmb_hour = ttk.Combobox(row, state="readonly", width=27,
                                     values=self.hour_labels(self.date_key))
        self.cmb_hour.pack(side="left", fill="x", expand=True)
        self.cmb_hour.bind("<<ComboboxSelected>>",
                           lambda e: self._set_hour(self.cmb_hour.current()))

        self.scale_step = ttk.Scale(box, from_=0, command=self._on_scale)
        self.scale_step.pack(fill="x", pady=(4, 2))
        row = ttk.Frame(box); row.pack()
        ttk.Button(row, text="◀", width=4, command=lambda: self._step(-1)).pack(side="left")
        ttk.Label(row, text=" ← → hour · ↑ ↓ date ", foreground="#666666").pack(side="left")
        ttk.Button(row, text="▶", width=4, command=lambda: self._step(+1)).pack(side="left")
        self._sync_time_widgets()

        # -- heliostat ---------------------------------------------------
        box = ttk.LabelFrame(p, text="heliostat", padding=6)
        box.pack(fill="x", pady=(0, 6))
        row = ttk.Frame(box); row.pack(fill="x")
        ttk.Label(row, text="id").pack(side="left")
        ent = ttk.Entry(row, textvariable=self.var_helio, width=8)
        ent.pack(side="left", padx=4)
        ent.bind("<Return>", lambda e: self._select_typed())
        ttk.Button(row, text="Go", width=4, command=self._select_typed).pack(side="left")
        ttk.Button(row, text="none", width=6,
                   command=lambda: self._select(None)).pack(side="left", padx=(4, 0))
        # "min"/"max" rather than "worst"/"best": whether large is bad depends on
        # the metric (big power_w is good, big r90_mm is not), and a label that
        # guesses would quietly select the opposite of what was meant.
        row = ttk.Frame(box); row.pack(fill="x", pady=(4, 0))
        ttk.Label(row, text="jump to").pack(side="left")
        ttk.Button(row, text="min", width=6, command=lambda: self._select_extreme("min")
                   ).pack(side="left", padx=(4, 0))
        ttk.Button(row, text="max", width=6, command=lambda: self._select_extreme("max")
                   ).pack(side="left", padx=4)
        ttk.Label(row, text="of colour metric", foreground="#666666").pack(side="left")

        # -- display -----------------------------------------------------
        box = ttk.LabelFrame(p, text="display", padding=6)
        box.pack(fill="x", pady=(0, 6))
        for label, var, values, cb in (
            ("colour by", self.var_colour, self.columns, self._on_colour),
            ("curve", self.var_metric, self.columns, self._on_metric),
        ):
            row = ttk.Frame(box); row.pack(fill="x", pady=1)
            ttk.Label(row, text=label, width=9).pack(side="left")
            c = ttk.Combobox(row, textvariable=var, values=values,
                             state="readonly", width=21)
            c.pack(side="left", fill="x", expand=True)
            c.bind("<<ComboboxSelected>>", cb)

        row = ttk.Frame(box); row.pack(fill="x", pady=(4, 1))
        ttk.Label(row, text="DNI W/m²", width=9).pack(side="left")
        e = ttk.Entry(row, textvariable=self.var_dni, width=9)
        e.pack(side="left")
        e.bind("<Return>", lambda ev: self._on_weights())
        ttk.Label(row, text="aperture r mm").pack(side="left", padx=(8, 2))
        e = ttk.Entry(row, textvariable=self.var_aperture, width=7)
        e.pack(side="left")
        e.bind("<Return>", lambda ev: self._on_weights())

        row = ttk.Frame(box); row.pack(fill="x", pady=(4, 1))
        ttk.Label(row, text="spot bins", width=9).pack(side="left")
        stored = self.cfg.receiver.grid_size
        sizes = sorted({stored // 8, stored // 4, stored // 2, stored,
                        stored * 2, stored * 4} - {0})
        c = ttk.Combobox(row, textvariable=self.var_bins, state="readonly", width=6,
                         values=[str(s) for s in sizes])
        c.pack(side="left")
        c.bind("<<ComboboxSelected>>", lambda ev: self._on_bins())
        self.lbl_bins = ttk.Label(row, text="", foreground="#666666")
        self.lbl_bins.pack(side="left", padx=(6, 0))

        row = ttk.Frame(box); row.pack(fill="x", pady=(4, 1))
        ttk.Label(row, text="spot view", width=9).pack(side="left")
        for text, value in (("image", "image"), ("encircled energy", "encircled")):
            ttk.Radiobutton(row, text=text, value=value, variable=self.var_spotview,
                            command=self._on_spotview).pack(side="left", padx=(0, 6))

        row = ttk.Frame(box); row.pack(fill="x", pady=(4, 1))
        ttk.Label(row, text="field", width=9).pack(side="left")
        ttk.Checkbutton(row, text="shadows", variable=self.var_show_shadow,
                        command=self._on_overlay).pack(side="left")
        ttk.Checkbutton(row, text="blocking", variable=self.var_show_block,
                        command=self._on_overlay).pack(side="left", padx=(6, 0))

        ttk.Checkbutton(box, text="apply shading × blocking (all views)",
                        variable=self.var_shading,
                        command=self._on_weights).pack(anchor="w", pady=(4, 0))
        ttk.Checkbutton(box, text="log colour scale on spots", variable=self.var_log,
                        command=self._refresh_all).pack(anchor="w")

        # -- quadoa ------------------------------------------------------
        box = ttk.LabelFrame(p, text="open in Quadoa", padding=6)
        box.pack(fill="x", pady=(0, 6))
        row = ttk.Frame(box); row.pack(fill="x")
        ttk.Button(row, text="Export selected", command=lambda: self._export(False)
                   ).pack(side="left", fill="x", expand=True)
        # Writing the occluders is pure text editing, so it needs no licence
        # seat -- which is why it is a separate button rather than an option on
        # the one above: it works while a sweep is running and that one does not.
        row = ttk.Frame(box); row.pack(fill="x", pady=(4, 0))
        ttk.Button(row, text="Export with shading + blocking geometry",
                   command=self._export_occluders).pack(fill="x", expand=True)
        row = ttk.Frame(box); row.pack(fill="x", pady=(4, 0))
        ttk.Label(row, text="or the").pack(side="left")
        ttk.Spinbox(row, from_=1, to=25, width=3, textvariable=self.var_nexport
                    ).pack(side="left", padx=3)
        ttk.Combobox(row, textvariable=self.var_extreme, values=["min", "max"],
                     state="readonly", width=4).pack(side="left")
        ttk.Button(row, text="by colour metric", command=lambda: self._export(True)
                   ).pack(side="left", fill="x", expand=True, padx=(3, 0))
        self.txt_log = tk.Text(box, height=7, width=36, font=("Consolas", 8),
                               background="#ffffff", relief="solid", borderwidth=1,
                               wrap="none", state="disabled")
        self.txt_log.pack(fill="x", pady=(6, 0))

        # -- readout -----------------------------------------------------
        box = ttk.LabelFrame(p, text="selected heliostat", padding=4)
        box.pack(fill="both", expand=True)
        self.tree_row = ttk.Treeview(box, columns=("v",), show="tree headings", height=16)
        self.tree_row.heading("#0", text="field")
        self.tree_row.heading("v", text="value")
        self.tree_row.column("#0", width=150, anchor="w")
        self.tree_row.column("v", width=130, anchor="e")
        self.tree_row.pack(fill="both", expand=True)

    def _build_table_tab(self) -> None:
        outer = ttk.Frame(self.book)
        self.book.add(outer, text="Table")
        # No figure here, so only the CSV half of the export bar -- but in the
        # same place, with the same wording, as on every plot tab.
        self._export_bar(outer, "Table", figure=False)
        frame = ttk.Frame(outer)
        frame.pack(fill="both", expand=True)
        cols = ["heliostat_id"] + [c for c in self.columns if c != "heliostat_id"]
        self.table_cols = cols[:14]
        self.tree = ttk.Treeview(frame, columns=self.table_cols, show="headings")
        for c in self.table_cols:
            self.tree.heading(c, text=c, command=lambda cc=c: self._sort_table(cc))
            self.tree.column(c, width=95, anchor="e", stretch=False)
        self.tree.column("heliostat_id", width=80)
        vsb = ttk.Scrollbar(frame, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        self.tree.bind("<<TreeviewSelect>>", self._on_table_select)
        self._sort_state = (None, False)

    # ------------------------------------------------------------------
    # Trace tab
    # ------------------------------------------------------------------
    def _build_trace_tab(self) -> None:
        """Set up a sweep: every option, the command it makes, and a Run button.

        Three columns, because the point is to see the whole specification at
        once: what to trace (left), what to trace it through and where to put it
        (middle), and the resulting command plus the live run (right).

        Deliberately independent of ``self.summary``: this tab is about the run
        that has not happened yet, so nothing here touches the store. (The window
        as a whole still needs a finished-or-running sweep to open -- see
        ``_load_run`` -- which is a limitation of the GUI, not of this tab.)
        """
        frame = ttk.Frame(self.book, padding=8)
        self.book.add(frame, text="Trace")
        frame.columnconfigure(2, weight=1)
        frame.rowconfigure(0, weight=1)

        col0 = ttk.Frame(frame); col0.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
        col1 = ttk.Frame(frame); col1.grid(row=0, column=1, sticky="nsew", padx=(0, 8))
        col2 = ttk.Frame(frame); col2.grid(row=0, column=2, sticky="nsew")
        col2.rowconfigure(2, weight=1)
        col2.columnconfigure(0, weight=1)

        cfg = self.cfg

        # -- variables ---------------------------------------------------
        self.var_t_newdate = tk.StringVar(value="")
        self.var_t_suggest = tk.StringVar(value="8")
        self.var_t_use_suggest = tk.BooleanVar(value=False)
        self.var_t_hour_step = tk.StringVar(value=_trace_num(cfg.sweep.hour_step))
        self.var_t_margin = tk.StringVar(value=_trace_num(cfg.sweep.sunrise_margin_min))
        self.var_t_all = tk.BooleanVar(value=True)
        self.var_t_rays = tk.StringVar(value=str(cfg.trace.rays_per_heliostat))
        self.var_t_rays_per_trace = tk.StringVar(value=str(cfg.trace.rays_per_trace))
        self.var_t_workers = tk.StringVar(value="1")
        self.var_t_unlock_workers = tk.BooleanVar(value=False)
        self.var_t_secondary = tk.StringVar(value=cfg.optics.secondary)
        self.var_t_focus = tk.StringVar(
            value="" if cfg.geometry.focus_height_mm is None
            else _trace_num(cfg.geometry.focus_height_mm))
        self.var_t_rim = tk.StringVar(
            value="" if cfg.geometry.secondary_rim_height_mm is None
            else _trace_num(cfg.geometry.secondary_rim_height_mm))
        self.var_t_nmirrors = tk.StringVar(value=str(cfg.optics.n_mirrors))
        self.var_t_flat = tk.BooleanVar(value=bool(cfg.optics.flat_mirrors))
        self.var_t_model = tk.StringVar(value=cfg.trace.model_file)
        self.var_t_occluders = tk.BooleanVar(value=False)
        self.var_t_output = tk.StringVar(value=self._trace_default_output())
        self.var_t_resume = tk.BooleanVar(value=False)
        self.var_t_state = tk.StringVar(value="idle")
        self.var_t_progress = tk.StringVar(value="")

        # -- dates -------------------------------------------------------
        box = ttk.LabelFrame(col0, text="dates", padding=6)
        box.pack(fill="x", pady=(0, 6))
        row = ttk.Frame(box); row.pack(fill="x")
        self.lst_t_dates = tk.Listbox(row, height=8, selectmode="extended",
                                      exportselection=False, font=("Consolas", 9))
        # exportselection=False or clicking anywhere else in the window clears it:
        # Tk hands the X selection to one widget at a time.
        sb = ttk.Scrollbar(row, orient="vertical", command=self.lst_t_dates.yview)
        self.lst_t_dates.configure(yscrollcommand=sb.set)
        self.lst_t_dates.pack(side="left", fill="x", expand=True)
        sb.pack(side="left", fill="y")
        for d in cfg.sweep.dates:
            self.lst_t_dates.insert("end", str(d))
        self.lst_t_dates.selection_set(0, "end")
        self.lst_t_dates.bind("<<ListboxSelect>>", lambda e: self._update_trace_command())

        row = ttk.Frame(box); row.pack(fill="x", pady=(4, 0))
        ttk.Button(row, text="all", width=5,
                   command=lambda: self._trace_select_dates(True)).pack(side="left")
        ttk.Button(row, text="none", width=6,
                   command=lambda: self._trace_select_dates(False)).pack(side="left", padx=3)
        ttk.Label(row, text="none = config.toml's dates", foreground="#666666"
                  ).pack(side="left", padx=(4, 0))

        row = ttk.Frame(box); row.pack(fill="x", pady=(4, 0))
        ttk.Label(row, text="add").pack(side="left")
        ent = ttk.Entry(row, textvariable=self.var_t_newdate, width=12)
        ent.pack(side="left", padx=3)
        ent.bind("<Return>", lambda e: self._trace_add_date())
        ttk.Button(row, text="+", width=3, command=self._trace_add_date).pack(side="left")
        ttk.Label(row, text="YYYY-MM-DD", foreground="#666666").pack(side="left", padx=(4, 0))

        row = ttk.Frame(box); row.pack(fill="x", pady=(4, 0))
        ttk.Checkbutton(row, text="instead, let energy.suggest_sweep_dates pick",
                        variable=self.var_t_use_suggest).pack(side="left")
        self.spin_t_suggest = ttk.Spinbox(row, from_=2, to=24, width=4,
                                          textvariable=self.var_t_suggest)
        self.spin_t_suggest.pack(side="left", padx=3)
        ttk.Label(box, text="declination-spaced dates, keeping config.toml's",
                  foreground="#666666").pack(anchor="w")

        # -- time grid ---------------------------------------------------
        box = ttk.LabelFrame(col0, text="time grid", padding=6)
        box.pack(fill="x", pady=(0, 6))
        row = ttk.Frame(box); row.pack(fill="x")
        ttk.Label(row, text="hour step", width=11).pack(side="left")
        ttk.Entry(row, textvariable=self.var_t_hour_step, width=7).pack(side="left")
        ttk.Label(row, text="h  MAXIMUM spacing", foreground="#666666"
                  ).pack(side="left", padx=(4, 0))
        row = ttk.Frame(box); row.pack(fill="x", pady=(3, 0))
        ttk.Label(row, text="sun margin", width=11).pack(side="left")
        ttk.Entry(row, textvariable=self.var_t_margin, width=7).pack(side="left")
        ttk.Label(row, text="min inside sunrise/sunset", foreground="#666666"
                  ).pack(side="left", padx=(4, 0))
        ttk.Label(box, text="the daylight window is divided into equal intervals no\n"
                            "wider than the step, so samples land on the true edges\n"
                            "rather than on a clock grid",
                  foreground="#666666").pack(anchor="w", pady=(3, 0))

        # -- heliostats --------------------------------------------------
        box = ttk.LabelFrame(col0, text="heliostats", padding=6)
        box.pack(fill="x", pady=(0, 6))
        ttk.Radiobutton(box, text="all 645 (--all-heliostats)", value=True,
                        variable=self.var_t_all).pack(anchor="w")
        ttk.Radiobutton(box, text=f"the downselect ({cfg.field.n_configs} from "
                                  f"{Path(cfg.field.downselect_file).name})",
                        value=False, variable=self.var_t_all).pack(anchor="w")

        # -- rays and workers --------------------------------------------
        box = ttk.LabelFrame(col0, text="rays and workers", padding=6)
        box.pack(fill="x", pady=(0, 6))
        row = ttk.Frame(box); row.pack(fill="x")
        ttk.Label(row, text="rays/heliostat", width=13).pack(side="left")
        ttk.Entry(row, textvariable=self.var_t_rays, width=9).pack(side="left")
        ttk.Label(row, text=f"config {cfg.trace.rays_per_heliostat:,}",
                  foreground="#666666").pack(side="left", padx=(4, 0))

        # The budget is split into chunks of this many rays, one traceRays call
        # each, so this entry -- not the one above -- is what sets the number of
        # iterations per heliostat.
        row = ttk.Frame(box); row.pack(fill="x", pady=(3, 0))
        ttk.Label(row, text="rays per call", width=13).pack(side="left")
        ttk.Entry(row, textvariable=self.var_t_rays_per_trace, width=9).pack(side="left")
        ttk.Label(row, text=f"config {cfg.trace.rays_per_trace:,}",
                  foreground="#666666").pack(side="left", padx=(4, 0))
        # Read-only and derived: the count of setRayDistributionCount1 +
        # traceRays + getRayPos round trips is the thing being chosen here, and
        # it is arithmetic on the two entries above rather than a third input.
        self.lbl_t_chunks = ttk.Label(box, text="", foreground="#336699",
                                      justify="left")
        self.lbl_t_chunks.pack(anchor="w", pady=(2, 0))
        _Tooltip(self.lbl_t_chunks,
                 "Each heliostat is ceil(rays/heliostat / rays per call) calls to\n"
                 "setRayDistributionCount1 + traceRays + getRayPos. The chunks sum\n"
                 "exactly to the ray budget; an uneven budget gets a short last one.\n\n"
                 "Whether the fixed per-call cost or the per-ray cost dominates has\n"
                 "NOT been measured -- scripts/probe_ray_cost.py is the experiment.")

        row = ttk.Frame(box); row.pack(fill="x", pady=(3, 0))
        ttk.Label(row, text="workers", width=13).pack(side="left")
        self.spin_t_workers = ttk.Spinbox(row, from_=1, to=1, width=4,
                                          textvariable=self.var_t_workers)
        self.spin_t_workers.pack(side="left")
        ttk.Checkbutton(row, text="allow > 1", variable=self.var_t_unlock_workers,
                        command=self._trace_unlock_workers).pack(side="left", padx=(6, 0))
        lbl = ttk.Label(box, text="capped at 1: ONE licence seat, and a failed second\n"
                                  "request leaks the first (measured)",
                        foreground="#a33")
        lbl.pack(anchor="w", pady=(3, 0))
        _Tooltip(lbl, WORKERS_WARNING)
        _Tooltip(self.spin_t_workers, WORKERS_WARNING)

        # -- secondary layout --------------------------------------------
        box = ttk.LabelFrame(col1, text="secondary layout", padding=6)
        box.pack(fill="x", pady=(0, 6))
        row = ttk.Frame(box); row.pack(fill="x")
        ttk.Label(row, text="secondary", width=12).pack(side="left")
        cmb = ttk.Combobox(row, textvariable=self.var_t_secondary, state="readonly",
                           width=13, values=list(SECONDARY_LAYOUTS))
        cmb.pack(side="left")
        cmb.bind("<<ComboboxSelected>>", lambda e: self._on_trace_layout())
        row = ttk.Frame(box); row.pack(fill="x", pady=(3, 0))
        ttk.Label(row, text="focus height", width=12).pack(side="left")
        self.ent_t_focus = ttk.Entry(row, textvariable=self.var_t_focus, width=10)
        self.ent_t_focus.pack(side="left")
        ttk.Label(row, text="mm  F1 on the axis", foreground="#666666"
                  ).pack(side="left", padx=(4, 0))
        row = ttk.Frame(box); row.pack(fill="x", pady=(3, 0))
        ttk.Label(row, text="rim height", width=12).pack(side="left")
        self.ent_t_rim = ttk.Entry(row, textvariable=self.var_t_rim, width=10)
        self.ent_t_rim.pack(side="left")
        ttk.Label(row, text="mm  hyperboloid rim", foreground="#666666"
                  ).pack(side="left", padx=(4, 0))
        row = ttk.Frame(box); row.pack(fill="x", pady=(3, 0))
        ttk.Label(row, text="n_mirrors", width=12).pack(side="left")
        ttk.Combobox(row, textvariable=self.var_t_nmirrors, state="readonly", width=4,
                     values=["1", "2"]).pack(side="left")
        ttk.Label(row, text="reflections -> throughput", foreground="#666666"
                  ).pack(side="left", padx=(4, 0))
        # Flat vs focused heliostats. In this box because it is the other axis of
        # the same comparison, but it is NOT a layout: it composes with all three
        # of them, which is why it is a checkbox beside the combobox rather than
        # a fourth entry inside it.
        ttk.Checkbutton(box, text="flat heliostats (no focusing)",
                        variable=self.var_t_flat,
                        command=self._update_trace_command).pack(anchor="w", pady=(4, 0))
        ttk.Label(box, text="Keeps the pointing, zeroes the mirror curvature.\n"
                            "Expect a much larger spot, much more spillage and\n"
                            "much less collected energy — that is the comparison,\n"
                            "not a fault.",
                  foreground="#666666", justify="left").pack(anchor="w")
        self.lbl_t_layout = ttk.Label(box, text="", foreground="#666666")
        self.lbl_t_layout.pack(anchor="w", pady=(3, 0))

        # -- model -------------------------------------------------------
        box = ttk.LabelFrame(col1, text="model", padding=6)
        box.pack(fill="x", pady=(0, 6))
        row = ttk.Frame(box); row.pack(fill="x")
        self.ent_t_model = ttk.Entry(row, textvariable=self.var_t_model, width=34)
        self.ent_t_model.pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="…", width=3, command=self._trace_pick_model).pack(side="left")
        ttk.Checkbutton(box, text="trace occlusion as real geometry (--occluders)",
                        variable=self.var_t_occluders,
                        command=self._on_trace_occluders).pack(anchor="w", pady=(3, 0))
        self.lbl_t_model = ttk.Label(box, text="", foreground="#666666")
        self.lbl_t_model.pack(anchor="w")

        # -- output ------------------------------------------------------
        box = ttk.LabelFrame(col1, text="output", padding=6)
        box.pack(fill="x", pady=(0, 6))
        row = ttk.Frame(box); row.pack(fill="x")
        ttk.Entry(row, textvariable=self.var_t_output, width=34).pack(
            side="left", fill="x", expand=True)
        ttk.Checkbutton(box, text="resume: keep the timesteps an existing run "
                                  "directory already has",
                        variable=self.var_t_resume).pack(anchor="w", pady=(3, 0))
        self.lbl_t_paths = ttk.Label(box, text="", foreground="#666666", justify="left")
        self.lbl_t_paths.pack(anchor="w")

        # -- estimate ----------------------------------------------------
        box = ttk.LabelFrame(col1, text="estimate", padding=6)
        box.pack(fill="x")
        self.lbl_t_estimate = ttk.Label(box, text="", justify="left")
        self.lbl_t_estimate.pack(anchor="w")

        # -- command -----------------------------------------------------
        box = ttk.LabelFrame(col2, text="command", padding=6)
        box.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        box.columnconfigure(0, weight=1)
        # Narrow on purpose: the column has grid weight, so it grows into
        # whatever is left over, while a wide request would squeeze the form.
        self.txt_t_cmd = tk.Text(box, height=7, width=46, wrap="word",
                                 font=("Consolas", 9), background="#ffffff",
                                 relief="solid", borderwidth=1)
        self.txt_t_cmd.grid(row=0, column=0, sticky="ew")
        row = ttk.Frame(box); row.grid(row=1, column=0, sticky="ew", pady=(4, 0))
        ttk.Button(row, text="Copy", width=8, command=self._trace_copy).pack(side="left")
        ttk.Label(row, text="read-only; runs from the repository root. Flags appear "
                            "only where\nthey differ from config.toml, except "
                            "--output/--rays/--workers.",
                  foreground="#666666", justify="left").pack(side="left", padx=(6, 0))

        # -- launch ------------------------------------------------------
        box = ttk.LabelFrame(col2, text="launch", padding=6)
        box.grid(row=1, column=0, sticky="ew", pady=(0, 6))
        row = ttk.Frame(box); row.pack(fill="x")
        self.btn_t_run = ttk.Button(row, text="Run", width=8, command=self._run_trace)
        self.btn_t_run.pack(side="left")
        self.btn_t_stop = ttk.Button(row, text="Stop", width=8, state="disabled",
                                     command=self._stop_trace)
        self.btn_t_stop.pack(side="left", padx=4)
        ttk.Label(row, textvariable=self.var_t_state, font=("", 9, "bold")
                  ).pack(side="left", padx=(6, 0))
        ttk.Label(box, justify="left", foreground="#666666",
                  text="The run is detached: closing this window leaves it going, and "
                       "leaves its\nlock directory in place. Run refuses while any "
                       "analysis_output/.*.lock exists,\nbecause that means something "
                       "already holds the one licence seat.").pack(anchor="w", pady=(4, 0))
        self.lbl_t_locks = ttk.Label(box, text="", foreground="#a33", justify="left")
        self.lbl_t_locks.pack(anchor="w")

        # -- monitor -----------------------------------------------------
        box = ttk.LabelFrame(col2, text="monitor", padding=6)
        box.grid(row=2, column=0, sticky="nsew")
        box.rowconfigure(2, weight=1)
        box.columnconfigure(0, weight=1)
        ttk.Label(box, textvariable=self.var_t_progress).grid(row=0, column=0, sticky="w")
        self.bar_t_progress = ttk.Progressbar(box, maximum=100.0)
        self.bar_t_progress.grid(row=1, column=0, sticky="ew", pady=(2, 4))
        self.txt_t_log = tk.Text(box, height=14, width=46, wrap="none",
                                 font=("Consolas", 8), background="#ffffff",
                                 relief="solid", borderwidth=1)
        vsb = ttk.Scrollbar(box, orient="vertical", command=self.txt_t_log.yview)
        hsb = ttk.Scrollbar(box, orient="horizontal", command=self.txt_t_log.xview)
        self.txt_t_log.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.txt_t_log.grid(row=2, column=0, sticky="nsew")
        vsb.grid(row=2, column=1, sticky="ns")
        # The sweep's progress lines run to ~95 characters; unwrapped and without
        # this the end of every one of them -- which is where the ETA is -- would
        # be off the edge of a pane that shares the window with the form.
        hsb.grid(row=3, column=0, sticky="ew")

        # Live preview: one callback on every variable, registered only now that
        # every widget exists, so building the form does not fire it half-built.
        for var in (self.var_t_suggest, self.var_t_use_suggest, self.var_t_hour_step,
                    self.var_t_margin, self.var_t_all, self.var_t_rays,
                    self.var_t_rays_per_trace,
                    self.var_t_workers, self.var_t_secondary, self.var_t_focus,
                    self.var_t_rim, self.var_t_nmirrors, self.var_t_model,
                    self.var_t_occluders, self.var_t_output, self.var_t_resume):
            var.trace_add("write", lambda *_: self._update_trace_command())

        self._on_trace_layout(update=False)
        self._on_trace_occluders(update=False)
        self._update_trace_command()
        self._poll_trace()

    # -- Trace tab: form -----------------------------------------------
    def _trace_default_output(self) -> str:
        """A run directory beside the one being viewed, with a name nothing uses.

        Beside it, not inside it: the sweep being explored is the neighbour of the
        one about to run, and the log and lock live one level up from both.
        """
        parent = self.store.root.parent
        rel = parent / "trace1"
        for n in range(1, 100):
            candidate = parent / f"trace{n}"
            if not candidate.exists() and not trace_log_path(candidate).exists():
                rel = candidate
                break
        try:
            return str(rel.relative_to(self.cfg.repo_root)).replace("\\", "/")
        except ValueError:
            return str(rel)

    def _trace_select_dates(self, everything: bool) -> None:
        if everything:
            self.lst_t_dates.selection_set(0, "end")
        else:
            self.lst_t_dates.selection_clear(0, "end")
        self._update_trace_command()

    def _trace_add_date(self) -> None:
        """Add a free-text ISO date to the list and select it."""
        text = self.var_t_newdate.get().strip()
        if not text:
            return
        try:
            date = _dt.date.fromisoformat(text)
        except ValueError:
            self._status(f"{text!r} is not an ISO date (YYYY-MM-DD)")
            return
        existing = list(self.lst_t_dates.get(0, "end"))
        if str(date) not in existing:
            self.lst_t_dates.insert("end", str(date))
            existing.append(str(date))
        self.lst_t_dates.selection_set(existing.index(str(date)))
        self.var_t_newdate.set("")
        self._update_trace_command()

    def _trace_unlock_workers(self) -> None:
        """Raise the worker cap only when explicitly unlocked -- see WORKERS_WARNING."""
        if self.var_t_unlock_workers.get():
            self.spin_t_workers.config(to=4)
        else:
            self.spin_t_workers.config(to=1)
            self.var_t_workers.set("1")
        self._update_trace_command()

    def _on_trace_layout(self, update: bool = True) -> None:
        """Enable the height entries the chosen layout actually needs.

        Also follows ``n_mirrors``: prime focus has one reflection and the other
        two have two, and config.toml only *warns* on a mismatch, so a form that
        left it alone would quietly produce a command line that misreports
        throughput.
        """
        layout = self.var_t_secondary.get()
        needs_focus = layout in ("prime_focus", "cassegrain")
        needs_rim = layout == "cassegrain"
        self.ent_t_focus.config(state="normal" if needs_focus else "disabled")
        self.ent_t_rim.config(state="normal" if needs_rim else "disabled")

        n = reflections_for(layout)
        if n is not None:
            self.var_t_nmirrors.set(str(n))

        note = {
            "axicon": "conical secondary; the aim point comes from each heliostat's\n"
                      "radius, so there is no single focus height",
            "prime_focus": "one reflection onto a detector at F1. REQUIRES focus height.\n"
                           "The shipped .optx has that height as a literal -- see README",
            "cassegrain": "hyperboloid relaying F1 to the receiver. REQUIRES focus and\n"
                          "rim heights; build the hyperboloid in Quadoa by hand",
        }.get(layout, "")
        self.lbl_t_layout.config(text=note)
        if update:
            self._update_trace_command()

    def _on_trace_occluders(self, update: bool = True) -> None:
        """--occluders brings its own model, so the model box goes read-only."""
        on = self.var_t_occluders.get()
        self.ent_t_model.config(state="disabled" if on else "normal")
        self.lbl_t_model.config(
            text=("--occluders traces models/heliostat_field_occluders.optx;\n"
                  "the CLI refuses --model-file alongside it"
                  if on else
                  "the .optx to trace, relative to the repository root"))
        if update:
            self._update_trace_command()

    def _trace_pick_model(self) -> None:
        path = filedialog.askopenfilename(
            title="Quadoa model", initialdir=str(self.cfg.repo_root / "models"),
            filetypes=[("Quadoa model", "*.optx"), ("all files", "*.*")])
        if not path:
            return
        try:
            path = str(Path(path).relative_to(self.cfg.repo_root)).replace("\\", "/")
        except ValueError:
            pass  # outside the repo: keep it absolute
        self.var_t_model.set(path)

    # -- Trace tab: the command ----------------------------------------
    def _trace_number(self, var, cast, name, problems, default=None):
        """One form field as a number, collecting its own error message."""
        text = var.get().strip()
        if not text:
            return default
        try:
            return cast(text)
        except ValueError:
            problems.append(f"{name}: {text!r} is not a number")
            return default

    def _trace_options(self) -> tuple[dict, list]:
        """The form as an options dict, plus everything wrong with it.

        Every value that matches config.toml is left out, so the command line
        stays a list of decisions rather than a restatement of the file. The
        exceptions are ``--output``, ``--rays`` and ``--workers``: where a run
        lands, what it costs and whether it fights for the licence seat are worth
        stating even when they happen to agree with the file.
        """
        cfg = self.cfg
        problems: list[str] = []
        opts: dict = {"config": self.config_path}

        selected = [self.lst_t_dates.get(i) for i in self.lst_t_dates.curselection()]
        if self.var_t_use_suggest.get():
            opts["suggest_dates"] = self._trace_number(
                self.var_t_suggest, int, "suggested dates", problems, 8)
        else:
            configured = [str(d) for d in cfg.sweep.dates]
            # An unchanged full selection is exactly what config.toml already
            # says, and spelling out twelve dates then adds noise, not
            # information. Any other selection is explicit.
            opts["dates"] = [] if selected == configured else selected

        opts["all_heliostats"] = bool(self.var_t_all.get())
        opts["rays"] = self._trace_number(self.var_t_rays, int, "rays", problems,
                                          cfg.trace.rays_per_heliostat)
        opts["workers"] = self._trace_number(self.var_t_workers, int, "workers",
                                            problems, 1)
        if opts["rays"] is not None and opts["rays"] <= 0:
            problems.append("rays must be positive")

        # The chunk size follows the ordinary rule -- emitted only where it
        # departs from config.toml -- because unlike the budget it does not
        # change what the run means, only how many round trips it takes to get
        # there. The CLI's interaction rule is mirrored exactly: with the flag
        # absent the chunk is clamped to the budget, and with it present a chunk
        # larger than the budget is refused rather than clamped.
        per_trace = self._trace_number(self.var_t_rays_per_trace, int,
                                       "rays per call", problems,
                                       cfg.trace.rays_per_trace)
        if per_trace is not None:
            if per_trace <= 0:
                problems.append("rays per call must be positive")
            elif per_trace != cfg.trace.rays_per_trace:
                if opts["rays"] and per_trace > opts["rays"]:
                    problems.append(
                        f"rays per call ({per_trace:,}) exceeds rays/heliostat "
                        f"({opts['rays']:,}); a chunk cannot be bigger than the "
                        f"budget. Equal means one call per heliostat")
                else:
                    opts["rays_per_trace"] = per_trace
        if opts["workers"] is not None and not 1 <= opts["workers"] <= 4:
            problems.append("workers must be 1-4 (the HASP key's limit)")

        step = self._trace_number(self.var_t_hour_step, float, "hour step", problems)
        if step is not None and step != cfg.sweep.hour_step:
            opts["hour_step"] = step
        margin = self._trace_number(self.var_t_margin, float, "sun margin", problems)
        if margin is not None and margin != cfg.sweep.sunrise_margin_min:
            opts["sunrise_margin_min"] = margin

        layout = self.var_t_secondary.get()
        if layout != cfg.optics.secondary:
            opts["secondary"] = layout
        focus = self._trace_number(self.var_t_focus, float, "focus height", problems)
        rim = self._trace_number(self.var_t_rim, float, "rim height", problems)
        if layout in ("prime_focus", "cassegrain"):
            if focus is None:
                problems.append(f"{layout} needs a focus height (F1 on the axis)")
            elif focus != cfg.geometry.focus_height_mm:
                opts["focus_height_mm"] = focus
        if layout == "cassegrain":
            if rim is None:
                problems.append("cassegrain needs a rim height (the shadow circle)")
            elif rim != cfg.geometry.secondary_rim_height_mm:
                opts["rim_height_mm"] = rim

        n_mirrors = self._trace_number(self.var_t_nmirrors, int, "n_mirrors", problems)
        if n_mirrors is not None and n_mirrors != cfg.optics.n_mirrors:
            opts["n_mirrors"] = n_mirrors

        # The tab's ordinary rule, applied to a boolean: say nothing where the
        # form agrees with config.toml, and say it explicitly in either direction
        # where it does not.
        flat = bool(self.var_t_flat.get())
        if flat != bool(cfg.optics.flat_mirrors):
            opts["flat_mirrors"] = flat

        opts["occluders"] = bool(self.var_t_occluders.get())
        model = self.var_t_model.get().strip()
        if model and model != cfg.trace.model_file:
            opts["model_file"] = model

        output = self.var_t_output.get().strip()
        if not output:
            problems.append("no output directory")
        opts["output"] = output
        opts["resume"] = bool(self.var_t_resume.get())
        return opts, problems

    def _trace_output_dir(self, opts) -> Path:
        """The run directory as an absolute path, resolved like the CLI does."""
        out = Path(opts.get("output") or "unnamed")
        return out if out.is_absolute() else self.cfg.repo_root / out

    def _update_trace_command(self, *_args) -> None:
        """Rebuild the preview, the paths, the estimate and the lock warning.

        Called on every keystroke, so everything in it must be cheap: the
        estimate is a few hundred sun positions, and the lock scan is a glob over
        one directory.
        """
        opts, problems = self._trace_options()
        argv = build_sweep_argv(opts)
        self._trace_argv = argv

        text = format_command(argv)
        if problems:
            text += "\n\nwill not run:\n" + "\n".join(f"  - {p}" for p in problems)
        self.txt_t_cmd.config(state="normal")
        self.txt_t_cmd.delete("1.0", "end")
        self.txt_t_cmd.insert("1.0", text)
        # Read-only, but only after the insert: a disabled Text rejects writes
        # from code as well as from the keyboard.
        self.txt_t_cmd.config(state="disabled")

        # Live, and from the same arithmetic the trace will use. Falls back to
        # config.toml's values so the label still reads while an entry is being
        # retyped and is momentarily blank.
        cfg = self.cfg
        self.lbl_t_chunks.config(text=describe_call_plan(
            opts.get("rays") or cfg.trace.rays_per_heliostat,
            opts.get("rays_per_trace") or cfg.trace.rays_per_trace))

        out = self._trace_output_dir(opts)
        self.lbl_t_paths.config(
            text=f"log  {trace_log_path(out)}\nlock {trace_lock_dir(out)}")
        self.lbl_t_estimate.config(text=self._trace_estimate(opts, problems))
        self._show_trace_locks()
        if not problems and (self._trace_proc is None or self._trace_proc.poll() is not None):
            self.btn_t_run.config(state="normal")
        else:
            self.btn_t_run.config(state="disabled")

    def _trace_estimate(self, opts, problems) -> str:
        """Traces, wall clock and disk, from the same arithmetic as `beamdown info`.

        Built on a *copy* of the config with the overrides applied, so working out
        what a candidate time grid costs cannot disturb the config the rest of the
        window is reading.
        """
        if problems:
            return "fix the problems above"
        # Keyed on exactly the inputs the estimate depends on, because this runs
        # on every keystroke anywhere in the form: rebuilding the time grid is
        # ~25 ms, which is enough to be felt while typing a directory name that
        # cannot possibly change the answer.
        key = (tuple(opts.get("dates") or ()), opts.get("suggest_dates"),
               opts.get("hour_step"), opts.get("sunrise_margin_min"),
               opts.get("all_heliostats"), opts.get("rays"))
        cached = self._trace_estimate_cache.get(key)
        if cached is not None:
            return cached

        import copy

        from . import solar as solar_mod
        from .config import apply_overrides

        cfg = copy.deepcopy(self.cfg)
        overrides = {}
        if opts.get("hour_step") is not None:
            overrides.setdefault("sweep", {})["hour_step"] = opts["hour_step"]
        if opts.get("sunrise_margin_min") is not None:
            overrides.setdefault("sweep", {})["sunrise_margin_min"] = opts["sunrise_margin_min"]
        apply_overrides(cfg, overrides)

        dates = None
        if opts.get("dates"):
            try:
                dates = [_dt.date.fromisoformat(d) for d in opts["dates"]]
            except ValueError:
                return "unparseable date in the list"
        elif opts.get("suggest_dates"):
            # suggest_sweep_dates needs the DNI/declination machinery; the count
            # is all the estimate needs and the grid per date is what varies.
            dates = list(cfg.sweep.dates)[:int(opts["suggest_dates"])] or None

        try:
            steps = solar_mod.build_time_grid(cfg, dates)
        except Exception as exc:                     # a nonsense date, say
            return f"cannot build the time grid: {type(exc).__name__}: {exc}"

        n_helio = 645 if opts.get("all_heliostats") else cfg.field.n_configs
        n_steps = len(steps)
        traces = n_helio * n_steps
        hours = traces * TRACE_SECONDS_PER_TRACE / 3600.0
        rays = int(opts.get("rays") or cfg.trace.rays_per_heliostat)
        raw_gb = traces * rays * 0.375 * 2 * 2 / 1e9   # ~37.5% of rays land, int16 x/y
        flux_gb = traces * cfg.receiver.grid_size ** 2 * 4 / 1e9
        note = "" if opts.get("suggest_dates") is None else (
            f"\n(estimated over {len(dates or ())} of config.toml's dates; "
            f"suggest_sweep_dates picks {opts['suggest_dates']} of its own)")
        text = (f"{n_helio} heliostats x {n_steps} timesteps = {traces:,} traces{note}\n"
                f"~{hours:.1f} h at {TRACE_SECONDS_PER_TRACE} s/trace, "
                f"{rays:,} rays each\n"
                f"~{raw_gb:.1f} GB raw rays + {flux_gb:.1f} GB flux maps")
        self._trace_estimate_cache[key] = text
        return text

    def _trace_copy(self) -> None:
        self.root.clipboard_clear()
        self.root.clipboard_append(format_command(self._trace_argv))
        self._status("command copied to the clipboard")

    # -- Trace tab: locks, launching, monitoring ------------------------
    def _trace_lock_roots(self) -> list:
        """Where to look for sweep locks.

        Both the repository's ``analysis_output`` and whatever directory this run
        would land in: the licence seat is global, so a lock held anywhere matters,
        and a run pointed somewhere else must still be able to see its own.
        An attribute rather than a constant so a test can point it at a tmp path.
        """
        roots = [self.cfg.repo_root / "analysis_output"]
        opts, _problems = self._trace_options()
        roots.append(self._trace_output_dir(opts).parent)
        return roots

    def _show_trace_locks(self) -> None:
        locks = scan_locks(self._trace_lock_roots())
        if not locks:
            self.lbl_t_locks.config(text="")
            return
        self.lbl_t_locks.config(
            text="lock held: " + ", ".join(f"{d['path'].name} (pid {d['pid']})"
                                           for d in locks))

    def _run_trace(self) -> None:
        """Launch the previewed command, or explain why not.

        Everything that can refuse does so before the process starts: a launched
        sweep takes the licence seat, and taking it from a run that is 20 hours in
        is not recoverable.
        """
        import subprocess

        if self._trace_proc is not None and self._trace_proc.poll() is None:
            self._status("a sweep launched from this window is still running")
            return

        opts, problems = self._trace_options()
        if problems:
            messagebox.showwarning("not starting", "\n".join(problems))
            return

        out = self._trace_output_dir(opts)
        refusal = launch_refusal(out, opts["resume"], scan_locks(self._trace_lock_roots()))
        if refusal:
            messagebox.showwarning("not starting", refusal)
            self._status(refusal.splitlines()[0])
            return

        lock_dir = trace_lock_dir(out)
        try:
            lock_dir.mkdir()          # atomic claim, exactly as the run scripts do
        except OSError as exc:
            messagebox.showwarning("not starting", f"cannot take {lock_dir}: {exc}")
            return

        argv = build_sweep_argv(opts)
        log_path = trace_log_path(out)
        try:
            log = open(log_path, "a", encoding="utf-8", errors="replace")
            log.write(f"\n[{_dt.datetime.now():%Y-%m-%d %H:%M:%S}] launched from the "
                      f"beamdown GUI\n{format_command(argv)}\n")
            log.flush()
            proc = subprocess.Popen(
                argv, cwd=str(self.cfg.repo_root), stdout=log,
                stderr=subprocess.STDOUT, **detached_popen_kwargs())
        except Exception as exc:
            # The lock must not outlive a launch that never happened.
            self._release_trace_lock(lock_dir)
            messagebox.showerror("could not launch", f"{type(exc).__name__}: {exc}")
            return
        finally:
            try:
                log.close()   # the child holds its own duplicate of the handle
            except Exception:
                pass

        (lock_dir / "pid").write_text(str(proc.pid), encoding="utf-8")
        self._trace_proc = proc
        self._trace_log_path = log_path
        self._trace_lock_dir = lock_dir
        self.var_t_state.set(f"running, pid {proc.pid}")
        self.btn_t_run.config(state="disabled")
        self.btn_t_stop.config(state="normal")
        self._log(f"sweep launched, pid {proc.pid} -> {log_path.name}")
        self._status(f"sweep running as pid {proc.pid}; watch the monitor pane")
        self._poll_trace(reschedule=False)

    def _stop_trace(self) -> None:
        """Kill the launched sweep and its children, and say what survives.

        There is no graceful stop to offer: the sweep is a detached console
        process with no signal handler, and on Windows nothing short of taskkill
        reaches it. What makes that acceptable is resume -- ``run_sweep`` skips
        any timestep already in the store (``store.has_timestep``), so only the
        timestep in flight is lost.
        """
        import subprocess

        proc = self._trace_proc
        if proc is None or proc.poll() is not None:
            self._status("no sweep from this window is running")
            return
        if not messagebox.askyesno(
                "stop the sweep?",
                f"Kill pid {proc.pid} and its workers?\n\n"
                f"Timesteps already written are kept -- relaunch with 'resume' "
                f"ticked and it continues from there. The timestep in flight is "
                f"lost, and the licence seat takes a moment to be released."):
            return
        try:
            if sys.platform == "win32":
                # /T for the tree: workers > 1 starts a real process pool.
                subprocess.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                               capture_output=True, check=False)
            else:
                import os
                import signal

                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        except Exception as exc:
            self._status(f"could not stop pid {proc.pid}: {exc}")
            return
        self._status(f"stopped pid {proc.pid}; resume will continue from the last "
                     f"completed timestep")
        self._poll_trace(reschedule=False)

    def _release_trace_lock(self, lock_dir=None) -> None:
        """Drop the lock directory, pid file and all (``rm -rf``, as the scripts)."""
        import shutil

        lock_dir = lock_dir or self._trace_lock_dir
        if lock_dir is None:
            return
        shutil.rmtree(lock_dir, ignore_errors=True)
        if lock_dir == self._trace_lock_dir:
            self._trace_lock_dir = None

    def _poll_trace(self, reschedule: bool = True) -> None:
        """Tail the log and update progress, about every two seconds.

        Its own ``after`` loop rather than a thread: reading the tail of a file is
        a millisecond and Tk variables may only be touched from this thread
        anyway. Runs whether or not a sweep is going, so the lock warning stays
        current when someone else starts one.
        """
        try:
            proc = self._trace_proc
            if self._trace_log_path is not None and self._trace_log_path.exists():
                text = tail_file(self._trace_log_path)
                if text != self.txt_t_log.get("1.0", "end-1c"):
                    self.txt_t_log.delete("1.0", "end")
                    self.txt_t_log.insert("1.0", text)
                    self.txt_t_log.see("end")
                info = parse_progress(text)
                if info.get("total"):
                    frac = 100.0 * info["done"] / info["total"]
                    self.bar_t_progress.config(value=frac)
                    eta = (f", eta {info['eta_min']:.0f} min"
                           if "eta_min" in info else "")
                    self.var_t_progress.set(
                        f"timestep {info['done']}/{info['total']}  ({frac:.0f}%){eta}")
                elif info.get("complete"):
                    self.var_t_progress.set("complete")

            if proc is not None:
                code = proc.poll()
                if code is None:
                    self.var_t_state.set(f"running, pid {proc.pid}")
                else:
                    self.var_t_state.set(
                        f"finished, exit {code}" if code == 0 else f"exited {code}")
                    self._release_trace_lock()
                    self._trace_proc = None
                    self.btn_t_stop.config(state="disabled")
                    self._update_trace_command()   # re-enables Run
            elif self._trace_lock_dir is None:
                self._show_trace_locks()
        except tk.TclError:
            return                                  # window going away
        finally:
            if reschedule:
                self.root.after(2000, self._poll_trace)

    def _on_close(self) -> None:
        """Closing the window leaves a launched sweep -- and its lock -- alone."""
        proc = self._trace_proc
        if proc is not None and proc.poll() is None:
            if not messagebox.askokcancel(
                    "a sweep is running",
                    f"pid {proc.pid} keeps going after this window closes, and its "
                    f"lock directory\n{self._trace_lock_dir}\nstays in place until "
                    f"it finishes -- so nothing else grabs the licence seat.\n\n"
                    f"Its log is {self._trace_log_path}. Close the window?"):
                return
        self.root.destroy()

    # ------------------------------------------------------------------
    # events
    # ------------------------------------------------------------------
    def _on_scale(self, value) -> None:
        """The slider runs over hours within the selected date, not over the run.

        Ignored while the widgets are being synced in code: ``Scale.set`` and
        ``Scale.config(to=...)`` both fire this command, so a date change that
        shortens the slider would otherwise re-enter and land on a stale hour.
        """
        if self._syncing:
            return
        hours = self.keys_by_date[self.date_key]
        i = int(np.clip(int(float(value)), 0, len(hours) - 1))
        if hours[i] != self.key:
            self._set_step(self.keys.index(hours[i]), move_scale=False)

    def _sync_time_widgets(self) -> None:
        """Point the date/hour selectors and the slider at the current key."""
        hours = self.keys_by_date[self.date_key]
        self._syncing = True
        try:
            self.cmb_date.current(self.dates.index(self.date_key))
            self.cmb_hour.config(values=self.hour_labels(self.date_key))
            self.cmb_hour.current(hours.index(self.key))
            self.scale_step.config(to=max(1, len(hours) - 1))
            self.scale_step.set(hours.index(self.key))
        finally:
            self._syncing = False

    def _set_step(self, i: int, move_scale: bool = True) -> None:
        self.step_i = int(np.clip(i, 0, len(self.keys) - 1))
        self._sync_time_widgets()
        self._refresh_all()

    def _set_date(self, i: int) -> None:
        """Change date, keeping the nearest hour so the view stays comparable."""
        date_key = self.dates[int(np.clip(i, 0, len(self.dates) - 1))]
        want = self._hour_of(self.key)
        best = min(self.keys_by_date[date_key], key=lambda k: abs(self._hour_of(k) - want))
        self._set_step(self.keys.index(best))

    def _set_hour(self, i: int) -> None:
        hours = self.keys_by_date[self.date_key]
        self._set_step(self.keys.index(hours[int(np.clip(i, 0, len(hours) - 1))]))

    @staticmethod
    def _hour_of(key: str) -> float:
        hhmm = key.split("_")[1]
        return int(hhmm[:2]) + int(hhmm[2:]) / 60.0

    def _arrow(self, fn, delta: int) -> None:
        """Arrow keys step time, unless something is being typed.

        The bindings are on the root window so the plot tabs respond without
        anything focused. That is wrong inside a text field: the Trace tab has
        entries for dates, ray counts and heights, where Left/Right belong to the
        caret -- and stepping the timestep there would also drop the caret and
        redraw a figure for no reason.
        """
        try:
            focus = self.root.focus_get()
        except (KeyError, tk.TclError):
            focus = None
        if isinstance(focus, (tk.Entry, ttk.Entry, tk.Text, tk.Listbox, ttk.Combobox)):
            return
        fn(delta)

    def _step(self, delta: int) -> None:
        """Step by hour, running on into the next or previous date at the ends."""
        self._set_step(self.step_i + delta)

    def _step_date(self, delta: int) -> None:
        self._set_date(self.dates.index(self.date_key) + delta)

    def _on_colour(self, _=None) -> None:
        self._refresh_all()

    def _on_metric(self, _=None) -> None:
        self._refresh_all()

    def _on_weights(self) -> None:
        """Weights, DNI or aperture changed -- every derived number is stale."""
        self._field_cache.clear()
        self._rows_cache.clear()
        self._refresh_all()  # _aperture_cache is keyed by radius, so it survives

    def _on_spotview(self) -> None:
        self._dirty.add("Spot")
        self._redraw_current()

    def _on_bins(self) -> None:
        """Bin count changes the picture only -- power is bin-independent."""
        self._field_cache.clear()
        self._dirty.add("Spot")
        self._redraw_current()
        self._update_bin_label()

    def _on_overlay(self) -> None:
        self._dirty.add("Field")
        self._redraw_current()

    def _on_dni_mode(self, _=None) -> None:
        """Annual-energy DNI model changed.

        The provider cache is keyed on mode (see ``_dni_provider``), so
        nothing needs clearing here -- a mode seen before is already built,
        and one seen for the first time is built on demand. ``_energy_key``
        already folds the mode in, so this new key is simply not in
        ``self._energy_cache`` yet; marking the tab dirty and redrawing is
        exactly what ``_ensure_energy`` needs to see that and kick off the
        background recompute, the same mechanism the shading checkbox uses.
        """
        self._dirty.add("Energy")
        self._redraw_current()

    def _update_bin_label(self) -> None:
        size = 2.0 * self.cfg.receiver.window_mm / self._bins()
        stored = int(self.cfg.receiver.grid_size)
        note = f"{size:.1f} mm"
        if self._bins() > stored:
            note += "  (from raw rays)"
        self.lbl_bins.config(text=note)

    def _select(self, hid: int, redraw_table: bool = True) -> None:
        """Change the selected heliostat and refresh every view.

        ``redraw_table`` exists because the table both *reports* the selection
        (by highlighting a row) and *sets* it (when a row is clicked). Rebuilding
        it in response to its own click would re-highlight, re-fire
        ``<<TreeviewSelect>>``, and recurse until the window locks up -- so a
        selection that came from the table leaves the table alone. It also keeps
        the scroll position, which a rebuild would throw away.
        """
        self.selected = None if hid is None else int(hid)
        self.var_helio.set("" if self.selected is None else str(self.selected))
        # Energy is field-total, not heliostat-specific, so a selection change
        # never needs to invalidate it -- but this assignment replaces the
        # whole dirty set rather than adding to it, so "Energy" has to be
        # named here or a pending first-draw flag set elsewhere (by
        # _refresh_all or the background compute finishing) is silently lost,
        # and the tab never draws the first time it is opened.
        self._dirty = {"Field", "Spot", "Through day", "Distribution", "Energy"}
        if redraw_table:
            self._dirty.add("Table")
        self._draw_readout()
        self._redraw_current()

    def _select_typed(self) -> None:
        """Empty box means no selection -- the field is the subject, not one mirror."""
        text = self.var_helio.get().strip()
        if not text:
            self._select(None)
            return
        try:
            hid = int(text)
        except ValueError:
            return
        if hid not in self.row_of:
            self._status(f"heliostat {hid} is not in this run")
            return
        self._select(hid)

    def _select_extreme(self, which: str) -> None:
        rows = self._rows_for_display()
        col = self.var_colour.get()
        if col not in rows.columns:
            return
        self._select(rows[col].idxmin() if which == "min" else rows[col].idxmax())

    def _on_field_click(self, event) -> None:
        if event.inaxes is None or event.xdata is None:
            return
        xy = getattr(self, "_field_xy", None)
        if xy is None or not len(xy):
            return
        d = np.hypot(xy[:, 0] - event.xdata, xy[:, 1] - event.ydata)
        self._select(self._field_ids[int(np.argmin(d))])

    def _on_table_select(self, _event) -> None:
        # Ignore the highlight the table sets on itself while being rebuilt;
        # only a genuine user click should move the selection.
        if self._syncing:
            return
        sel = self.tree.selection()
        if not sel:
            return
        hid = int(self.tree.item(sel[0], "values")[0])
        if hid != self.selected:
            self._select(hid, redraw_table=False)

    def _sort_table(self, col: str) -> None:
        prev_col, prev_desc = self._sort_state
        self._sort_state = (col, not prev_desc if col == prev_col else False)
        self._draw_table()

    def _open_run(self) -> None:
        path = filedialog.askdirectory(title="sweep output directory",
                                       initialdir=str(self.store.root.parent))
        if not path:
            return
        from .store import RunStore

        try:
            object.__setattr__(self.cfg.storage, "root", path)
            self.store = RunStore(Path(path), cfg=self.cfg, mode="r")
            self._load_run()
        except Exception as exc:
            messagebox.showerror("could not open", f"{type(exc).__name__}: {exc}")
            return
        self.step_i = min(self.step_i, len(self.keys) - 1)
        if self.selected not in self.row_of:
            self.selected = self.ids[0]
            self.var_helio.set(str(self.selected))
        self.lbl_run.config(text=self.store.root.name)
        self.cmb_date.config(values=self.date_labels())
        self._sync_time_widgets()
        self.var_bins.set(str(self.cfg.receiver.grid_size))
        for var in (self.var_colour, self.var_metric):
            if var.get() not in self.columns:
                var.set(self.columns[0])
        self._refresh_all()

    def _reload(self) -> None:
        """Pick up timesteps written since the run was opened."""
        before = len(self.keys)
        self._load_run()
        self.step_i = min(self.step_i, len(self.keys) - 1)
        self.cmb_date.config(values=self.date_labels())
        self._sync_time_widgets()
        self._refresh_all()
        self._status(f"reloaded: {len(self.keys)} timesteps (+{len(self.keys)-before})")

    # ------------------------------------------------------------------
    # drawing
    # ------------------------------------------------------------------
    def _aperture(self):
        """Entrance-aperture radius in mm, or None for the full receiver window.

        A radius, not a half-width, so it means the same thing as the aperture
        swept by ``beamdown compare`` and as ``spot_metrics``' spillage. It used
        to be a view-only crop, which is why changing it moved no number.
        """
        try:
            v = float(self.var_aperture.get())
        except (ValueError, tk.TclError):
            return None
        return v if v > 0 else None

    def _refresh_all(self) -> None:
        self._dirty = {"Field", "Spot", "Through day", "Distribution", "Energy", "Table"}
        self._draw_readout()
        self._redraw_current()
        self._update_bin_label()
        n = len(self.summary[self.summary.timestep == self.key])
        self.lbl_run_info.config(
            text=f"{len(self.ids)} heliostats · {len(self.keys)} timesteps · {n} rows here"
                 f"\n{self.run_optics_label()}")

    def run_optics_label(self) -> str:
        """What this run's optics were, read from the manifest, not from config.

        The two differ constantly -- config.toml is the file as it is *now*, and
        a stored run was produced by whatever the command line said then. Flat
        and focused runs of the same layout differ by a large factor in collected
        energy and are otherwise indistinguishable in the store, so a loaded run
        has to say which it was.

        Falls back to config.toml only for runs written before the manifest
        carried these keys, and says so rather than asserting it.
        """
        manifest = self.store.manifest
        layout = manifest.get("secondary")
        flat = manifest.get("flat_mirrors")
        stale = layout is None or flat is None
        if layout is None:
            layout = self.cfg.optics.secondary
        if flat is None:
            flat = bool(self.cfg.optics.flat_mirrors)
        shape = "flat heliostats" if flat else "focused heliostats"
        return f"{layout} · {shape}" + (" (from config; run predates it)" if stale else "")

    def _redraw_current(self) -> None:
        try:
            name = self.book.tab(self.book.select(), "text")
        except tk.TclError:
            return
        # The Design tab is not a view of the loaded run -- it evaluates
        # geometry from config and the field file alone -- so it is deliberately
        # outside the run-driven _dirty bookkeeping, which every _refresh_all
        # would otherwise reset out from under it. It draws once, on first
        # sight, and after that only when its own controls move.
        if name == "Design":
            if not self._design_drawn:
                self._design_refresh()
            return
        if name not in self._dirty:
            return
        self._dirty.discard(name)
        try:
            {"Field": self._draw_field, "Spot": self._draw_spot,
             "Through day": self._draw_curve, "Distribution": self._draw_hist,
             "Energy": self._draw_energy, "Table": self._draw_table}[name]()
        except Exception as exc:
            self._status(f"{name}: {type(exc).__name__}: {exc}")
            traceback.print_exc()

    def _style(self, ax) -> None:
        """Minimal frame, from the shared paper style -- see beamdown.plot_style."""
        plot_style.style_axes(ax)

    # ------------------------------------------------------------------
    # export: the figure as drawn, and the numbers behind it
    # ------------------------------------------------------------------
    #
    # Two actions per tab, one implementation. The dialogs are thin wrappers
    # over ``_export_figure_to`` / ``_export_data_to``, which take a path -- so
    # the tests exercise the real export into a temp directory without a file
    # dialog, and what the tests check is what the buttons do.

    @staticmethod
    def _safe(text: str) -> str:
        """A filename fragment: no separators, no spaces, no surprises."""
        keep = [c if (c.isalnum() or c in "-_.") else "_" for c in str(text)]
        return "".join(keep).strip("_") or "untitled"

    def _export_stem(self, name: str) -> str:
        """A default filename that says what the file contains.

        Run name, view, timestep and whatever the view is keyed on -- never
        "figure1". Two exports from different tabs, different timesteps or
        different metrics cannot collide, so a directory of these stays
        readable a month later.
        """
        if name == "Design":
            res = self._design_result or {}
            layout = res.get("layout", "design")
            if layout == "axicon":
                tail = f"tip{res.get('tip_mm', 0)/1000:g}m_ang{res.get('angle_deg', 0):g}"
            elif layout == "cassegrain":
                tail = f"rim{res.get('rim_z_mm', 0)/1000:g}m_f1{res.get('f1_mm', 0)/1000:g}m"
            else:
                tail = f"f1{res.get('f1_mm', 0)/1000:g}m"
            return self._safe(f"design_{layout}_{tail}")

        run = self.store.root.name if self.store is not None else "run"
        parts = [run, name.lower().replace(" ", "")]
        if name in ("Field", "Spot", "Distribution", "Table"):
            parts.append(self.key)
        if name == "Field" or name == "Distribution":
            parts.append(self.var_colour.get())
        elif name == "Spot":
            parts.append(self.var_spotview.get())
            parts.append(f"h{self.selected}" if self.selected is not None else "field")
        elif name == "Through day":
            parts.append(self.var_metric.get())
        elif name == "Energy":
            # _dni_mode() is None when the selector is on "(config default)",
            # whose literal text makes an ugly filename -- name the config
            # rather than the placeholder.
            parts.append("dni_" + (self._dni_mode() or "config"))
        return self._safe("_".join(str(p) for p in parts))

    def _export_figure_to(self, name: str, path) -> list:
        """Write the named tab's figure as PNG (600 dpi) + PDF. Returns the paths."""
        fig = self.figures.get(name)
        if fig is None:
            raise KeyError(f"{name} has no figure to save")
        return plot_style.save_figure(fig, path)

    def _export_data_to(self, name: str, path) -> list:
        """Write the numbers behind the named tab. Returns the paths written.

        One CSV per logical table. Most views have exactly one; the Design tab
        has two (the readout, and the cross-section curve), which are written
        beside each other with a suffix rather than crammed into one file with
        a discriminator column no spreadsheet would thank you for.
        """
        base = Path(path)
        if base.suffix.lower() == ".csv":
            base = base.with_suffix("")
        if base.parent and str(base.parent) not in ("", "."):
            base.parent.mkdir(parents=True, exist_ok=True)

        written = []
        for suffix, frame in self._export_frames(name):
            out = base.with_name(base.name + suffix).with_suffix(".csv")
            frame.to_csv(out, index=False)
            written.append(out)
        return written

    def _export_frames(self, name: str) -> list:
        """``[(filename suffix, DataFrame)]`` -- the processed data behind a view.

        Deliberately the *displayed* numbers, not the raw store: the aperture,
        the DNI override and the shading toggle are all applied, so the CSV
        matches the picture beside it. Reading the sweep's untouched output is
        what ``store.summary()`` is for.
        """
        if name in ("Field", "Table"):
            rows = self._rows_for_display().reset_index()
            if name == "Table":
                col, desc = self._sort_state
                if col and col in rows.columns:
                    rows = rows.sort_values(col, ascending=not desc)
            return [("", rows)]

        if name == "Spot":
            if self.var_spotview.get() == "encircled":
                return [("", self._encircled_frame())]
            return [("", self._flux_frame())]

        if name == "Through day":
            return [("", self._through_day_frame())]

        if name == "Distribution":
            return [("", self._histogram_frame())]

        if name == "Energy":
            return [("", self._energy_frame())]

        if name == "Design":
            res = dict(self._design_result or {})
            notes = res.pop("notes", [])
            rows = [{"quantity": k, "value": v} for k, v in res.items()]
            rows += [{"quantity": f"note_{i}", "value": n}
                     for i, n in enumerate(notes)]
            summary = pd.DataFrame(rows, columns=["quantity", "value"])
            section = self._design_section
            frames = [("", summary)]
            if section is not None and len(section):
                frames.append(("_section", section))
            return frames

        raise KeyError(f"no data export defined for {name!r}")

    # -- the per-view tables --------------------------------------------
    def _flux_frame(self):
        """The receiver flux map(s) as a long table: x_mm, y_mm, flux_w_m2."""
        w = self.cfg.receiver.window_mm
        bins = self._bins()
        edges = np.linspace(-w, w, bins + 1)
        mid = 0.5 * (edges[:-1] + edges[1:])
        # imshow was given origin="lower" and extent [-w, w, -w, w], so
        # flux[j, i] sits at (x = mid[i], y = mid[j]) -- meshgrid in that order.
        X, Y = np.meshgrid(mid, mid)
        parts = []
        for label, flux in self._spot_pair():
            if flux is None:
                continue
            parts.append(pd.DataFrame({
                "panel": label,
                "x_mm": X.ravel(),
                "y_mm": Y.ravel(),
                "flux_w_m2": np.asarray(flux, float).ravel(),
            }))
        return (pd.concat(parts, ignore_index=True) if parts
                else pd.DataFrame(columns=["panel", "x_mm", "y_mm", "flux_w_m2"]))

    def _encircled_frame(self):
        """Cumulative power against aperture radius -- the curve as plotted."""
        from .metrics import bin_radius

        rr = bin_radius(self.cfg, self._bins()).ravel()
        order = np.argsort(rr)
        grid = np.linspace(0.0, self.cfg.receiver.window_mm, 240)
        area = self.bin_area_m2
        parts = []
        for label, flux in self._spot_pair():
            if flux is None:
                continue
            cum = np.concatenate(([0.0], np.cumsum(flux.ravel()[order] * area)))
            power = cum[np.searchsorted(rr[order], grid, side="right")]
            parts.append(pd.DataFrame({
                "panel": label,
                "radius_mm": grid,
                "enclosed_power_w": power,
            }))
        return (pd.concat(parts, ignore_index=True) if parts
                else pd.DataFrame(columns=["panel", "radius_mm", "enclosed_power_w"]))

    def _through_day_frame(self):
        """One row per timestep: the field curve, and the selected heliostat's."""
        metric = self.var_metric.get()
        adjusted = pd.concat([self._rows_for_display(k).assign(timestep=k)
                              for k in self.keys]).reset_index()
        out = (adjusted.groupby(["date", "hour", "timestep"], as_index=False)
               .agg(field_power_w=("power_w", "sum")))
        out["field_power_kw"] = out["field_power_w"] / 1000.0
        if self.selected is not None and metric in adjusted.columns:
            mine = (adjusted[adjusted.heliostat_id == self.selected]
                    [["timestep", metric]]
                    .rename(columns={metric: f"heliostat_{metric}"}))
            out = out.merge(mine, on="timestep", how="left")
            out.insert(0, "heliostat_id", self.selected)
        return out.sort_values(["date", "hour"]).reset_index(drop=True)

    def _histogram_frame(self):
        """The bars as drawn: bin edges and how many heliostats fell in each."""
        col = self.var_colour.get()
        vals = self._rows_for_display()[col].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        nbins = min(50, max(8, len(vals) // 8))
        counts, edges = np.histogram(vals, bins=nbins)
        return pd.DataFrame({
            "column": col,
            "bin_left": edges[:-1],
            "bin_right": edges[1:],
            "bin_centre": 0.5 * (edges[:-1] + edges[1:]),
            "heliostats": counts,
        })

    def _energy_frame(self):
        """The annual curve: one row per day, plus the fitted sinusoid on it."""
        result = self._energy_cache.get(self._energy_key())
        if result is None:
            raise RuntimeError("annual energy is still computing -- try again "
                               "once the Energy tab has drawn")
        daily = result["annual"]["daily"].copy()
        traced = set(self.dates_as_dates)
        daily["traced"] = daily["date"].isin(traced)
        daily["energy_mwh"] = daily["energy_kwh"] / 1000.0
        sine = result.get("sine")
        if sine is not None:
            doy = pd.to_datetime(daily["date"]).dt.dayofyear.to_numpy(float)
            daily["fitted_energy_mwh"] = sine["predict"](doy) / 1000.0
        return daily.reset_index(drop=True)

    # -- the dialogs ----------------------------------------------------
    def _save_figure_dialog(self, name: str) -> None:
        stem = self._export_stem(name)
        path = filedialog.asksaveasfilename(
            title=f"Save the {name} figure — {plot_style.describe()}",
            initialfile=stem + ".png", defaultextension=".png",
            filetypes=[("PNG and PDF together", "*.png"), ("all files", "*.*")])
        if not path:
            return
        try:
            written = self._export_figure_to(name, path)
        except Exception as exc:
            messagebox.showerror("could not save", f"{type(exc).__name__}: {exc}")
            return
        self._status("wrote " + ", ".join(str(p) for p in written))

    def _save_data_dialog(self, name: str) -> None:
        stem = self._export_stem(name)
        path = filedialog.asksaveasfilename(
            title=f"Save the numbers behind the {name} view",
            initialfile=stem + ".csv", defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("all files", "*.*")])
        if not path:
            return
        try:
            written = self._export_data_to(name, path)
        except Exception as exc:
            messagebox.showerror("could not save", f"{type(exc).__name__}: {exc}")
            return
        self._status("wrote " + ", ".join(str(p) for p in written))

    def _draw_field(self) -> None:
        """The field seen from directly above, mirror by mirror.

        Each heliostat is its actual rectangle projected down, not a dot, so the
        pointing is legible: at low sun the mirrors stand up and draw as slivers.
        The metric is carried on a pale grey ramp, leaving saturated colour free
        for the two loss overlays -- on viridis, a translucent wash over the dark
        purple end is invisible and over the yellow end is unreadable.
        """
        from matplotlib.cm import ScalarMappable
        from matplotlib.collections import PolyCollection
        from matplotlib.colors import Normalize

        fig = self.figures["Field"]
        fig.clear()
        ax = fig.add_subplot(111)
        rows = self._rows_for_display().reset_index()
        col = self.var_colour.get()
        vals = rows[col].to_numpy(float)
        ids, outlines, shadows, cone = self._field_polygons()

        finite = vals[np.isfinite(vals)]
        norm = Normalize(vmin=float(finite.min()) if finite.size else 0.0,
                         vmax=float(finite.max()) if finite.size else 1.0)
        cmap = _metric_colormap()

        ax.add_collection(PolyCollection(
            outlines, array=vals, cmap=cmap, norm=norm,
            edgecolors="#9aa4ad", linewidths=0.3, zorder=2))

        # Two different things, deliberately drawn differently.
        #
        # The pale wash is the literal shadow pattern on the ground, and it is
        # thrown clear of the mirror that cast it -- 29 m at 9.7 deg from a 5 m
        # pedestal. It shows how the field packs, not which mirror lost what: the
        # mirror actually shaded sits at pedestal height, where the offset is
        # nothing, so it is the *nearest* up-sun neighbour, not the one 29 m away.
        #
        # The saturated fill on each mirror is the loss the sweep computed for
        # that mirror. That is the number, so it goes where the number belongs.
        by_id = rows.set_index("heliostat_id")

        def loss_fill(column, rgb):
            lost = 1.0 - by_id[column].reindex(ids).fillna(1.0).to_numpy(float)
            rgba = np.zeros((len(ids), 4))
            rgba[:, :3] = rgb
            rgba[:, 3] = np.clip(lost * 1.15, 0.0, 0.9)
            return PolyCollection(outlines, facecolors=rgba, edgecolors="none")

        # -- where the secondary actually is, as opposed to where its shadow is --
        #
        # Two different circles that a reader would otherwise conflate. This one
        # is nailed to the axis at (0, 0) and never moves; the dashed one below
        # is the silhouette and walks across the field with the sun. Dotted, a
        # different colour, and separately labelled in the legend for exactly
        # that reason. Prime focus has no body and draws nothing.
        footprint = self.secondary_footprint()
        if footprint is not None:
            ax.plot(footprint[:, 0], footprint[:, 1], linestyle=":",
                    color=BODY_EDGE, linewidth=1.6, zorder=5,
                    label=f"{self.cfg.optics.secondary} above the field "
                          f"(r {self.cfg.geometry.axicon_aperture_radius_mm/1000:.0f} m, "
                          f"fixed on the axis)")

        if self.var_show_shadow.get():
            ax.add_collection(PolyCollection(
                shadows, facecolors=GROUND_COLOUR, edgecolors="none",
                alpha=0.16, zorder=1))
            # The secondary's shadow, drawn at mirror height, so what it covers
            # is what it shades. Outlined rather than washed, because at mid
            # elevations it lands squarely on the inner field and is the point.
            if len(cone):
                ax.add_collection(PolyCollection(
                    [cone], facecolors=GROUND_COLOUR, alpha=0.22,
                    edgecolors=CONE_EDGE, linewidths=1.4, linestyle="--", zorder=2,
                    label="its shadow at this instant (moves with the sun)"))
            if "eta_shade" in by_id.columns:
                c = loss_fill("eta_shade", SHADOW_COLOUR)
                c.set_zorder(3)
                ax.add_collection(c)
        if self.var_show_block.get() and "eta_block" in by_id.columns:
            c = loss_fill("eta_block", BLOCK_COLOUR)
            c.set_zorder(4)
            ax.add_collection(c)

        sel = rows[rows.heliostat_id == self.selected]
        if len(sel):
            ax.scatter(sel.x_m, sel.y_m, s=190, facecolors="none",
                       edgecolors="#d6604d", linewidths=2.0, zorder=6)
            ax.annotate(f"#{self.selected}", (float(sel.x_m.iloc[0]), float(sel.y_m.iloc[0])),
                        textcoords="offset points", xytext=(10, 8),
                        fontsize=9, color="#d6604d", weight="bold")

        pad = max(self.cfg.field.mirror_width_mm / 1000.0, 3.0)
        span = np.concatenate([outlines.reshape(-1, 2), shadows.reshape(-1, 2)])
        if len(cone) and self.var_show_shadow.get():
            span = np.concatenate([span, cone])
        ax.set_xlim(span[:, 0].min() - pad, span[:, 0].max() + pad)
        ax.set_ylim(span[:, 1].min() - pad, span[:, 1].max() + pad)
        ax.set_aspect("equal")
        ax.set_xlabel("x (m)", fontsize=9)
        ax.set_ylabel("y (m)", fontsize=9)

        overlays = []
        if self.var_show_shadow.get() and "eta_shade" in by_id.columns:
            overlays.append(f"blue = shaded ({1 - by_id.eta_shade.mean():.0%} lost), "
                            f"pale wash = shadows on the ground")
            # Label from the layout, not the word "axicon": with a Cassegrain the
            # dashed outline is a disc, and with prime focus there is no outline
            # at all and no loss to report.
            if "eta_secondary" in by_id.columns and len(cone):
                overlays[-1] += (f", dashed = {self.cfg.optics.secondary} "
                                 f"({1 - by_id.eta_secondary.mean():.1%})")
        if self.var_show_block.get() and "eta_block" in by_id.columns:
            overlays.append(f"red = blocked ({1 - by_id.eta_block.mean():.0%} lost)")
        ax.set_title(f"{self.step_label(self.key)}\ncoloured by {col}"
                     + (f"   ·   {'  ·  '.join(overlays)}" if overlays else "")
                     + f"\n{self._weights_label()}", fontsize=9.5)
        self._style(ax)
        # Only the two secondary circles are labelled, and only when at least one
        # of them is on screen -- an always-on legend would take a corner of the
        # plot to say nothing. It exists so the fixed footprint and the moving
        # shadow can never be read as the same object.
        handles, _labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(loc="upper right", fontsize=7.5, framealpha=0.9,
                      handlelength=2.6, borderpad=0.5).set_zorder(7)
        ax.grid(True, color="#eeeeee", linewidth=0.5)
        ax.set_axisbelow(True)
        cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                          fraction=0.043, pad=0.02)
        cb.set_label(col, fontsize=8.5)
        cb.ax.tick_params(labelsize=7.5)

        self._field_xy = rows[["x_m", "y_m"]].to_numpy()
        self._field_ids = rows.heliostat_id.to_numpy()
        plot_style.finish(fig)
        self.canvases["Field"].draw_idle()

    def _weights_label(self) -> str:
        state = "shading×blocking on" if self.var_shading.get() else "weights off"
        radius = self._aperture()
        aperture = f"aperture r{radius:.0f} mm" if radius else "full window"
        return f"DNI {self._dni():.0f} W/m²   ·   {state}   ·   {aperture}"

    def _spot_pair(self):
        """(label, flux map) per panel.

        With no heliostat selected the field gets the whole figure rather than
        half of it beside an empty axis -- "no selection" is a legitimate way to
        look at a sweep, not a missing value to apologise for.
        """
        field = (f"all {len(self.ids)} heliostats", self._field_flux())
        if self.selected is None:
            return [field]
        return [(f"heliostat {self.selected}", self._heliostat_flux()), field]

    def _draw_spot(self) -> None:
        if self.var_spotview.get() == "encircled":
            self._draw_encircled()
            return

        from matplotlib.colors import LogNorm
        from matplotlib.patches import Circle

        from .metrics import radial_mask
        from .plots import flux_colormap

        fig = self.figures["Spot"]
        fig.clear()
        panels = self._spot_pair()
        axes = np.atleast_1d(fig.subplots(1, len(panels)))
        w = self.cfg.receiver.window_mm
        radius = self._aperture()
        bins = self._bins()
        mask = radial_mask(self.cfg, radius, bins) if radius else None
        area = self.bin_area_m2
        cmap = flux_colormap()

        for ax, (title, flux) in zip(axes, panels):
            if flux is None:
                ax.text(0.5, 0.5, "not available", ha="center", transform=ax.transAxes)
                continue
            kw = {}
            if self.var_log.get() and flux.max() > 0:
                floor = max(flux[flux > 0].min(), flux.max() / 1e4)
                kw["norm"] = LogNorm(vmin=floor, vmax=flux.max())
            else:
                kw["vmin"] = 0.0
            im = ax.imshow(flux / 1000.0 if not kw.get("norm") else flux,
                           origin="lower", cmap=cmap, extent=[-w, w, -w, w],
                           aspect="equal", interpolation="nearest", **kw)

            # The aperture is where the numbers come from, so it is drawn.
            if radius:
                ax.add_patch(Circle((0, 0), radius, fill=False, linewidth=1.3,
                                    edgecolor="#59c6f3", linestyle="--", zorder=4))
                ax.set_xlim(-radius, radius)
                ax.set_ylim(-radius, radius)

            # A scalar weight scales every pixel *and* the autoscaled colour bar,
            # so the picture cannot show it -- these numbers are the only place a
            # single heliostat's shading is visible.
            total = flux.sum() * area
            captured = (flux[mask].sum() * area) if mask is not None else total
            head = f"{title}\n{captured/1000:.1f} kW"
            if radius:
                head += f" inside r{radius:.0f}"
                if total > 0:
                    head += f", spill {1 - captured/total:.1%}"
            # Peak is resolution-dependent by nature, so the bin size is named
            # beside it -- a peak quoted without one is not comparable.
            head += (f"\npeak {flux.max()/1000.0:.1f} kW/m² "
                     f"at {2*self.cfg.receiver.window_mm/bins:.1f} mm bins")

            ax.set_title(head, fontsize=9.5)
            ax.set_xlabel("x (mm)", fontsize=8.5)
            ax.set_ylabel("y (mm)", fontsize=8.5)
            self._style(ax)
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
            cb.set_label("W/m²" if kw.get("norm") else "kW/m²", fontsize=8)
            cb.ax.tick_params(labelsize=7.5)

        fig.suptitle(f"{self.step_label(self.key)}   ·   {self._weights_label()}",
                     fontsize=9)
        plot_style.finish(fig)
        self.canvases["Spot"].draw_idle()

    def _draw_encircled(self) -> None:
        """Cumulative power vs aperture radius, beside the spot it came from.

        Reads the same flux maps the image view draws, so the curve at the
        aperture radius equals the captured power printed on the image -- if the
        two ever disagree, one of them is wrong.
        """
        from .metrics import bin_radius
        from .plots import encircled_energy_axes

        fig = self.figures["Spot"]
        fig.clear()
        panels = self._spot_pair()
        axes = np.atleast_1d(fig.subplots(1, len(panels)))
        radius = self._aperture()
        area = self.bin_area_m2
        rr = bin_radius(self.cfg, self._bins()).ravel()
        order = np.argsort(rr)
        grid = np.linspace(0.0, self.cfg.receiver.window_mm, 240)

        colours = ("#d6604d", "#2166ac") if len(panels) == 2 else ("#2166ac",)
        for ax, (title, flux), colour in zip(axes, panels, colours):
            if flux is None:
                ax.text(0.5, 0.5, "not available", ha="center", transform=ax.transAxes)
                continue
            cum = np.concatenate(([0.0], np.cumsum(flux.ravel()[order] * area)))
            power = cum[np.searchsorted(rr[order], grid, side="right")]
            encircled_energy_axes(ax, grid, power, aperture_mm=radius, colour=colour)
            ax.set_title(f"{title}\n{cum[-1]/1000:.1f} kW total", fontsize=9.5)
            ax.set_xlabel("radius about the receiver axis (mm)", fontsize=8.5)
            ax.set_ylabel("enclosed power (kW)", fontsize=8.5)

        fig.suptitle(f"{self.step_label(self.key)}   ·   {self._weights_label()}   ·   "
                     f"encircled energy about the axis, not the centroid", fontsize=9)
        plot_style.finish(fig)
        self.canvases["Spot"].draw_idle()

    def _draw_curve(self) -> None:
        fig = self.figures["Through day"]
        fig.clear()
        metric = self.var_metric.get()
        # Every timestep re-derived under the current weights, so the curve
        # agrees with the spot and the table rather than showing whatever the
        # sweep happened to bake in.
        adjusted = pd.concat([self._rows_for_display(k).assign(timestep=k)
                              for k in self.keys]).reset_index()
        current_date = adjusted[adjusted.timestep == self.key].iloc[0]["date"]

        # Separate axes rather than a second y-scale: field total is ~600x the
        # single heliostat, and overlaying them on one axis would hide both. With
        # nothing selected the upper axis has nothing to say, so it is not drawn.
        if self.selected is None:
            ax1 = None
            ax2 = fig.subplots(1, 1)
        else:
            ax1, ax2 = fig.subplots(2, 1, sharex=True)
            mine = adjusted[adjusted.heliostat_id == self.selected]
            for date, grp in mine.groupby("date"):
                grp = grp.sort_values("hour")
                now = date == current_date
                ax1.plot(grp.hour, grp[metric], "-o",
                         linewidth=2.2 if now else 1.0, markersize=4.5 if now else 2.5,
                         color="#d6604d" if now else "#cccccc",
                         zorder=3 if now else 1, label=str(date) if now else None)

        field = (adjusted.groupby(["date", "hour"], as_index=False)
                 .agg(power_w=("power_w", "sum")))
        for date, grp in field.groupby("date"):
            grp = grp.sort_values("hour")
            now = date == current_date
            ax2.plot(grp.hour, grp.power_w / 1000.0, "-o",
                     linewidth=2.2 if now else 1.0, markersize=4.5 if now else 2.5,
                     color="#4393c3" if now else "#cccccc", zorder=3 if now else 1)

        hour = float(adjusted[adjusted.timestep == self.key].iloc[0]["hour"])
        for ax in ([ax1, ax2] if ax1 is not None else [ax2]):
            ax.axvline(hour, color="#888888", linewidth=0.8, linestyle="--", zorder=0)
            self._style(ax)
            ax.grid(True, axis="y", color="#eeeeee", linewidth=0.5)
            ax.set_axisbelow(True)
        if ax1 is not None:
            ax1.set_ylabel(metric, fontsize=8.5)
            ax1.set_title(f"heliostat {self.selected} — {metric}   "
                          f"(highlighted {current_date}; grey = other dates)",
                          fontsize=9.5)
        # The heliostat count rides in the ylabel rather than in a title of its
        # own: it is the same one line of information, and a title between two
        # stacked panels costs a band of white space across the whole figure.
        ax2.set_ylabel(f"field total, {len(self.ids)} heliostats (kW)", fontsize=8.5)
        ax2.set_xlabel("local hour", fontsize=8.5)
        plot_style.finish(fig)
        self.canvases["Through day"].draw_idle()

    def _draw_hist(self) -> None:
        fig = self.figures["Distribution"]
        fig.clear()
        ax = fig.add_subplot(111)
        rows = self._rows_for_display()
        col = self.var_colour.get()
        vals = rows[col].to_numpy(float)
        vals = vals[np.isfinite(vals)]

        ax.hist(vals, bins=min(50, max(8, len(vals) // 8)),
                color="#9ecae1", edgecolor="white", linewidth=0.6)
        r = self._selected_row()
        if r is not None and np.isfinite(float(r[col])):
            v = float(r[col])
            ax.axvline(v, color="#d6604d", linewidth=2.0)
            rank = int((vals < v).sum()) + 1
            ax.annotate(f"#{self.selected}: {_fmt(v)}  (rank {rank} of {len(vals)})",
                        xy=(v, ax.get_ylim()[1] * 0.92), fontsize=9, color="#d6604d",
                        ha="left" if v < np.median(vals) else "right",
                        xytext=(6 if v < np.median(vals) else -6, 0),
                        textcoords="offset points")
        ax.set_xlabel(col, fontsize=9)
        ax.set_ylabel("heliostats", fontsize=9)
        ax.set_title(f"{col} across the field — {self.step_label(self.key)}\n"
                     f"mean {vals.mean():.4g}   min {vals.min():.4g}   "
                     f"max {vals.max():.4g}", fontsize=9.5)
        self._style(ax)
        ax.grid(True, axis="y", color="#eeeeee", linewidth=0.5)
        ax.set_axisbelow(True)
        plot_style.finish(fig)
        self.canvases["Distribution"].draw_idle()

    def _draw_energy(self) -> None:
        """Total energy for the selected day, the annual curve, and its audit trail.

        ``annual_energy`` is computed once per shading state in a background
        thread (see ``_ensure_energy``); this draws whatever is cached right
        now and, if nothing is cached yet, kicks off the compute and shows a
        placeholder. The job's completion re-adds "Energy" to ``self._dirty``
        and asks for a redraw (see ``_poll_jobs``), so the real figure appears
        without the user having to do anything once the background thread
        lands.
        """
        result = self._energy_cache.get(self._energy_key())
        fig = self.figures["Energy"]
        fig.clear()

        if result is None:
            self._ensure_energy()
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5,
                   "computing annual energy in the background\n"
                   "(walks every hour of the year on first view of this tab)",
                   ha="center", va="center", fontsize=10, color="#888888",
                   transform=ax.transAxes)
            ax.set_axis_off()
            self.canvases["Energy"].draw_idle()
            self.lbl_energy_headline.config(text="computing…")
            self.lbl_energy_annual.config(text="")
            self.lbl_energy_check.config(text="")
            return

        annual, sine, checks = result["annual"], result["sine"], result["checks"]
        date = self._current_date()
        daily = annual["daily"]
        row = daily[daily["date"] == date]

        if len(row):
            mwh_today = float(row["energy_kwh"].iloc[0]) / 1000.0
            self.lbl_energy_headline.config(text=f"{date}:  {mwh_today:.2f} MWh collected")
        else:
            self.lbl_energy_headline.config(text=f"{date}:  outside the computed year")

        # Deliberately not self._weights_label(): that reports the per-view "DNI
        # W/m2" and aperture overrides used by the Field/Spot/Table tabs, which
        # do not apply here -- the annual integral uses the real DNI provider
        # below and has no aperture at all. Naming that provider instead of the
        # spot override avoids implying the annual total came from a constant
        # 1000 W/m2 trace when it did not.
        shading_state = "shading×blocking on" if self.var_shading.get() else "weights off"
        self.lbl_energy_annual.config(
            text=f"Annual total {annual['annual_energy_mwh']:,.1f} MWh   ·   "
                 f"annual optical efficiency {annual['annual_optical_efficiency']:.3f}   ·   "
                 f"{annual['extrapolated_fraction']:.1%} of daylight hours extrapolated "
                 f"beyond the traced declination hull   ·   {shading_state}   ·   "
                 f"DNI model: {self._dni_provider().describe()}")

        check_row = checks[checks["date"] == date] if len(checks) else checks
        if len(check_row):
            c = check_row.iloc[0]
            self.lbl_energy_check.config(
                text=f"cross-check, {date}: direct trapezoid over traced samples "
                     f"{c['traced_energy_kwh']/1000:.2f} MWh   vs   interpolated-model route "
                     f"{c['interpolated_energy_kwh']/1000:.2f} MWh   "
                     f"(residual {c['residual_frac']:+.1%})")
        else:
            self.lbl_energy_check.config(
                text=f"{date} was not ray-traced -- no direct cross-check available")

        from .plots import annual_energy_axes

        ax = fig.add_subplot(111)
        annual_energy_axes(ax, daily, sine_fit=sine, traced_dates=self.dates_as_dates,
                           hourly=annual["hourly"])
        if len(row):
            ax.axvline(pd.Timestamp(date), color="#d6604d", linewidth=1.4, zorder=5)
        plot_style.finish(fig)
        self.canvases["Energy"].draw_idle()

    def _draw_table(self) -> None:
        self.tree.delete(*self.tree.get_children())
        rows = self._rows_for_display().reset_index()
        col, desc = self._sort_state
        if col and col in rows.columns:
            rows = rows.sort_values(col, ascending=not desc)
        for c in self.table_cols:
            mark = "" if c != col else ("  ▼" if desc else "  ▲")
            self.tree.heading(c, text=c + mark)

        # Formatting column-wise beats iterrows(): 645 rows x 14 columns is
        # 9,000 cell lookups, and per-row Series indexing makes that visible.
        text = {c: [_fmt(v) for v in rows[c].to_numpy()] for c in self.table_cols}
        ids = rows.heliostat_id.to_numpy()

        self._syncing = True
        try:
            focus = None
            for i in range(len(rows)):
                iid = self.tree.insert("", "end",
                                       values=[text[c][i] for c in self.table_cols])
                if int(ids[i]) == self.selected:
                    focus = iid
            if focus:
                self.tree.selection_set(focus)
                self.tree.see(focus)
        finally:
            self._syncing = False

    def _draw_readout(self) -> None:
        self.tree_row.delete(*self.tree_row.get_children())
        r = self._selected_row()
        if r is None:
            return
        for name in ["heliostat_id", "x_m", "y_m", "radius_m", "-",
                     "solar_az_deg", "solar_el_deg", "-",
                     "rot_az_deg", "rot_el_deg", "aoi_deg", "cosine_efficiency", "-",
                     "eta_shade", "eta_block", "eta_occlusion",
                     "transmission", "spillage", "-",
                     "rays_landed", "rays_outside_window", "power_w",
                     "power_in_aperture_w", "aperture_spillage", "-",
                     "peak_flux_w_m2", "centroid_x_mm", "centroid_y_mm",
                     "rms_radius_mm", "r50_mm", "r90_mm"]:
            if name == "-":
                self.tree_row.insert("", "end", text="", values=("",))
            elif name in r.index:
                self.tree_row.insert("", "end", text=name, values=(_fmt(r[name]),))

    # ------------------------------------------------------------------
    # Quadoa export
    # ------------------------------------------------------------------
    def _export(self, by_metric: bool) -> None:
        from . import inspect_model as IM

        if by_metric:
            try:
                n = max(1, int(self.var_nexport.get()))
            except ValueError:
                n = 1
            ids = IM.pick_heliostats(self.summary, timestep=self.key,
                                     by=self.var_colour.get(), n=n,
                                     worst=self.var_extreme.get() == "min")
        elif self.selected is None:
            self._status("no heliostat selected — pick one, or export by metric")
            return
        else:
            ids = [self.selected]

        self._status(f"exporting {ids} … needs a free Quadoa license seat")
        self._log(f"export {ids} at {self.key}")
        cfg, key, summary = self.cfg, self.key, self.summary

        def work():
            try:
                report = IM.export_for_inspection(cfg, ids, timestep=key, summary=summary)
                self._jobs.put(("export_ok", report))
            except Exception as exc:
                self._jobs.put(("export_fail", exc))

        threading.Thread(target=work, daemon=True).start()

    def _export_occluders(self) -> None:
        """Write the selected heliostat with its occluders as real geometry.

        The same model the sweep traces, with this heliostat's slots filled in,
        so what opens in Quadoa is the geometry the numbers came from rather than
        a reconstruction of it. No licence seat: it is a text edit.
        """
        if self.selected is None:
            self._status("no heliostat selected")
            return
        from . import build_occluder_model as B

        try:
            report = B.build_from_slot_model(self.cfg, self.summary,
                                             self.selected, self.key)
        except Exception as exc:
            self._log(f"FAILED: {type(exc).__name__}: {exc}")
            self._status(f"{type(exc).__name__}: {exc}")
            return
        self._log(report.describe())
        self._status(f"wrote {report.path}")

    def _poll_jobs(self) -> None:
        try:
            while True:
                kind, payload = self._jobs.get_nowait()
                if kind == "export_ok":
                    self._log(payload.describe())
                    self._status(f"wrote {payload.path}")
                elif kind == "export_fail":
                    name = type(payload).__name__
                    if "License" in name or "dll load failed" in str(payload).lower():
                        msg = ("No free Quadoa license seat. A running sweep or an open "
                               "Quadoa window holds it; close one and retry.")
                    else:
                        msg = f"{name}: {payload}"
                    self._log("FAILED: " + msg)
                    self._status(msg)
                elif kind == "energy_ok":
                    self._energy_cache[payload["key"]] = payload
                    if self._energy_job_key == payload["key"]:
                        self._energy_job_key = None
                    self._status("annual energy ready")
                    self._dirty.add("Energy")
                    self._redraw_current()
                elif kind == "design_build":
                    self._design_proc_busy = False
                    self.btn_design_figure.config(state="normal")
                    self.btn_design_full.config(
                        state=("normal"
                               if (self._design_result.get("layout") == "cassegrain"
                                   and self._design_result.get("feasible"))
                               else "disabled"))
                    self._design_log(payload["text"])
                    # The builder's last stdout line is its own verdict ("PASS
                    # -- next: ...", "refusing to overwrite ... -- pass
                    # --force"), so the status line quotes it rather than
                    # paraphrasing it.
                    self._status(payload.get("verdict")
                                 or f"model build exited {payload['rc']}")
                elif kind == "energy_fail":
                    if self._energy_job_key == payload["key"]:
                        self._energy_job_key = None
                    exc = payload["exc"]
                    self._status(f"annual energy failed: {type(exc).__name__}: {exc}")
                    traceback.print_exc()
        except queue.Empty:
            pass
        self.root.after(150, self._poll_jobs)

    def _log(self, text: str) -> None:
        self.txt_log.config(state="normal")
        self.txt_log.insert("end", text.rstrip() + "\n")
        self.txt_log.see("end")
        self.txt_log.config(state="disabled")

    def _status(self, text: str) -> None:
        self.var_status.set(text)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="beamdown.gui")
    ap.add_argument("--config", default=None)
    ap.add_argument("--output", default=None, help="sweep output directory")
    args = ap.parse_args(argv)

    from .config import load_config
    from .store import RunStore

    cfg = load_config(args.config)
    if args.output:
        object.__setattr__(cfg.storage, "root", args.output)

    store = RunStore(cfg.output_root, cfg=cfg, mode="r")

    root = tk.Tk()
    try:
        root.call("ttk::style", "theme", "use", "vista")
    except tk.TclError:
        pass
    BeamdownGUI(root, cfg, store, config_path=args.config)
    root.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
