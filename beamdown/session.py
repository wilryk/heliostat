"""A single Quadoa session, wrapped.

This is the only module that imports ``quadoa``, and it does so lazily, so the
rest of the package works without a license attached.

Two things this handles that the original scripts did not:

**The DLL search path.** ``import quadoa`` fails with a DLL-initialisation error
unless ``C:\\Program Files\\Quadoa`` is added via :func:`os.add_dll_directory`
first. That was the cause of the import failures.

**Flaky licensing.** The USB HASP key intermittently refuses a session even when
nothing else holds it -- observed failing three consecutive attempts and then
recovering unprompted. Every session creation retries with backoff.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass

import numpy as np

_QUADOA = None

# HASP seat exhaustion surfaces as a DLL-initialisation failure on import (and
# as a modal H0038 dialog). It is not transient, so it must never be retried.
_LICENSE_MARKERS = (
    "dll load failed",
    "h0038",
    "too many concurrent",
    "license",
    "hasp",
)


class LicenseUnavailable(RuntimeError):
    """No Quadoa license seat was free."""


def is_license_error(exc: BaseException) -> bool:
    text = f"{type(exc).__name__}: {exc}".lower()
    return any(m in text for m in _LICENSE_MARKERS)


def _import_quadoa(quadoa_folder: str):
    """Import quadoa once, after fixing the DLL search path."""
    global _QUADOA
    if _QUADOA is None:
        if os.path.isdir(quadoa_folder):
            os.add_dll_directory(quadoa_folder)
        import quadoa as _q

        _QUADOA = _q
    return _QUADOA


def running_quadoa_processes() -> list[dict]:
    """Quadoa-related processes currently running.

    A hung process holding a license seat is the most likely explanation for
    only 3 of 4 seats being available, so this is checked before blaming the
    pool size.
    """
    import subprocess

    try:
        out = subprocess.run(
            ["tasklist", "/FO", "CSV", "/NH"],
            capture_output=True, text=True, timeout=30,
        ).stdout
    except Exception:
        return []

    found = []
    for line in out.splitlines():
        parts = [p.strip('"') for p in line.split('","')]
        if not parts:
            continue
        name = parts[0].strip('"').lower()
        if any(k in name for k in ("quadoa", "hasp", "python")):
            found.append({"name": parts[0].strip('"'),
                          "pid": parts[1] if len(parts) > 1 else "?"})
    return found


def describe_api(quadoa_folder: str, filter_substrings=None) -> str:
    """Dump QuadoaCore's methods and docstrings.

    The API is only documented in these docstrings, so this is how questions
    like "does setConfig alone apply parameter changes?" get answered.
    """
    import inspect

    import gc

    quadoa = _import_quadoa(quadoa_folder)
    core = quadoa.QuadoaCore()
    try:
        lines = []
        for name, member in inspect.getmembers(core, predicate=inspect.ismethod):
            if name.startswith("_"):
                continue
            if filter_substrings and not any(s.lower() in name.lower() for s in filter_substrings):
                continue
            doc = inspect.getdoc(member) or "(no docstring)"
            lines.append(f"{name}\n{'-' * len(name)}\n{doc}\n")
        return "\n".join(lines)
    finally:
        del core
        gc.collect()


@dataclass
class TraceResult:
    """Receiver-plane rays from one heliostat at one instant."""

    xy_mm: np.ndarray      # (N, 2) float32, NaN rays already dropped
    rays_emitted: int
    rays_landed: int

    @property
    def transmission(self) -> float:
        return self.rays_landed / self.rays_emitted if self.rays_emitted else 0.0


class QuadoaSession:
    """One persistent QuadoaCore, reused across many heliostats.

    Creating a session is slow and license-limited, so the sweep creates a small
    fixed number of these and reuses them for every heliostat and timestep.
    """

    def __init__(self, cfg, seq: int | None = None, surface: int | None = None):
        self.cfg = cfg
        self.seq = cfg.trace.analysis_seq if seq is None else seq
        self.surface = cfg.trace.analysis_surface if surface is None else surface
        self._surface_explicit = surface is not None
        self._has_axicon = None
        self.core = None
        self._current_config = None
        self._open()
        self._resolve_surface()

    @property
    def has_axicon_slot(self) -> bool:
        """Whether the loaded model carries the secondary's own shadow plane.

        Probed once by reading the parameter back, rather than assumed from the
        file name: writing ``ax0_x`` to a model without it is silently ignored,
        which would leave the axicon's shadow untraced while everything claimed
        otherwise.
        """
        if self._has_axicon is None:
            try:
                v = self.core.getMulticonfParam("ax0_x", 0)
                self._has_axicon = v is not None and not np.isnan(float(v))
            except Exception:
                self._has_axicon = False
        return self._has_axicon

    def _resolve_surface(self) -> None:
        """Take the detector index from the sequence, not from config.

        ``analysis_surface`` in config.toml is an index *within the sequence*, so
        it is only valid for a model with that exact surface count. Adding
        occluders to the path moved the receiver from index 3 to 17, and reading
        index 3 then samples a parked obscuration -- which still returns rays, so
        nothing errors and every number is quietly wrong. The sequence knows
        where its own detector is; ask it.
        """
        if self._surface_explicit or self.core is None:
            return
        image = int(self.core.getSequenceImageSurface(self.seq))
        if image != self.surface:
            self.surface = image

    # -- lifecycle --------------------------------------------------------
    def _open(self) -> None:
        cfg = self.cfg
        last_error = None

        for attempt in range(max(1, cfg.trace.max_retries)):
            try:
                quadoa = _import_quadoa(cfg.trace.quadoa_folder)
                core = quadoa.QuadoaCore()
                for glass in ("SCHOTT.glas", "Experimental Misc.glas"):
                    path = os.path.join(cfg.trace.quadoa_folder, "glass", glass)
                    if os.path.exists(path):
                        core.loadMaterialFile(path)

                model_path = str(cfg.model_path)
                if not os.path.exists(model_path):
                    raise FileNotFoundError(f"Model not found: {model_path}")
                core.loadModelFile(model_path)
                core.applyChangesAndInitModel()

                self.core = core
                return
            except Exception as exc:
                last_error = exc
                # A seat-exhaustion failure is not transient, and every retry
                # pops another modal H0038 dialog. Fail immediately instead.
                if is_license_error(exc):
                    raise LicenseUnavailable(
                        "No Quadoa license seat available. Close other Quadoa "
                        "sessions (or lower trace.n_workers) and retry. "
                        f"Underlying error: {exc}"
                    ) from None
                time.sleep(1.5 * (attempt + 1))

        raise RuntimeError(
            f"Could not open a Quadoa session after {cfg.trace.max_retries} "
            f"attempts. Last error: {last_error}"
        )

    def close(self) -> None:
        """Release the license seat.

        Dropping the last Python reference is what frees the seat, so the
        reference is deleted and a collection forced. Merely assigning ``None``
        left seats held long enough to make later sessions look 'flaky'.
        """
        import gc

        if self.core is not None:
            try:
                del self.core
            except Exception:
                pass
        self.core = None
        gc.collect()

    def __enter__(self) -> "QuadoaSession":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    # -- parameters -------------------------------------------------------
    def set_global_geometry(self) -> None:
        """Write the geometry values that are shared by every configuration.

        These are ``single_param`` in the model, so the config index is
        irrelevant -- but they must match what the shape solver assumed, or
        pointing and model disagree silently.

        *Which* parameters those are is the layout's business, not this method's:
        see :meth:`SecondaryStrategy.global_params`. For the axicon the set is
        ``sec_height``, ``rec_offset``, ``axi_angle``, written in that order,
        exactly as this method used to do unconditionally. A prime-focus or
        Cassegrain model has no ``axi_angle`` parameter to write.
        """
        from .secondary import get_strategy

        strategy = get_strategy(self.cfg)
        for name, value in strategy.global_params(self.cfg.geometry).items():
            self.core.setMulticonfParam(name, 0, value)

    def set_sun(self, solar_az_deg: float, solar_el_deg: float) -> None:
        """Sun direction. Also a ``single_param``: one value for the model."""
        self.core.setMulticonfParam("solaz", 0, float(solar_az_deg))
        self.core.setMulticonfParam("solze", 0, float(90.0 - solar_el_deg))

    def set_heliostat(self, x_mm, y_mm, solution, config_num: int = 0) -> None:
        """Write one heliostat's position, pointing and shape into a config."""
        c = self.core
        c.setMulticonfParam("posx", config_num, float(x_mm))
        c.setMulticonfParam("posy", config_num, float(y_mm))
        c.setMulticonfParam("rot_az", config_num, float(solution.rot_az_deg))
        c.setMulticonfParam("rot_el", config_num, float(solution.rot_el_deg))
        c.setMulticonfParam("c3", config_num, float(solution.c3))
        c.setMulticonfParam("c4", config_num, float(solution.c4))
        c.setMulticonfParam("c5", config_num, float(solution.c5))

    def set_occluders(self, plan, n_shade: int, n_block: int) -> None:
        """Move the occluder slots onto this heliostat's neighbours.

        Every slot gets a real neighbour -- see :func:`occluder_slots._rank` for
        why there is no "off" position. A slot beyond the occluding ones holds a
        neighbour that blocks nothing, which is exactly neutral.

        These are ``single_param`` entries: one value for the whole model, not
        one per configuration. The bulk sweep writes configuration 0 only, so
        that costs nothing and keeps the multiconfig block from growing by 56
        columns' worth of parameters.
        """
        c = self.core
        ax = getattr(plan, "axicon_xy", None)
        if ax is not None and self.has_axicon_slot:
            c.setMulticonfParam("ax0_x", 0, float(ax[0]))
            c.setMulticonfParam("ax0_y", 0, float(ax[1]))

        for prefix, slots, count in (("so", plan.shading, n_shade),
                                     ("bo", plan.blocking, n_block)):
            if not slots:
                raise ValueError(
                    f"heliostat {plan.heliostat_id} has no {prefix} slots to "
                    f"write; every slot must hold a real neighbour")
            for k in range(count):
                s = slots[min(k, len(slots) - 1)]
                name = f"{prefix}{k}"
                c.setMulticonfParam(f"{name}_x", 0, float(s.x_mm))
                c.setMulticonfParam(f"{name}_y", 0, float(s.y_mm))
                c.setMulticonfParam(f"{name}_az", 0, float(s.rot_az_deg))
                c.setMulticonfParam(f"{name}_el", 0, float(s.rot_el_deg))

    # Whether a full model re-init is needed after writing multiconfig params.
    # Expected to be False -- selecting the configuration should be enough --
    # but it is verified empirically by :meth:`check_reinit_needed` rather than
    # assumed, because getting it wrong is expensive in opposite directions:
    # unnecessary re-inits cost time on every one of ~29,000 traces, while a
    # missing one silently produces identical spots for every heliostat.
    reinit_after_params: bool = False

    def activate(self, config_num: int = 0, reinit: bool | None = None) -> None:
        self.core.setConfig(config_num)
        self._current_config = config_num
        if self.reinit_after_params if reinit is None else reinit:
            self.core.applyChangesAndInitModel()

    def check_reinit_needed(self, strategy, rays: int = 120000, verbose: bool = True) -> bool:
        """Is ``applyChangesAndInitModel`` needed for a parameter change to take?

        The question is not whether two traces of the same configuration differ
        -- they always do, because the source is Monte-Carlo -- but whether
        *switching heliostats* without a re-init actually moves the spot.

        So: trace A, then switch to a well-separated B without re-init, then
        trace B with re-init. If B-without-reinit still looks like A, the
        parameter write did not take and the re-init is mandatory. Monte-Carlo
        noise is measured directly, by tracing A twice, and used as the yardstick
        rather than assumed.
        """
        cfg = self.cfg
        az, el = 77.9, 50.1
        self.set_global_geometry()
        self.set_sun(az, el)

        a = (22300.0, -60000.0)
        b = (-48167.7, 61589.5)
        sol_a = strategy.solve(*a, az, el, cfg.geometry)
        sol_b = strategy.solve(*b, az, el, cfg.geometry)

        def run(x, y, sol, reinit):
            self.set_heliostat(x, y, sol, cfg.trace.bulk_config)
            self.activate(cfg.trace.bulk_config, reinit=reinit)
            r = self.trace(rays=rays)
            return r.xy_mm.mean(axis=0), r.rays_landed

        # Estimate the Monte-Carlo noise floor from repeated traces of the same
        # configuration, rather than from a single pair.
        a_runs = [run(*a, sol_a, True) for _ in range(3)]
        a_cent = [c for c, _ in a_runs]
        noise = max(
            float(np.linalg.norm(a_cent[i] - a_cent[j]))
            for i in range(len(a_cent)) for j in range(i + 1, len(a_cent))
        )

        c_b_no, n_b_no = run(*b, sol_b, False)   # switch WITHOUT reinit
        c_b_yes, n_b_yes = run(*b, sol_b, True)  # switch WITH reinit

        moved = float(np.linalg.norm(c_b_no - a_cent[0]))
        agrees = float(np.linalg.norm(c_b_no - c_b_yes))
        signal = float(np.linalg.norm(c_b_yes - a_cent[0]))

        floor = max(noise, 1e-6)
        # Decisive criterion: does the no-reinit result land where the reinit
        # result lands, to within the measured noise? The "moved" check only
        # confirms the two probe heliostats are distinguishable at all.
        distinguishable = moved > 2.0 * floor
        consistent = agrees < 1.5 * floor
        took_effect = distinguishable and consistent
        needed = not took_effect

        if verbose:
            for k, (c, n) in enumerate(a_runs, 1):
                print(f"  A trace {k}        : centroid ({c[0]:+8.2f}, {c[1]:+8.2f})  {n} rays")
            print(f"  B, no reinit     : centroid ({c_b_no[0]:+8.2f}, {c_b_no[1]:+8.2f})  {n_b_no} rays")
            print(f"  B, with reinit   : centroid ({c_b_yes[0]:+8.2f}, {c_b_yes[1]:+8.2f})  {n_b_yes} rays")
            print(f"  Monte-Carlo noise floor    : {noise:7.3f} mm  (3 traces of A)")
            print(f"  A -> B with reinit (signal): {signal:7.3f} mm")
            print(f"  A -> B no reinit  (moved)  : {moved:7.3f} mm"
                  f"   {'> 2x noise: heliostats distinguishable' if distinguishable else '<= 2x noise: INCONCLUSIVE'}")
            print(f"  B no-reinit vs B reinit    : {agrees:7.3f} mm"
                  f"   {'< 1.5x noise: agree' if consistent else '>= 1.5x noise: DIFFER'}")
            if not distinguishable:
                print("  -> INCONCLUSIVE: pick probe heliostats whose spots separate further,")
                print("     or raise `rays` to lower the noise floor.")
            print(f"  -> re-init is {'REQUIRED' if needed else 'NOT needed (skipping it)'}")

        self.reinit_after_params = needed
        return needed

    # -- tracing ----------------------------------------------------------
    def trace(self, rays: int | None = None, chunk_sizes=None) -> TraceResult:
        """Trace the active configuration and return receiver-plane x/y in mm.

        The ray budget is split into chunks because a single enormous
        ``traceRays`` call allocates the whole ray set at once. How many chunks
        that is comes from ``[trace] rays_per_trace`` (``--rays-per-trace`` on
        the command line), and it is a real cost: each chunk is one
        ``setRayDistributionCount1`` + ``traceRays`` + ``getRayPos`` round trip.
        Which of the fixed per-call cost and the marginal per-ray cost dominates
        has never been measured -- ``scripts/probe_ray_cost.py`` is the
        experiment that would settle it.
        """
        cfg = self.cfg
        if chunk_sizes is None:
            # config.chunk_plan, not local arithmetic: the chunks must sum
            # EXACTLY to the budget, and the same split is quoted by
            # `beamdown info` and by the GUI's derived call-count label.
            from .config import chunk_plan

            chunk_sizes = (
                cfg.trace.chunk_sizes if rays is None
                else chunk_plan(rays, cfg.trace.rays_per_trace)
            )

        collected = []
        emitted = 0
        for n in chunk_sizes:
            self.core.setRayDistributionCount1(self.seq, int(n))
            self.core.traceRays(self.seq, 0, 0)
            pos = np.array(self.core.getRayPos(self.seq, 0, 0, self.surface), copy=True)
            if pos.size == 0:
                continue
            emitted += pos.shape[1]
            good = ~np.isnan(pos).any(axis=0)
            collected.append(pos[:2, good].T.astype(np.float32))

        xy = (
            np.concatenate(collected, axis=0)
            if collected
            else np.empty((0, 2), dtype=np.float32)
        )
        return TraceResult(xy_mm=xy, rays_emitted=emitted, rays_landed=int(xy.shape[0]))

    def trace_heliostat(self, x_mm, y_mm, solution, config_num: int = 0) -> TraceResult:
        """Set up and trace one heliostat in a single call."""
        self.set_heliostat(x_mm, y_mm, solution, config_num)
        self.activate(config_num)
        return self.trace()

    # -- validation -------------------------------------------------------
    def self_test(self, strategy, verbose: bool = True) -> dict:
        """Confirm parameter writes actually reach the trace.

        Guards against the failure mode that would quietly ruin an entire sweep:
        ``setMulticonfParam`` succeeding but the traced result never changing,
        leaving 645 identical spots. Traces two well-separated heliostats and
        checks that the receiver footprints genuinely differ.
        """
        cfg = self.cfg
        self.set_global_geometry()
        self.set_sun(77.9, 50.1)

        probes = [(22300.0, -60000.0), (-48167.0, 61589.0)]
        out = []
        for x, y in probes:
            sol = strategy.solve(x, y, 77.9, 50.1, cfg.geometry)
            res = self.trace_heliostat(x, y, sol, cfg.trace.bulk_config)
            centroid = res.xy_mm.mean(axis=0) if res.rays_landed else np.array([np.nan, np.nan])
            out.append(
                {
                    "x_mm": x, "y_mm": y,
                    "rot_az": sol.rot_az_deg, "rot_el": sol.rot_el_deg,
                    "rays_landed": res.rays_landed,
                    "transmission": res.transmission,
                    "centroid_mm": centroid,
                }
            )
            if verbose:
                print(
                    f"  heliostat ({x/1000:+7.1f}, {y/1000:+7.1f}) m -> "
                    f"{res.rays_landed:7d} rays ({res.transmission:5.1%}), "
                    f"centroid ({centroid[0]:+8.2f}, {centroid[1]:+8.2f}) mm"
                )

        separation = float(np.linalg.norm(out[0]["centroid_mm"] - out[1]["centroid_mm"]))
        distinct = separation > 1.0 and all(o["rays_landed"] > 0 for o in out)
        if verbose:
            print(f"  centroid separation: {separation:.2f} mm -> "
                  f"{'PASS: parameters reach the trace' if distinct else 'FAIL: results are identical'}")
        return {"probes": out, "separation_mm": separation, "passed": distinct}
