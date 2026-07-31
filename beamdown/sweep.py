"""The sweep driver: trace every heliostat, one at a time, at every timestep.

Structure of one timestep
-------------------------
1. Solve pointing and shape for all heliostats (pure numpy, milliseconds).
2. Compute shading and blocking efficiencies (analytic, seconds for the field).
3. Trace each heliostat individually through Quadoa (the expensive part).
4. Bin, quantise, and write.

Steps 1, 2 and 4 are cheap; only step 3 needs a license, so it is the only part
that runs in the worker pool.

Worker pool
-----------
Quadoa sessions are expensive to create and limited by the USB HASP key, so a
small fixed number (``trace.n_workers``, 1-4) are created once as pool
initialisers and reused for every heliostat of every timestep. With
``n_workers = 1`` the pool is bypassed entirely, which keeps debugging simple
and avoids Windows spawn overhead.

Resumability
------------
A timestep is written atomically-ish and ``--resume`` skips any timestep already
present, so an interrupted multi-hour sweep continues rather than restarting.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import field as field_mod
from . import shading as shading_mod
from . import solar as solar_mod
from .secondary import get_strategy
from .store import INT16_MAX, RunStore, TimestepResult

_WORKER: dict = {}


# --------------------------------------------------------------------------
# Worker-side
# --------------------------------------------------------------------------

def _init_worker(config_path: str, lock=None, overrides: dict | None = None) -> None:
    """Create one persistent Quadoa session in this process.

    Session creation is serialised across workers by ``lock``. Two workers
    requesting a license seat at the same instant can collide even when seats
    are free, and because seat exhaustion is deliberately not retried (each
    retry pops a modal dialog), a collision would otherwise permanently retire
    a worker. Acquiring seats one at a time removes the race entirely.
    """
    from .config import apply_overrides, load_config
    from .session import QuadoaSession

    if isinstance(overrides, str):
        # Legacy initargs: before overrides existed, the third element was the
        # occluder model path alone. A driver started under the OLD code keeps
        # passing that tuple for its whole run, and Pool respawns a dead worker
        # with the ORIGINAL initargs against whatever sweep.py now says -- so a
        # respawn during a long in-flight sweep lands exactly here. Translate
        # rather than crash the run's replacement worker.
        overrides = {"trace": {"model_file": overrides}}

    cfg = load_config(config_path)
    # The worker builds its own cfg from disk, so EVERY override the driver
    # applied has to be replayed here -- a value set only on the driver's copy
    # does not reach the trace. This started as the --occluders model swap alone
    # (without it an occluder sweep silently traced the plain model) and is now
    # the whole override set, because the same hole applied to --rays: the driver
    # reported and recorded the overridden ray budget while the worker read
    # config.toml's and emitted that instead.
    #
    # Not re-validated here: the driver validated the same combination at startup
    # (beamdown.cli._load_with_overrides), and a warning raised once per worker
    # process would say the same thing four times.
    apply_overrides(cfg, overrides)
    _WORKER["cfg"] = cfg
    # Built from the whole cfg, after apply_overrides, so [optics] flat_mirrors
    # -- which only ever reaches this process through that replay -- is baked
    # into the strategy the worker actually solves with.
    _WORKER["strategy"] = get_strategy(cfg)

    if lock is not None:
        lock.acquire()
    try:
        session = QuadoaSession(cfg)
        session.set_global_geometry()
        _WORKER["session"] = session
        _WORKER["error"] = None
    except Exception as exc:
        # Never raise from a pool initialiser: multiprocessing.Pool responds by
        # respawning the worker forever, which would loop the license dialog.
        # Record the failure and let the health check retire this worker.
        _WORKER["session"] = None
        _WORKER["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        if lock is not None:
            time.sleep(0.4)  # let the seat register before the next worker asks
            lock.release()


def _ping(_=None):
    """Report whether this worker holds a usable session."""
    return (os.getpid(), _WORKER.get("session") is not None, _WORKER.get("error"))


def _open_session_patiently(config_path: str, wait_s: float, progress=print,
                            interval_s: float = 30.0,
                            overrides: dict | None = None) -> None:
    """Open the single worker session, waiting for a seat if necessary.

    Seats are released lazily -- a session closed seconds ago can still be held
    at the license-manager level -- so a sweep that may run for hours should not
    abort because a seat is thirty seconds from freeing. Unlike the worker-pool
    case, this issues one seat request at a time at a slow cadence, so it cannot
    produce a storm of license dialogs.
    """
    deadline = time.time() + max(0.0, wait_s)
    announced = False
    while True:
        _init_worker(config_path, None, overrides)
        if _WORKER.get("session") is not None:
            if announced:
                progress("  license seat acquired")
            return

        error = _WORKER.get("error", "")
        if "LicenseUnavailable" not in str(error) or time.time() >= deadline:
            raise RuntimeError(f"Could not open a Quadoa session: {error}")

        if not announced:
            progress(f"  no license seat free; waiting up to "
                     f"{max(0.0, wait_s)/60:.0f} min (checking every "
                     f"{interval_s:.0f}s). Close other Quadoa sessions to speed this up.")
            announced = True
        time.sleep(interval_s)


def make_pool(cfg, config_path: str, workers: int, progress=print,
              overrides: dict | None = None):
    """Start a worker pool, degrading if the license cannot supply every seat.

    Rather than probing the license ceiling up front -- which pops a modal
    dialog for every seat request that fails -- the pool is started at the
    requested size and each worker is asked whether it actually got a session.
    If some did not, the pool is rebuilt at the size that did work. A long sweep
    then proceeds at reduced parallelism instead of dying at startup.
    """
    import multiprocessing as mp

    if mp.current_process().name != "MainProcess":
        raise RuntimeError(
            "run_sweep() with workers > 1 was reached from a child process. On "
            "Windows the 'spawn' start method re-imports the calling module, so a "
            "script that calls run_sweep() at module level will recursively spawn "
            "pools until it grinds to a halt. Wrap the call:\n\n"
            "    if __name__ == '__main__':\n"
            "        run_sweep(cfg, ...)\n"
        )

    ctx = mp.get_context("spawn")
    for n in range(workers, 0, -1):
        lock = ctx.Manager().Lock() if n > 1 else None
        pool = ctx.Pool(n, initializer=_init_worker, initargs=(config_path, lock, overrides))
        # Enough pings to reach every worker with high probability.
        health = pool.map(_ping, range(n * 4), chunksize=1)
        healthy = {pid for pid, ok, _ in health if ok}
        errors = {err for _, ok, err in health if not ok and err}

        if len(healthy) >= n:
            progress(f"  worker pool: {n} session(s) ready")
            return pool, n

        progress(f"  worker pool: asked for {n}, {len(healthy)} got a license seat")
        for err in list(errors)[:1]:
            progress(f"    {err.splitlines()[0]}")
        pool.close()
        pool.join()
        time.sleep(1.0)  # let seats release before retrying smaller

    raise RuntimeError("Could not start even one Quadoa worker session")


def _trace_one(job):
    """Trace a single heliostat. Runs in a worker process."""
    hid, x_mm, y_mm, az, el, plan = job
    cfg = _WORKER["cfg"]
    session = _WORKER["session"]
    strategy = _WORKER["strategy"]
    if session is None:
        raise RuntimeError(f"Worker has no Quadoa session: {_WORKER.get('error')}")

    solution = strategy.solve(x_mm, y_mm, az, el, cfg.geometry)
    session.set_sun(az, el)
    if plan is not None:
        from .occluder_slots import N_BLOCK, N_SHADE

        session.set_occluders(plan, N_SHADE, N_BLOCK)
    result = session.trace_heliostat(x_mm, y_mm, solution, cfg.trace.bulk_config)

    # Drop rays outside the storable window *before* binning and quantising, so
    # the raw store and the binned counts describe exactly the same ray set.
    window = cfg.receiver.window_mm
    if result.rays_landed:
        inside = RunStore.inside_window(result.xy_mm, window)
        xy = result.xy_mm[inside]
        n_outside = int((~inside).sum())
        edges = cfg.receiver.edges
        counts, _, _ = np.histogram2d(xy[:, 1], xy[:, 0], bins=[edges, edges])
    else:
        xy = result.xy_mm
        n_outside = 0
        counts = np.zeros((cfg.receiver.grid_size, cfg.receiver.grid_size))

    return (
        hid,
        counts.astype(np.uint32),
        RunStore.quantise(xy, window),
        result.rays_emitted,
        int(xy.shape[0]),
        solution,
        n_outside,
    )


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def run_sweep(
    cfg,
    dates=None,
    heliostat_indices=None,
    resume: bool = True,
    workers: int | None = None,
    config_path: str | None = None,
    steps=None,
    license_wait_s: float = 600.0,
    occluders: bool = False,
    overrides: dict | None = None,
    progress=print,
) -> RunStore:
    """Run the full sweep and return the populated store.

    ``steps`` overrides the computed time grid, which is useful for smoke tests
    and for re-running a specific subset of timesteps.

    ``overrides`` is ``{section: {field: value}}`` as built by
    :func:`beamdown.cli._override_map` -- the config values this run was told to
    use instead of what config.toml says. They are applied to ``cfg`` here (a
    no-op when the caller already did it) *and* handed to every worker, which
    loads config.toml for itself and would otherwise trace the file's values
    while the driver reported the overridden ones. They are also written into the
    manifest, so a stored run records what it was actually asked to do.
    """
    from .config import apply_overrides

    workers = workers or cfg.trace.n_workers
    config_path = config_path or str(cfg.repo_root / "config.toml")
    # Idempotent -- the CLI has already applied these -- so that a caller which
    # passes overrides without applying them cannot end up with workers tracing
    # one config and the driver reporting another.
    overrides = {section: dict(values) for section, values in (overrides or {}).items()}
    apply_overrides(cfg, overrides)
    traced_secondary = False

    if occluders:
        # A different model: the one carrying occluder slots in the ray path.
        # Set on cfg so the worker's own load_config sees it too.
        model = cfg.repo_root / "models" / "heliostat_field_occluders.optx"
        if not model.exists():
            raise FileNotFoundError(
                f"{model} does not exist. Build it with "
                f"scripts/build_occluder_field_model.py before sweeping with "
                f"--occluders."
            )
        object.__setattr__(cfg.trace, "model_file", str(model))
        overrides.setdefault("trace", {})["model_file"] = str(model)
        # Does this model also carry the secondary's own shadow plane? Read it
        # from the file rather than assuming, so an older model without the
        # slot keeps its scalar treatment instead of losing the axicon
        # silently.
        traced_secondary = "ax0_x" in model.read_text(encoding="utf-8")

    full_field = field_mod.load_field(cfg)
    fld = full_field if heliostat_indices is None else full_field.subset(heliostat_indices)
    strategy = get_strategy(cfg)
    if steps is None:
        steps = solar_mod.build_time_grid(cfg, dates)

    store = RunStore(cfg.output_root, cfg=cfg, mode="w")
    store.write_manifest({
        "n_heliostats": len(fld),
        "heliostat_ids": fld.ids.tolist(),
        "timesteps": [s.key for s in steps],
        "dates": [str(d) for d in (dates or cfg.sweep.dates)],
        "workers": workers,
        # Whether neighbour shading and blocking are already in the stored ray
        # counts. Anything that weights those counts afterwards must read this:
        # applying eta_shade x eta_block to an occluder run charges the same
        # loss twice, which reads as a uniform few-percent deficit at every
        # aperture radius and looks deceptively like a real result.
        "occluders": bool(occluders),
        "traced_secondary": bool(traced_secondary),
        # What this run was told to use instead of config.toml. config.toml is
        # copied into the run directory too, but that copy is the file as it was,
        # not as the command line amended it -- and the file legitimately changes
        # between runs, so without this a stored run cannot say which ray budget
        # or which secondary layout produced it.
        "overrides": overrides,
    })

    keep_raw = cfg.storage.raw_rays != "none"
    min_el = min(s.solar_el_deg for s in steps)
    radius = shading_mod.search_radius_for(
        min_el, cfg.field.mirror_height_mm, cfg.field.mirror_width_mm
    )
    # Occlusion is computed over the WHOLE field, never over the traced subset.
    # A heliostat is shaded by the mirrors that are physically next to it,
    # whether or not those happen to be in this run. Taking neighbours from the
    # subset made a 25-heliostat downselect -- chosen by farthest-point sampling,
    # so deliberately spread out -- report almost no occlusion at all.
    neighbours = field_mod.neighbour_pairs(full_field, radius)
    secondary = shading_mod.secondary_body(cfg)
    # Row of each traced heliostat within the full field.
    full_row = {int(h): k for k, h in enumerate(full_field.ids)}
    take = np.array([full_row[int(h)] for h in fld.ids])

    progress(
        f"Sweep: {len(fld)} heliostats x {len(steps)} timesteps = "
        f"{len(fld) * len(steps):,} traces"
    )
    progress(
        # describe(), not name: it is the one string that says both which layout
        # is running AND whether the heliostats are flat, and a flat run must not
        # be mistakable for a focused one in the log.
        f"  secondary={strategy.describe()}  workers={workers}  "
        f"rays/heliostat={cfg.trace.rays_per_heliostat:,}  raw_rays={cfg.storage.raw_rays}"
    )
    progress(f"  neighbour search radius {radius/1000:.1f} m "
             f"(min sun elevation {min_el:.1f} deg)")
    if overrides:
        flat = ", ".join(f"{section}.{key}={value}"
                         for section, values in sorted(overrides.items())
                         for key, value in sorted(values.items()))
        progress(f"  overriding config.toml: {flat}")
    if occluders:
        progress(f"  occluders traced as geometry; secondary shadow "
                 f"{'traced too' if traced_secondary else 'still a scalar'}")

    pool = None
    if workers > 1:
        pool, workers = make_pool(cfg, config_path, workers, progress,
                                  overrides=overrides)
    else:
        _open_session_patiently(config_path, license_wait_s, progress,
                                overrides=overrides)

    t_start = time.time()
    try:
        for step_no, step in enumerate(steps, 1):
            if resume and store.has_timestep(step.key):
                progress(f"[{step_no}/{len(steps)}] {step.label} -- already done, skipping")
                continue

            t0 = time.time()
            all_solutions = [
                strategy.solve(float(full_field.x_mm[i]), float(full_field.y_mm[i]),
                               step.solar_az_deg, step.solar_el_deg, cfg.geometry)
                for i in range(len(full_field))
            ]
            geoms, aims = shading_mod.build_geometries(full_field, all_solutions, cfg)
            all_shade, all_block, all_sec = shading_mod.shading_blocking(
                geoms, aims, step.solar_az_deg, step.solar_el_deg, neighbours,
                secondary=secondary,
            )
            solutions = [all_solutions[k] for k in take]
            eta_shade, eta_block = all_shade[take], all_block[take]
            eta_secondary = all_sec[take]

            plans = None
            if occluders:
                from . import occluder_slots as slots_mod

                all_plans = slots_mod.plan_field(
                    geoms, aims, full_field.ids, neighbours,
                    shading_mod.sun_vector(step.solar_az_deg, step.solar_el_deg),
                    body=secondary,
                    # prime_focus has no body over the field and therefore no
                    # shadow plane to place; anything else does.
                    has_secondary=cfg.optics.secondary != "prime_focus",
                )
                plans = [all_plans[k] for k in take]
                over = slots_mod.overflow_report(plans)
                if over["overflowed"]:
                    progress(f"  WARNING: {over['overflowed']} heliostat(s) had more "
                             f"occluders than slots; worst dropped "
                             f"{over['worst_dropped_fraction']:.4f} of the aperture")

            jobs = [
                (int(fld.ids[i]), float(fld.x_mm[i]), float(fld.y_mm[i]),
                 step.solar_az_deg, step.solar_el_deg,
                 plans[i] if plans is not None else None)
                for i in range(len(fld))
            ]
            results = (
                pool.map(_trace_one, jobs, chunksize=4) if pool
                else [_trace_one(j) for j in jobs]
            )

            store.write_timestep(
                _assemble(cfg, fld, step, results, solutions,
                          eta_shade, eta_block, eta_secondary, keep_raw,
                          occluders=occluders,
                          traced_secondary=traced_secondary)
            )

            dt = time.time() - t0
            done = step_no
            eta_total = (time.time() - t_start) / done * (len(steps) - done)
            progress(
                f"[{step_no}/{len(steps)}] {step.label}  el={step.solar_el_deg:5.1f}  "
                f"{dt:6.1f}s  ({dt/max(1,len(fld))*1000:5.1f} ms/heliostat)  "
                f"eta {eta_total/60:5.1f} min"
            )
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    progress(f"Sweep complete in {(time.time() - t_start)/60:.1f} min -> {store.root}")
    return store


def _assemble(cfg, fld, step, results, solutions, eta_shade, eta_block,
              eta_secondary, keep_raw, occluders: bool = False,
              traced_secondary: bool = False):
    """Turn worker output into a :class:`TimestepResult`."""
    from . import metrics as metrics_mod

    order = {int(hid): k for k, (hid, *_rest) in enumerate(results)}
    grid = cfg.receiver.grid_size

    counts = np.zeros((len(fld), grid, grid), dtype=np.uint32)
    rays_chunks: list[np.ndarray] = []
    index = np.zeros((len(fld), 3), dtype=np.int64)
    rows = []
    cursor = 0

    for i in range(len(fld)):
        hid = int(fld.ids[i])
        hid_, cnt, quant, emitted, landed, solution, n_outside = results[order[hid]]
        counts[i] = cnt

        if keep_raw:
            rays_chunks.append(quant)
            index[i] = (hid, cursor, quant.shape[0])
            cursor += quant.shape[0]

        # With the occluders in the ray path, neighbour shading and blocking are
        # already in the ray counts -- applying them again as a scalar would
        # square the loss. The secondary is the exception: it shades the mirror
        # from sunlight, but it sits in the *outgoing* leg of the sequence, so
        # its shadow is not traced and stays a scalar. That is a fair
        # approximation for it alone, because its shadow is a hard-edged 30 m
        # disc: a heliostat is essentially all in it or all out, unlike a
        # neighbour's shadow, which carves a patch out of one mirror.
        if occluders:
            eff = 1.0 if traced_secondary else float(eta_secondary[i])
        else:
            eff = float(eta_shade[i] * eta_block[i])
        xy = quant.astype(np.float32) * (cfg.receiver.window_mm / INT16_MAX)
        row = metrics_mod.spot_metrics(xy, emitted, cfg, efficiency=eff)
        row.update({
            "date": str(step.date),
            "hour": step.hour,
            "timestep": step.key,
            "heliostat_id": hid,
            "x_m": float(fld.x_m[i]),
            "y_m": float(fld.y_m[i]),
            "radius_m": float(fld.radius_mm[i] / 1000.0),
            "solar_az_deg": step.solar_az_deg,
            "solar_el_deg": step.solar_el_deg,
            "rot_az_deg": solution.rot_az_deg,
            "rot_el_deg": solution.rot_el_deg,
            "aoi_deg": solution.aoi_deg,
            "cosine_efficiency": solution.cosine_efficiency,
            # eta_shade is the whole shading loss -- neighbours and the secondary
            # unioned. eta_secondary is the secondary alone, carried only so its
            # share can be read off; multiplying the two would double-count the
            # mirror that both of them shade.
            "eta_shade": float(eta_shade[i]),
            "eta_secondary": float(eta_secondary[i]),
            "eta_block": float(eta_block[i]),
            "rays_outside_window": int(n_outside),
        })
        rows.append(row)

    lead = ["date", "hour", "timestep", "heliostat_id", "x_m", "y_m", "radius_m"]
    frame = pd.DataFrame(rows)
    frame = frame[lead + [c for c in frame.columns if c not in lead]]

    return TimestepResult(
        key=step.key,
        date=str(step.date),
        hour=step.hour,
        solar_az_deg=step.solar_az_deg,
        solar_el_deg=step.solar_el_deg,
        heliostat_ids=fld.ids,
        rays_emitted=int(cfg.trace.rays_per_heliostat),
        counts=counts,
        rays=np.concatenate(rays_chunks, axis=0) if rays_chunks else None,
        index=index if keep_raw else None,
        rows=frame,
    )


def check_worker_capacity(cfg, config_path: str | None = None, max_workers: int = 4,
                          confirm: bool = False):
    """Find how many concurrent Quadoa sessions the license actually allows.

    **This probe deliberately exhausts the license**, and the HASP runtime shows
    a modal "H0038 / too many concurrent users" dialog for each seat request
    that fails. It therefore stops at the first failure and must be opted into
    with ``confirm=True``.

    Run it from a clean Python process with no other Quadoa sessions open --
    including the interactive Quadoa application -- or the answer is just
    "how many seats were left over", not the real capacity.
    """
    import multiprocessing as mp

    if not confirm:
        print("  skipped: pass confirm=True to probe (this pops modal license dialogs)")
        return None

    from .session import running_quadoa_processes

    others = [p for p in running_quadoa_processes() if "quadoa" in p["name"].lower()]
    if others:
        print(f"  WARNING: {len(others)} Quadoa process(es) already running; "
              f"they hold seats and will understate capacity:")
        for p in others:
            print(f"    {p['name']} (pid {p['pid']})")

    config_path = config_path or str(cfg.repo_root / "config.toml")
    ctx = mp.get_context("spawn")
    usable = 0

    for n in range(1, max_workers + 1):
        pool = None
        try:
            pool = ctx.Pool(n, initializer=_init_worker, initargs=(config_path,))
            jobs = [(i, 22300.0, -60000.0, 77.9, 50.1) for i in range(n)]
            t0 = time.time()
            out = pool.map(_trace_one, jobs)
            dt = time.time() - t0
            landed = [r[4] for r in out]
            if all(l > 0 for l in landed):
                usable = n
                print(f"  {n} worker(s): OK   {dt:5.2f}s, rays landed {landed}")
            else:
                print(f"  {n} worker(s): BAD  some traces returned no rays: {landed}")
                break
        except Exception as exc:
            print(f"  {n} worker(s): FAILED -- {type(exc).__name__}")
            print(f"    stopping here; {usable} concurrent session(s) confirmed working")
            break
        finally:
            if pool is not None:
                pool.close()
                pool.join()
            time.sleep(1.0)  # let seats be released before the next probe

    print(f"\n  Usable concurrent sessions: {usable}")
    if usable:
        print(f"  Set trace.n_workers = {usable} in config.toml")
    return usable
