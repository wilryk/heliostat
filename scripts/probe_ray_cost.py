"""What does a ray actually cost? The first honest measurement.

    python scripts/probe_ray_cost.py --dry-run          # no licence, no trace
    python scripts/probe_ray_cost.py 2>&1 | tee analysis_output/probe_ray_cost.log

**DO NOT RUN THIS WHILE A SWEEP HOLDS THE LICENCE SEAT.** There is ONE USB HASP
seat. A failed request for a second has been measured to leak the first, leaving
the running sweep's next session waiting with zero seats and a modal H0038 dialog
for a human to dismiss. The script refuses to start if any
``analysis_output/.*.lock`` directory exists; ``--force`` overrides that, and you
should not need it. Check first::

    ls -d analysis_output/.*.lock

**Why this is the first honest measurement.** Until the override-propagation fix,
``--rays`` never reached the sweep workers: a worker calls ``load_config`` itself,
from disk, and the driver's overridden copy stayed in the driver. Every stored run
in ``analysis_output/`` therefore traced ``config.toml``'s 120,000 rays no matter
what the command line asked for -- including full8, whose manifest records 60,000.
So there is **no** empirical data anywhere in this repository about how trace time
scales with ray count, and any claim of a measured ray-count saving predates this
script and is false. This script is how that gets fixed.

**The question.** One heliostat trace is not one operation. It is:

  (a) a fixed per-heliostat cost -- seven ``setMulticonfParam`` writes plus
      ``setConfig`` (plus 56 more writes with ``--occluders``);
  (b) a fixed per-``traceRays``-call cost, paid ``n_chunks`` times, where
      ``n_chunks = ceil(rays_per_heliostat / rays_per_trace)``;
  (c) a marginal per-ray cost -- the trace itself, the ``getRayPos``
      marshalling, and the numpy copy + NaN filter.

Only (c) falls when you ask for fewer rays. If (a) + (b) dominate then halving
the ray budget does not halve the run, and the hypothesis worth testing -- the
user's -- is that 12,000 rays costs nearly what 6,000 does.

**How it measures.** Each of the three is timed **separately** with
``time.perf_counter``, never inferred by subtraction:

  * (a) directly: the clock brackets ``set_heliostat`` + ``activate``.
  * (b) and (c) from a regression over the **individual per-call timings**
    (``t_call ~ b + c * chunk_rays``), not over whole-heliostat totals -- so the
    intercept is a per-call cost measured across calls of genuinely different
    sizes, rather than a residual left over after guessing something else.
  * the trace side (``setRayDistributionCount1`` + ``traceRays``) and the fetch
    side (``getRayPos`` + numpy) are clocked apart, so "rays are expensive" can
    be attributed to the tracer or to the marshalling.

The fitted model is then cross-checked against the measured per-cell totals; a
large residual means the model is wrong and its extrapolations should not be
believed.

**How it avoids the usual ways a benchmark lies.**

  * One session on the 25-heliostat downselect at ONE timestep: minutes, not
    hours, and no licence churn.
  * The first cell is repeated up front as a discarded warm-up, because the first
    trace pays one-off model initialisation.
  * **Medians, not means.** This machine has 4 physical cores, Quadoa's tracer
    already runs at ~4x parallelism inside one session, and a previous benchmark
    here was ruined by concurrent numpy work. A single descheduled call moves a
    mean and does not move a median.
  * Every cell traces the same heliostats, so geometry differences cancel between
    cells rather than being confounded with the grid.

**Speed is never reported alone.** Per-heliostat spot noise scales as
1/sqrt(rays), and field power sums 645 independent traces, so the decision table
prints the expected relative noise on field power beside every wall-clock number
-- against the ~0.46% energy-integration residual, which is the error that
actually limits the annual number.
"""

from __future__ import annotations

import argparse
import math
import statistics
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

# -- established figures, not invented here --------------------------------
#
# Per-heliostat Monte-Carlo noise floor at 120,000 rays, measured in this repo by
# repeated traces of the SAME configuration (session.check_reinit_needed and
# scripts/verify_occluder_trace.py both quote it). Relative, on a spot centroid /
# landed-count basis.
NOISE_FLOOR_120K = 0.0067
NOISE_FLOOR_RAYS = 120000
# Field power sums this many independent traces, so the field-level relative
# noise is the per-heliostat figure divided by sqrt(N).
N_FIELD = 645
# Timesteps in the full annual grid (12 monthly dates, uniform daylight window)
# -- analysis_output/full8's manifest lists exactly this many.
N_TIMESTEPS = 161
# The residual between the two independent annual-energy integrations
# (scripts/report_energy.py). This is the error the annual number actually
# carries, and the yardstick any Monte-Carlo noise should be judged against.
ENERGY_RESIDUAL = 0.0046

# Aggregate seconds per trace over full7/full8 at 120,000 rays. Used ONLY to
# estimate how long this probe itself will take, never as a result.
NOMINAL_S_PER_TRACE_120K = 0.55


def field_noise(rays: int) -> float:
    """Expected relative Monte-Carlo noise on FIELD power at this ray budget.

    Per-heliostat noise scales as 1/sqrt(rays) from the measured floor; the field
    sums N_FIELD independent traces, so their noise adds in quadrature and the
    relative figure falls by sqrt(N_FIELD).
    """
    per_heliostat = NOISE_FLOOR_120K * math.sqrt(NOISE_FLOOR_RAYS / float(rays))
    return per_heliostat / math.sqrt(N_FIELD)


def chunk_sizes_for(rays: int, n_chunks: int) -> list[int]:
    """Split ``rays`` into exactly ``n_chunks`` calls, summing exactly to it.

    Uses the package's own splitter so the probe measures the split the sweep
    would actually make, remainder and all.
    """
    from beamdown.config import chunk_plan

    per = -(-int(rays) // int(n_chunks))          # ceil
    return chunk_plan(int(rays), per)


# --------------------------------------------------------------------------
# The fit
# --------------------------------------------------------------------------

def fit_line(x, y) -> tuple[float, float]:
    """Least-squares ``y = a + b*x``. Returns ``(a, b)``.

    Two points are enough and are the common case here (a grid of a few distinct
    chunk sizes), so this stays plain lstsq rather than anything with a
    covariance estimate the sample size could not support.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or np.allclose(x, x[0]):
        # Cannot separate intercept from slope: report the mean as the intercept
        # and no slope, which the caller will see as an obviously degenerate fit.
        return float(np.median(y)), 0.0
    design = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return float(coef[0]), float(coef[1])


class CostModel:
    """``T_heliostat = A + B * n_chunks + C * rays``, fitted from measured parts."""

    def __init__(self, a_setup: float, b_call: float, c_ray: float,
                 c_ray_trace: float = float("nan"),
                 c_ray_fetch: float = float("nan"),
                 b_call_trace: float = float("nan"),
                 b_call_fetch: float = float("nan")):
        self.A = a_setup          # s per heliostat, MEASURED directly
        self.B = b_call           # s per traceRays call, fitted intercept
        self.C = c_ray            # s per ray, fitted slope
        self.C_trace = c_ray_trace
        self.C_fetch = c_ray_fetch
        self.B_trace = b_call_trace
        self.B_fetch = b_call_fetch

    def predict(self, rays: int, n_chunks: int) -> float:
        return self.A + self.B * n_chunks + self.C * rays

    def sweep_hours(self, rays: int, n_chunks: int,
                    n_heliostats: int = N_FIELD,
                    n_timesteps: int = N_TIMESTEPS) -> float:
        return self.predict(rays, n_chunks) * n_heliostats * n_timesteps / 3600.0

    def halving_saving(self, rays: int, n_chunks: int = 1) -> float:
        """Fraction of wall clock saved by halving the ray budget from ``rays``."""
        full = self.predict(rays, n_chunks)
        half = self.predict(rays // 2, n_chunks)
        return (full - half) / full if full > 0 else 0.0

    def overhead_dominated_below(self, n_chunks: int = 1) -> float:
        """Ray budget below which halving the rays saves less than 10%.

        ``saving(r) = (C r / 2) / (A + B n + C r) = 0.10``  =>
        ``0.4 C r = 0.1 (A + B n)``  =>  ``r* = 0.25 (A + B n) / C``.
        Below r* the run is overhead-dominated and buying precision with rays is
        nearly free.
        """
        if self.C <= 0:
            return float("inf")
        return 0.25 * (self.A + self.B * n_chunks) / self.C


def fit_cost_model(rows: list[dict]) -> CostModel:
    """Fit from the raw timing rows, measuring what can be measured.

    ``A`` is not fitted at all -- it is the median of the directly clocked
    setup times. ``B`` and ``C`` come from a regression over the **individual
    per-call timings** keyed by that call's own ray count, which is why the
    intercept means "cost of making a call" rather than "whatever was left".
    """
    live = [r for r in rows if not r["warmup"]]
    if not live:
        raise ValueError("no non-warmup rows to fit")

    a_setup = float(statistics.median(r["t_setup_s"] for r in live))

    # Pool every individual chunk timing by its ray count, then take the median
    # per distinct size before fitting: an unequally sampled grid would otherwise
    # let the most-repeated chunk size dominate the least squares, and the median
    # is what survives a descheduled call.
    by_size_trace: dict[int, list[float]] = {}
    by_size_fetch: dict[int, list[float]] = {}
    for r in live:
        for n, t_tr, t_fe in zip(r["chunk_rays"], r["t_trace_chunks"],
                                 r["t_fetch_chunks"]):
            by_size_trace.setdefault(int(n), []).append(float(t_tr))
            by_size_fetch.setdefault(int(n), []).append(float(t_fe))

    sizes = sorted(by_size_trace)
    med_trace = [statistics.median(by_size_trace[n]) for n in sizes]
    med_fetch = [statistics.median(by_size_fetch[n]) for n in sizes]
    med_total = [t + f for t, f in zip(med_trace, med_fetch)]

    b_tr, c_tr = fit_line(sizes, med_trace)
    b_fe, c_fe = fit_line(sizes, med_fetch)
    b_all, c_all = fit_line(sizes, med_total)
    return CostModel(a_setup, b_all, c_all, c_tr, c_fe, b_tr, b_fe)


def cell_medians(rows: list[dict]) -> dict[tuple[int, int], dict]:
    """Median timings per (rays, n_chunks) cell, warm-up excluded."""
    cells: dict[tuple[int, int], list[dict]] = {}
    for r in rows:
        if r["warmup"]:
            continue
        cells.setdefault((r["rays"], r["n_chunks"]), []).append(r)

    out = {}
    for key, group in sorted(cells.items()):
        out[key] = {
            "n": len(group),
            "setup": statistics.median(r["t_setup_s"] for r in group),
            "trace": statistics.median(r["t_trace_s"] for r in group),
            "fetch": statistics.median(r["t_fetch_s"] for r in group),
            "total": statistics.median(r["t_total_s"] for r in group),
            "landed": statistics.median(r["rays_landed"] for r in group),
        }
    return out


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def report(rows: list[dict], model: CostModel, args) -> None:
    cells = cell_medians(rows)

    print()
    print("  measured medians per cell "
          f"({args.per_cell} heliostats x {args.repeats} rep(s), warm-up discarded)")
    print("  " + "-" * 84)
    print(f"  {'rays':>8} {'calls':>6} {'chunk':>8} {'setup':>9} {'trace':>9} "
          f"{'fetch':>9} {'total':>9} {'model':>9} {'resid':>8}")
    print(f"  {'':>8} {'':>6} {'rays':>8} {'ms':>9} {'ms':>9} {'ms':>9} "
          f"{'ms':>9} {'ms':>9} {'%':>8}")
    print("  " + "-" * 84)
    resid = []
    for (rays, n_chunks), m in sorted(cells.items()):
        predicted = model.predict(rays, n_chunks)
        pct = 100.0 * (m["total"] - predicted) / predicted if predicted else float("nan")
        resid.append(abs(pct))
        print(f"  {rays:>8,} {n_chunks:>6} {rays // n_chunks:>8,} "
              f"{1e3 * m['setup']:>9.1f} {1e3 * m['trace']:>9.1f} "
              f"{1e3 * m['fetch']:>9.1f} {1e3 * m['total']:>9.1f} "
              f"{1e3 * predicted:>9.1f} {pct:>+8.1f}")
    print("  " + "-" * 84)

    print()
    print("  cost model   T_heliostat = A + B*n_chunks + C*rays")
    print(f"    A = {1e3 * model.A:8.2f} ms per heliostat      "
          f"(MEASURED: set_heliostat + activate, 7 setMulticonfParam + setConfig"
          f"{' + 56 occluder writes' if args.occluders else ''})")
    print(f"    B = {1e3 * model.B:8.2f} ms per traceRays call "
          f"(fitted intercept over per-call timings; "
          f"trace side {1e3 * model.B_trace:.2f}, fetch side {1e3 * model.B_fetch:.2f})")
    print(f"    C = {1e6 * model.C:8.3f} us per ray            "
          f"(fitted slope; trace side {1e6 * model.C_trace:.3f}, "
          f"fetch side {1e6 * model.C_fetch:.3f})")
    print(f"        i.e. {1e3 * model.C * 120000:.1f} ms of pure ray cost at "
          f"120,000 rays, {1e3 * model.C * 12000:.1f} ms at 12,000")
    if resid:
        worst = max(resid)
        print(f"    model vs measured cell medians: worst residual {worst:+.1f}%"
              + ("  -- the linear model holds" if worst < 10 else
                 "  -- LARGE: do not trust the extrapolations below"))

    fixed_at_120k = model.A + model.B * 2
    variable_at_120k = model.C * 120000
    share = fixed_at_120k / (fixed_at_120k + variable_at_120k)
    print()
    print(f"  at the config.toml default (120,000 rays in 2 calls of 60,000): "
          f"{100 * share:.0f}% of a heliostat's time is fixed cost, "
          f"{100 * (1 - share):.0f}% is per-ray")
    r_star = model.overhead_dominated_below(n_chunks=1)
    print(f"  overhead-dominated below {r_star:,.0f} rays/heliostat "
          f"(halving the budget there saves <10% of wall clock)")
    if args.occluders:
        print("  A includes the 56 occluder-slot writes, as an --occluders sweep pays")
    else:
        print("  NOTE: measured WITHOUT --occluders, so A understates an occluder")
        print("        sweep (full8-style) by 56 extra setMulticonfParam writes per")
        print("        heliostat. Re-run with --occluders to predict one of those.")

    # -- the decision table -------------------------------------------------
    print()
    print(f"  decision table -- full sweep, {args.n_heliostats} heliostats x "
          f"{args.n_timesteps} timesteps = "
          f"{args.n_heliostats * args.n_timesteps:,} traces")
    print("  " + "-" * 92)
    print(f"  {'rays':>8} {'calls':>6} {'s/heliostat':>12} {'sweep h':>9} "
          f"{'vs 120k':>9} {'halving':>9} {'field noise':>12} {'vs 0.46%':>10}")
    print(f"  {'':>8} {'':>6} {'predicted':>12} {'':>9} {'saves':>9} "
          f"{'saves':>9} {'on power':>12} {'residual':>10}")
    print("  " + "-" * 92)
    base = model.sweep_hours(120000, chunks_for_policy(120000, args),
                             args.n_heliostats, args.n_timesteps)
    for rays in sorted(set(args.rays) | {120000}):
        n_chunks = chunks_for_policy(rays, args)
        hours = model.sweep_hours(rays, n_chunks, args.n_heliostats, args.n_timesteps)
        noise = field_noise(rays)
        print(f"  {rays:>8,} {n_chunks:>6} {model.predict(rays, n_chunks):>12.3f} "
              f"{hours:>9.1f} {100 * (1 - hours / base) if base else 0:>8.0f}% "
              f"{100 * model.halving_saving(rays, n_chunks):>8.0f}% "
              f"{100 * noise:>11.3f}% {noise / ENERGY_RESIDUAL:>9.2f}x")
    print("  " + "-" * 92)
    print(f"  field noise = {NOISE_FLOOR_120K} / sqrt({N_FIELD}) scaled by "
          f"sqrt({NOISE_FLOOR_RAYS:,}/rays): the measured per-heliostat floor, "
          f"averaged over")
    print(f"  {N_FIELD} independent traces. 'vs residual' compares it with the "
          f"{100 * ENERGY_RESIDUAL:.2f}% energy-integration residual, which is the")
    print("  error the annual number actually carries -- a ray budget whose noise "
          "is well under 1.00x")
    print("  is not what limits the answer, and buying more rays there buys nothing.")
    print()
    print("  Read the two halves together. The left half says what speed costs; "
          "the right half says")
    print("  what precision costs. Neither is a recommendation on its own.")


def chunks_for_policy(rays: int, args) -> int:
    """How many calls a real sweep would make at this budget."""
    if args.policy_chunk <= 0:
        return 1
    return max(1, math.ceil(rays / args.policy_chunk))


def write_csv(rows: list[dict], path: Path) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["rep", "cell", "warmup", "rays", "n_chunks", "chunk_rays",
              "heliostat_id", "rays_landed", "t_setup_s", "t_trace_s",
              "t_fetch_s", "t_total_s", "t_trace_chunks", "t_fetch_chunks"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            out = dict(r)
            # Per-call detail survives into the CSV as a ';'-joined list, so the
            # regression can be re-run or re-weighted later without re-tracing.
            out["chunk_rays"] = ";".join(str(int(n)) for n in r["chunk_rays"])
            out["t_trace_chunks"] = ";".join(f"{t:.6f}" for t in r["t_trace_chunks"])
            out["t_fetch_chunks"] = ";".join(f"{t:.6f}" for t in r["t_fetch_chunks"])
            w.writerow({k: out.get(k) for k in fields})


# --------------------------------------------------------------------------
# Measurement
# --------------------------------------------------------------------------

def build_grid(args) -> list[tuple[int, int]]:
    """The (rays, n_chunks) cells to measure, absurd ones dropped."""
    grid = []
    for rays in sorted(set(args.rays)):
        for n_chunks in sorted(set(args.chunks)):
            sizes = chunk_sizes_for(rays, n_chunks)
            if len(sizes) != n_chunks:
                continue                      # cannot be split that many ways
            if min(sizes) < args.min_chunk:
                continue                      # chunks too small to be meaningful
            grid.append((rays, n_chunks))
    return grid


def licence_locks(repo_root: Path) -> list[Path]:
    """Sweep lock directories, which mean the one HASP seat is spoken for."""
    out = Path(repo_root, "analysis_output")
    return sorted(p for p in out.glob(".*.lock") if p.is_dir()) if out.is_dir() else []


def measure(args) -> list[dict]:
    """Trace the grid on one session and return the raw timing rows."""
    from beamdown import field as F
    from beamdown import solar as solar_mod
    from beamdown.config import load_config
    from beamdown.secondary import get_strategy
    from beamdown.session import QuadoaSession

    cfg = load_config(args.config)

    plans = None
    OS_N_SHADE = OS_N_BLOCK = 0        # bound properly below when --occluders
    if args.occluders:
        model = cfg.repo_root / "models" / "heliostat_field_occluders.optx"
        if not model.exists():
            raise SystemExit(f"{model} does not exist; build it or drop --occluders")
        object.__setattr__(cfg.trace, "model_file", str(model))

    full_field = F.load_field(cfg)
    idx, prov = F.load_or_build_downselect(cfg, full_field)
    fld = full_field.subset(idx)
    strategy = get_strategy(cfg)

    import datetime as _dt

    date = _dt.date.fromisoformat(args.date)
    steps = solar_mod.build_time_grid(cfg, [date])
    step = steps[args.step_index if args.step_index is not None else len(steps) // 2]

    # Every cell traces the SAME heliostats, so per-heliostat geometry
    # differences cancel between cells instead of being confounded with the grid.
    n_take = min(args.per_cell, len(fld))
    picks = list(range(n_take))
    solutions = [strategy.solve(float(fld.x_mm[i]), float(fld.y_mm[i]),
                                step.solar_az_deg, step.solar_el_deg, cfg.geometry)
                 for i in picks]

    if args.occluders:
        from beamdown import occluder_slots as OS
        from beamdown import shading as S

        radius = S.search_radius_for(step.solar_el_deg, cfg.field.mirror_height_mm,
                                     cfg.field.mirror_width_mm)
        neighbours = F.neighbour_pairs(full_field, radius)
        all_sols = [strategy.solve(float(full_field.x_mm[i]), float(full_field.y_mm[i]),
                                   step.solar_az_deg, step.solar_el_deg, cfg.geometry)
                    for i in range(len(full_field))]
        geoms, aims = S.build_geometries(full_field, all_sols, cfg)
        field_plans = OS.plan_field(geoms, aims, full_field.ids, neighbours,
                                    S.sun_vector(step.solar_az_deg, step.solar_el_deg),
                                    body=S.secondary_body(cfg))
        row_of = {int(h): k for k, h in enumerate(full_field.ids)}
        plans = [field_plans[row_of[int(fld.ids[i])]] for i in picks]
        OS_N_SHADE, OS_N_BLOCK = OS.N_SHADE, OS.N_BLOCK

    grid = build_grid(args)
    print(f"  heliostats: {n_take} of {len(fld)} ({prov})")
    print(f"  timestep:   {step.key}  az {step.solar_az_deg:.1f} "
          f"el {step.solar_el_deg:.1f}")
    print(f"  model:      {Path(cfg.trace.model_file).name}")
    print(f"  grid:       {len(grid)} cells x {n_take} heliostats "
          f"x {args.repeats} rep(s) = {len(grid) * n_take * args.repeats} traces")
    print()

    rows: list[dict] = []
    session = QuadoaSession(cfg)
    try:
        seq, surface = session.seq, session.surface
        core = session.core
        session.set_global_geometry()
        session.set_sun(step.solar_az_deg, step.solar_el_deg)

        def trace_cell(rays: int, n_chunks: int, rep: int, cell: int,
                       warmup: bool) -> None:
            sizes = chunk_sizes_for(rays, n_chunks)
            for k, i in enumerate(picks):
                # (a) fixed per-heliostat cost, clocked on its own -- INCLUDING
                # the occluder-slot writes when there are any, because a
                # full8-style sweep pays those 56 setMulticonfParam calls once
                # per heliostat exactly like the seven pointing ones. Timing
                # them outside this bracket would have made A look identical
                # with and without --occluders, which is the whole reason the
                # flag exists.
                t0 = time.perf_counter()
                if plans is not None:
                    session.set_occluders(plans[k], OS_N_SHADE, OS_N_BLOCK)
                session.set_heliostat(float(fld.x_mm[i]), float(fld.y_mm[i]),
                                      solutions[k], cfg.trace.bulk_config)
                session.activate(cfg.trace.bulk_config)
                t_setup = time.perf_counter() - t0

                # (b) and (c), per call, clocked apart. This mirrors
                # session.trace() exactly rather than calling it, because
                # trace() reports one number for the whole loop and the whole
                # point here is the split.
                t_trace_chunks, t_fetch_chunks = [], []
                landed = 0
                emitted = 0
                for n in sizes:
                    t0 = time.perf_counter()
                    core.setRayDistributionCount1(seq, int(n))
                    core.traceRays(seq, 0, 0)
                    t_trace_chunks.append(time.perf_counter() - t0)

                    t0 = time.perf_counter()
                    pos = np.array(core.getRayPos(seq, 0, 0, surface), copy=True)
                    if pos.size:
                        # pos is (2|3, N) BEFORE the NaN filter, so N is what the
                        # sequence actually emitted. Recording it is not
                        # bookkeeping: setRayDistributionCount1 is a literal ray
                        # count on some sequences and a PER-AXIS GRID DENSITY on
                        # others in this very model, and a grid sequence turns a
                        # request for 12,000 into ~0.785*12000^2 rays. Every
                        # per-ray cost here is divided by the requested count, so
                        # on a grid sequence the whole fit would be wrong by
                        # four orders of magnitude while still looking tidy.
                        # _check_ray_semantics below asserts emitted == requested.
                        emitted += int(pos.shape[1])
                        good = ~np.isnan(pos).any(axis=0)
                        xy = pos[:2, good].T.astype(np.float32)
                        landed += int(xy.shape[0])
                    t_fetch_chunks.append(time.perf_counter() - t0)

                t_trace = sum(t_trace_chunks)
                t_fetch = sum(t_fetch_chunks)
                rows.append(dict(
                    rep=rep, cell=cell, warmup=warmup, rays=rays,
                    n_chunks=n_chunks, chunk_rays=list(sizes),
                    heliostat_id=int(fld.ids[i]), rays_landed=landed,
                    rays_emitted=emitted,
                    t_setup_s=t_setup, t_trace_s=t_trace, t_fetch_s=t_fetch,
                    t_total_s=t_setup + t_trace + t_fetch,
                    t_trace_chunks=t_trace_chunks, t_fetch_chunks=t_fetch_chunks,
                ))
                if emitted != rays:
                    ratio = emitted / rays if rays else float("nan")
                    raise SystemExit(
                        f"\nABORTING: sequence {seq} emitted {emitted:,} rays when "
                        f"{rays:,} were requested ({ratio:.4g}x).\n"
                        "setRayDistributionCount1 is a literal ray count on some "
                        "sequences of this model and a per-axis GRID DENSITY on "
                        "others; a grid sequence emits about 0.785*n^2 rays.\n"
                        "Every per-ray number this script fits divides by the "
                        "REQUESTED count, so continuing would produce a tidy, "
                        "confident and completely wrong cost model.\n"
                        f"Fix the sequence (config.toml [trace] analysis_seq, "
                        f"currently {seq}) or convert the request to a density "
                        "before trusting any output."
                    )

        # Warm-up: the very first trace pays one-off model initialisation, and
        # folding that into a real cell would inflate whichever cell happened to
        # go first. Run the first cell once and throw it away.
        warm_rays, warm_chunks = grid[0]
        print(f"  warm-up (discarded): {warm_rays:,} rays in {warm_chunks} call(s)",
              flush=True)
        trace_cell(warm_rays, warm_chunks, rep=0, cell=-1, warmup=True)

        t_started = time.perf_counter()
        for rep in range(1, args.repeats + 1):
            for cell, (rays, n_chunks) in enumerate(grid):
                t0 = time.perf_counter()
                trace_cell(rays, n_chunks, rep=rep, cell=cell, warmup=False)
                dt = time.perf_counter() - t0
                print(f"  rep {rep}  cell {cell + 1}/{len(grid)}  "
                      f"{rays:>7,} rays / {n_chunks} call(s) -> "
                      f"{1e3 * dt / n_take:7.1f} ms per heliostat", flush=True)
        print(f"\n  traced {len(rows)} heliostats in "
              f"{time.perf_counter() - t_started:.1f} s")
    finally:
        session.close()
    return rows


def synthetic_rows(args, a=0.120, b=0.045, c=4.0e-6, seed=0) -> list[dict]:
    """Fake timings from a KNOWN A, B, C, for --dry-run.

    Exists so the reporting and the fit can be exercised -- and shown to recover
    the truth -- without a licence seat. Nothing here is a measurement, and
    --dry-run says so on every line it prints.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for rep in range(1, args.repeats + 1):
        for cell, (rays, n_chunks) in enumerate(build_grid(args)):
            sizes = chunk_sizes_for(rays, n_chunks)
            for k in range(args.per_cell):
                t_trace = [0.7 * b + 0.75 * c * n + rng.normal(0, 0.002) for n in sizes]
                t_fetch = [0.3 * b + 0.25 * c * n + rng.normal(0, 0.002) for n in sizes]
                t_setup = a + rng.normal(0, 0.004)
                rows.append(dict(
                    rep=rep, cell=cell, warmup=False, rays=rays, n_chunks=n_chunks,
                    chunk_rays=list(sizes), heliostat_id=k, rays_landed=int(0.375 * rays),
                    t_setup_s=t_setup, t_trace_s=sum(t_trace), t_fetch_s=sum(t_fetch),
                    t_total_s=t_setup + sum(t_trace) + sum(t_fetch),
                    t_trace_chunks=t_trace, t_fetch_chunks=t_fetch))
    return rows


# --------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Measure the fixed/per-call/per-ray split of a heliostat trace.")
    ap.add_argument("--config", default=None, help="path to config.toml")
    ap.add_argument("--rays", nargs="+", type=int,
                    default=[6000, 12000, 30000, 60000, 120000],
                    help="ray budgets per heliostat to sweep")
    ap.add_argument("--chunks", nargs="+", type=int, default=[1, 2, 4],
                    help="traceRays calls per heliostat to sweep")
    ap.add_argument("--min-chunk", type=int, default=1000,
                    help="drop cells whose chunks would be smaller than this")
    ap.add_argument("--per-cell", type=int, default=5,
                    help="heliostats traced in every cell (the same ones each time)")
    ap.add_argument("--repeats", type=int, default=1,
                    help="passes over the whole grid; medians are taken over all")
    ap.add_argument("--date", default="2026-03-20", help="ISO date for the one timestep")
    ap.add_argument("--step-index", type=int, default=None,
                    help="which step of that day's grid (default: the middle one)")
    ap.add_argument("--occluders", action="store_true",
                    help="trace the occluder model and write the 56 slot parameters, "
                         "as a full8-style sweep does -- this is part of the fixed "
                         "per-heliostat cost A")
    ap.add_argument("--policy-chunk", type=int, default=60000,
                    help="chunk size the decision table assumes a real sweep would "
                         "use, i.e. config.toml's rays_per_trace; 0 means one call")
    ap.add_argument("--n-heliostats", type=int, default=N_FIELD)
    ap.add_argument("--n-timesteps", type=int, default=N_TIMESTEPS)
    ap.add_argument("--csv", default="analysis_output/probe_ray_cost.csv",
                    help="where the raw per-trace timings go")
    ap.add_argument("--dry-run", action="store_true",
                    help="no Quadoa, no licence: print the grid and exercise the "
                         "fit on synthetic timings with a known A, B, C")
    ap.add_argument("--force", action="store_true",
                    help="run even though a sweep lock exists. Do not use this.")
    args = ap.parse_args(argv)

    repo_root = Path(__file__).resolve().parent.parent
    grid = build_grid(args)
    if not grid:
        print("no cells survive the grid filters; loosen --min-chunk")
        return 1

    print("probe_ray_cost: where does a heliostat trace's time actually go?")
    print("  cells (rays/calls): "
          + ", ".join(f"{r:,}/{n}" for r, n in grid))

    if args.dry_run:
        print()
        print("  *** DRY RUN: no Quadoa session, no licence seat, NO MEASUREMENT. ***")
        print("  The numbers below come from synthetic timings with a known")
        print("  A = 120.00 ms, B = 45.00 ms, C = 4.000 us, and exist only to show")
        print("  that the fit recovers them and that the report renders.")
        rows = synthetic_rows(args)
        model = fit_cost_model(rows)
        report(rows, model, args)
        print()
        print(f"  fit recovered A = {1e3 * model.A:.2f} ms (true 120.00), "
              f"B = {1e3 * model.B:.2f} ms (true 45.00), "
              f"C = {1e6 * model.C:.3f} us (true 4.000)")
        return 0

    locks = licence_locks(repo_root)
    if locks and not args.force:
        print()
        for lock in locks:
            pid = (lock / "pid")
            who = pid.read_text(encoding="utf-8").strip() if pid.exists() else "(no pid)"
            print(f"  REFUSING: {lock} exists (pid {who})")
        print("  A sweep holds the one HASP licence seat. Asking for a second has")
        print("  been measured to leak the first, leaving that sweep stuck with a")
        print("  modal H0038 dialog. Wait for it to finish. (--force overrides.)")
        return 2

    estimate_s = sum(
        args.per_cell * args.repeats * NOMINAL_S_PER_TRACE_120K * rays / 120000.0
        for rays, _ in grid)
    print(f"  rough estimate {estimate_s / 60:.0f} min, assuming time is "
          f"proportional to rays --")
    print("  which is precisely the hypothesis under test, so treat it as a floor.")
    print()

    rows = measure(args)
    model = fit_cost_model(rows)
    report(rows, model, args)

    csv_path = Path(args.csv)
    if not csv_path.is_absolute():
        csv_path = repo_root / csv_path
    write_csv(rows, csv_path)
    print(f"\n  wrote {len(rows)} raw timing rows to {csv_path}")
    print("  (per-call detail is kept, so the regression can be redone without "
          "re-tracing)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
