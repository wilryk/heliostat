"""Command line interface.

    python -m beamdown info                     # config, time grid, cost estimate
    python -m beamdown sweep [options]          # run the ray-trace sweep
    python -m beamdown figures [options]        # build figures from a finished sweep
    python -m beamdown rank [options]           # rank heliostats
    python -m beamdown compare RUN_A RUN_B      # is one sweep better than the other?
    python -m beamdown gui [options]            # interactive desktop explorer
    python -m beamdown inspect [options]        # export an .optx to open in Quadoa
    python -m beamdown fetch-dni --source pvgis # download a DNI series
"""

from __future__ import annotations

import argparse
import datetime as _dt
import sys
import warnings
from pathlib import Path

# Flags that change which secondary layout runs, and therefore which
# ``[geometry]`` keys are required and which ``n_mirrors`` is consistent. Their
# presence is what makes ``_load_with_overrides`` re-run the layout validation.
_LAYOUT_FLAGS = ("secondary", "n_mirrors", "focus_height_mm", "rim_height_mm")


def _override_map(cfg, args) -> dict:
    """The config values this command line overrides, as ``{section: {field: v}}``.

    Returned as data rather than applied on the spot because the same set has to
    be replayed inside every sweep worker: a worker calls ``load_config`` itself,
    from disk (:func:`beamdown.sweep._init_worker`), so an override applied only
    to the driver's copy never reaches the trace.

    ``getattr`` with a default throughout because these flags exist on some
    subcommands and not others, and ``_load`` is shared by all of them.
    """
    over: dict[str, dict] = {}

    def put(section: str, key: str, value) -> None:
        if value is not None:
            over.setdefault(section, {})[key] = value

    rays = getattr(args, "rays", None)
    per_trace = getattr(args, "rays_per_trace", None)
    put("trace", "rays_per_heliostat", rays)
    if per_trace is not None:
        # Both given: honour both literally. The chunk is the thing the user is
        # actually choosing -- it decides how many setRayDistributionCount1 +
        # traceRays + getRayPos round trips each heliostat costs -- so it is not
        # quietly clamped here. validate_trace refuses a chunk bigger than the
        # budget; see its docstring for the rule.
        put("trace", "rays_per_trace", per_trace)
    elif rays:
        # --rays alone: the historical clamp. rays_per_trace is the per-traceRays
        # chunk size and config.toml's value is a ceiling, so a chunk larger than
        # the whole budget would claim a split that never happens. TraceSpec
        # clamps this in __post_init__, which an override placed afterwards
        # bypasses.
        put("trace", "rays_per_trace", min(rays, cfg.trace.rays_per_trace))
    put("trace", "model_file", getattr(args, "model_file", None))

    put("optics", "secondary", getattr(args, "secondary", None))
    put("optics", "n_mirrors", getattr(args, "n_mirrors", None))
    # Tri-state on purpose: --flat-mirrors is True, --focused-mirrors is False,
    # and neither given is None, which ``put`` skips so config.toml stands. A
    # plain store_true would default to False and therefore override the file
    # with "focused" on every command line that never mentioned mirrors.
    #
    # It has to travel in this dict, not just on the driver's cfg: workers call
    # load_config from disk in _init_worker and replay these, so a flat run whose
    # flag stopped at the driver would trace 645 FOCUSED heliostats while the
    # console and the manifest both said flat.
    put("optics", "flat_mirrors", getattr(args, "flat_mirrors", None))
    put("geometry", "focus_height_mm", getattr(args, "focus_height_mm", None))
    put("geometry", "secondary_rim_height_mm", getattr(args, "rim_height_mm", None))

    put("sweep", "hour_step", getattr(args, "hour_step", None))
    put("sweep", "sunrise_margin_min", getattr(args, "sunrise_margin_min", None))
    return over


def _load_with_overrides(args):
    """Load config.toml, apply every override flag, and return ``(cfg, overrides)``.

    The overrides come back alongside the config because the sweep needs to hand
    them to its worker processes; see :func:`_override_map`.
    """
    from . import config as C

    touches_layout = any(getattr(args, name, None) is not None for name in _LAYOUT_FLAGS)

    if touches_layout:
        # load_config validates the layout described by the *file*, which is not
        # the layout about to run. A self-consistent
        # "--secondary prime_focus --n-mirrors 1" would still be judged against
        # whatever config.toml says (warning about n_mirrors, or refusing over a
        # [geometry] key this command line is itself supplying), so the file's
        # verdict is suppressed here and the check re-run below with the
        # overrides in place. _validate_layout is the only thing load_config
        # warns from, so nothing else is being swallowed.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cfg = C.load_config(args.config)
    else:
        cfg = C.load_config(args.config)

    if getattr(args, "output", None):
        object.__setattr__(cfg.storage, "root", args.output)

    overrides = _override_map(cfg, args)
    C.apply_overrides(cfg, overrides)

    try:
        # Unconditional: unlike the layout check this is cheap, and the only way
        # to reach a bad combination is a command line, which is exactly what is
        # being read here. Fail now, before a licence seat is taken.
        C.validate_trace(cfg)
    except ValueError as exc:
        raise SystemExit(f"beamdown: {exc}") from None

    if touches_layout:
        try:
            C.validate_layout(cfg)
        except ValueError as exc:
            # Fail now, before a licence seat is taken, with the config module's
            # own message -- not hours into a sweep, and above all not silently
            # with 645 heliostats aimed at a point that was never configured.
            raise SystemExit(f"beamdown: {exc}") from None
    return cfg, overrides


def _load(args):
    cfg, _overrides = _load_with_overrides(args)
    return cfg


def _resolve_dates(cfg, args):
    if getattr(args, "dates", None):
        return tuple(_dt.date.fromisoformat(d) for d in args.dates)
    if getattr(args, "suggest_dates", None):
        from . import energy

        return tuple(energy.suggest_sweep_dates(
            cfg, args.suggest_dates, branch="ascending", must_include=cfg.sweep.dates
        ))
    return cfg.sweep.dates


def _resolve_heliostats(cfg, args):
    from . import field as F

    fld = F.load_field(cfg)
    if getattr(args, "all_heliostats", False):
        return None, fld, "all 645"
    idx, prov = F.load_or_build_downselect(cfg, fld)
    return idx, fld, prov


# --------------------------------------------------------------------------

def cmd_info(args) -> int:
    from . import solar, energy, field as F

    cfg = _load(args)
    print(solar.describe_time_grid(cfg, _resolve_dates(cfg, args)))

    fld = F.load_field(cfg)
    steps = solar.build_time_grid(cfg, _resolve_dates(cfg, args))
    idx, _, prov = _resolve_heliostats(cfg, args)
    n_helio = len(fld) if idx is None else len(idx)

    print(f"\nheliostats: {n_helio} ({prov})")
    print(f"rays/heliostat: {cfg.trace.rays_per_heliostat:,}   workers: {cfg.trace.n_workers}")
    # The number of traceRays calls, not just the ray budget: one heliostat is
    # ceil(rays_per_heliostat / rays_per_trace) round trips, and how the fixed
    # per-call cost trades against the per-ray cost is unmeasured -- see
    # scripts/probe_ray_cost.py.
    sizes = cfg.trace.chunk_sizes
    print(f"rays/traceRays call: {cfg.trace.rays_per_trace:,}   -> "
          f"{cfg.trace.n_chunks} call(s) per heliostat "
          f"({' + '.join(f'{n:,}' for n in sizes)} = {sum(sizes):,})")
    print(f"traces: {n_helio * len(steps):,}")
    print(f"estimated time @0.55 s/trace: {n_helio * len(steps) * 0.55 / 3600:.1f} h")
    raw_gb = n_helio * len(steps) * 45000 * 2 * 2 / 1e9
    flux_gb = n_helio * len(steps) * cfg.receiver.grid_size ** 2 * 4 / 1e9
    print(f"estimated storage: {raw_gb:.1f} GB raw + {flux_gb:.1f} GB flux maps")

    print("\ndeclination coverage:")
    print(energy.declination_coverage(cfg).to_string(index=False))
    return 0


def cmd_sweep(args) -> int:
    from . import sweep as S

    if args.occluders and args.model_file:
        # run_sweep picks models/heliostat_field_occluders.optx itself for
        # --occluders and would silently overwrite anything given here, so say so
        # rather than tracing a model the command line did not ask for.
        raise SystemExit(
            "beamdown sweep: --occluders selects models/heliostat_field_occluders.optx "
            "itself, so it cannot be combined with --model-file. Drop one of them."
        )

    cfg, overrides = _load_with_overrides(args)
    idx, _, prov = _resolve_heliostats(cfg, args)
    dates = _resolve_dates(cfg, args)
    print(f"heliostats: {prov}")

    S.run_sweep(
        cfg,
        dates=dates,
        heliostat_indices=idx,
        workers=args.workers or cfg.trace.n_workers,
        resume=not args.no_resume,
        config_path=args.config,
        occluders=getattr(args, "occluders", False),
        overrides=overrides,
    )
    return 0


def cmd_figures(args) -> int:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from . import plots, dni as D
    from .store import RunStore

    cfg = _load(args)
    store = RunStore(cfg.output_root, cfg=cfg, mode="r")
    summary = store.summary()
    keys = store.timestep_keys()
    outdir = Path(args.figdir) if args.figdir else cfg.output_root / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"{len(keys)} timesteps, {len(summary)} summary rows -> {outdir}")

    # Per-timestep shading x blocking, applied as weights on the stored counts.
    eff = {}
    if not args.no_shading:
        for key in keys:
            sub = summary[summary.timestep == key].sort_values("heliostat_id")
            eff[key] = (sub.eta_shade * sub.eta_block).to_numpy()

    provider = D.load_dni_provider(cfg)

    made = []
    fig = plots.through_day_panels(
        store, cfg, keys=keys, efficiency_by_key=eff or None,
        dni_provider=provider, crop_mm=args.crop,
        save_path=outdir / "through_day_panels.png",
    )
    plt.close(fig); made.append("through_day_panels.png")

    fig = plots.power_through_day(summary, save_path=outdir / "power_through_day.png")
    plt.close(fig); made.append("power_through_day.png")

    fig = plots.efficiency_breakdown(summary, save_path=outdir / "efficiency_breakdown.png")
    plt.close(fig); made.append("efficiency_breakdown.png")

    mid = keys[len(keys) // 2]
    ids = sorted(summary.heliostat_id.unique())
    fig = plots.single_vs_field(store, cfg, mid, heliostat_row=0, heliostat_id=ids[0],
                                crop_mm=args.crop,
                                save_path=outdir / "single_vs_field.png")
    plt.close(fig); made.append("single_vs_field.png")

    for col in ("power_w", "eta_shade", "eta_block", "transmission"):
        if col in summary.columns:
            fig = plots.field_scatter(summary, col, save_path=outdir / f"field_{col}.png")
            plt.close(fig); made.append(f"field_{col}.png")

    for name in made:
        print(f"  wrote {name}")
    return 0


def cmd_rank(args) -> int:
    from . import metrics
    from .store import RunStore

    cfg = _load(args)
    store = RunStore(cfg.output_root, cfg=cfg, mode="r")
    ranked = metrics.rank_heliostats(store.summary(), by=args.by,
                                     ascending=not args.best_first)
    print(ranked.head(args.n).to_string(index=False))
    return 0


def cmd_occluders(args) -> int:
    """Emit the occluding rectangles for one heliostat, to rebuild in Quadoa."""
    from . import occluders as O
    from .store import RunStore

    cfg = _load(args)
    if args.output:
        object.__setattr__(cfg.storage, "root", args.output)
    summary = RunStore(cfg.output_root, cfg=cfg, mode="r").summary()
    found, totals = O.occluders_for(cfg, summary, args.heliostat, args.timestep)
    text = O.describe(found, totals, cfg, args.heliostat, args.timestep)
    print(text)
    if args.write:
        Path(args.write).write_text(text, encoding="utf-8")
        print(f"\n  written to {args.write}")
    return 0


def cmd_rescale(args) -> int:
    """Recompute shading/blocking on a finished run. No trace, no license."""
    from . import rescale as R

    cfg = _load(args)
    report = R.recompute_weights(args.run, cfg, apply=args.apply)
    print(report.describe())
    if not args.apply:
        print("\n  re-run with --apply to write it (the original is copied aside first)")
    return 0


def cmd_compare(args) -> int:
    from . import compare as C
    from . import field as F
    from .store import RunStore

    cfg = _load(args)
    radii = ([float(r) for r in args.radii] if args.radii else C.DEFAULT_RADII_MM)

    stores = []
    for root in (args.run_a, args.run_b):
        object.__setattr__(cfg.storage, "root", root)
        stores.append(RunStore(cfg.output_root, cfg=cfg, mode="r"))

    labels = args.labels or [Path(args.run_a).name, Path(args.run_b).name]
    report = C.compare_runs(*stores, cfg, label_a=labels[0], label_b=labels[1],
                            radii_mm=radii, reference_mm=args.aperture,
                            dni_w_m2=args.dni, use_shading=not args.no_shading)
    print(report.describe())

    if args.attribute:
        print("\nwhere the change landed, banded by how far the correction moved "
              "each heliostat")
        table = C.attribute_by_foreshortening(report, cfg, F.load_field(cfg))
        print(table.to_string(index=False))
        print(f"\n  monotonic in delta_w: {table.attrs['monotonic']}")
        print(f"  correlation with delta_w:  predicted |d astig| r = "
              f"{table.attrs['r_d_astig']:+.3f}   1/L^2 alone r = "
              f"{table.attrs['r_inv_L2']:+.3f}")
        print("  The lowest band is the control: the correction barely moves those")
        print("  heliostats, so their delta_w should sit at zero. A gain there means")
        print("  something other than the shape correction moved.")

    if args.csv:
        report.per_heliostat.to_csv(args.csv, index=False)
        print(f"\nwrote per-heliostat detail to {args.csv}")
    return 0


def cmd_inspect(args) -> int:
    from . import inspect_model as IM
    from .store import RunStore

    cfg = _load(args)
    store = RunStore(cfg.output_root, cfg=cfg, mode="r")
    summary = store.summary()

    timestep = args.timestep
    if timestep is None:
        keys = store.timestep_keys()
        timestep = keys[len(keys) // 2]
        print(f"no --timestep given, using {timestep}")
    elif timestep not in set(summary.timestep):
        raise SystemExit(
            f"timestep {timestep!r} not in the store. Available:\n  "
            + "\n  ".join(store.timestep_keys())
        )

    if args.heliostat:
        ids = [int(h) for h in args.heliostat]
    else:
        ids = IM.pick_heliostats(summary, timestep=None if args.over_run else timestep,
                                 by=args.by, n=args.n, worst=not args.best)
        which = "best" if args.best else "worst"
        scope = "averaged over the run" if args.over_run else f"at {timestep}"
        print(f"selected the {which} {len(ids)} by {args.by} ({scope}): {ids}")

    report = IM.export_for_inspection(cfg, ids, timestep=timestep, summary=summary,
                                      out_path=args.out)
    print()
    print(report.describe())
    return 0


def cmd_gui(args) -> int:
    from .gui import main as gui_main

    argv = []
    if args.config:
        argv += ["--config", args.config]
    if args.output:
        argv += ["--output", args.output]
    return gui_main(argv)


def cmd_fetch_dni(args) -> int:
    from . import dni as D

    cfg = _load(args)
    path = D.fetch(args.source, cfg, year=args.year)
    print(f"wrote {path}")
    print('Set [dni] mode = "table" in config.toml to use it.')
    return 0


def build_parser() -> argparse.ArgumentParser:
    """The whole CLI, built without running anything.

    Split out of :func:`main` so a test can parse an argument vector -- and the
    GUI's Trace tab can be checked against the real flag set -- without
    dispatching to a command that would take a licence seat.
    """
    # Imported here rather than at module scope to keep the lazy-import style of
    # the rest of this file; config only pulls in tomli.
    from .config import SECONDARY_LAYOUTS

    ap = argparse.ArgumentParser(prog="beamdown")
    ap.add_argument("--config", default=None, help="path to config.toml")
    sub = ap.add_subparsers(dest="command", required=True)

    def common(p, traces=False):
        p.add_argument("--output", default=None, help="override storage.root")
        p.add_argument("--dates", nargs="+", default=None, help="ISO dates")
        p.add_argument("--suggest-dates", type=int, default=None,
                       metavar="N", help="use N declination-spaced dates")
        p.add_argument("--all-heliostats", action="store_true",
                       help="all 645 instead of the downselect")
        if traces:
            p.add_argument("--rays", type=int, default=None,
                           metavar="N",
                           help="override [trace] rays_per_heliostat: the total "
                                "ray budget for one heliostat at one instant")
            p.add_argument("--rays-per-trace", type=int, default=None,
                           metavar="N",
                           help="override [trace] rays_per_trace: rays per "
                                "traceRays call, so each heliostat costs "
                                "ceil(rays / rays-per-trace) round trips of "
                                "setRayDistributionCount1 + traceRays + "
                                "getRayPos. With --rays and WITHOUT this, the "
                                "chunk stays config.toml's value clamped to the "
                                "budget; with BOTH, both are honoured literally "
                                "and a chunk larger than the budget is an error. "
                                "--rays-per-trace equal to --rays is one call "
                                "per heliostat")
            p.add_argument("--workers", type=int, default=None)
            p.add_argument("--no-resume", action="store_true")

    p = sub.add_parser("info", help="show config, time grid and cost estimate")
    common(p); p.set_defaults(func=cmd_info)

    p = sub.add_parser("sweep", help="run the ray-trace sweep")
    common(p, traces=True)
    p.add_argument("--occluders", action="store_true",
                   help="trace neighbour shading and blocking as real geometry "
                        "instead of applying them as scalars afterwards")
    # -- config overrides -------------------------------------------------
    #
    # Every one of these exists so a run can be specified entirely on the command
    # line. Editing config.toml to select a layout is shared mutable state: a
    # sweep already running re-reads the file at the end of the run for its
    # report, so a mid-run edit corrupts that, and two people cannot set up two
    # different runs at once. They all end up in _override_map, are applied to
    # this process's config AND replayed inside each worker, and are recorded in
    # the run's manifest.
    p.add_argument("--secondary", default=None, choices=list(SECONDARY_LAYOUTS),
                   help="override [optics] secondary. prime_focus and cassegrain "
                        "additionally need --focus-height-mm (cassegrain also "
                        "--rim-height-mm) unless config.toml already sets them")
    p.add_argument("--focus-height-mm", type=float, default=None,
                   help="override [geometry] focus_height_mm -- the single on-axis "
                        "point F1 the whole field aims at (prime_focus, cassegrain)")
    p.add_argument("--rim-height-mm", type=float, default=None,
                   help="override [geometry] secondary_rim_height_mm -- height of "
                        "the Cassegrain hyperboloid rim whose shadow falls on the field")
    p.add_argument("--n-mirrors", type=int, default=None, choices=[1, 2],
                   help="override [optics] n_mirrors: reflections in the path, so "
                        "throughput = reflectivity**n. 1 for prime_focus, 2 for "
                        "axicon/cassegrain. Applied BEFORE the layout check, so a "
                        "self-consistent command line does not warn")
    # Flat vs focused heliostats: a second comparison axis, orthogonal to
    # --secondary, so it composes with all three layouts. Both directions exist
    # because the flag has to be able to DISAGREE with config.toml either way --
    # the GUI's Trace tab only emits a flag where it differs from the file.
    p.add_argument("--flat-mirrors", dest="flat_mirrors", action="store_true",
                   default=None,
                   help="override [optics] flat_mirrors: FLAT heliostats -- keep "
                        "the pointing, force the heliostat Zernike z3/z4/z5 to "
                        "zero so each mirror is its bare radius = inf plane. "
                        "Composes with every --secondary. Expect a much larger "
                        "spot, much more spillage and much less collected "
                        "energy than the focused run: that is the comparison")
    p.add_argument("--focused-mirrors", dest="flat_mirrors", action="store_false",
                   default=None,
                   help="the opposite of --flat-mirrors, for when config.toml "
                        "sets flat_mirrors = true and this run should not be flat")
    p.add_argument("--model-file", default=None,
                   help="override [trace] model_file; different layouts need "
                        "different .optx files. Not compatible with --occluders, "
                        "which selects its own model")
    p.add_argument("--hour-step", type=float, default=None,
                   help="override [sweep] hour_step: the MAXIMUM spacing between "
                        "samples, not a clock grid (see solar.build_time_grid)")
    p.add_argument("--sunrise-margin-min", type=float, default=None,
                   help="override [sweep] sunrise_margin_min: minutes of margin "
                        "inside sunrise and sunset")
    p.set_defaults(func=cmd_sweep)

    p = sub.add_parser("figures", help="build figures from a finished sweep")
    common(p)
    p.add_argument("--figdir", default=None)
    p.add_argument("--crop", type=float, default=None, help="crop panels to +/- mm")
    p.add_argument("--no-shading", action="store_true")
    p.set_defaults(func=cmd_figures)

    p = sub.add_parser("rank", help="rank heliostats")
    common(p)
    p.add_argument("--by", default="power_w")
    p.add_argument("-n", type=int, default=15)
    p.add_argument("--best-first", action="store_true")
    p.set_defaults(func=cmd_rank)

    p = sub.add_parser("occluders",
                       help="list the shading/blocking rectangles for one heliostat")
    common(p)
    p.add_argument("--heliostat", type=int, required=True)
    p.add_argument("--timestep", required=True, help="e.g. 20260320_1600")
    p.add_argument("--write", default=None, help="also save the listing to a file")
    p.set_defaults(func=cmd_occluders)

    p = sub.add_parser("rescale",
                       help="recompute shading/blocking on a finished run")
    common(p)
    p.add_argument("run", help="sweep output directory")
    p.add_argument("--apply", action="store_true",
                   help="write the result; without this it only reports")
    p.set_defaults(func=cmd_rescale)

    p = sub.add_parser("compare", help="compare two sweeps of the same field")
    common(p)
    p.add_argument("run_a", help="baseline run directory")
    p.add_argument("run_b", help="run directory to test against it")
    p.add_argument("--labels", nargs=2, default=None, help="names for the two runs")
    p.add_argument("--radii", nargs="+", type=float, default=None,
                   help="aperture radii in mm to sweep")
    p.add_argument("--aperture", type=float, default=700.0,
                   help="reference aperture radius for the per-timestep table")
    p.add_argument("--dni", type=float, default=1000.0)
    p.add_argument("--no-shading", action="store_true")
    p.add_argument("--attribute", action="store_true",
                   help="band the change by 1/L^2 to check where it came from")
    p.add_argument("--csv", default=None, help="write per-heliostat detail here")
    p.set_defaults(func=cmd_compare)

    p = sub.add_parser("inspect", help="export an .optx set up for specific heliostats")
    common(p)
    p.add_argument("--heliostat", nargs="+", type=int, default=None,
                   help="heliostat id(s); omit to select by metric")
    p.add_argument("--timestep", default=None, help="e.g. 20260922_0800")
    p.add_argument("--by", default="power_w", help="metric used when selecting")
    p.add_argument("-n", type=int, default=1, help="how many to export")
    p.add_argument("--best", action="store_true", help="select best rather than worst")
    p.add_argument("--over-run", action="store_true",
                   help="rank over the whole run instead of at this timestep")
    p.add_argument("--out", default=None, help="output .optx path")
    p.set_defaults(func=cmd_inspect)

    p = sub.add_parser("gui", help="interactive desktop explorer")
    common(p); p.set_defaults(func=cmd_gui)

    p = sub.add_parser("fetch-dni", help="download an hourly DNI series")
    common(p)
    p.add_argument("--source", default="pvgis", choices=["pvgis", "nasa"])
    p.add_argument("--year", type=int, default=None)
    p.set_defaults(func=cmd_fetch_dni)

    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
