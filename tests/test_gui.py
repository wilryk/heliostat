"""GUI tests: every view renders, the table does not recurse, and the Trace tab
builds the command it claims to.

    python tests/test_gui.py [run_directory]

Runs against a real sweep (default ``analysis_output/demo25``) but needs no
Quadoa license -- the GUI only reads stored counts.

**Nothing here launches a sweep.** The Trace tab's launch path is exercised with
``subprocess.Popen`` monkeypatched, its output and lock directories pointed at a
temporary directory, and its lock scan pointed away from ``analysis_output`` --
because a real launch would ask for the single HASP licence seat, and a failed
request for it has been measured to leak the seat a running sweep holds.

The recursion test drives everything through ``root.update()`` rather than
calling the draw functions directly. That distinction is the whole point: the
table highlights the selected row with ``selection_set``, Tk delivers
``<<TreeviewSelect>>`` only while the event loop runs, and a handler that
rebuilt the table in response would loop until the window locked up. A test that
calls ``_draw_table()`` directly passes even when the bug is present.
"""

from __future__ import annotations

import sys
import time
import tkinter as tk
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

RUN = sys.argv[1] if len(sys.argv) > 1 else "analysis_output/demo25"


def build():
    from beamdown.config import load_config
    from beamdown.gui import BeamdownGUI
    from beamdown.store import RunStore

    cfg = load_config(str(REPO / "config.toml"))
    object.__setattr__(cfg.storage, "root", RUN)
    store = RunStore(cfg.output_root, cfg=cfg, mode="r")
    root = tk.Tk()
    root.geometry("1400x850")
    return root, BeamdownGUI(root, cfg, store), cfg


def pump(root, n=25):
    for _ in range(n):
        try:
            root.update()
        except tk.TclError:
            return
        time.sleep(0.005)


def check_command_builder() -> bool:
    """The Trace tab's pure logic: options dict -> argv, and the refusals.

    No Tk and no filesystem beyond a tmp directory. This is the part that has to
    be right even if the GUI never renders: the command preview the user reads
    and the argv the Run button executes come from this one function, so a bug
    here is a run that does something other than what it says.
    """
    import tempfile

    from beamdown import gui as G

    ok = True

    def check(label, condition):
        nonlocal ok
        ok &= bool(condition)
        print(f"    {'OK  ' if condition else 'FAIL'} {label}")

    # -- every option maps to the flag the CLI actually defines ----------
    argv = G.build_sweep_argv({
        "config": "other.toml",
        "output": "analysis_output/t1",
        "dates": ["2026-03-20", "2026-06-21"],
        "all_heliostats": True,
        "rays": 60000,
        "rays_per_trace": 20000,
        "workers": 1,
        "hour_step": 0.5,
        "sunrise_margin_min": 15.0,
        "secondary": "cassegrain",
        "focus_height_mm": 24000.0,
        "rim_height_mm": 20000.0,
        "n_mirrors": 2,
        "model_file": "models/cassegrain.optx",
        "occluders": False,
        "resume": True,
    }, python="py")
    print(f"  {G.format_command(argv)}")
    check("--config precedes the subcommand",
          argv.index("--config") < argv.index("sweep"))
    for flag, value in (("--output", "analysis_output/t1"), ("--rays", "60000"),
                        ("--rays-per-trace", "20000"),
                        ("--workers", "1"), ("--hour-step", "0.5"),
                        ("--sunrise-margin-min", "15"),
                        ("--secondary", "cassegrain"),
                        ("--focus-height-mm", "24000"),
                        ("--rim-height-mm", "20000"), ("--n-mirrors", "2"),
                        ("--model-file", "models/cassegrain.optx")):
        check(f"{flag} {value}", flag in argv and argv[argv.index(flag) + 1] == value)
    check("--dates takes both dates",
          argv[argv.index("--dates") + 1:argv.index("--dates") + 3]
          == ["2026-03-20", "2026-06-21"])
    check("--all-heliostats present", "--all-heliostats" in argv)
    check("resume on omits --no-resume", "--no-resume" not in argv)
    # A real parse is the only proof that these are flags beamdown accepts.
    from beamdown.cli import build_parser

    parsed = build_parser().parse_args(argv[argv.index("-m") + 2:])
    check("argparse accepts the whole line", parsed.command == "sweep")
    check("parsed values round-trip",
          parsed.rays == 60000 and parsed.secondary == "cassegrain"
          and parsed.focus_height_mm == 24000.0 and parsed.n_mirrors == 2
          and parsed.hour_step == 0.5 and parsed.dates == ["2026-03-20", "2026-06-21"])
    check("--rays-per-trace round-trips as an int", parsed.rays_per_trace == 20000)

    # -- the derived call count, which is what --rays-per-trace really buys ----
    #
    # The label the Trace tab shows is this function, and it splits with
    # config.chunk_plan -- the same call the trace makes -- so it cannot describe
    # a split that will not happen. The sum is the load-bearing part: a dropped
    # remainder would emit fewer rays than the run claims.
    from beamdown.config import chunk_plan

    check("one call when the chunk is the whole budget",
          chunk_plan(60000, 60000) == [60000]
          and "1 traceRays call per heliostat" in G.describe_call_plan(60000, 60000))
    check("an even split is N equal chunks",
          chunk_plan(120000, 60000) == [60000, 60000]
          and "2 traceRays calls per heliostat" in G.describe_call_plan(120000, 60000))
    check("an uneven budget gets a short last chunk, still summing exactly",
          chunk_plan(100000, 30000) == [30000, 30000, 30000, 10000]
          and "4 traceRays calls" in G.describe_call_plan(100000, 30000)
          and "= 100,000 rays" in G.describe_call_plan(100000, 30000))
    check("a chunk larger than the budget collapses to one call",
          chunk_plan(6000, 30000) == [6000])
    for total, per in ((120000, 60000), (100000, 30000), (6000, 12000),
                       (12000, 12000), (7, 3), (60000, 12000)):
        check(f"chunk_plan({total}, {per}) sums to {total}",
              sum(chunk_plan(total, per)) == total)

    # -- defaults, and the flags that are always written out -------------
    bare = G.build_sweep_argv({"output": "analysis_output/t2"}, python="py")
    check("workers defaults to 1 even when unset",
          bare[bare.index("--workers") + 1] == "1")
    check("no dates -> no --dates and no --suggest-dates",
          "--dates" not in bare and "--suggest-dates" not in bare)
    check("resume off (the default) emits --no-resume", "--no-resume" in bare)
    check("nothing else is invented",
          not any(a.startswith("--secondary") or a.startswith("--model-file")
                  for a in bare))
    # Unlike --rays, the chunk size is NOT one of the always-written flags: it is
    # a departure from config.toml like any other, so silence means "the file's".
    check("no --rays-per-trace unless it differs from config.toml",
          "--rays-per-trace" not in bare)

    suggested = G.build_sweep_argv({"output": "o", "dates": ["2026-01-01"],
                                    "suggest_dates": 8}, python="py")
    check("--dates wins over --suggest-dates, never both",
          "--dates" in suggested and "--suggest-dates" not in suggested)

    occ = G.build_sweep_argv({"output": "o", "occluders": True,
                              "model_file": "models/x.optx"}, python="py")
    check("--occluders suppresses --model-file (the CLI refuses both)",
          "--occluders" in occ and "--model-file" not in occ)

    # -- flat heliostats: a tri-state, because "not mentioned" is a real state --
    #
    # The tab only says what departs from config.toml, so the dict has to be able
    # to express agreement (key absent) separately from "explicitly focused"
    # (key False) -- otherwise a run against a flat config.toml could not be
    # asked for focused mirrors at all.
    check("flat_mirrors absent emits neither flag",
          not any(a.endswith("-mirrors") for a in bare))
    flat = G.build_sweep_argv({"output": "o", "flat_mirrors": True}, python="py")
    check("flat_mirrors True emits --flat-mirrors",
          "--flat-mirrors" in flat and "--focused-mirrors" not in flat)
    focused = G.build_sweep_argv({"output": "o", "flat_mirrors": False}, python="py")
    check("flat_mirrors False emits --focused-mirrors",
          "--focused-mirrors" in focused and "--flat-mirrors" not in focused)
    from beamdown.cli import build_parser as _bp

    check("argparse accepts both and they land on one tri-state dest",
          _bp().parse_args(["sweep", "--flat-mirrors"]).flat_mirrors is True
          and _bp().parse_args(["sweep", "--focused-mirrors"]).flat_mirrors is False
          and _bp().parse_args(["sweep"]).flat_mirrors is None)

    # -- the refusals ---------------------------------------------------
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        # Deliberately NOT analysis_output: a lock there belongs to a real run.
        lock = root / ".full9.lock"
        lock.mkdir()
        (lock / "pid").write_text("31337", encoding="utf-8")
        locks = G.scan_locks([root])
        check("scan_locks finds the lock and its pid",
              len(locks) == 1 and locks[0]["pid"] == "31337")
        refusal = G.launch_refusal(root / "new", True, locks)
        check("a lock refuses the launch, resume or not",
              refusal is not None and ".full9.lock" in refusal and "31337" in refusal
              and "seat" in refusal)

        check("no lock and a free name is allowed",
              G.launch_refusal(root / "new", False, []) is None)
        (root / "existing").mkdir()
        collision = G.launch_refusal(root / "existing", False, [])
        check("an existing output directory refuses without resume",
              collision is not None and "resume" in collision)
        check("...and is allowed with resume",
              G.launch_refusal(root / "existing", True, []) is None)
        check("a missing parent directory refuses",
              G.launch_refusal(root / "nope" / "run", True, []) is not None)

        lock_without_pid = root / ".bare.lock"
        lock_without_pid.mkdir()
        check("a lock with no pid file still reports",
              any(d["pid"] == "(no pid file)" for d in G.scan_locks([root])))

    # -- log parsing ----------------------------------------------------
    log = ("Sweep: 645 heliostats x 44 timesteps\n"
           "[1/44] 2026-03-20 08:00  el= 20.9   505.9s  (784.3 ms/heliostat)  "
           "eta 362.6 min\n"
           "[2/44] 2026-03-20 09:00 -- already done, skipping\n")
    info = G.parse_progress(log)
    check("progress reads the latest [i/N]",
          info.get("done") == 2 and info.get("total") == 44)
    check("ETA survives a skipped timestep with no ETA of its own",
          abs(info.get("eta_min", 0) - 362.6) < 1e-9)
    check("an empty log yields nothing rather than a wrong zero",
          G.parse_progress("starting") == {})
    check("completion is noticed",
          G.parse_progress(log + "Sweep complete in 5.0 min -> x").get("complete"))
    return ok


def check_cli_overrides() -> bool:
    """Each layout flag reaches the config, and a bad combination fails early.

    Lives in the GUI test because it is the Trace tab's contract: the tab may
    only offer what can be said on a command line, since editing config.toml to
    select options is shared state that a running sweep re-reads.
    """
    from beamdown.cli import build_parser
    from beamdown import cli as CLI

    ok = True

    def check(label, condition):
        nonlocal ok
        ok &= bool(condition)
        print(f"    {'OK  ' if condition else 'FAIL'} {label}")

    ap = build_parser()

    try:
        ap.parse_args(["sweep", "--help"])
        code = None
    except SystemExit as exc:
        code = exc.code
    check("sweep --help exits 0", code == 0)

    args = ap.parse_args(["sweep", "--rays", "1000", "--secondary", "cassegrain",
                          "--focus-height-mm", "24000", "--rim-height-mm", "20000",
                          "--n-mirrors", "2", "--hour-step", "2",
                          "--sunrise-margin-min", "20",
                          "--model-file", "models/other.optx"])
    cfg, overrides = CLI._load_with_overrides(args)
    check("secondary override lands", cfg.optics.secondary == "cassegrain")
    check("focus/rim heights land",
          cfg.geometry.focus_height_mm == 24000.0
          and cfg.geometry.secondary_rim_height_mm == 20000.0)
    check("ray budget lands, chunk clamped to it",
          cfg.trace.rays_per_heliostat == 1000 and cfg.trace.rays_per_trace == 1000)
    check("time grid overrides land",
          cfg.sweep.hour_step == 2.0 and cfg.sweep.sunrise_margin_min == 20.0)
    check("model file lands", cfg.trace.model_file == "models/other.optx")
    # The same mapping is what run_sweep replays inside each worker, which loads
    # config.toml for itself -- an override missing from here would silently not
    # reach the trace.
    check("overrides are reported for the workers",
          overrides["optics"]["secondary"] == "cassegrain"
          and overrides["trace"]["rays_per_heliostat"] == 1000
          and overrides["trace"]["model_file"] == "models/other.optx"
          and overrides["geometry"]["focus_height_mm"] == 24000.0
          and overrides["sweep"]["hour_step"] == 2.0)

    from beamdown.config import apply_overrides, load_config

    fresh = load_config(str(REPO / "config.toml"))
    apply_overrides(fresh, overrides)
    check("apply_overrides reproduces the driver's config from the dict alone",
          fresh.optics.secondary == "cassegrain"
          and fresh.trace.rays_per_heliostat == 1000
          and fresh.geometry.focus_height_mm == 24000.0)

    # -- --flat-mirrors must reach the WORKERS, not just the driver -----------
    #
    # This is the bug class this repository has actually shipped: a worker calls
    # load_config from disk in sweep._init_worker and replays the override dict,
    # so a flag that lands only on the driver's cfg never reaches the trace. A
    # flat run would then print "flat", write "flat" in its manifest, and trace
    # 645 focused heliostats.
    from beamdown.secondary import FlatHeliostats, get_strategy

    flat_args = ap.parse_args(["sweep", "--flat-mirrors"])
    cfg_flat, over_flat = CLI._load_with_overrides(flat_args)
    check("--flat-mirrors lands on the driver's config",
          cfg_flat.optics.flat_mirrors is True)
    check("and is in the override map the workers replay",
          over_flat["optics"]["flat_mirrors"] is True)

    worker_flat = load_config(str(REPO / "config.toml"))
    check("a worker's own config.toml is focused before the replay",
          worker_flat.optics.flat_mirrors is False)
    apply_overrides(worker_flat, over_flat)
    check("a worker-style reload comes out flat",
          worker_flat.optics.flat_mirrors is True)
    # The end of the chain: the object the worker actually solves with.
    worker_strategy = get_strategy(worker_flat)
    sol = worker_strategy.solve(45000.0, 12000.0, 135.0, 35.0, worker_flat.geometry)
    check("the strategy a worker builds from that config is flat",
          isinstance(worker_strategy, FlatHeliostats)
          and (sol.c3, sol.c4, sol.c5) == (0.0, 0.0, 0.0))

    # Silence is not "focused": a command line that never mentions mirrors must
    # leave the key out entirely, or it would override config.toml on every run.
    quiet = ap.parse_args(["sweep", "--rays", "1000"])
    _cfg_quiet, over_quiet = CLI._load_with_overrides(quiet)
    check("no mirror flag puts nothing in the override map",
          "flat_mirrors" not in over_quiet.get("optics", {}))

    focused_args = ap.parse_args(["sweep", "--focused-mirrors"])
    cfg_focused, over_focused = CLI._load_with_overrides(focused_args)
    check("--focused-mirrors is carried explicitly, as False not as absence",
          cfg_focused.optics.flat_mirrors is False
          and over_focused["optics"]["flat_mirrors"] is False)

    # Orthogonality: flat composes with every layout, and does not disturb the
    # layout machinery it sits beside.
    combo = ap.parse_args(["sweep", "--secondary", "prime_focus", "--n-mirrors", "1",
                           "--focus-height-mm", "27000", "--flat-mirrors"])
    cfg_combo, over_combo = CLI._load_with_overrides(combo)
    check("--flat-mirrors composes with a full layout override",
          cfg_combo.optics.secondary == "prime_focus"
          and cfg_combo.optics.flat_mirrors is True
          and over_combo["optics"]["flat_mirrors"] is True
          and over_combo["optics"]["secondary"] == "prime_focus")

    # -- --fixed-shapes: the same trip, for the frozen mirror figure -----------
    #
    # Same bug class as --flat-mirrors above, with a worse failure mode: a table
    # path that stopped at the driver would leave every worker re-figuring its
    # mirror at every timestep while the console, the manifest and the run name
    # all said the figure was fixed -- and the answer would land between the two
    # cases, which is exactly where a real fixed-figure answer belongs.
    import tempfile

    from beamdown.secondary import FixedShapeError, FixedShapeHeliostats

    tmp_shapes = Path(tempfile.mkdtemp()) / "fixed_shapes_fixture.csv"
    # Two heliostats is enough: one to hit, one to prove the lookup is by
    # position and not by row order.
    tmp_shapes.write_text(
        "# fixture: not a real figure, just three distinguishable numbers\n"
        "heliostat,x_mm,y_mm,c3,c4,c5\n"
        "0,45000.0,12000.0,1.5e-05,-8.5e-06,3.1e-07\n"
        "1,-30000.0,70000.0,2.0e-05,-1.0e-05,4.0e-07\n"
    )

    fixed_args = ap.parse_args(["sweep", "--fixed-shapes", str(tmp_shapes)])
    cfg_fixed, over_fixed = CLI._load_with_overrides(fixed_args)
    check("--fixed-shapes lands on the driver's config",
          cfg_fixed.optics.fixed_shapes == str(tmp_shapes))
    check("and is in the override map the workers replay",
          over_fixed["optics"]["fixed_shapes"] == str(tmp_shapes))
    check("no figure flag puts nothing in the override map",
          "fixed_shapes" not in over_quiet.get("optics", {}))

    worker_fixed = load_config(str(REPO / "config.toml"))
    check("a worker's own config.toml has no fixed figure before the replay",
          worker_fixed.optics.fixed_shapes == "")
    apply_overrides(worker_fixed, over_fixed)
    fixed_strategy = get_strategy(worker_fixed)
    hit = fixed_strategy.solve(45000.0, 12000.0, 135.0, 35.0, worker_fixed.geometry)
    focused_hit = get_strategy(worker_fixed.optics.secondary).solve(
        45000.0, 12000.0, 135.0, 35.0, worker_fixed.geometry)
    check("the strategy a worker builds from that config carries the table",
          isinstance(fixed_strategy, FixedShapeHeliostats)
          and (hit.c3, hit.c4, hit.c5) == (1.5e-05, -8.5e-06, 3.1e-07))
    # The whole point of the option: pointing is the focused answer, bit for bit,
    # and only the figure moved.
    check("pointing still tracks exactly as the focused solve does",
          (hit.rot_az_deg, hit.rot_el_deg) == (focused_hit.rot_az_deg,
                                               focused_hit.rot_el_deg)
          and (hit.c3, hit.c4, hit.c5) != (focused_hit.c3, focused_hit.c4,
                                           focused_hit.c5))

    missing = False
    try:
        fixed_strategy.solve(1.0, 2.0, 135.0, 35.0, worker_fixed.geometry)
    except FixedShapeError:
        missing = True
    check("a heliostat absent from the table is a hard error, not a fall-back",
          missing)

    # Both flags write c3/c4/c5, so argparse refuses the pair rather than
    # letting one of them silently do nothing.
    try:
        ap.parse_args(["sweep", "--flat-mirrors", "--fixed-shapes", str(tmp_shapes)])
        exclusive = False
    except SystemExit:
        exclusive = True
    check("--flat-mirrors and --fixed-shapes are mutually exclusive", exclusive)
    # ...and the config.toml-says-flat case, which argparse cannot see.
    flat_and_fixed = load_config(str(REPO / "config.toml"))
    apply_overrides(flat_and_fixed, {"optics": {"flat_mirrors": True,
                                                "fixed_shapes": str(tmp_shapes)}})
    try:
        get_strategy(flat_and_fixed)
        refused = False
    except ValueError:
        refused = True
    check("a flat config plus a fixed table is refused by get_strategy", refused)
    # --focused-mirrors is the documented way out of that, so it must compose.
    both_ok = ap.parse_args(["sweep", "--focused-mirrors",
                             "--fixed-shapes", str(tmp_shapes)])
    cfg_both, _ = CLI._load_with_overrides(both_ok)
    check("--focused-mirrors composes with --fixed-shapes",
          cfg_both.optics.flat_mirrors is False
          and isinstance(get_strategy(cfg_both), FixedShapeHeliostats))

    # The manifest key: written only when the run has a fixed figure, so an
    # ABSENT key keeps meaning "re-figured every timestep" for every run written
    # before the option existed.
    from beamdown.store import RunStore

    manifests = {}
    for label, source in (("fixed", cfg_both), ("plain", cfg_focused)):
        root = Path(tempfile.mkdtemp()) / label
        RunStore(root, cfg=source, mode="w").write_manifest()
        manifests[label] = RunStore(root, cfg=source, mode="r").manifest
    check("a fixed-figure run records fixed_shapes in its manifest",
          manifests["fixed"].get("fixed_shapes") == str(tmp_shapes))
    check("an ordinary run leaves the key out entirely",
          "fixed_shapes" not in manifests["plain"])

    # -- --rays-per-trace: the number of traceRays calls per heliostat ---------
    #
    # This is the flag that controls the ITERATION count, which --rays does not:
    # a heliostat is ceil(rays / rays_per_trace) round trips of
    # setRayDistributionCount1 + traceRays + getRayPos, and until this flag
    # existed that count could only be changed by editing config.toml.
    both = ap.parse_args(["sweep", "--rays", "60000", "--rays-per-trace", "12000"])
    cfg3, over3 = CLI._load_with_overrides(both)
    check("both flags are honoured literally, neither clamps the other",
          cfg3.trace.rays_per_heliostat == 60000 and cfg3.trace.rays_per_trace == 12000)
    check("the derived call count follows",
          cfg3.trace.n_chunks == 5 and sum(cfg3.trace.chunk_sizes) == 60000)
    check("the chunk override travels to the workers",
          over3["trace"]["rays_per_trace"] == 12000
          and over3["trace"]["rays_per_heliostat"] == 60000)

    # A worker builds its config from disk and replays the dict; without the
    # chunk in there it would trace config.toml's 60,000-ray chunks -- one call
    # per heliostat instead of five -- while the driver reported otherwise.
    worker_cfg = load_config(str(REPO / "config.toml"))
    apply_overrides(worker_cfg, over3)
    check("a worker-style reload reproduces the chunking exactly",
          worker_cfg.trace.rays_per_trace == 12000
          and worker_cfg.trace.chunk_sizes == [12000] * 5)

    # --rays alone keeps the historical clamp: config.toml's chunk capped at the
    # budget, so a 1,000-ray smoke run does not claim 60,000-ray chunks.
    check("--rays alone still clamps the chunk to the budget",
          cfg.trace.rays_per_trace == 1000 and cfg.trace.n_chunks == 1)

    # --rays-per-trace alone is measured against config.toml's budget.
    alone = ap.parse_args(["sweep", "--rays-per-trace", "30000"])
    cfg4, _ = CLI._load_with_overrides(alone)
    check("--rays-per-trace alone re-chunks config.toml's budget",
          cfg4.trace.rays_per_trace == 30000
          and sum(cfg4.trace.chunk_sizes) == cfg4.trace.rays_per_heliostat)

    # A budget that does not divide evenly must still sum EXACTLY: a dropped
    # remainder emits fewer rays than the run says it did, and every watt
    # reported from it is scaled by the shortfall.
    uneven = ap.parse_args(["sweep", "--rays", "100000", "--rays-per-trace", "30000"])
    cfg5, _ = CLI._load_with_overrides(uneven)
    check("a non-dividing pair keeps a short final chunk that closes the sum",
          cfg5.trace.chunk_sizes == [30000, 30000, 30000, 10000]
          and sum(cfg5.trace.chunk_sizes) == 100000 and cfg5.trace.n_chunks == 4)

    # Contradictory flags are refused, not silently clamped: "--rays 6000
    # --rays-per-trace 30000" asked for two different things explicitly.
    contradiction = ap.parse_args(["sweep", "--rays", "6000",
                                   "--rays-per-trace", "30000"])
    try:
        CLI._load_with_overrides(contradiction)
        refused, message = False, ""
    except SystemExit as exc:
        refused, message = True, str(exc)
    check("rays_per_trace > rays fails at startup, naming both numbers",
          refused and "6,000" in message and "30,000" in message
          and "rays_per_trace" in message)

    zero = ap.parse_args(["sweep", "--rays-per-trace", "0"])
    try:
        CLI._load_with_overrides(zero)
        refused_zero = False
    except SystemExit:
        refused_zero = True
    check("--rays-per-trace 0 fails rather than dividing by zero", refused_zero)

    # prime_focus with no focus height is exactly the mistake that would aim 645
    # heliostats at a point nobody chose. It must fail here, not mid-sweep.
    bad = ap.parse_args(["sweep", "--secondary", "prime_focus", "--n-mirrors", "1"])
    try:
        CLI._load_with_overrides(bad)
        failed = False
        message = ""
    except SystemExit as exc:
        failed, message = True, str(exc)
    check("prime_focus without --focus-height-mm fails at startup",
          failed and "focus_height_mm" in message)

    # A self-consistent layout + n_mirrors must not warn, even though config.toml
    # describes a different layout: the override is applied before validation.
    import warnings

    good = ap.parse_args(["sweep", "--secondary", "prime_focus", "--n-mirrors", "1",
                          "--focus-height-mm", "27000"])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cfg2, _ = CLI._load_with_overrides(good)
    check("a consistent prime_focus command line warns about nothing",
          cfg2.optics.n_mirrors == 1 and not caught)

    # And an inconsistent one still does, because the check really did re-run.
    mismatch = ap.parse_args(["sweep", "--secondary", "prime_focus",
                              "--focus-height-mm", "27000"])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        CLI._load_with_overrides(mismatch)
    check("prime_focus with config.toml's n_mirrors = 2 still warns",
          any("n_mirrors" in str(w.message) for w in caught))
    return ok


def check_exports(root, g, cfg) -> bool:
    """Every tab's two exports, run for real into a temp directory.

    The dialogs are skipped, not faked: ``_save_figure_dialog`` is a
    ``filedialog`` call wrapping ``_export_figure_to(name, path)``, and it is
    that inner function -- the one the button actually does its work in -- which
    is exercised here. Nothing is written outside ``tempfile``; in particular
    nothing goes near ``analysis_output/``.

    What is checked for each tab: a PNG and a PDF land with plausible sizes, the
    PNG really is 600 dpi (read back out of the file rather than trusted), the
    CSV parses with pandas and is not empty, and the default filename says what
    the file contains rather than "figure1".
    """
    import tempfile

    import pandas as pd

    from beamdown import plot_style

    ok = True

    def check(label, condition):
        nonlocal ok
        ok &= bool(condition)
        print(f"    {'OK  ' if condition else 'FAIL'} {label}")

    # The Energy tab exports the cached annual walk; make sure it is there
    # before asking for it, exactly as the tab itself does.
    if not g._energy_cache:
        g._ensure_energy()
        for _ in range(600):
            if g._energy_cache:
                break
            pump(root, 1)

    # -- the style is actually installed --------------------------------
    import matplotlib

    rc = matplotlib.rcParams
    check(f"paper style applied: white figure, {rc['lines.linewidth']:g} pt data "
          f"lines, constrained layout, savefig {rc['savefig.dpi']:g} dpi",
          rc["figure.facecolor"] == "white" and rc["savefig.facecolor"] == "white"
          and float(rc["lines.linewidth"]) >= 2.0
          and bool(rc["figure.constrained_layout.use"])
          and float(rc["savefig.dpi"]) == plot_style.EXPORT_DPI)
    check("every GUI figure carries a layout engine, so none uses the "
          "default margins",
          all(f.get_layout_engine() is not None for f in g.figures.values()))

    tabs = ["Field", "Spot", "Through day", "Distribution", "Energy", "Design",
            "Table"]
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for name in tabs:
            stem = g._export_stem(name)
            # The name has to identify the content: run or layout, the view, and
            # whatever the view is keyed on. "figure1" would pass a len() check,
            # so this asserts the actual parts are in it.
            keyed = (name == "Design" or stem.startswith(Path(RUN).name))
            check(f"{name:<13s} default name {stem!r}",
                  keyed and " " not in stem and len(stem) > 8)

            written = []
            if name in g.figures:
                written = g._export_figure_to(name, tmp / stem)
                exts = sorted(p.suffix for p in written)
                png = [p for p in written if p.suffix == ".png"][0]
                pdf = [p for p in written if p.suffix == ".pdf"][0]
                # 600 dpi is read back out of the PNG's own pHYs chunk (pixels
                # per metre), not taken on trust from the call we just made.
                from PIL import Image

                with Image.open(png) as im:
                    dpi_x = im.info.get("dpi", (0, 0))[0]
                    size = im.size
                check(f"{name:<13s} PNG+PDF {exts}, {size[0]}x{size[1]} px at "
                      f"{dpi_x:.0f} dpi, {png.stat().st_size//1024} kB / "
                      f"{pdf.stat().st_size//1024} kB",
                      exts == [".pdf", ".png"] and abs(dpi_x - 600) < 1
                      and png.stat().st_size > 5000 and pdf.stat().st_size > 1000)

            csvs = g._export_data_to(name, tmp / stem)
            frames = [pd.read_csv(c) for c in csvs]
            shapes = ", ".join(f"{c.name} {f.shape[0]}x{f.shape[1]}"
                               for c, f in zip(csvs, frames))
            check(f"{name:<13s} CSV parses: {shapes}",
                  bool(csvs) and all(len(f) > 0 and len(f.columns) > 1
                                     for f in frames))

        # -- the flux CSV must BE the flux map, not a re-derivation -------
        #
        # The whole point of "save the processed data" is that it is the same
        # numbers as the picture beside it. Summing the exported grid over its
        # bin area has to reproduce the power the tab reports.
        g.var_spotview.set("image")
        g._draw_spot()
        flux = g._flux_frame()
        field_panel = flux[flux.panel.str.startswith("all ")]
        exported_w = float(field_panel.flux_w_m2.sum()) * g.bin_area_m2
        shown_w = float(g._field_flux().sum()) * g.bin_area_m2
        print(f"    flux CSV total {exported_w/1e3:,.1f} kW vs the map's "
              f"{shown_w/1e3:,.1f} kW")
        check("the exported flux grid is the drawn flux map",
              abs(exported_w - shown_w) < 1e-6 * max(shown_w, 1.0)
              and len(field_panel) == g._bins() ** 2)

        # -- and the encircled CSV must be the encircled curve -------------
        #
        # Not "ends at the field total": the curve runs to the window radius,
        # and the square grid's corner bins lie outside that circle, so the last
        # point is legitimately a hair under the total. The check that means
        # something is the physical one -- power inside radius r, off the
        # exported curve, must equal power inside the same circular mask of the
        # drawn map, which is the identity main() already pins for the plot.
        import numpy as _np

        from beamdown.metrics import radial_mask

        g.var_spotview.set("encircled")
        g._draw_spot()
        enc = g._encircled_frame()
        field_curve = enc[enc.panel.str.startswith("all ")]
        r_test = 700.0
        from_csv = float(_np.interp(r_test, field_curve.radius_mm,
                                    field_curve.enclosed_power_w))
        from_map = float(g._field_flux()[radial_mask(cfg, r_test, g._bins())].sum()
                         * g.bin_area_m2)
        last = float(field_curve.enclosed_power_w.iloc[-1])
        print(f"    encircled CSV at r{r_test:.0f} {from_csv/1e3:,.1f} kW vs the "
              f"map inside the same circle {from_map/1e3:,.1f} kW; curve ends at "
              f"{last/1e3:,.1f} kW of {shown_w/1e3:,.1f} kW total (corner bins "
              f"lie outside the window radius)")
        check("the exported encircled curve is the drawn map's enclosed power",
              abs(from_csv - from_map) < 2e-3 * max(from_map, 1.0)
              and last <= shown_w * (1 + 1e-12)
              and (shown_w - last) < 1e-2 * shown_w)
        g.var_spotview.set("image")

        # -- plot_style's own contract ------------------------------------
        #
        # save_figure is what every button and every script goes through, so its
        # extension handling is load-bearing: a path that came from a "Save
        # as..." dialog arrives with .png already on it and must not become
        # "name.png.pdf".
        written = plot_style.save_figure(g.figures["Distribution"],
                                         tmp / "already.png")
        check(f"save_figure strips an existing extension: "
              f"{[p.name for p in written]}",
              sorted(p.name for p in written) == ["already.pdf", "already.png"])

    return ok


def check_design_tab(root, g) -> bool:
    """The Design tab: it evaluates geometry, never the loaded run.

    **Nothing here builds a model.** The export buttons shell out to the
    ``scripts/build_*_model.py`` programs, which write into ``models/``; what is
    checked is the argv they would be handed, exactly as ``check_command_builder``
    does for the Trace tab, so the test stays side-effect free.

    The tab is also the one view that must work with no run at all -- it reads
    the field file and config through ``beamdown.design_eval`` and nothing else
    -- so the first check removes the summary outright and re-evaluates.
    """
    from beamdown import design_eval as DE

    ok = True

    def check(label, condition):
        nonlocal ok
        ok &= bool(condition)
        print(f"    {'OK  ' if condition else 'FAIL'} {label}")

    # -- the tab must not need a loaded run -----------------------------
    summary, store = g.summary, g.store
    g.summary, g.store = None, None
    try:
        g._design_refresh()
        independent = bool(g._design_result) and g._design_drawn
    except Exception as exc:
        independent = False
        print(f"    Design tab touched the run: {type(exc).__name__}: {exc}")
    finally:
        g.summary, g.store = summary, store
    check("the tab evaluates and draws with no run loaded", independent)

    # -- the built axicon is the reference, and it is computed, not typed ----
    built = DE.built_axicon()
    check(f"built axicon indexes to exactly 1.0 "
          f"(cap {built['max_sagittal_correction']:.4e} /mm, "
          f"r90 {built['r90_mm']:.1f} mm)",
          abs(DE.eval_axicon(DE.BUILT_TIP_MM, DE.BUILT_ANGLE_DEG)["energy_index"]
              - 1.0) < 1e-12
          and abs(built["max_sagittal_correction"] - 7.115e-06) < 1e-9)

    # -- every layout evaluates, and the GUI shows what design_eval returned --
    for layout, params, want in (
            ("axicon", {"tip": 27.0, "angle": 20.0},
             DE.eval_axicon(27000.0, 20.0)),
            ("cassegrain", {"rim": 30.0, "f1": 36.0},
             DE.eval_cassegrain(30000.0, 36000.0)),
            ("prime_focus", {"pf": 36.0}, DE.eval_prime_focus(36000.0))):
        g.var_design_layout.set(layout)
        for key, value in params.items():
            g._design_vars[key].set(value)
        g._design_layout_changed()
        g._design_refresh()
        pump(root, 3)
        got = g._design_result
        check(f"{layout:<12s} energy {got.get('energy_index', float('nan')):.4f}x, "
              f"occlusion {got.get('occlusion', 0)*100:.1f}%, "
              f"r90 {got.get('r90_mm', 0):.0f} mm -- matches design_eval",
              got.get("layout") == layout and got.get("feasible")
              and got["energy_index"] == want["energy_index"]
              and got["r90_mm"] == want["r90_mm"])
        shown = [lbl.cget("text") for lbl in g._design_lines if lbl.cget("text")]
        check(f"{layout:<12s} readout is written out ({len(shown)} lines, "
              f"first: {shown[0] if shown else '<empty>'})", len(shown) >= 4)

    # -- an impossible design refuses, in words, rather than drawing a lie ---
    g.var_design_layout.set("cassegrain")
    g._design_vars["rim"].set(30.0)
    g._design_vars["f1"].set(41.0)
    g._design_layout_changed()
    g._design_refresh()
    bad = g._design_result
    check(f"F1 above the aperture-fill limit is refused: "
          f"{(bad.get('notes') or ['<silent>'])[0][:60]}…",
          not bad.get("feasible") and bool(bad.get("notes")))
    check("the full-field cassegrain button is disabled while infeasible",
          str(g.btn_design_full.cget("state")) == "disabled")

    # -- the export argv, without exporting anything ------------------------
    #
    # Same contract as the Trace tab's command builder: what the button hands to
    # the builder is checked, the builder is never run. --force never appears --
    # an existing target must draw the script's own refusal.
    seen = {}
    real_run = g._design_run
    g._design_run = lambda cmds, chain_on_output=None, chain_suffix="_dish": \
        seen.update(cmds=cmds, chain=chain_on_output, suffix=chain_suffix)
    try:
        g.var_design_date.set("2026-02-20")
        g.var_design_hour.set("9.454")

        g.var_design_layout.set("axicon")
        g._design_vars["tip"].set(29.0)
        g._design_vars["angle"].set(18.5)
        g._design_layout_changed()
        g._design_refresh()
        g._design_export_figure()
        argv = seen["cmds"][0]
        print(f"    {' '.join(argv[1:])}")
        check("axicon export carries the cone overrides",
              argv[1].endswith("build_figure_model.py")
              and argv[argv.index("--tip-height-mm") + 1] == "29000.0"
              and argv[argv.index("--axicon-angle-deg") + 1] == "18.5"
              and argv[argv.index("--date") + 1] == "2026-02-20"
              and argv[argv.index("--hour") + 1] == "9.454")
        check("axicon export chains nothing", seen["chain"] is None)
        check("no --force is ever passed", "--force" not in argv)

        g.var_design_layout.set("cassegrain")
        g._design_vars["rim"].set(30.0)
        g._design_vars["f1"].set(36.0)
        g._design_layout_changed()
        g._design_refresh()
        g._design_export_figure()
        argv, chain = seen["cmds"][0], seen["chain"]
        print(f"    {' '.join(argv[1:])}")
        print(f"    then {' '.join(chain[1:])} --base <written> --out <+{seen['suffix']}>")
        check("cassegrain export sets the layout, F1 and the rim",
              argv[argv.index("--secondary") + 1] == "cassegrain"
              and argv[argv.index("--focus-height-mm") + 1] == "36000.0"
              and argv[argv.index("--rim-height-mm") + 1] == "30000.0")
        check("cassegrain export chains the dish builder at the same geometry",
              chain is not None and chain[1].endswith("build_cassegrain_model.py")
              and chain[chain.index("--rim-z-mm") + 1] == "30000.0"
              and chain[chain.index("--f1-mm") + 1] == "36000.0"
              and "--force" not in chain and seen["suffix"] == "_dish30")

        g.var_design_layout.set("prime_focus")
        g._design_vars["pf"].set(38.0)
        g._design_layout_changed()
        g._design_refresh()
        g._design_export_figure()
        argv = seen["cmds"][0]
        print(f"    {' '.join(argv[1:])}")
        check("prime-focus export sets the layout and F1, and chains nothing",
              argv[argv.index("--secondary") + 1] == "prime_focus"
              and argv[argv.index("--focus-height-mm") + 1] == "38000.0"
              and seen["chain"] is None)
    finally:
        g._design_run = real_run

    # -- the builders really do accept those flags --------------------------
    #
    # The argv above is only worth checking if the scripts define those flags --
    # the same worry check_command_builder answers by parsing its line with
    # beamdown's own parser. Reading the source rather than importing the
    # scripts, because importing them runs their module-level sys.path surgery.
    for script, flags in (
            ("build_figure_model.py",
             ("--tip-height-mm", "--axicon-angle-deg", "--rim-height-mm",
              "--secondary", "--focus-height-mm")),
            ("build_cassegrain_model.py", ("--rim-z-mm", "--f1-mm", "--force"))):
        src = (REPO / "scripts" / script).read_text(encoding="utf-8")
        check(f"{script} defines {', '.join(flags)}",
              all(f'"{f}"' in src for f in flags))

    return ok


def check_trace_tab(root, g, cfg) -> bool:
    """The Trace tab: form -> command, and the launch path with Popen faked.

    The launch is never real. ``subprocess.Popen`` is replaced, the output and
    lock directories are a tmp path, and the lock scan is pointed at that tmp
    path rather than ``analysis_output`` -- a genuine launch would ask for the one
    HASP licence seat, and a failed request for it leaks the seat whatever is
    already running holds.
    """
    import subprocess
    import tempfile

    from beamdown import gui as G

    ok = True

    def check(label, condition):
        nonlocal ok
        ok &= bool(condition)
        print(f"    {'OK  ' if condition else 'FAIL'} {label}")

    def preview():
        return g.txt_t_cmd.get("1.0", "end-1c")

    # -- the tab must not need a loaded run -----------------------------
    #
    # It describes a sweep that has not happened, so nothing in it may read the
    # store. Removing the summary is the bluntest possible proof.
    summary = g.summary
    g.summary = None
    try:
        g._update_trace_command()
        independent = True
    except Exception as exc:
        independent = False
        print(f"    Trace tab touched the run: {type(exc).__name__}: {exc}")
    finally:
        g.summary = summary
    check("the tab rebuilds with no summary loaded", independent)

    # -- form -> command ------------------------------------------------
    check("default command is a valid sweep line",
          "-m beamdown sweep" in preview() and "--workers 1" in preview())

    g.lst_t_dates.selection_clear(0, "end")
    g.lst_t_dates.selection_set(0)
    g.lst_t_dates.selection_set(2)
    g.var_t_rays.set("1000")
    g.var_t_hour_step.set("2")
    g._update_trace_command()
    dates = [g.lst_t_dates.get(0), g.lst_t_dates.get(2)]
    check("a partial date selection is spelled out",
          all(d in preview() for d in dates) and "--dates" in preview())
    check("hour step differing from config appears", "--hour-step 2" in preview())
    check("rays follow the form", "--rays 1000" in preview())

    # -- rays per traceRays call, and the read-only count it implies ----------
    #
    # The entry is an input; the label is derived. The label is what the user is
    # actually choosing -- how many setRayDistributionCount1 + traceRays +
    # getRayPos round trips each heliostat costs -- so it must track the pair
    # live and must never disagree with the command being previewed.
    def chunk_label():
        return g.lbl_t_chunks.cget("text")

    g.var_t_rays.set("120000")
    g.var_t_rays_per_trace.set(str(cfg.trace.rays_per_trace))
    g._update_trace_command()
    check("at config.toml's chunk the flag is omitted, the label still derived",
          "--rays-per-trace" not in preview()
          and "traceRays call" in chunk_label())

    g.var_t_rays_per_trace.set("12000")
    g._update_trace_command()
    check("a chunk differing from config appears as a flag",
          "--rays-per-trace 12000" in preview())
    check("and the derived label counts the calls",
          "10 traceRays calls per heliostat" in chunk_label()
          and "120,000 rays" in chunk_label())

    g.var_t_rays.set("100000")
    g.var_t_rays_per_trace.set("30000")
    g._update_trace_command()
    check("an uneven budget is reported as a short final chunk",
          "4 traceRays calls per heliostat" in chunk_label()
          and "3 x 30,000 + 10,000" in chunk_label()
          and "= 100,000 rays" in chunk_label())

    g.var_t_rays_per_trace.set("100000")
    g._update_trace_command()
    check("one call per heliostat is sayable",
          "1 traceRays call per heliostat" in chunk_label()
          and "--rays-per-trace 100000" in preview())

    g.var_t_rays.set("6000")
    g._update_trace_command()
    check("a chunk larger than the budget refuses, as the CLI would",
          "will not run" in preview() and "rays per call" in preview()
          and "--rays-per-trace" not in preview())

    g.var_t_rays_per_trace.set("0")
    g._update_trace_command()
    check("a zero chunk refuses too", "must be positive" in preview())

    g.var_t_rays.set("1000")
    g.var_t_rays_per_trace.set(str(cfg.trace.rays_per_trace))
    g._update_trace_command()

    g.lst_t_dates.selection_set(0, "end")
    g._update_trace_command()
    check("selecting every configured date drops --dates as redundant",
          "--dates" not in preview())

    # -- free-text dates, and the suggest-dates alternative -------------
    g.var_t_newdate.set("2026-08-08")
    g._trace_add_date()
    check("a free-text date is added and selected",
          "2026-08-08" in list(g.lst_t_dates.get(0, "end"))
          and "2026-08-08" in preview())
    g.var_t_newdate.set("not-a-date")
    before = g.lst_t_dates.size()
    g._trace_add_date()
    check("a malformed date is rejected, not added",
          g.lst_t_dates.size() == before)

    g.var_t_use_suggest.set(True)
    g.var_t_suggest.set("6")
    g._update_trace_command()
    check("suggest-dates replaces the explicit list, never both",
          "--suggest-dates 6" in preview() and "--dates" not in preview())
    g.var_t_use_suggest.set(False)
    check("the estimate quotes traces and hours",
          "traces" in g.lbl_t_estimate.cget("text")
          and "s/trace" in g.lbl_t_estimate.cget("text"))

    # -- the layout choice injects the overrides the layout needs -------
    g.var_t_secondary.set("cassegrain")
    g._on_trace_layout()
    check("cassegrain enables both height entries",
          str(g.ent_t_focus.cget("state")) == "normal"
          and str(g.ent_t_rim.cget("state")) == "normal")
    check("cassegrain with no heights refuses, naming them",
          "will not run" in preview() and "focus height" in preview())
    g.var_t_focus.set("24000")
    g.var_t_rim.set("20000")
    g._update_trace_command()
    check("filled in, the layout flags appear",
          "--secondary cassegrain" in preview()
          and "--focus-height-mm 24000" in preview()
          and "--rim-height-mm 20000" in preview())

    g.var_t_secondary.set("prime_focus")
    g._on_trace_layout()
    check("prime_focus follows n_mirrors to 1 by itself",
          g.var_t_nmirrors.get() == "1" and "--n-mirrors 1" in preview())
    check("prime_focus disables the rim height it has no use for",
          str(g.ent_t_rim.cget("state")) == "disabled")
    g.var_t_secondary.set(cfg.optics.secondary)
    g._on_trace_layout()
    check("back to the configured layout, no --secondary at all",
          "--secondary" not in preview() and "--n-mirrors" not in preview())

    # -- flat heliostats: the fourth comparison axis, orthogonal to the layout --
    check("the checkbox starts where config.toml is",
          g.var_t_flat.get() == bool(cfg.optics.flat_mirrors)
          and not any(f in preview() for f in ("--flat-mirrors", "--focused-mirrors")))
    g.var_t_flat.set(not bool(cfg.optics.flat_mirrors))
    g._update_trace_command()
    want_flag = "--focused-mirrors" if cfg.optics.flat_mirrors else "--flat-mirrors"
    check(f"ticking it against config emits {want_flag}", want_flag in preview())
    opts, _problems = g._trace_options()
    check("and the options dict carries it as a real boolean",
          opts.get("flat_mirrors") is (not bool(cfg.optics.flat_mirrors)))
    # It must survive a layout change: the two are independent knobs.
    g.var_t_secondary.set("cassegrain")
    g.var_t_focus.set("24000")
    g.var_t_rim.set("20000")
    g._on_trace_layout()
    check("it composes with a layout override rather than being reset by it",
          want_flag in preview() and "--secondary cassegrain" in preview())
    g.var_t_secondary.set(cfg.optics.secondary)
    g._on_trace_layout()
    g.var_t_flat.set(bool(cfg.optics.flat_mirrors))
    g._update_trace_command()
    check("back at config.toml's value it says nothing again",
          not any(f in preview() for f in ("--flat-mirrors", "--focused-mirrors"))
          and "flat_mirrors" not in g._trace_options()[0])

    # -- a fixed-figure run's manifest key is inert here ----------------
    #
    # The Trace tab does not offer --fixed-shapes (it is a per-run table built
    # offline), but the GUI must still LOAD a run that used it. Every manifest
    # read in this module is a .get() by name, so an unknown key can only hurt if
    # something enumerates -- checked here rather than assumed.
    g.store.manifest["fixed_shapes"] = "data/fixed_shapes_fixture.csv"
    line = g.run_optics_label()
    g.store.manifest.pop("fixed_shapes")
    check("a manifest carrying fixed_shapes still reads its optics line",
          "heliostats" in line)

    # -- workers are capped at 1 until explicitly unlocked --------------
    check("workers default to 1", g.var_t_workers.get() == "1"
          and int(g.spin_t_workers.cget("to")) == 1)
    g.var_t_unlock_workers.set(True)
    g._trace_unlock_workers()
    check("unlocking raises the cap to the HASP key's 4",
          int(g.spin_t_workers.cget("to")) == 4)
    g.var_t_workers.set("3")
    g._update_trace_command()
    check("and the command says so", "--workers 3" in preview())
    g.var_t_unlock_workers.set(False)
    g._trace_unlock_workers()
    check("re-locking snaps back to 1",
          g.var_t_workers.get() == "1" and "--workers 1" in preview())

    # -- occluders and the model file are mutually exclusive ------------
    g.var_t_occluders.set(True)
    g._on_trace_occluders()
    check("--occluders appears and the model box goes read-only",
          "--occluders" in preview()
          and str(g.ent_t_model.cget("state")) == "disabled")
    g.var_t_occluders.set(False)
    g._on_trace_occluders()
    g.var_t_model.set("models/heliostat_field_model.optx")
    g._update_trace_command()
    check("a model file otherwise appears",
          "--model-file models/heliostat_field_model.optx" in preview())
    g.var_t_model.set(cfg.trace.model_file)

    # -- launching, with Popen faked -----------------------------------
    class FakePopen:
        calls: list = []

        def __init__(self, argv, **kw):
            FakePopen.calls.append((list(argv), kw))
            self.pid = 4242
            self.code = None

        def poll(self):
            return self.code

    warned: list = []
    real_popen = subprocess.Popen
    real_warn = G.messagebox.showwarning
    real_error = G.messagebox.showerror
    real_ask = G.messagebox.askyesno
    subprocess.Popen = FakePopen
    G.messagebox.showwarning = lambda title, msg, **k: warned.append(msg)
    G.messagebox.showerror = lambda title, msg, **k: warned.append(msg)
    G.messagebox.askyesno = lambda *a, **k: True

    try:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            # Point everything at the tmp directory, including the lock scan --
            # otherwise the real analysis_output/.full8.lock (a 24-hour run) would
            # refuse every launch below, and the launch path would go untested.
            g._trace_lock_roots = lambda: [tmp]
            g.var_t_output.set(str(tmp / "run1").replace("\\", "/"))
            g.var_t_rays.set("1000")
            g._update_trace_command()

            g._run_trace()
            check("Popen was called exactly once", len(FakePopen.calls) == 1)
            argv, kw = FakePopen.calls[-1]
            check("argv is exactly what the preview showed",
                  argv == g._trace_argv
                  and G.format_command(argv) == preview().split("\n")[0])
            check("it runs from the repository root",
                  Path(kw["cwd"]) == Path(cfg.repo_root))
            check("stderr is folded into stdout",
                  kw["stderr"] == subprocess.STDOUT)
            log_path = G.trace_log_path(tmp / "run1")
            check("stdout goes to the per-run log beside the output directory",
                  getattr(kw["stdout"], "name", "") == str(log_path)
                  and log_path.exists())
            check("the log records the command that was launched",
                  "-m beamdown sweep" in log_path.read_text(encoding="utf-8"))
            if sys.platform == "win32":
                flags = kw.get("creationflags", 0)
                check("detached, in its own process group, so it outlives the GUI",
                      bool(flags & subprocess.DETACHED_PROCESS)
                      and bool(flags & subprocess.CREATE_NEW_PROCESS_GROUP))
            else:
                check("start_new_session so it outlives the GUI",
                      kw.get("start_new_session") is True)

            lock = G.trace_lock_dir(tmp / "run1")
            check("the run's lock is taken, holding the child's pid",
                  lock.is_dir() and (lock / "pid").read_text().strip() == "4242")
            check("Run is disabled and Stop enabled while it lives",
                  str(g.btn_t_run.cget("state")) == "disabled"
                  and str(g.btn_t_stop.cget("state")) == "normal")

            # A second launch while the first lives must not start anything.
            g._run_trace()
            check("no second sweep while one is running", len(FakePopen.calls) == 1)

            # The monitor reads the log and the progress out of it.
            with open(log_path, "a", encoding="utf-8") as fh:
                fh.write("[2/44] 2026-03-20 09:00  el=35.5  1.0s  eta 12.5 min\n")
            g._poll_trace(reschedule=False)
            check("progress and ETA are parsed out of the log",
                  "2/44" in g.var_t_progress.get()
                  and "12" in g.var_t_progress.get())
            check("the log tail is shown",
                  "eta 12.5 min" in g.txt_t_log.get("1.0", "end-1c"))

            # The child exits; the lock must go and Run must come back.
            g._trace_proc.code = 0
            g._poll_trace(reschedule=False)
            check("a finished run releases its lock",
                  not lock.exists() and g._trace_proc is None)
            check("Run is available again",
                  str(g.btn_t_run.cget("state")) == "normal"
                  and "finished" in g.var_t_state.get())

            # -- refusals, at the button ---------------------------------
            #
            # The faked launch wrote no run directory, so stand one in for the
            # sweep that would have: an existing output directory is what the
            # collision guard is about.
            (tmp / "run1").mkdir()
            warned.clear()
            g._update_trace_command()
            g._run_trace()          # run1 exists now, resume is off
            check("relaunching over the finished run is refused",
                  len(FakePopen.calls) == 1 and warned
                  and "resume" in warned[-1])
            g.var_t_resume.set(True)
            g._update_trace_command()
            g._run_trace()
            check("with resume ticked it launches", len(FakePopen.calls) == 2)
            g._trace_proc.code = 0
            g._poll_trace(reschedule=False)

            warned.clear()
            other = tmp / ".someoneelse.lock"
            other.mkdir()
            (other / "pid").write_text("1675", encoding="utf-8")
            g.var_t_output.set(str(tmp / "run2").replace("\\", "/"))
            g.var_t_resume.set(False)
            g._update_trace_command()
            check("a foreign lock is shown on the tab",
                  "someoneelse" in g.lbl_t_locks.cget("text")
                  and "1675" in g.lbl_t_locks.cget("text"))
            g._run_trace()
            check("and refuses the launch, quoting the lock and its pid",
                  len(FakePopen.calls) == 2 and warned
                  and ".someoneelse.lock" in warned[-1] and "1675" in warned[-1])
            check("no lock was taken for the refused run",
                  not G.trace_lock_dir(tmp / "run2").exists())
    finally:
        subprocess.Popen = real_popen
        G.messagebox.showwarning = real_warn
        G.messagebox.showerror = real_error
        G.messagebox.askyesno = real_ask
        # Never leave a fake process looking live: _on_close would then ask
        # about it, and root.destroy() below must not block.
        g._trace_proc = None
        g._trace_log_path = None
        g._trace_lock_dir = None
    return ok


def main() -> int:
    import numpy as np

    root, g, cfg = build()
    ok = True
    print(f"{RUN}: {len(g.ids)} heliostats, {len(g.keys)} timesteps")

    # -- every view draws ------------------------------------------------
    for name, fn in (("Field", g._draw_field), ("Spot", g._draw_spot),
                     ("Through day", g._draw_curve), ("Distribution", g._draw_hist),
                     ("Table", g._draw_table)):
        t0 = time.perf_counter()
        fn()
        print(f"  {name:14s} {time.perf_counter()-t0:5.2f}s")
    g._draw_readout()

    # -- what the GUI shows must equal what the summary says -------------
    g.var_shading.set(True)
    g._field_cache.clear()
    power_map = g._field_flux().sum() * cfg.receiver.bin_area_m2
    power_tab = g._step_rows().power_w.sum()
    rel = abs(power_map - power_tab) / power_tab
    print(f"\n  field flux {power_map:,.1f} W vs summary {power_tab:,.1f} W "
          f"-> rel err {rel:.2e}")
    ok &= rel < 1e-9

    # -- the weights toggle has to reach every view, not just the spot ---
    #
    # power_w was written with shading x blocking already folded in, so a view
    # reading the summary directly shows weighted numbers whatever the checkbox
    # says. That is what made the weights look as though they were applied to
    # the field spot but not to a single heliostat's.
    g.var_shading.set(True)
    g._on_weights()
    on = g._rows_for_display()
    # Which columns carry the weight depends on the run: a sweep that traced its
    # occluders has shading and blocking in the ray counts already, and only the
    # secondary's shadow left as a scalar.
    eta = g._eta_series(on)
    g.var_shading.set(False)
    g._on_weights()
    off = g._rows_for_display()

    # off x eta must reproduce on, including where eta is exactly zero -- the
    # axicon shades some heliostats completely, and a naive divide-back-out
    # turns that into a nan rather than the unshaded power.
    worst = float(np.abs(off.power_w.to_numpy(float) * eta
                         - on.power_w.to_numpy(float)).max())
    n_zero = int((eta == 0).sum())
    print(f"\n  weights off  : unweighted x eta reproduces the weighted power, "
          f"worst {worst:.2e} W ({n_zero} heliostat(s) at eta = 0)")
    ok &= worst < 1e-6 and np.isfinite(off.power_w.to_numpy(float)).all()
    ok &= float(off.shading_blocking_efficiency.max()) == 1.0
    # The raw summary must be left alone -- the adjustment is a display layer.
    ok &= abs(g._step_rows().power_w.sum() - power_tab) < 1e-6

    # -- the aperture has to change the numbers, not just the axis limits -
    g.var_shading.set(True)
    captured = {}
    for radius in ("300", "700", ""):
        g.var_aperture.set(radius)
        g._on_weights()
        captured[radius] = float(g._rows_for_display().power_in_aperture_w.sum())
    print(f"  aperture     : r300 {captured['300']/1e3:,.1f} kW  "
          f"r700 {captured['700']/1e3:,.1f} kW  full {captured['']/1e3:,.1f} kW")
    ok &= captured["300"] < captured["700"] < captured[""] * (1 + 1e-12)
    # With no aperture set, every landed ray counts, so this is just power_w.
    ok &= abs(captured[""] - power_tab) / power_tab < 1e-12

    # -- the encircled-energy curve must agree with the image it replaces -
    from beamdown.metrics import bin_radius, radial_mask

    g.var_aperture.set("700")
    g._on_weights()
    area = cfg.receiver.bin_area_m2
    rr = bin_radius(cfg).ravel()
    order = np.argsort(rr)
    worst = 0.0
    for _label, flux in g._spot_pair():
        cum = np.concatenate(([0.0], np.cumsum(flux.ravel()[order] * area)))
        curve = cum[np.searchsorted(rr[order], 700.0, side="right")]
        image = float(flux[radial_mask(cfg, 700.0)].sum() * area)
        worst = max(worst, abs(curve - image) / max(image, 1e-9))
    print(f"  encircled    : curve at r700 matches the image, worst {worst:.2e}")
    ok &= worst < 1e-9

    for view in ("image", "encircled"):
        g.var_spotview.set(view)
        g._draw_spot()

    # -- rebinning changes the picture, never the power ------------------
    #
    # Coarsening block-sums the stored bins; refining goes back to the raw rays
    # and must reapply each heliostat's weight, which a single histogram over the
    # concatenated ray file would silently drop.
    stored = int(cfg.receiver.grid_size)
    totals, peaks = {}, {}
    for bins in (stored // 4, stored, stored * 2):
        g.var_bins.set(str(bins))
        g._field_cache.clear()
        flux = g._field_flux()
        totals[bins] = float(flux.sum() * g.bin_area_m2)
        peaks[bins] = float(flux.max())
        ok &= flux.shape == (bins, bins)
    spread = max(totals.values()) - min(totals.values())
    print(f"\n  rebin        : totals {[f'{v/1e3:,.1f} kW' for v in totals.values()]}, "
          f"spread {spread:.2e} W")
    ok &= spread / max(totals.values()) < 1e-12
    # Peak flux, unlike power, *should* rise as the bins get smaller.
    ok &= peaks[stored // 4] < peaks[stored] <= peaks[stored * 2]
    g.var_bins.set(str(stored))
    g._field_cache.clear()

    # -- the field view is a real projection, and the shadows are real ---
    low = min(g.keys, key=lambda k: float(
        g.summary[g.summary.timestep == k].iloc[0]["solar_el_deg"]))
    g._set_step(g.keys.index(low))
    el = float(g._step_rows().solar_el_deg.iloc[0])
    _ids, outlines, shadows, cone = g._field_polygons()
    offset = np.hypot(*(shadows.mean(1) - outlines.mean(1)).T)
    want = getattr(cfg.field, "draw_pedestal_height_mm", 5000.0) / 1000.0 / np.tan(
        np.deg2rad(el))
    print(f"  field        : sun el {el:.1f}°, shadow offset {offset.mean():.2f} m "
          f"vs pedestal/tan(el) = {want:.2f} m")
    ok &= abs(offset.mean() - want) < 1e-6 and float(np.ptp(offset)) < 1e-6
    # A steeply tilted mirror projects to less than its own area from above.
    e1 = outlines[:, 1] - outlines[:, 0]
    e2 = outlines[:, 3] - outlines[:, 0]
    area = np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0])
    full = cfg.field.mirror_area_m2
    ok &= area.max() <= full * (1 + 1e-9) and area.min() < full * 0.9

    # The secondary's shadow is drawn at mirror height, so at this elevation it
    # must land beyond the aperture radius and stay a bounded polygon.
    centre = cone.mean(0)
    throw = cfg.geometry.secondary_height_mm / 1000.0 / np.tan(np.deg2rad(el))
    print(f"  secondary    : shadow centroid {np.hypot(*centre):.1f} m from the axis "
          f"(vertex throw {throw:.1f} m), {len(cone)} hull vertices")
    ok &= len(cone) >= 3 and np.hypot(*centre) > cfg.geometry.axicon_aperture_radius_mm / 1000.0

    # -- the secondary's PHYSICAL footprint, which is not its shadow ------
    #
    # The dashed silhouette above walks across the field with the sun. This one
    # is where the body actually hangs: a circle of axicon_aperture_radius_mm on
    # the axis, fixed. They are the same size, so the plan view has to keep them
    # distinguishable -- different linestyle, different colour, and a legend
    # entry each.
    want_r = cfg.geometry.axicon_aperture_radius_mm / 1000.0
    foot = g.secondary_footprint()
    radii = np.hypot(foot[:, 0], foot[:, 1])
    # A closed polyline: the last point repeats the first, so the centroid is
    # taken over the distinct ones.
    centre_off = float(np.hypot(*foot[:-1].mean(0)))
    print(f"  footprint    : {len(foot)} points at r {radii.min():.3f}-{radii.max():.3f} m "
          f"(want {want_r:.1f} m), centroid {centre_off:.2e} m from the axis")
    ok &= (abs(radii - want_r).max() < 1e-9
           and centre_off < 1e-9
           and np.allclose(foot[0], foot[-1])          # closed, so plot() joins up
           and foot is g.secondary_footprint())        # cached, not rebuilt

    # Drawn for the two layouts that have a body, absent for the one that does
    # not -- and prime_focus must not crash on the way past.
    # Both are restored in the finally: the Trace tab's rule is "emit a flag only
    # where it differs from config.toml", so a rim height left behind here would
    # silently delete --rim-height-mm from a later test's command preview.
    real_layout = cfg.optics.secondary
    real_rim = cfg.geometry.secondary_rim_height_mm
    drawn = {}
    try:
        for layout, height in (("axicon", None), ("cassegrain", 20000.0),
                               ("prime_focus", None)):
            object.__setattr__(cfg.optics, "secondary", layout)
            if height is not None:
                object.__setattr__(cfg.geometry, "secondary_rim_height_mm", height)
            g._footprint_cache = None
            g._poly_cache.clear()
            g._draw_field()
            ax = g.figures["Field"].axes[0]
            labels = [t.get_text() for t in ax.get_legend().get_texts()] \
                if ax.get_legend() else []
            drawn[layout] = (g.secondary_footprint() is not None,
                             any("above the field" in t for t in labels),
                             any("shadow" in t for t in labels))
    finally:
        object.__setattr__(cfg.optics, "secondary", real_layout)
        object.__setattr__(cfg.geometry, "secondary_rim_height_mm", real_rim)
        g._footprint_cache = None
        g._poly_cache.clear()
    for layout, (has_body, in_legend, shadow_in_legend) in drawn.items():
        print(f"  field plan   : {layout:12s} footprint drawn={has_body}  "
              f"legend: body={in_legend} shadow={shadow_in_legend}")
    ok &= drawn["axicon"] == (True, True, True)
    ok &= drawn["cassegrain"] == (True, True, True)
    ok &= drawn["prime_focus"] == (False, False, False)

    g._draw_field()

    # -- no heliostat selected is a legitimate state ---------------------
    g._select(None)
    for fn in (g._draw_field, g._draw_spot, g._draw_curve, g._draw_hist,
               g._draw_table, g._draw_readout):
        fn()
    ok &= g._selected_row() is None and g._heliostat_flux() is None
    ok &= len(g._spot_pair()) == 1
    g.var_spotview.set("encircled")
    g._draw_spot()
    g.var_spotview.set("image")
    print(f"  deselect     : all views drew with no selection, "
          f"spot shows {len(g._spot_pair())} panel")
    g._select(g.ids[0])
    ok &= len(g._spot_pair()) == 2

    # -- date and hour are separate axes of the same index ---------------
    g._set_date(0)
    g._set_hour(0)
    first = g.key
    g._step(+1)
    ok &= g._hour_of(g.key) > g._hour_of(first) and g.date_key == first.split("_")[0]
    g._step_date(+1)
    print(f"  time         : {first} -> +1 hour -> +1 date -> {g.key} "
          f"({len(g.dates)} dates)")
    ok &= g.date_key == g.dates[1]
    # Changing date keeps the nearest hour, so the two views stay comparable.
    ok &= abs(g._hour_of(g.key) - g._hour_of(first) - 1) < 1.01
    g._set_step(len(g.keys) // 2)

    # -- the table must not recurse through its own selection event ------
    calls = {"n": 0, "depth": 0, "max": 0}
    real = g._draw_table

    def counted():
        calls["n"] += 1
        calls["depth"] += 1
        calls["max"] = max(calls["max"], calls["depth"])
        if calls["depth"] > 12:
            raise RecursionError("runaway _draw_table recursion")
        try:
            real()
        finally:
            calls["depth"] -= 1

    g._draw_table = counted
    g._dirty.add("Table")
    tab = [i for i in range(g.book.index("end"))
           if g.book.tab(i, "text") == "Table"][0]
    g.book.select(tab)
    pump(root, 30)
    print(f"\n  open Table tab   : {calls['n']} draw(s), depth {calls['max']}, "
          f"{len(g.tree.get_children())} rows")
    ok &= calls["max"] <= 1

    rows = g.tree.get_children()
    target = rows[len(rows) // 2]
    want = int(g.tree.item(target, "values")[0])
    calls["n"] = 0
    g.tree.selection_set(target)
    pump(root, 30)
    print(f"  click a row      : {calls['n']} rebuild(s), selected {g.selected} "
          f"(clicked {want})")
    ok &= g.selected == want and calls["n"] == 0

    calls["n"] = 0
    t0 = time.perf_counter()
    for _ in range(4):
        g._step(+1)
        pump(root, 10)
    print(f"  step time x4     : {calls['n']} draw(s), depth {calls['max']}, "
          f"{time.perf_counter()-t0:.2f}s")
    ok &= calls["max"] <= 1

    g.var_colour.set("r90_mm")
    g._select_extreme("max")
    pump(root, 20)
    sel = g.tree.selection()
    shown = int(g.tree.item(sel[0], "values")[0]) if sel else None
    print(f"  external select  : selected {g.selected}, table shows {shown}")
    ok &= shown == g.selected

    # -- Energy tab: draws in the background and caches on redraw --------
    #
    # annual_energy walks 8760 hours, so the tab must not block the UI thread
    # and must not repeat the walk on a second draw. Counting calls to the
    # real function (rather than trusting a timer) is the only way to prove
    # the cache, not just the wall clock, is what made the second draw fast.
    import datetime as _dt

    import pandas as pd

    from beamdown import dni as D
    from beamdown import energy as E

    calls = {"n": 0}
    real_annual = E.annual_energy

    def counted_annual(*a, **kw):
        calls["n"] += 1
        return real_annual(*a, **kw)

    E.annual_energy = counted_annual
    try:
        energy_tab = [i for i in range(g.book.index("end"))
                      if g.book.tab(i, "text") == "Energy"][0]
        g.book.select(energy_tab)
        pump(root, 5)  # first draw: cache miss, kicks off the background thread
        for _ in range(400):
            if g._energy_cache:
                break
            pump(root, 1)
        pump(root, 10)
        first_calls = calls["n"]
        print(f"\n  Energy tab       : {'ready' if g._energy_cache else 'NOT READY'} "
              f"after background compute, annual_energy called {first_calls} time(s)")
        print(f"    headline: {g.lbl_energy_headline.cget('text')}")
        print(f"    annual  : {g.lbl_energy_annual.cget('text')}")
        print(f"    check   : {g.lbl_energy_check.cget('text')}")
        ok &= first_calls == 1 and bool(g._energy_cache)

        g._dirty.add("Energy")
        g._draw_energy()
        print(f"  Energy redraw    : annual_energy called {calls['n']} time(s) total "
              f"(still {first_calls} expected -- a cached redraw must not recompute)")
        ok &= calls["n"] == first_calls
    finally:
        E.annual_energy = real_annual

    # -- DNI mode selector: switching must recompute, not reuse the cache ---
    #
    # _energy_key folds the selected DNI mode in for exactly this reason: a
    # live switch has to invalidate the annual-energy cache, or the Energy tab
    # keeps showing the previous model's MWh after the control changes. (a)
    # proves the cached result for the new mode actually differs from the old
    # one; (b) is a sanity check that the difference is the right *shape* --
    # roughly proportional to the change in annual DNI totals, not some
    # unrelated number -- using a generous tolerance because the DNI-weighted
    # average optical efficiency is not quite identical between models (a flat
    # 1000 W/m2 weights every daylight hour equally; the real diurnal/seasonal
    # curve does not).
    from beamdown.gui import DNI_MODE_DEFAULT

    default_mode = g.var_dni_mode.get()
    ok &= default_mode == DNI_MODE_DEFAULT
    default_key = g._energy_key()
    default_result = g._energy_cache.get(default_key)
    ok &= default_result is not None

    g.var_dni_mode.set("constant")
    g._on_dni_mode()
    for _ in range(400):
        if g._energy_key() in g._energy_cache:
            break
        pump(root, 1)
    pump(root, 10)
    constant_key = g._energy_key()
    constant_result = g._energy_cache.get(constant_key)
    ok &= constant_result is not None and constant_key != default_key

    if default_result is not None and constant_result is not None:
        mwh_default = default_result["annual"]["annual_energy_mwh"]
        mwh_constant = constant_result["annual"]["annual_energy_mwh"]
        dni_default = default_result["annual"]["annual_dni_kwh_m2"]
        dni_constant = constant_result["annual"]["annual_dni_kwh_m2"]

        print(f"\n  DNI mode switch  : {default_mode} {mwh_default:,.1f} MWh "
              f"(DNI {dni_default:,.1f} kWh/m2)  vs  constant {mwh_constant:,.1f} MWh "
              f"(DNI {dni_constant:,.1f} kWh/m2)")
        ok &= mwh_constant != mwh_default  # (a) cache actually changed

        ratio_energy = mwh_constant / mwh_default
        ratio_dni = dni_constant / dni_default
        print(f"    ratio  energy {ratio_energy:.4f}  vs  DNI totals {ratio_dni:.4f}")
        ok &= abs(ratio_energy / ratio_dni - 1.0) < 0.15  # (b)

    # Leave the selector as found -- the default-mode result is already
    # cached, so this does not trigger another background compute.
    g.var_dni_mode.set(default_mode)
    g._on_dni_mode()

    # -- traced_day_energy vs annual_energy's daily row, on real data -----
    #
    # The two routes share nothing after power_w: one interpolates an
    # efficiency surface onto solar.hours_of_year's fixed 24-hour grid, the
    # other trapezoids the actual traced samples anchored at true sunrise and
    # sunset. The tolerance is generous because the anchor itself is an
    # approximation (a straight ramp from the first/last sample down to zero
    # at the horizon), which is worst on short-day, low-elevation dates.
    provider = D.load_dni_provider(cfg)
    annual = (g._energy_cache[g._energy_key()]["annual"] if g._energy_cache
              else E.annual_energy(g.summary, cfg, provider))
    checks = E.cross_check_daily_energy(g.summary, cfg, provider, annual=annual)
    print(f"\n  daily energy cross-check ({provider.describe()}):")
    for _, row in checks.iterrows():
        print(f"    {row['date']}  traced {row['traced_energy_kwh']/1000:8.2f} MWh   "
              f"interpolated {row['interpolated_energy_kwh']/1000:8.2f} MWh   "
              f"residual {row['residual_frac']:+7.2%}")
    ok &= len(checks) == len(g.dates)
    ok &= bool((checks["residual_frac"].abs() < 0.15).all())

    # -- sine fit recovers a known synthetic signal ------------------------
    #
    # Pins the least-squares maths on its own, independent of whatever the
    # real DNI series looks like: a signal built exactly from the model the
    # fit assumes must come back out to machine precision.
    doy = np.arange(1, 366, dtype=float)
    w = 2.0 * np.pi / 365.25
    A_true, C_true, D_true = 500.0, 120.0, -80.0
    y_true = A_true + C_true * np.sin(w * doy) + D_true * np.cos(w * doy)
    synth = pd.DataFrame({
        "date": [_dt.date(2026, 1, 1) + _dt.timedelta(days=int(d - 1)) for d in doy],
        "energy_kwh": y_true,
    })
    fit = E.fit_annual_sine(synth)
    print(f"\n  sine fit self-test: A {fit['A']:.6f} (want {A_true}), "
          f"C {fit['C']:.6f} (want {C_true}), D {fit['D']:.6f} (want {D_true}), "
          f"R2 {fit['r_squared']:.9f}")
    ok &= (abs(fit["A"] - A_true) < 1e-6 and abs(fit["C"] - C_true) < 1e-6
           and abs(fit["D"] - D_true) < 1e-6 and fit["r_squared"] > 1 - 1e-9)

    # -- energy scales linearly with DNI -----------------------------------
    #
    # DNI is documented (dni.py) as a pure post-processing multiplier on a
    # trace normalised to exactly 1000 W/m^2; this is the one-line proof.
    date0 = g.dates_as_dates[0]
    low = E.traced_day_energy(g.summary, cfg, D.ConstantDNI(400.0), date0)
    high = E.traced_day_energy(g.summary, cfg, D.ConstantDNI(800.0), date0)
    ratio = high["energy_kwh"] / low["energy_kwh"]
    print(f"\n  DNI linearity    : {date0} at 400 W/m2 -> {low['energy_kwh']:.1f} kWh, "
          f"800 W/m2 -> {high['energy_kwh']:.1f} kWh, ratio {ratio:.9f}")
    ok &= abs(ratio - 2.0) < 1e-9

    # -- Trace tab: setting up a sweep -----------------------------------
    print("\n  command builder (pure):")
    ok &= check_command_builder()
    print("\n  CLI overrides the Trace tab depends on:")
    ok &= check_cli_overrides()
    print("\n  Design tab (licence-free, run-independent):")
    ok &= check_design_tab(root, g)
    print("\n  Figure and data export (into a temp dir, never analysis_output):")
    ok &= check_exports(root, g, cfg)
    print("\n  Trace tab:")
    ok &= check_trace_tab(root, g, cfg)

    root.destroy()
    print("\nPASS" if ok else "\nFAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
