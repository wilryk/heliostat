"""Does a built figure model really carry 25 usable configurations?

NEEDS A LICENCE SEAT.  Everything checkable without one is already checked by
``scripts/build_figure_model.py --check``: that the XML is well formed, that the
column count is 25, that all 175 per-config values and 5 globals round-trip
through the file, and that nothing outside ``<multiconfig>`` moved.  This is the
other half -- what only Quadoa can answer -- and it CANNOT run while a sweep
holds the key.  Run it once, after building a figure model, and read every
line::

    python scripts/verify_figure_model.py models/figure_model_25cfg_20260621_1200.optx

It refuses to start while any ``analysis_output/.*.lock`` exists.  That is not
politeness: there is ONE seat, a failed acquisition raises a modal H0038 dialog
on the owner's desktop, and this script would be taking the seat away from a
multi-hour run.

Six checks, in the order a failure is cheapest to diagnose:

1. **Quadoa agrees there are 25 configurations.**  ``getNrConfigs() == 25``.
   The column count was grown by text surgery on the ``<multiconfig>`` header
   plus one ``<variable name="val_i">`` per ``<param>``; the Python API has no
   method that creates a configuration, so this is the only way to confirm the
   surgery is legible to the loader rather than merely well-formed XML.  If this
   reads 24, the header was rewritten but a ``<param>`` was left short (or the
   reverse) and Quadoa fell back to the smaller count.

2. **Per-config and global parameters are told apart, by Quadoa, not by us.**
   The ``.optx`` says ``rot_az`` is a ``<param>`` and ``solaz`` is a
   ``<single_param>``; the API confirms it with the NaN trick --
   ``getMulticonfParam(name, i)`` for ``i > 0`` returns a real number for a
   per-config parameter and NaN for a global one.  This is the check that would
   have caught the whole class of bug this exercise is about: a per-config write
   aimed at a global parameter appears to succeed, silently applies to the
   entire model, and leaves 25 identical configurations.  Asserted for all seven
   ``<param>`` and all five ``<single_param>``.

3. **The stored values are the values Quadoa loaded.**  Every one of the
   25 x 7 per-config values and 5 globals is read back through
   ``getMulticonfParam`` and compared to what the file holds.  Check 1 proves
   the columns exist; this proves they were populated rather than defaulted.
   It also catches the specific historical failure: the shipped
   ``models/figure_model_25cfg.optx`` has config 0 holding a stale heliostat and
   configs 1..24 all zero, which passes check 1 and fails here loudly.

4. **Ray-count semantics on the sequence being used.**  README records
   ``setRayDistributionCount1`` as a LITERAL ray count on sequences 0 and 3 and
   a per-axis GRID DENSITY on 1 and 2.  ``figure_seq = 1``, so a figure model
   asked for 60,000 would emit billions.  Probe both the analysis sequence and
   the figure sequence with 300 and count what comes out -- the measurement,
   not the inherited claim, because it is measured per sequence and this file
   is not the file it was measured on.

5. **The noise floor.**  Trace the SAME configuration three times.  The source
   is Monte-Carlo, so ray counts and centroids differ run to run with nothing
   changed; comparing two single traces proves nothing.  The floor is the
   largest centroid excursion across those three.

6. **Switching configuration moves the spot, well above that floor.**  Pick the
   two most widely separated of the 25 heliostats from the file, ``setConfig``
   between them, and require the centroid shift to exceed the floor by the
   margin in ``MARGIN``.  This is the check that fails when
   ``applyChangesAndInitModel`` is genuinely needed after ``setConfig`` and is
   not being called -- the failure mode that gives 25 plausible, identical
   spots.  ``QuadoaSession.reinit_after_params`` is False on the claim that
   selecting a configuration is enough; on a 25-column model that claim has
   never been tested, so if this check fails, re-run with ``--reinit`` before
   concluding the writes are broken.

Nothing here writes to the model file or to ``analysis_output/``.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import numpy as np  # noqa: E402

# How many times the config-switch signal must exceed the Monte-Carlo floor.
MARGIN = 10.0

# Probe count for the ray-count-semantics check. Small on purpose: if the
# sequence turns out to be a grid density, 300 still only costs ~70,000 rays.
PROBE_COUNT = 300

EXPECT_CONFIGS = 25


def find_locks() -> list[Path]:
    """Every ``analysis_output/.*.lock`` directory. Non-empty means: do not run."""
    root = REPO / "analysis_output"
    if not root.exists():
        return []
    return sorted(p for p in root.glob(".*.lock"))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("model", type=Path,
                   help="the built figure model, e.g. "
                        "models/figure_model_25cfg_20260621_1200.optx")
    p.add_argument("--config", type=Path, default=None, help="config.toml to use")
    p.add_argument("--rays", type=int, default=20000,
                   help="rays per trace for the noise-floor and switch checks")
    p.add_argument("--reinit", action="store_true",
                   help="call applyChangesAndInitModel after setConfig (check 6)")
    a = p.parse_args(argv)

    # -- the gate, before anything imports quadoa --------------------------
    locks = find_locks()
    if locks:
        print("REFUSING TO RUN -- a run holds the licence seat:")
        for lk in locks:
            print(f"  {lk}")
        print("\nThere is ONE seat. Taking it would fail the running sweep, and a")
        print("failed acquisition pops a modal H0038 dialog. Wait for the lock to")
        print("clear, then re-run. Do not delete a lock whose python is alive.")
        return 3

    if not a.model.exists():
        print(f"no such model: {a.model}")
        print("build one first:  python scripts/build_figure_model.py "
              "--date YYYY-MM-DD --hour H --check")
        return 2

    from build_figure_model import PARAMS, SINGLE_PARAMS, read_multiconfig
    from beamdown.config import load_config
    from beamdown.session import QuadoaSession

    cfg = load_config(a.config)
    text = a.model.read_text(encoding="utf-8")
    cols, want_p, want_g = read_multiconfig(text)
    print(f"model: {a.model}")
    print(f"  file says: columns={cols}, {len(want_p)} params, "
          f"{len(want_g)} single_params\n")

    results: list[tuple[bool, str]] = []

    def check(ok, msg):
        results.append((bool(ok), msg))
        print(f"  [{'ok' if ok else 'FAIL'}] {msg}")

    # The analysis sequence is used for the tracing checks: its
    # setRayDistributionCount1 is a literal ray count, which is what makes a
    # ray budget mean what it says. figure_seq is only probed, never traced at
    # volume -- see check 4.
    session = QuadoaSession(cfg, seq=cfg.trace.analysis_seq)
    try:
        core = session.core
        core.loadModelFile(str(a.model))
        core.applyChangesAndInitModel()
        session._resolve_surface()

        # -- 1. Quadoa sees 25 configurations -----------------------------
        print("1. configuration count")
        n = int(core.getNrConfigs())
        check(n == EXPECT_CONFIGS,
              f"getNrConfigs() = {n}, expected {EXPECT_CONFIGS} -- the "
              f"<multiconfig columns> surgery is legible to the loader")
        check(n == cols, f"loader agrees with the file's own header ({cols})")

        # -- 2. per-config vs global, by Quadoa's own answer ---------------
        print("\n2. per-config vs global (the NaN trick)")
        probe_cfg = min(1, n - 1)
        per_ok = [name for name in PARAMS
                  if not np.isnan(float(core.getMulticonfParam(name, probe_cfg)))]
        glob_ok = [name for name in SINGLE_PARAMS
                   if np.isnan(float(core.getMulticonfParam(name, probe_cfg)))]
        check(per_ok == PARAMS,
              f"all {len(PARAMS)} <param> read a real value at config "
              f"{probe_cfg}: {per_ok}")
        check(glob_ok == SINGLE_PARAMS,
              f"all {len(SINGLE_PARAMS)} <single_param> read NaN at config "
              f"{probe_cfg} (one value for the whole model): {glob_ok}")

        # -- 3. the loaded values are the stored values --------------------
        print("\n3. stored values survived the load")
        bad = []
        for name in PARAMS:
            for i in range(n):
                got = float(core.getMulticonfParam(name, i))
                if not np.isclose(got, want_p[name][i], rtol=0, atol=1e-9):
                    bad.append(f"{name}[{i}] {got!r} != {want_p[name][i]!r}")
        check(not bad,
              f"all {n * len(PARAMS)} per-config values match the file"
              + (f" -- {len(bad)} differ, first: {bad[0]}" if bad else ""))
        gbad = []
        for name in SINGLE_PARAMS:
            got = float(core.getMulticonfParam(name, 0))
            if not np.isclose(got, want_g[name], rtol=0, atol=1e-9):
                gbad.append(f"{name} {got!r} != {want_g[name]!r}")
        check(not gbad, f"all {len(SINGLE_PARAMS)} globals match the file"
                        + (f" -- {gbad}" if gbad else ""))

        xy = np.array([want_p["posx"], want_p["posy"]]).T
        distinct = len({tuple(r) for r in xy})
        check(distinct == n,
              f"{distinct}/{n} configurations hold distinct positions -- the "
              f"all-zeros failure mode is absent")

        # -- 4. ray-count semantics, measured on THIS file -----------------
        print("\n4. setRayDistributionCount1 semantics (measured, not assumed)")
        for label, seq in (("analysis_seq", cfg.trace.analysis_seq),
                           ("figure_seq", cfg.trace.figure_seq)):
            img = int(core.getSequenceImageSurface(seq))
            core.setRayDistributionCount1(seq, PROBE_COUNT)
            core.traceRays(seq, 0, 0)
            got = np.array(core.getRayPos(seq, 0, 0, img), copy=True).shape[1]
            kind = "LITERAL" if got == PROBE_COUNT else f"GRID ({got / PROBE_COUNT:.0f}x)"
            expect = "LITERAL" if seq in (0, 3) else "GRID"
            check(kind.split()[0] == expect,
                  f"{label}={seq} (image surface {img}): asked {PROBE_COUNT}, "
                  f"got {got} -> {kind}; README says {expect}")

        # -- 5. the Monte-Carlo noise floor --------------------------------
        print("\n5. noise floor: same configuration, three traces")
        far = int(np.argmax(np.hypot(xy[:, 0], xy[:, 1])))
        near = int(np.argmin(np.hypot(xy[:, 0], xy[:, 1])))
        session.activate(near, reinit=a.reinit)
        cents = []
        for k in range(3):
            res = session.trace(rays=a.rays)
            if not res.rays_landed:
                check(False, f"trace {k + 1} of config {near} landed nothing")
                break
            c = res.xy_mm.mean(axis=0)
            cents.append(c)
            print(f"      trace {k + 1}: {res.rays_landed} rays, centroid "
                  f"({c[0]:+.2f}, {c[1]:+.2f}) mm")
        if len(cents) < 3:
            raise SystemExit(1)
        cents = np.array(cents)
        floor = float(max(np.hypot(*(cents[i] - cents[j]))
                          for i in range(3) for j in range(i + 1, 3)))
        check(floor >= 0.0,
              f"noise floor = {floor:.3f} mm (largest centroid excursion across "
              f"3 identical traces of config {near})")

        # -- 6. switching configuration moves the spot ---------------------
        print(f"\n6. config {near} -> {far} (the two most separated heliostats)")
        base_c = cents.mean(axis=0)
        session.activate(far, reinit=a.reinit)
        res = session.trace(rays=a.rays)
        if not res.rays_landed:
            check(False, f"config {far} landed nothing")
        else:
            c = res.xy_mm.mean(axis=0)
            shift = float(np.hypot(*(c - base_c)))
            sep = float(np.hypot(*(xy[far] - xy[near])) / 1000.0)
            print(f"      config {far}: {res.rays_landed} rays, centroid "
                  f"({c[0]:+.2f}, {c[1]:+.2f}) mm")
            print(f"      heliostats are {sep:.1f} m apart in the field")
            check(shift > MARGIN * floor,
                  f"centroid moved {shift:.2f} mm = {shift / max(floor, 1e-9):.0f}x "
                  f"the noise floor (need >{MARGIN:.0f}x) -- setConfig reaches "
                  f"the trace"
                  + ("" if a.reinit else "; reinit_after_params=False holds on a "
                                         "25-column model"))
    finally:
        session.close()

    bad_n = sum(not ok for ok, _ in results)
    print(f"\n  {len(results) - bad_n}/{len(results)} checks passed")
    print("  " + ("PASS -- the figure model carries 25 real, distinct, traceable "
                  "configurations"
                  if not bad_n else
                  "FAIL -- do not use this figure model for a paper figure"))
    return 1 if bad_n else 0


if __name__ == "__main__":
    sys.exit(main())
