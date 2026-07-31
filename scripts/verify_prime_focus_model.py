"""Does ``models/heliostat_field_prime_focus.optx`` actually trace prime focus?

NEEDS A LICENCE SEAT. Everything the builder could check without one is already
checked by ``scripts/build_prime_focus_model.py``; this is the other half, and it
cannot run while a sweep holds the key. Run it once, before the first
prime-focus sweep, and read every line::

    python scripts/verify_prime_focus_model.py

Six checks, in the order that a failure is cheapest to diagnose:

1. **The parameter exists.** ``pf_height`` reads back as ``focus_height_mm``
   after ``set_global_geometry``. Writing a parameter a model does not have is
   silently ignored by Quadoa, so this is the one check that catches "the wrong
   ``.optx`` is loaded" -- against the stock model the write vanishes, the
   detector stays at its literal 27000, and every number below is plausible and
   wrong.
2. **The detector is where the sequence looks.**
   ``getSequenceImageSurface(3)`` must resolve to the ``prime_focus`` surface.
   ``config.toml``'s ``analysis_surface = 3`` is an index into a four-surface
   sequence and this one has three; ``QuadoaSession._resolve_surface`` is
   supposed to fix that by asking the sequence, and this is where that is
   confirmed rather than assumed.
3. **Ray-count semantics.** README records ``setRayDistributionCount1`` as a
   LITERAL ray count on sequence 3 and a GRID DENSITY on sequences 1 and 2, all
   measured on the axicon model. Rewriting sequence 3's surface list could in
   principle move it into the other regime, and a grid-density sequence asked
   for 60,000 emits hundreds of millions. Same probe as the README's: request
   300 and count what comes out.
4. **Writes reach the trace.** ``session.self_test`` with the prime-focus
   strategy: two well-separated heliostats must produce two different spots.
5. **The spot lands where the geometry says.** A heliostat at ``(x, y, 0)``
   aimed at ``(0, 0, 47000)`` puts its spot centre ON THE AXIS, because the
   detector plane *is* the aim plane -- so the predicted centroid is ``(0, 0)``
   for every heliostat in the field, with no ray tracing needed to predict it.
   That is a much sharper test than it looks: it fails if the detector is at the
   wrong height (the cone has moved on by then), if the pointing solve disagrees
   with the model's coordinate convention, or if ``rot_az``/``rot_el`` are
   reaching the wrong surface.
6. **The Zernikes reach this file's ``helio_surf``.** The same heliostat traced
   flat must throw a much larger spot. ``c3``/``c4``/``c5`` are written to a
   surface that sequence 3 still lists after the rewrite, but "still listed" and
   "still bent by the parameters" are different claims; a flat mirror at 30-120 m
   throws metres where a focused one throws centimetres, so the difference is
   unmissable if the writes land and invisible if they do not.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

MODEL = "models/heliostat_field_prime_focus.optx"
PF_HEIGHT = 47000.0

# Where the probe heliostat sits. Far enough out that a flat mirror's spot is
# unmistakably bigger than a focused one, and off both axes so a sign error in
# either coordinate shows up instead of cancelling.
PROBE = (22300.0, -60000.0)
PROBE_SUN = (77.9, 50.1)      # azimuth, elevation (deg)

# Tolerance on the centroid. The prediction is exact -- the aim point is on the
# axis and the detector plane is the aim plane -- so this is Monte-Carlo noise on
# the mean of ~45,000 rays spread over a spot of order 100 mm, plus whatever
# asymmetry the sun's 0.0024 rad cone leaves. 50 mm is generous by a wide margin
# and still an order of magnitude under any real misplacement: putting the
# detector at 27000 instead of 47000 does not shift the centroid, it changes the
# SIZE, which is what `r90` below is for.
CENTROID_TOL_MM = 50.0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", default=MODEL)
    ap.add_argument("--rays", type=int, default=120000)
    ap.add_argument("--focus-height-mm", type=float, default=PF_HEIGHT)
    args = ap.parse_args(argv)

    from beamdown.config import load_config, validate_layout
    from beamdown.secondary import get_strategy
    from beamdown.session import QuadoaSession

    cfg = load_config(None)
    model = cfg.repo_root / args.model
    if not model.exists():
        print(f"FAIL: {model} does not exist -- run "
              f"scripts/build_prime_focus_model.py first")
        return 1

    # The config on disk is the axicon's. Override to the prime-focus run rather
    # than requiring config.toml to be edited, so this is runnable at any time
    # (including while another layout's sweep is queued) and so the numbers it
    # prints cannot depend on what config.toml happened to say.
    object.__setattr__(cfg.trace, "model_file", str(model))
    object.__setattr__(cfg.trace, "rays_per_heliostat", args.rays)
    object.__setattr__(cfg.geometry, "focus_height_mm", float(args.focus_height_mm))
    object.__setattr__(cfg.optics, "secondary", "prime_focus")
    object.__setattr__(cfg.optics, "n_mirrors", 1)
    validate_layout(cfg)

    strategy = get_strategy(cfg)
    print(f"  model    : {model.name}")
    print(f"  layout   : {strategy.describe()}")
    print(f"  F1       : (0, 0, {cfg.geometry.focus_height_mm:.0f}) mm")
    print(f"  sequence : {cfg.trace.analysis_seq} (config.toml analysis_seq)")

    results: list[tuple[bool, str]] = []

    def check(ok: bool, msg: str) -> None:
        results.append((bool(ok), msg))
        print(f"  [{'ok' if ok else 'FAIL'}] {msg}")

    session = QuadoaSession(cfg)
    try:
        core = session.core

        # -- 1. pf_height is real and carries the config's value --------------
        print("\n1. the pf_height parameter")
        session.set_global_geometry()
        written = strategy.global_params(cfg.geometry)
        print(f"  global_params -> {written}")
        read = core.getMulticonfParam("pf_height", 0)
        got = float("nan") if read is None else float(read)
        check(np.isfinite(got),
              f"pf_height exists in the loaded model (read back {got!r})")
        check(np.isfinite(got) and abs(got - cfg.geometry.focus_height_mm) < 1e-6,
              f"pf_height == focus_height_mm: {got} vs "
              f"{cfg.geometry.focus_height_mm}")

        # -- 2. the sequence's own detector ------------------------------------
        print("\n2. the detector sequence 3 reports")
        image = int(core.getSequenceImageSurface(cfg.trace.analysis_seq))
        print(f"  getSequenceImageSurface({cfg.trace.analysis_seq}) = {image}; "
              f"config analysis_surface = {cfg.trace.analysis_surface}; "
              f"session resolved to {session.surface}")
        check(session.surface == image,
              f"session traces the sequence's own image surface ({image}), not "
              f"config's {cfg.trace.analysis_surface}")
        # The name, not just the index. Read it from the file: the API exposes
        # the index, and an index equal to 2 would also be true of a model that
        # still ended on `secondary`.
        import re
        seq = list(re.finditer(r"<sequence\b.*?</sequence>",
                               model.read_text(encoding="utf-8"), re.S))
        listed = re.findall(r'<surf id="([^"]*)"',
                            seq[cfg.trace.analysis_seq].group(0))
        at = listed[image] if image < len(listed) else "(out of range)"
        check(at == "prime_focus",
              f"surface {image} of sequence {cfg.trace.analysis_seq} is {at!r} "
              f"({' -> '.join(listed)})")

        # -- 3. ray-count semantics -------------------------------------------
        print("\n3. setRayDistributionCount1 semantics (README's probe)")
        probe_n = 300
        core.setRayDistributionCount1(cfg.trace.analysis_seq, probe_n)
        core.traceRays(cfg.trace.analysis_seq, 0, 0)
        pos = np.array(core.getRayPos(cfg.trace.analysis_seq, 0, 0,
                                      session.surface), copy=True)
        emitted = 0 if pos.size == 0 else pos.shape[1]
        print(f"  requested {probe_n} -> emitted {emitted}")
        check(emitted == probe_n,
              f"LITERAL ray count on sequence {cfg.trace.analysis_seq} "
              f"({emitted} for {probe_n}; a grid density would give thousands "
              f"-- and the whole ray budget would be off by that factor)")

        # -- 4. parameter writes reach the trace -------------------------------
        print("\n4. self_test: do writes reach the trace?")
        st = session.self_test(strategy)
        check(st["passed"],
              f"two heliostats give two spots, {st['separation_mm']:.1f} mm apart")

        # -- 5. the centroid is on the axis -----------------------------------
        print("\n5. centroid against the analytic prediction")
        az, el = PROBE_SUN
        session.set_sun(az, el)
        sol = strategy.solve(*PROBE, az, el, cfg.geometry)
        res = session.trace_heliostat(*PROBE, sol, cfg.trace.bulk_config)
        if not res.rays_landed:
            check(False, "no rays landed on the detector at all")
        else:
            c = res.xy_mm.mean(axis=0)
            r = np.hypot(res.xy_mm[:, 0] - c[0], res.xy_mm[:, 1] - c[1])
            r90 = float(np.percentile(r, 90))
            off = float(np.hypot(*c))
            print(f"  heliostat ({PROBE[0]/1000:+.1f}, {PROBE[1]/1000:+.1f}) m, "
                  f"sun az {az} el {el}")
            print(f"  aim ({sol.extras['aim_x_mm']:.0f}, "
                  f"{sol.extras['aim_y_mm']:.0f}, {sol.extras['aim_z_mm']:.0f}) mm, "
                  f"slant {sol.focal_dist_mm/1000:.1f} m")
            print(f"  {res.rays_landed} rays ({res.transmission:.1%}), centroid "
                  f"({c[0]:+.1f}, {c[1]:+.1f}) mm, r90 {r90:.1f} mm")
            check(off < CENTROID_TOL_MM,
                  f"centroid is on axis: |({c[0]:+.1f}, {c[1]:+.1f})| = "
                  f"{off:.1f} mm < {CENTROID_TOL_MM:.0f} mm "
                  f"(predicted (0, 0): the detector plane IS the aim plane)")
            focused_r90 = r90

            # -- 6. flat mirrors throw a much bigger spot ---------------------
            print("\n6. flat heliostats, same geometry")
            flat = get_strategy(cfg, flat=True)
            fsol = flat.solve(*PROBE, az, el, cfg.geometry)
            check(fsol.c3 == fsol.c4 == fsol.c5 == 0.0
                  and (sol.c3, sol.c4, sol.c5) != (0.0, 0.0, 0.0),
                  f"flat solve zeroes the Zernikes: focused c3/c4/c5 = "
                  f"({sol.c3:.4g}, {sol.c4:.4g}, {sol.c5:.4g}) -> "
                  f"({fsol.c3:.4g}, {fsol.c4:.4g}, {fsol.c5:.4g}); pointing "
                  f"unchanged ({fsol.rot_az_deg:.4f}, {fsol.rot_el_deg:.4f})")
            fres = session.trace_heliostat(*PROBE, fsol, cfg.trace.bulk_config)
            if not fres.rays_landed:
                check(False, "flat trace landed nothing")
            else:
                fc = fres.xy_mm.mean(axis=0)
                fr = np.hypot(fres.xy_mm[:, 0] - fc[0], fres.xy_mm[:, 1] - fc[1])
                fr90 = float(np.percentile(fr, 90))
                print(f"  {fres.rays_landed} rays, centroid "
                      f"({fc[0]:+.1f}, {fc[1]:+.1f}) mm, r90 {fr90:.1f} mm")
                check(fr90 > 3.0 * focused_r90,
                      f"flat spot is {fr90/focused_r90:.1f}x the focused one "
                      f"(r90 {fr90:.0f} vs {focused_r90:.0f} mm) -- the c3/c4/c5 "
                      f"writes reach this file's helio_surf")
                # A flat spot much larger than the 2500 mm draw radius is also
                # the direct evidence that float_ap does not clip. If this one
                # fails while the ratio above passes, look at the float_ap
                # before believing anything about spillage.
                check(fr90 > 2500.0 or float(np.max(fr)) > 2500.0,
                      f"rays land beyond the 2500 mm float_ap draw radius "
                      f"(max {float(np.max(fr)):.0f} mm) -- the detector is "
                      f"unbounded, spillage stays a post-processing step")
    finally:
        session.close()

    bad = sum(not ok for ok, _ in results)
    print(f"\n  {len(results) - bad}/{len(results)} checks passed")
    print("  " + ("PASS -- the prime-focus model traces what the Python side thinks"
                  if not bad else
                  "FAIL -- do not run a prime-focus sweep until this is green"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
