"""Does tracing the occluders reproduce the analytic shading and blocking?

The scalar model says a heliostat delivers ``rays_landed x eta_shade x
eta_block``. The occluder model removes those rays in the trace instead. If both
are right, the traced ray count divided by the unoccluded one must equal
``eta_shade x eta_block`` -- for every heliostat, not just on average.

Run against a finished scalar sweep (full5) so the unoccluded counts and the
analytic efficiencies both come from the same recorded run:

    python scripts/verify_occluder_trace.py --run analysis_output/full5
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="analysis_output/full5")
    ap.add_argument("--rays", type=int, default=120000)
    ap.add_argument("--timesteps", nargs="+",
                    default=["20260320_1800", "20260922_0800", "20260320_1300"])
    ap.add_argument("--per-step", type=int, default=6,
                    help="heliostats per timestep, spanning the loss range")
    args = ap.parse_args(argv)

    from beamdown import field as F
    from beamdown import occluder_slots as OS
    from beamdown import shading as S
    from beamdown.config import load_config
    from beamdown.secondary import get_strategy
    from beamdown.session import QuadoaSession
    from beamdown.store import RunStore

    cfg = load_config(None)
    model = cfg.repo_root / "models" / "heliostat_field_occluders.optx"
    if not model.exists():
        print(f"FAIL: {model} does not exist yet")
        return 1
    object.__setattr__(cfg.trace, "model_file", str(model))
    object.__setattr__(cfg.trace, "rays_per_heliostat", args.rays)
    traced_secondary = "ax0_x" in model.read_text(encoding="utf-8")
    print(f"  model: {model.name}   secondary shadow "
          f"{'traced' if traced_secondary else 'not traced'}")

    object.__setattr__(cfg.storage, "root", args.run)
    summary = RunStore(cfg.output_root, cfg=cfg, mode="r").summary()

    fld = F.load_field(cfg)
    strategy = get_strategy(cfg)
    neighbours = F.neighbour_pairs(
        fld, S.search_radius_for(float(summary.solar_el_deg.min()),
                                 cfg.field.mirror_height_mm, cfg.field.mirror_width_mm))

    rows = []
    session = QuadoaSession(cfg)
    try:
        session.set_global_geometry()
        for key in args.timesteps:
            step = summary[summary.timestep == key]
            if not len(step):
                print(f"  skipping {key}: not in {args.run}")
                continue
            az = float(step.solar_az_deg.iloc[0])
            el = float(step.solar_el_deg.iloc[0])
            sols = [strategy.solve(float(fld.x_mm[i]), float(fld.y_mm[i]), az, el,
                                   cfg.geometry) for i in range(len(fld))]
            geoms, aims = S.build_geometries(fld, sols, cfg)
            to_sun = S.sun_vector(az, el)
            cone = S.secondary_body(cfg) if traced_secondary else None
            plans = OS.plan_field(geoms, aims, fld.ids, neighbours, to_sun,
                                  body=cone)

            # Compare against NEIGHBOUR shading only. The summary's eta_shade is
            # the union of neighbours and the secondary, and the secondary is
            # deliberately not in the traced path -- it shades the mirror from
            # sunlight, but the sequence meets it on the way out. Comparing
            # against the union would score the axicon's shadow as a trace
            # failure: one heliostat here reads eta_shade 0.0000 purely because
            # the axicon covers it.
            eta_nb, eta_bl, _ = S.shading_blocking(geoms, aims, az, el, neighbours,
                                                   secondary=None)
            # The union, not the product: a shaded patch sends no beam, so it
            # cannot also be blocked, and multiplying removes it twice. The
            # secondary joins the same union once the model traces it, and must
            # stay out of the prediction while it is still a scalar -- comparing
            # against the wrong one scores a correct trace as a failure.
            eta_joint = S.occlusion_efficiency(geoms, aims, az, el, neighbours,
                                               secondary=cone)

            # Span the loss range rather than sampling at random: the cases that
            # matter are the heavily occluded ones, and a random draw at a
            # high-sun timestep would be all eta = 1.
            by_id = step.set_index("heliostat_id")
            nb_eff = pd.Series(eta_joint, index=list(fld.ids))
            loss = (1.0 - nb_eff.reindex(by_id.index)).sort_values()
            pick = loss.index[np.linspace(0, len(loss) - 1, args.per_step).astype(int)]

            session.set_sun(az, el)
            for hid in pick:
                i = int(np.flatnonzero(fld.ids == hid)[0])
                r = by_id.loc[hid]
                session.set_occluders(plans[i], OS.N_SHADE, OS.N_BLOCK)
                t0 = time.perf_counter()
                res = session.trace_heliostat(float(fld.x_mm[i]), float(fld.y_mm[i]),
                                              sols[i], cfg.trace.bulk_config)
                dt = time.perf_counter() - t0
                predicted = float(eta_joint[i])
                observed = res.rays_landed / float(r.rays_landed)
                rows.append(dict(timestep=key, heliostat_id=int(hid), el=el,
                                 eta_shade=float(eta_nb[i]),
                                 eta_block=float(eta_bl[i]),
                                 predicted=predicted, observed=observed,
                                 diff=observed - predicted,
                                 n_shade=len(plans[i].shading),
                                 n_block=len(plans[i].blocking),
                                 seconds=dt))
                print(f"  {key} h{int(hid):<4d} el {el:5.1f}  "
                      f"slots {len(plans[i].shading)}s/{len(plans[i].blocking)}b  "
                      f"predicted {predicted:.4f}  traced {observed:.4f}  "
                      f"diff {observed - predicted:+.4f}   {dt:.2f}s", flush=True)
    finally:
        session.close()

    t = pd.DataFrame(rows)
    if t.empty:
        print("FAIL: nothing traced")
        return 1

    # Monte-Carlo noise on a ratio of two independent counts of ~45,000.
    noise = float(np.sqrt(2.0 / 45000.0))
    print("\n  " + "-" * 66)
    print(f"  {len(t)} heliostats traced")
    print(f"  mean  |traced - predicted| = {t['diff'].abs().mean():.5f}")
    print(f"  worst |traced - predicted| = {t['diff'].abs().max():.5f}")
    print(f"  bias  (traced - predicted) = {t['diff'].mean():+.5f}")
    print(f"  Monte-Carlo noise floor    = {noise:.5f}")
    print(f"  mean seconds per trace     = {t.seconds.mean():.2f}")

    ok = bool(t["diff"].abs().max() < 6 * noise and abs(t["diff"].mean()) < 2 * noise)
    print("\n  " + ("PASS -- the traced occlusion matches the analytic scalars"
                    if ok else
                    "FAIL -- the trace and the scalars disagree by more than noise"))
    t.to_csv("analysis_output/verify_occluder_trace.csv", index=False)
    print("  wrote analysis_output/verify_occluder_trace.csv")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
