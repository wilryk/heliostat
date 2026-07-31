"""End-to-end smoke test: a small real sweep, then read the results back.

Uses the 25 downselected heliostats and 2 timesteps so it finishes in about a
minute, but exercises the entire path: solve -> shade -> trace -> bin ->
quantise -> store -> summary -> flux -> metrics -> annual energy.

    python tests/smoke_sweep.py [--workers N] [--rays N]
"""

from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--rays", type=int, default=None)
    ap.add_argument("--steps", type=int, default=2, help="number of timesteps")
    ap.add_argument("--keep", action="store_true", help="keep the output directory")
    args = ap.parse_args()

    from beamdown.config import load_config
    from beamdown import field as F, dni as D, energy, metrics, sweep as S
    from beamdown.store import RunStore

    cfg = load_config()
    out = REPO / "_smoke_output"
    if out.exists():
        shutil.rmtree(out)
    object.__setattr__(cfg.storage, "root", str(out))

    if args.rays:
        # Workers reload config from disk, so an in-memory override would not
        # reach them; the ray budget has to come from config.toml.
        print("note: --rays is ignored here; set trace.rays_per_heliostat in config.toml")

    fld = F.load_field(cfg)
    idx, prov = F.load_or_build_downselect(cfg, fld)
    print(f"heliostats: {idx.size} ({prov})")

    from beamdown import solar

    dates = cfg.sweep.dates[:1]          # one date
    steps_all = solar.build_time_grid(cfg, dates)
    mid = len(steps_all) // 2
    steps = steps_all[mid - 1:mid + args.steps - 1]
    print(f"date {dates[0]} has {len(steps_all)} timesteps; using {len(steps)}: "
          + ", ".join(s.label for s in steps))
    print(f"rays/heliostat: {cfg.trace.rays_per_heliostat:,}  chunks {cfg.trace.chunk_sizes}")

    t0 = time.time()
    store = S.run_sweep(
        cfg,
        heliostat_indices=idx,
        steps=steps,
        workers=args.workers or cfg.trace.n_workers,
        resume=False,
    )
    elapsed = time.time() - t0
    keys_expected = len(steps)

    # ---- read back -------------------------------------------------------
    print("\n--- reading results back ---")
    reader = RunStore(cfg.output_root, cfg=cfg, mode="r")
    keys = reader.timestep_keys()
    print(f"timesteps written: {len(keys)}")

    summary = reader.summary()
    print(f"summary rows: {len(summary)}  columns: {len(summary.columns)}")

    key = keys[len(keys) // 2]
    counts = reader.read_counts(key)
    print(f"\n{key}: counts array {counts.shape}, total rays {int(np.asarray(counts).sum()):,}")

    rays = reader.read_rays(key, heliostat_id=int(idx[0]))
    print(f"  raw rays for heliostat {idx[0]}: {rays.shape}, "
          f"x [{rays[:,0].min():.1f}, {rays[:,0].max():.1f}] mm")

    # Linearity: the field map must equal the sum of per-heliostat maps.
    field_map = reader.field_flux(key)
    per_helio_sum = np.asarray(counts).sum(axis=0)
    from beamdown.store import scale_factor
    expected = per_helio_sum * scale_factor(cfg, cfg.trace.rays_per_heliostat) / cfg.receiver.bin_area_m2
    lin_err = float(np.abs(field_map - expected).max())
    print(f"  linearity check (field == sum of parts): max err {lin_err:.3e} W/m2")

    # Energy conservation between the summary table and the ray counts.
    step_rows = summary[summary.timestep == key]
    power_from_rows = step_rows.power_w.sum()
    landed = step_rows.rays_landed.sum()
    from beamdown.store import occlusion_weight_columns
    eff = np.ones(len(step_rows))
    for col in occlusion_weight_columns(reader.manifest, step_rows.columns):
        eff = eff * step_rows[col].to_numpy(float)
    power_from_counts = float(
        (np.asarray(counts).sum(axis=(1, 2)) * eff).sum()
        * scale_factor(cfg, cfg.trace.rays_per_heliostat)
    )
    print(f"  power from summary rows : {power_from_rows:12.1f} W")
    print(f"  power from stored counts: {power_from_counts:12.1f} W")
    print(f"  relative difference     : {abs(power_from_rows-power_from_counts)/max(power_from_rows,1e-9):.2e}")

    print(f"\n  rays landed this step: {landed:,} of "
          f"{cfg.trace.rays_per_heliostat * len(idx):,} emitted "
          f"({landed/(cfg.trace.rays_per_heliostat*len(idx)):.1%})")
    print(f"  shading eta : {step_rows.eta_shade.min():.4f} .. {step_rows.eta_shade.max():.4f}")
    print(f"  blocking eta: {step_rows.eta_block.min():.4f} .. {step_rows.eta_block.max():.4f}")
    print(f"  peak flux   : {field_map.max()/1000:.2f} kW/m2")

    # ---- ranking and annual energy --------------------------------------
    print("\n--- worst heliostats by delivered power ---")
    ranked = metrics.rank_heliostats(summary, by="power_w")
    print(ranked.head(5)[["rank", "heliostat_id", "x_m", "y_m", "radius_m",
                          "power_w_sum", "mean_transmission"]].to_string(index=False))

    print("\n--- annual energy (from these timesteps only; indicative) ---")
    provider = D.load_dni_provider(cfg)
    try:
        result = energy.annual_energy(summary, cfg, provider, n_heliostats=len(idx))
        print(f"  DNI model            : {provider.describe()}")
        print(f"  traced timesteps     : {result['traced_timesteps']}")
        print(f"  distinct declinations: {result['traced_declinations']}")
        print(f"  mirror area          : {result['mirror_area_m2']:.1f} m2")
        print(f"  annual energy        : {result['annual_energy_mwh']:.2f} MWh")
        print(f"  annual optical eff.  : {result['annual_optical_efficiency']:.4f}")
        print(f"  extrapolated fraction: {result['extrapolated_fraction']:.1%}"
              "   <- high is expected with one date")
    except Exception as exc:
        print(f"  annual energy needs >=3 declinations; got: {type(exc).__name__}: {exc}")

    n_traces = len(idx) * len(keys)
    print(f"\nsweep: {n_traces} traces in {elapsed:.1f}s "
          f"({elapsed/max(1,n_traces)*1000:.0f} ms/trace)")
    full = 645 * 109 * (elapsed / max(1, n_traces))
    print(f"extrapolated full sweep (645 x 109 timesteps): {full/3600:.1f} h")

    ok = lin_err < 1e-6 and len(summary) == len(idx) * len(keys)
    print(f"\n{'PASS' if ok else 'FAIL'}")

    if not args.keep:
        shutil.rmtree(out, ignore_errors=True)
    else:
        print(f"output kept at {out}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
