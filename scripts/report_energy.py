"""Annual and per-day collected energy for a finished sweep.

Everything here is post-processing over a stored summary: no license, seconds to
run. DNI is a pure multiplier on the trace (which is normalised to exactly
1000 W/m^2), so the DNI model can be changed and this re-run without re-tracing.

    python scripts/report_energy.py --run analysis_output/full8
    python scripts/report_energy.py --run analysis_output/full8 --dni-mode table

The per-day number is computed two independent ways -- a trapezoid over the
traced samples, and the interpolated efficiency surface walked over 24 hourly
steps -- because the residual between them measures whether the time grid is
sampled densely enough. A few tenths of a percent means it is; several percent
would mean the morning and evening wings are being cut by straight lines.
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# The five pairs of configured dates that land within 0.2 deg of the same solar
# declination. Same sun geometry, different calendar date and different DNI, so
# their optical efficiency has to agree -- a free end-to-end check on pointing,
# occluder planning and the trace itself.
DECLINATION_PAIRS = [
    ("2026-01-21", "2026-11-21"),
    ("2026-02-20", "2026-10-21"),
    ("2026-03-20", "2026-09-22"),
    ("2026-04-21", "2026-08-21"),
    ("2026-05-21", "2026-07-21"),
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="analysis_output/full8")
    ap.add_argument("--dni-mode", default=None,
                    help="override [dni] mode: constant | table | monthly")
    args = ap.parse_args(argv)

    from beamdown import dni as D
    from beamdown import energy as E
    from beamdown.config import load_config
    from beamdown.store import RunStore

    cfg = load_config(None)
    object.__setattr__(cfg.storage, "root", args.run)
    summary = RunStore(cfg.output_root, cfg=cfg, mode="r").summary()
    prov = (D.provider_for(cfg, args.dni_mode) if args.dni_mode
            else D.load_dni_provider(cfg))

    dupes = int(summary.duplicated(["timestep", "heliostat_id"]).sum())
    print(f"{args.run}")
    print(f"  rows {len(summary):,}   timesteps {summary.timestep.nunique()}   "
          f"heliostats {summary.heliostat_id.nunique()}   duplicates {dupes}")
    print(f"  sun elevation {summary.solar_el_deg.min():.2f} to "
          f"{summary.solar_el_deg.max():.2f} deg")
    print(f"  DNI model: {prov.describe()}")

    res = E.annual_energy(summary, cfg, prov)
    print(f"\n  annual energy             {res['annual_energy_mwh']:>12,.1f} MWh")
    print(f"  annual DNI                {res['annual_dni_kwh_m2']:>12,.1f} kWh/m2")
    print(f"  mirror area               {res['mirror_area_m2']:>12,.0f} m2")
    print(f"  annual optical efficiency {res['annual_optical_efficiency']:>12.4f}")
    print(f"  extrapolated fraction     {res['extrapolated_fraction']:>12.4f}")
    print(f"  traced timesteps          {res['traced_timesteps']:>12d}")
    print(f"  distinct declinations     {res['traced_declinations']:>12d}")

    fit = E.fit_annual_sine(res["daily"])
    print(f"\n  fitted sinusoid over the {len(res['daily'])} modelled days:")
    print(f"    mean       {fit['mean'] / 1000.0:>9,.2f} MWh/day")
    print(f"    amplitude  {fit['amplitude'] / 1000.0:>9,.2f} MWh/day")
    print(f"    peaks on   day-of-year {fit['phase_day_of_year']:.0f}")
    print(f"    R^2        {fit['r_squared']:>9.4f}")

    print("\n  per-traced-day energy, two independent routes:")
    daily = res["daily"].set_index("date")
    # Traced dates come from the summary, not from cfg.sweep.dates: the config
    # lists the 12 dates we now want, but an older run may hold fewer, and every
    # configured date appears in `daily` regardless because the model covers all
    # 365 days. Only a date actually in the summary has samples to integrate.
    traced = sorted({dt.datetime.strptime(k.split("_")[0], "%Y%m%d").date()
                     for k in summary.timestep.unique()})
    worst = 0.0
    for d in traced:
        if d not in daily.index:
            continue
        direct = E.traced_day_energy(summary, cfg, prov, d)["energy_kwh"] / 1000.0
        interp = float(daily.loc[d, "energy_kwh"]) / 1000.0
        resid = 100.0 * (direct - interp) / interp if interp else float("nan")
        worst = max(worst, abs(resid))
        print(f"    {d}   trapezoid {direct:7.2f} MWh   "
              f"model {interp:7.2f} MWh   residual {resid:+6.2f}%")
    print(f"    worst |residual| = {worst:.2f}%")

    eta = (E.optical_efficiency(summary, cfg)
           .groupby("date")["eta_optical"].mean())
    print("\n  near-duplicate declination pairs -- optical efficiency must agree:")
    for a, b in DECLINATION_PAIRS:
        da, db = dt.date.fromisoformat(a), dt.date.fromisoformat(b)
        if da in eta.index and db in eta.index:
            print(f"    {a} {eta[da]:.4f}   vs   {b} {eta[db]:.4f}   "
                  f"diff {eta[da] - eta[db]:+.4f}")

    out = Path(args.run) / "daily_energy.csv"
    res["daily"].to_csv(out, index=False)
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
