"""Validate the fetched DNI table before switching config.toml to table mode.

Run after ``python -m beamdown fetch-dni --source pvgis``:

    python -m scripts.check_dni

Checks, in order:

1. Path agreement between what ``dni.fetch`` writes and what
   ``dni.load_dni_provider`` reads.
2. Basic sanity: row/day counts, DNI range, annual total, monthly means
   (wet/dry season signal expected for this site).
3. Timezone alignment: DNI must be ~0 whenever computed solar elevation is
   negative, and the DNI-weighted centre of each day must sit close to solar
   noon. A systematic multi-hour offset here would mean the UTC->local
   conversion in ``dni._tidy`` is wrong.
4. Interpolation spot checks at the fractional-hour sample times the sweep
   will actually use.
5. How many (month, day) buckets end up with a hop count other than 24
   because the +/-3h local shift can walk a UTC row across a month boundary.

This script only reads data; it never touches config.toml, dni.py's
production code path (aside from calling it), or Quadoa.
"""

from __future__ import annotations

import dataclasses
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, ".")

from beamdown import config as C
from beamdown import dni as D
from beamdown import solar as S


def section(title: str) -> None:
    print()
    print(f"== {title} ==")


def main() -> int:
    cfg = C.load_config()
    site = cfg.site

    # ------------------------------------------------------------------
    section("1. path agreement")
    fetch_target = cfg.path(f"dni_{'pvgis'}.csv")
    load_target = cfg.path(cfg.dni.table_file)
    print(f"fetch() writes : {fetch_target}")
    print(f"loader reads   : {load_target}")
    if fetch_target != load_target:
        print("MISMATCH -- loader will not find what fetch() wrote.")
        return 1
    if not fetch_target.exists():
        print(f"{fetch_target} does not exist -- run the fetch first.")
        return 1
    print("OK -- same path.")

    # Exercise the real load_dni_provider() code path without touching
    # config.toml: build a copy of cfg with [dni] mode="table".
    table_cfg = dataclasses.replace(cfg, dni=dataclasses.replace(cfg.dni, mode="table"))
    provider = D.load_dni_provider(table_cfg)
    assert isinstance(provider, D.TableDNI), type(provider)
    print(f"loaded via load_dni_provider(): {provider.describe()}")

    frame = pd.read_csv(fetch_target)

    # ------------------------------------------------------------------
    section("2. basic sanity")
    n_rows = len(frame)
    day_pairs = frame.groupby(["month", "day"]).size()
    print(f"rows: {n_rows} (expect 8760)")
    print(f"distinct (month, day) pairs: {day_pairs.shape[0]} (expect 365)")

    dni = frame["dni_w_m2"].to_numpy(float)
    dni_clip = np.clip(dni, 0.0, None)  # PVGIS emits tiny negative noise at night
    print(f"dni range: {dni.min():.2f} .. {dni.max():.2f} W/m2")
    print(f"dni mean (all hours, incl. night): {dni_clip.mean():.2f} W/m2")

    annual_kwh_m2 = dni_clip.sum() * 1.0 / 1000.0  # hourly samples -> Wh/m2 -> kWh/m2
    print(f"annual total: {annual_kwh_m2:.1f} kWh/m2/yr (expect ~1700-2200 for this site)")
    if not (1700.0 <= annual_kwh_m2 <= 2200.0):
        print("*** OUT OF EXPECTED RANGE -- investigate before trusting this table ***")

    monthly = frame.assign(dni_clip=dni_clip).groupby("month")["dni_clip"].mean()
    print("\nmonthly mean DNI (W/m2, all hours incl. night):")
    for m in range(1, 13):
        bar = "#" * int(monthly.get(m, 0.0) / 5)
        print(f"  {m:2d}: {monthly.get(m, 0.0):6.1f}  {bar}")
    wet = monthly[[11, 12, 1, 2, 3]].mean()
    dry = monthly[[5, 6, 7, 8, 9]].mean()
    print(f"wet season (Nov-Mar) mean: {wet:.1f} W/m2")
    print(f"dry season (May-Sep) mean: {dry:.1f} W/m2")
    if dry <= wet:
        print("*** expected dry season > wet season for Mato Grosso -- signal looks flat/wrong ***")
    else:
        print(f"OK -- dry season is {dry - wet:.1f} W/m2 higher, a real seasonal signal.")

    # ------------------------------------------------------------------
    section("3. timezone alignment")
    year = 2026  # arbitrary non-leap year; table is matched on (month, day) only

    # First, a handful of representative days printed for inspection. A single
    # day's DNI-weighted centroid is a noisy estimator of solar noon -- a TMY
    # is built from real historical days, so an individual day's cloud pattern
    # (a clear morning under an afternoon storm, say) can shift its own peak by
    # an hour or two with no timezone error involved. The printout is for
    # eyeballing plausibility; the pass/fail call below uses the full year.
    display_days = [(1, 15), (3, 20), (6, 21), (7, 15), (9, 22), (12, 21)]
    for month, day in display_days:
        grp = frame[(frame.month == month) & (frame.day == day)].sort_values("hour")
        hours_d = grp["hour"].to_numpy(float)
        vals_d = np.clip(grp["dni_w_m2"].to_numpy(float), 0.0, None)
        rise, set_ = S.sunrise_sunset(site.latitude, site.longitude, site.timezone,
                                       year, month, day)
        solar_noon = 0.5 * (rise + set_)
        centroid = float(np.sum(hours_d * vals_d) / vals_d.sum()) if vals_d.sum() > 0 else float("nan")
        peak_hour = float(hours_d[np.argmax(vals_d)])
        print(f"  {month:2d}/{day:<2d}  solar noon {solar_noon:5.2f}h  "
              f"DNI peak {peak_hour:5.2f}h  DNI centroid {centroid:5.2f}h  "
              f"offset {centroid - solar_noon:+5.2f}h  (sunrise {rise:.2f} sunset {set_:.2f})")

    # Now the two checks that actually decide pass/fail, done at the table's
    # own hourly resolution (not an interpolated sub-grid -- interpolating
    # across the sunrise/sunset hour necessarily "bleeds" a fraction of the
    # neighbouring hour's value for up to about an hour either side of the
    # terminator, which is a resolution limitation, not a timezone error; that
    # is covered separately in section 4).
    #
    # (a) DNI must be exactly ~0 at every table hour whose computed elevation
    #     is unambiguously negative (more than a degree below the horizon, to
    #     stay clear of refraction/twilight edge cases).
    # (b) Averaged over the whole year (so day-to-day weather noise cancels),
    #     the DNI-weighted centre of each day should sit close to solar noon.
    #     A real ~3h systematic offset -- the size of this site's own UTC
    #     offset -- would mean the conversion in `_tidy` was missed or
    #     double-applied. Small (well under an hour) systematic skew is
    #     plausible on its own merits: this is a wet-season convective climate
    #     where afternoon cloud buildup is a real, physical, non-bug effect.
    night_violations = []
    offsets = []
    for (m, d), grp in frame.groupby(["month", "day"]):
        if (m, d) == (2, 29):
            continue  # see section 5 -- not a real calendar day for this TMY
        hours_d = grp["hour"].to_numpy(float)
        vals_d = np.clip(grp["dni_w_m2"].to_numpy(float), 0.0, None)
        _, elev = S.sun_position(site.latitude, site.longitude, site.timezone, year, m, d, hours_d)

        clearly_night = elev < -1.0
        bad = vals_d[clearly_night]
        if (bad > 1.0).any():
            night_violations.append((m, d, float(bad.max())))

        rise, set_ = S.sunrise_sunset(site.latitude, site.longitude, site.timezone, year, m, d)
        solar_noon = 0.5 * (rise + set_)
        if vals_d.sum() > 0:
            offsets.append(float(np.sum(hours_d * vals_d) / vals_d.sum()) - solar_noon)

    offsets = np.array(offsets)
    print(f"\ndays with DNI > 1 W/m2 while elevation < -1deg: {len(night_violations)} / {len(offsets)}")
    if night_violations:
        print(f"  worst: {night_violations[:5]}")

    print(f"annual mean offset (DNI centroid - solar noon), n={len(offsets)} days: "
          f"{offsets.mean():+.3f}h  (std {offsets.std():.3f}h, mean |offset| {np.abs(offsets).mean():.3f}h)")

    tz_pass = len(night_violations) == 0 and abs(offsets.mean()) < 1.5
    if abs(offsets.mean()) > 2.0:
        print("*** offset looks like a systematic ~3h shift -- UTC->local conversion is likely "
              "wrong (missed or double-applied) ***")
    print(f"\nTIMEZONE ALIGNMENT: {'PASS' if tz_pass else 'FAIL'}")

    print("\nPVGIS 'Irradiance Time Offset (h)' (0.0817h ~ 5 min) is not applied by "
          "_fetch_pvgis_tmy -- it only reads time(UTC) as given. At the ~1h table "
          "resolution and the fractional-hour sampling checked below, a 5-minute "
          "offset is well under the interpolation grid spacing and immaterial; noted "
          "here rather than silently ignored.")

    # ------------------------------------------------------------------
    section("4. interpolation spot checks")
    # Pick a day with a clean sunrise/sunset and show fractional-hour behaviour,
    # including the non-clock sample times the sweep will actually use.
    month, day = 6, 21
    date = pd.Timestamp(year, month, day).date()
    rise, set_ = S.sunrise_sunset(site.latitude, site.longitude, site.timezone, year, month, day)
    grp = frame[(frame.month == month) & (frame.day == day)].sort_values("hour")
    print(f"date {date}, sunrise {rise:.2f}h, sunset {set_:.2f}h")
    print("table hours bracketing sunrise/sunset (the raw hourly values being interpolated):")
    for h in (6.0, 7.0, 8.0, 17.0, 18.0, 19.0):
        row = grp[grp.hour == h]
        if len(row):
            print(f"  {h:5.1f}h -> {row['dni_w_m2'].iloc[0]:7.2f} W/m2")
    sample_hours = sorted({rise - 0.5, rise, rise + 0.1, 6 + 42 / 60, 7 + 41 / 60, 8 + 40 / 60,
                            12.0, set_ - 0.1, set_, set_ + 0.5})
    print("interpolated queries (sweep will use non-clock times like these):")
    for h in sample_hours:
        hh = int(h)
        mm = int(round((h - hh) * 60))
        print(f"  {hh:02d}:{mm:02d} ({h:6.3f}h) -> {provider.dni(date, h):7.2f} W/m2")
    print("note: a query inside the [06:00, 07:00) table bin can be nonzero even a little "
          "before the computed instantaneous sunrise -- the table's 07:00 sample already "
          "reflects a partly-lit hour, and linear interpolation spreads that value back "
          "across the bin. That is a real limitation of hourly-resolution data, not a bug; "
          "it is bounded to at most one table-hour width either side of sunrise/sunset.")

    # ------------------------------------------------------------------
    section("5. month-boundary hop count from the local shift")
    bad = day_pairs[day_pairs != 24]
    print(f"(month, day) buckets with != 24 entries: {len(bad)}")
    if len(bad):
        print(bad.to_string())
        print("\nSee dni.py docstring note below for why -- PVGIS TMY splices twelve "
              "different representative years together (one per calendar month); the "
              "+/-3h local shift at the Feb/Mar splice can walk into a Feb 29 that only "
              "exists in the March donor year's calendar, if that donor year happens to "
              "be a leap year. This produces a near-empty (2, 29) bucket and a short-by-3 "
              "(2, 28) bucket. It does not affect energy: the missing hours fall at local "
              "21:00-23:59, well after sunset here, and no sweep date lands on Feb 29.")
    else:
        print("none -- every calendar day has exactly 24 local hours.")

    return 0 if tz_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
