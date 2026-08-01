"""Annual angle-of-incidence statistics for the axicon field. NO licence needed.

AOI here is the heliostat mirror's angle of incidence: half the angle between
the direction to the sun and the direction to the aim point (the mirror normal
bisects the two, ``beamdown.secondary.mirror.heliostat_orientation``).  For the
axicon the aim point is SUN-INDEPENDENT -- each heliostat aims at the receiver's
mirror image in the cone flank, a fixed point on the far side of the axis
(``receiver_correction``: radial ``x_r = -drop*sin(2a)``, height
``y_r = +drop*cos(2a)`` above the tip).  So the whole year's AOI field is pure
geometry: one fixed unit vector per heliostat against one sun vector per
timestep.

The time grid is the sweep's own (``beamdown.solar.build_time_grid``, sunrise
+margin to sunset-margin, cfg's hour_step) over all 365 days of 2026, so the
averages sample the year exactly the way the energy runs do.  Three annual
means are reported:

  unweighted   trapezoid over daylight time, all heliostats equal
  DNI          additionally weighted by the configured DNI provider
  DNI*cos(AOI) weighted by what each mirror actually collects -- the mean AOI
               "of the collected energy" (no occlusion/spillage, stated).

It also reports where each heliostat's aim ray crosses the optical axis
(projecting past the cone), because that is the number a prime-focus F1 must
reproduce for comparable blocking/shadowing -- and then scans the traced
declination dates for the instant whose field-mean AOI sits closest to the
annual mean, which is the instant to build a figure model at.

Usage::

    python scripts/aoi_stats.py
    python scripts/aoi_stats.py --pf-height-mm 36000   # also show prime-focus AOI
"""

from __future__ import annotations

import argparse
import datetime as _dt
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import numpy as np  # noqa: E402

# The 7 distinct declinations every sweep traces (see run_prime_focus.sh).
TRACED_DATES = ["2026-12-21", "2026-01-21", "2026-02-20", "2026-03-20",
                "2026-04-21", "2026-05-21", "2026-06-21"]


def sun_vectors(az_deg: np.ndarray, el_deg: np.ndarray) -> np.ndarray:
    """(T, 3) unit vectors toward the sun, mirror.py's exact convention."""
    az = np.deg2rad(az_deg)
    el = np.deg2rad(el_deg)
    return np.stack([
        np.cos(el) * np.cos(np.pi / 2 - az),
        np.cos(el) * np.sin(np.pi / 2 - az),
        np.sin(el),
    ], axis=1)


def day_weights(hours: np.ndarray) -> np.ndarray:
    """Trapezoid weights (in hours) for one day's uniformly spaced samples."""
    n = len(hours)
    if n == 1:
        return np.array([1.0])
    h = (hours[-1] - hours[0]) / (n - 1)
    w = np.full(n, h)
    w[0] = w[-1] = h / 2.0
    return w


def annual_means(aoi: np.ndarray, w_time: np.ndarray, dni: np.ndarray) -> dict:
    """The three weighted annual mean AOIs, over (heliostat, timestep)."""
    per_t = aoi.mean(axis=0)                       # plain field mean per instant
    out = {"unweighted": float(np.average(per_t, weights=w_time)),
           "dni": float(np.average(per_t, weights=w_time * dni))}
    w_full = (w_time * dni)[None, :] * np.cos(np.deg2rad(aoi))
    out["dni_cos"] = float(np.average(aoi, weights=w_full))
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--pf-height-mm", type=float, default=None,
                   help="also report AOI stats for a prime-focus field aimed at "
                        "(0,0,THIS) -- e.g. 36000 for 9 m above the axicon tip")
    p.add_argument("--config", type=Path, default=None)
    a = p.parse_args(argv)

    from beamdown import field as F
    from beamdown.config import load_config
    from beamdown.dni import provider_for
    from beamdown.secondary.axicon import receiver_correction
    from beamdown.solar import build_time_grid

    cfg = load_config(a.config)
    g = cfg.geometry
    if cfg.optics.secondary != "axicon":
        print(f"config layout is {cfg.optics.secondary!r}; this script reads the "
              f"axicon geometry keys regardless -- they are what it is about.")

    full = F.load_field(cfg)
    x = np.asarray(full.x_mm, dtype=float)
    y = np.asarray(full.y_mm, dtype=float)
    R = np.hypot(x, y)

    # Fixed aim point per heliostat: the receiver's image in the cone flank.
    drop = g.secondary_height_mm - g.receiver_height_mm
    x_r, y_r, _, _ = receiver_correction(R, g.secondary_height_mm, drop,
                                         g.axicon_angle_deg)
    aim = np.stack([x / R * x_r, y / R * x_r,
                    np.full_like(R, g.secondary_height_mm + y_r)], axis=1)
    mirror = np.stack([x, y, np.zeros_like(x)], axis=1)
    to_aim = aim - mirror
    to_aim /= np.linalg.norm(to_aim, axis=1, keepdims=True)

    # Where the aim ray crosses the axis (radius R down to 0 through x_r < 0).
    z_axis = aim[:, 2] * R / (R - x_r)
    far = int(np.argmax(R))
    tip = g.secondary_height_mm
    print(f"aim-ray axis crossing, z above the axicon tip ({tip:,.0f} mm):")
    print(f"  aim image point: radial {float(np.unique(np.round(x_r, 6))[0]):,.1f} mm, "
          f"z {tip + float(np.unique(np.round(y_r, 6))[0]):,.1f} mm (same for all)")
    print(f"  axis crossing:  min {z_axis.min() - tip:+,.0f}   "
          f"mean {z_axis.mean() - tip:+,.0f}   max {z_axis.max() - tip:+,.0f} mm")
    print(f"  farthest heliostat (R = {R[far] / 1000:.1f} m): crosses at "
          f"z = {z_axis[far]:,.0f} mm = tip {z_axis[far] - tip:+,.0f} mm")
    print()

    # The year, sampled the way the sweeps sample it.
    year_dates = [_dt.date(2026, 1, 1) + _dt.timedelta(days=i) for i in range(365)]
    steps = build_time_grid(cfg, year_dates)
    az = np.array([s.solar_az_deg for s in steps])
    el = np.array([s.solar_el_deg for s in steps])
    sun = sun_vectors(az, el)

    w_time = np.empty(len(steps))
    i = 0
    while i < len(steps):
        j = i
        while j < len(steps) and steps[j].date == steps[i].date:
            j += 1
        w_time[i:j] = day_weights(np.array([s.hour for s in steps[i:j]]))
        i = j

    provider = provider_for(cfg)
    dni = np.array([provider.dni(s.date, s.hour) for s in steps])
    print(f"time grid: {len(steps)} timesteps over 365 days "
          f"(cfg hour_step {cfg.sweep.hour_step:g} h); DNI {provider.describe()}")
    print()

    def report(label: str, target: np.ndarray) -> dict:
        aoi = 0.5 * np.degrees(np.arccos(np.clip(target @ sun.T, -1.0, 1.0)))
        m = annual_means(aoi, w_time, dni)
        print(f"{label}: annual mean AOI over {len(R)} heliostats")
        print(f"  unweighted (daylight time)   {m['unweighted']:7.3f} deg")
        print(f"  DNI-weighted                 {m['dni']:7.3f} deg")
        print(f"  DNI*cos(AOI)-weighted        {m['dni_cos']:7.3f} deg   "
              f"(mean AOI of the collected energy; no occlusion/spillage)")
        print()
        return {"aoi": aoi, **m}

    ax = report("AXICON", to_aim)

    if a.pf_height_mm is not None:
        f1 = np.array([0.0, 0.0, float(a.pf_height_mm)])
        to_f1 = f1[None, :] - mirror
        to_f1 /= np.linalg.norm(to_f1, axis=1, keepdims=True)
        report(f"PRIME FOCUS @ z={a.pf_height_mm:,.0f} mm", to_f1)

    # Which traced-date instant looks most like the annual average?
    target = ax["dni"]
    print(f"instants closest to the axicon DNI-weighted mean ({target:.3f} deg),")
    print(f"per traced declination date (field-mean AOI at each grid instant):")
    dates_set = {_dt.date.fromisoformat(d) for d in TRACED_DATES}
    per_t = ax["aoi"].mean(axis=0)
    best_rows = []
    for d in sorted(dates_set):
        idx = [k for k, s in enumerate(steps) if s.date == d]
        k = min(idx, key=lambda k: abs(per_t[k] - target))
        s = steps[k]
        best_rows.append((abs(per_t[k] - target), d, s.hour, per_t[k],
                          s.solar_az_deg, s.solar_el_deg))
    for delta, d, hour, val, saz, sel in sorted(best_rows):
        hh, mm = int(hour), int(round((hour % 1) * 60))
        print(f"  {d}  {hh:02d}:{mm:02d}  field-mean AOI {val:7.3f} deg "
              f"(|delta| {delta:5.3f})   sun az {saz:7.2f}  el {sel:6.2f}")
    kbest = min(range(len(steps)), key=lambda k: abs(per_t[k] - target))
    s = steps[kbest]
    print(f"  (unrestricted best over the whole year: {s.date} {s.hour:.3f} h, "
          f"field-mean {per_t[kbest]:.3f} deg)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
