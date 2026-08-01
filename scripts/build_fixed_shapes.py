"""Per-heliostat FIXED mirror figures for the fixed-shape scenarios. NO licence.

The focused sweeps re-solve c3/c4/c5 at every instant -- an idealised mirror
that re-figures itself as the sun moves.  A real mirror is ground once.  This
script computes, for every heliostat, the three "ground once" candidates:

  spherical   correct radius of curvature for the heliostat's slant range
              (R = 2 * focal distance), zero astigmatism.  Analytic, no
              averaging: ``heliostat_shape(rot, R, R)`` -> c3 = c5 = 0.
  mean_cos    the year-round average of the instantaneous solve, weighted by
              trapezoid daylight time x cos(AOI) -- so the large-AOI hours,
              whose astigmatism is hugely different AND whose collected power
              is lowest, do not dominate the average.
  median      the component-wise median of the instantaneous solve over the
              same grid.  Contender to mean_cos; a 4-day test run of each
              decides which represents the fixed astigmatic mirror.

Everything physical is BORROWED, never reimplemented: the instantaneous
coefficients come from the configured strategy's own ``solve`` (so they are
exactly what a sweep would write, Quadoa conventions included), the spherical
figure from ``mirror.heliostat_shape`` + ``mirror.to_quadoa_zernike`` in the
``rad_s == rad_t`` limit, and the time grid from ``solar.build_time_grid``
over all 365 days of 2026 (the sweeps' own sampling).

The POINTING is not frozen -- heliostats track in every scenario.  Only the
c3/c4/c5 shape is.

Output: ``data/fixed_shapes_<tag>_<mode>.csv`` with columns
``heliostat,x_mm,y_mm,c3,c4,c5`` and ``#`` metadata lines (read back with
``pandas.read_csv(comment="#")``).  Consumed by ``beamdown sweep
--fixed-shapes <csv>`` once that option lands.

Usage::

    python scripts/build_fixed_shapes.py --secondary prime_focus --focus-height-mm 36000
    python scripts/build_fixed_shapes.py                      # config's layout (axicon)
"""

from __future__ import annotations

import argparse
import datetime as _dt
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import numpy as np  # noqa: E402

MODES = ("spherical", "mean_cos", "median")


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    """Median of ``values`` under ``weights`` (interpolated at 0.5)."""
    order = np.argsort(values)
    v, w = values[order], weights[order]
    cum = np.cumsum(w) - 0.5 * w
    return float(np.interp(0.5 * w.sum(), cum, v))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--secondary", default=None,
                   choices=["axicon", "prime_focus", "cassegrain"],
                   help="override config's layout (with --focus-height-mm)")
    p.add_argument("--focus-height-mm", type=float, default=None)
    p.add_argument("--day-stride", type=int, default=1,
                   help="sample every Nth day of 2026 (1 = all 365)")
    p.add_argument("--out-dir", type=Path, default=REPO / "data")
    p.add_argument("--config", type=Path, default=None)
    a = p.parse_args(argv)

    from beamdown import field as F
    from beamdown.config import load_config, validate_layout
    from beamdown.secondary import get_strategy
    from beamdown.secondary.mirror import heliostat_shape, to_quadoa_zernike
    from beamdown.solar import build_time_grid

    cfg = load_config(a.config)
    if a.secondary is not None:
        object.__setattr__(cfg.optics, "secondary", a.secondary)
        object.__setattr__(cfg.optics, "n_mirrors",
                           1 if a.secondary == "prime_focus" else 2)
        if a.focus_height_mm is not None:
            object.__setattr__(cfg.geometry, "focus_height_mm",
                               float(a.focus_height_mm))
        validate_layout(cfg)
    strategy = get_strategy(cfg)

    if cfg.optics.secondary in ("prime_focus", "cassegrain"):
        tag = (f"{'pf' if cfg.optics.secondary == 'prime_focus' else 'cass'}"
               f"{cfg.geometry.focus_height_mm:.0f}")
    else:
        tag = "axicon"

    full = F.load_field(cfg)
    x = np.asarray(full.x_mm, dtype=float)
    y = np.asarray(full.y_mm, dtype=float)
    n = len(x)

    year = [_dt.date(2026, 1, 1) + _dt.timedelta(days=i)
            for i in range(0, 365, a.day_stride)]
    steps = build_time_grid(cfg, year)
    t_count = len(steps)

    # Trapezoid daylight-time weight per timestep, day by day.
    w_time = np.empty(t_count)
    i = 0
    while i < t_count:
        j = i
        while j < t_count and steps[j].date == steps[i].date:
            j += 1
        hours = np.array([s.hour for s in steps[i:j]])
        if len(hours) == 1:
            w_time[i:j] = 1.0
        else:
            h = (hours[-1] - hours[0]) / (len(hours) - 1)
            w_time[i:j] = h
            w_time[i] = w_time[j - 1] = h / 2.0
        i = j

    print(f"layout {strategy.describe()!r}, tag {tag!r}: {n} heliostats x "
          f"{t_count} timesteps ({len(year)} days, stride {a.day_stride})")

    c3 = np.empty((n, t_count))
    c4 = np.empty((n, t_count))
    c5 = np.empty((n, t_count))
    aoi = np.empty((n, t_count))
    t0 = time.monotonic()
    for k, s in enumerate(steps):
        for i in range(n):
            sol = strategy.solve(x[i], y[i], s.solar_az_deg, s.solar_el_deg,
                                 cfg.geometry)
            c3[i, k], c4[i, k], c5[i, k] = sol.c3, sol.c4, sol.c5
            aoi[i, k] = sol.aoi_deg
        if (k + 1) % 500 == 0 or k + 1 == t_count:
            dt = time.monotonic() - t0
            print(f"  [{k + 1}/{t_count}] {dt:.0f}s elapsed, "
                  f"eta {dt / (k + 1) * (t_count - k - 1):.0f}s")

    w = w_time[None, :] * np.cos(np.deg2rad(aoi))       # (n, T)

    tables: dict[str, np.ndarray] = {}
    tables["mean_cos"] = np.stack([
        np.average(c, weights=w, axis=1) for c in (c3, c4, c5)], axis=1)
    tables["median"] = np.stack([
        np.array([weighted_median(c[i], w_time) for i in range(n)])
        for c in (c3, c4, c5)], axis=1)

    # Spherical: the aoi -> 0 limit of the same pipeline. rad_s == rad_t == 2f
    # makes the astigmatic terms exactly zero and c4 the pure-focus term.
    sph = np.empty((n, 3))
    for i in range(n):
        sol = strategy.solve(x[i], y[i], 0.0, 45.0, cfg.geometry)  # any sun
        rad = 2.0 * sol.focal_dist_mm
        _c0, b3, b4, b5 = heliostat_shape(0.0, rad, rad)
        sph[i] = to_quadoa_zernike(b3, b4, b5)
    tables["spherical"] = sph

    a.out_dir.mkdir(parents=True, exist_ok=True)
    stamp = _dt.date.today().isoformat()
    for mode in MODES:
        t = tables[mode]
        out = a.out_dir / f"fixed_shapes_{tag}_{mode}.csv"
        with open(out, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(f"# fixed mirror figure, mode={mode}, layout={tag}\n")
            fh.write(f"# generated {stamp} by scripts/build_fixed_shapes.py: "
                     f"{t_count} timesteps over {len(year)} days of 2026, "
                     f"day stride {a.day_stride}\n")
            fh.write("# mean_cos weight = trapezoid daylight hours x cos(AOI); "
                     "median = time-weighted component-wise median\n")
            fh.write("heliostat,x_mm,y_mm,c3,c4,c5\n")
            for i in range(n):
                # float() first: repr of a numpy-2 scalar is "np.float64(...)",
                # which the sweep's loader rejects by design.
                fh.write(f"{i},{float(x[i])!r},{float(y[i])!r},{float(t[i, 0])!r},"
                         f"{float(t[i, 1])!r},{float(t[i, 2])!r}\n")
        print(f"wrote {out}  ({n} heliostats)")

    print("\nfield means of each candidate figure (Quadoa Zernike units):")
    print(f"  {'mode':<10} {'c3':>12} {'c4':>12} {'c5':>12}")
    for mode in MODES:
        t = tables[mode]
        print(f"  {mode:<10} {t[:, 0].mean():>12.4g} {t[:, 1].mean():>12.4g} "
              f"{t[:, 2].mean():>12.4g}")
    d = tables["mean_cos"] - tables["median"]
    ref = np.abs(tables["mean_cos"]).mean(axis=0)
    print(f"\n  mean_cos vs median, mean |delta| per component: "
          f"{np.abs(d).mean(axis=0)} (scale of mean_cos: {ref})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
