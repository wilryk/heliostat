"""Figures.

Design notes
------------
**One shared colour scale across every panel of a through-the-day figure.** The
legacy code auto-scaled ``vmax`` per plot, which makes a time sequence actively
misleading -- a weak morning spot and a strong noon spot render identically.
Every multi-panel function here computes one scale over all panels and states it
on a single shared colourbar.

**Sequential, single-hue ramp for magnitude.** The truncated ``magma`` from
``quadoa_tools/QuadoaIrradiance.py`` is kept so new figures sit alongside earlier
ones, and it is perceptually uniform, which a rainbow map is not. Built through
the modern colormap API -- ``cm.get_cmap`` is removed in matplotlib 3.11.

**Axes recede, data dominates.** Thin spines, muted ticks, no chartjunk.
"""

from __future__ import annotations

import datetime as _dt

import numpy as np
import pandas as pd

_FLUX_CMAP = None


def flux_colormap():
    """Truncated magma, matching the legacy irradiance figures."""
    global _FLUX_CMAP
    if _FLUX_CMAP is None:
        import matplotlib as mpl
        from matplotlib.colors import ListedColormap

        base = mpl.colormaps["magma"].resampled(128)
        colors = base(np.linspace(0, 1, 128))
        _FLUX_CMAP = ListedColormap(np.concatenate((colors[30:31], colors[45:127]), axis=0))
    return _FLUX_CMAP


def _style_axis(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.6)
        ax.spines[side].set_color("#888888")
    ax.tick_params(colors="#555555", labelsize=8, width=0.6)


def _parse_key_date(key: str) -> _dt.date:
    return _dt.datetime.strptime(key.split("_")[0], "%Y%m%d").date()


def _parse_key_hour(key: str) -> float:
    hhmm = key.split("_")[1]
    return int(hhmm[:2]) + int(hhmm[2:]) / 60.0


def through_day_panels(
    store,
    cfg,
    keys=None,
    efficiency_by_key=None,
    dni_provider=None,
    crop_mm=None,
    percentile=99.9,
    figsize_per_panel=1.5,
    title="Receiver flux through the day",
    save_path=None,
):
    """The article figure: one row per date, one column per hour.

    ``efficiency_by_key`` maps a timestep key to a per-heliostat efficiency
    array (shading x blocking); omit it for the unshaded case.
    """
    import matplotlib.pyplot as plt

    keys = list(keys) if keys is not None else store.timestep_keys()
    if not keys:
        raise ValueError("No timesteps in the store")

    dates = sorted({_parse_key_date(k) for k in keys})
    hours = sorted({_parse_key_hour(k) for k in keys})
    by_cell = {(_parse_key_date(k), _parse_key_hour(k)): k for k in keys}

    # Compute every map once, then derive a single shared scale.
    maps: dict[str, np.ndarray] = {}
    for k in keys:
        dni = dni_provider.dni(_parse_key_date(k), _parse_key_hour(k)) if dni_provider else 1000.0
        eff = efficiency_by_key.get(k) if efficiency_by_key else None
        maps[k] = store.field_flux(k, cfg=cfg, dni_w_m2=dni, efficiency=eff)

    stacked = np.concatenate([m.ravel() for m in maps.values()])
    positive = stacked[stacked > 0]
    vmax = float(np.percentile(positive, percentile)) if positive.size else 1.0

    # ``extent`` describes what the array covers, which is always the full
    # window. Setting it to the crop box instead would rescale the map into that
    # box -- the same picture with smaller axis numbers -- rather than crop it.
    # Cropping is done below with set_xlim/set_ylim, as in single_vs_field.
    window = cfg.receiver.window_mm
    extent = [-window, window, -window, window]

    nrows, ncols = len(dates), len(hours)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize_per_panel * ncols + 1.6, figsize_per_panel * nrows + 1.0),
        squeeze=False,
    )

    cmap = flux_colormap()
    img = None
    for r, date in enumerate(dates):
        for c, hour in enumerate(hours):
            ax = axes[r][c]
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_linewidth(0.5); s.set_color("#cccccc")

            key = by_cell.get((date, hour))
            if key is None:
                ax.set_facecolor("#f7f7f7")
                ax.text(0.5, 0.5, "–", ha="center", va="center",
                        color="#bbbbbb", fontsize=9, transform=ax.transAxes)
            else:
                img = ax.imshow(maps[key] / 1000.0, origin="lower", cmap=cmap,
                                vmin=0.0, vmax=vmax / 1000.0, extent=extent,
                                aspect="equal", interpolation="nearest")
                peak = maps[key].max() / 1000.0
                ax.text(0.03, 0.97, f"{peak:.0f}", transform=ax.transAxes,
                        ha="left", va="top", fontsize=6.5, color="white", alpha=0.85)
            if crop_mm:
                ax.set_xlim(-crop_mm, crop_mm); ax.set_ylim(-crop_mm, crop_mm)

            if r == 0:
                h = int(hour)
                m = int(round((hour - h) * 60))
                ax.set_title(f"{h:02d}:{m:02d}", fontsize=8.5, color="#333333", pad=4)
            if c == 0:
                ax.set_ylabel(f"{date:%b %d}", fontsize=8.5, color="#333333",
                              rotation=0, ha="right", va="center", labelpad=26)

    fig.suptitle(title, fontsize=12, y=0.985)
    if img is not None:
        cbar = fig.colorbar(img, ax=[a for row in axes for a in row],
                            fraction=0.018, pad=0.012)
        cbar.set_label("Receiver flux (kW/m²)", fontsize=9)
        cbar.ax.tick_params(labelsize=8, width=0.6, colors="#555555")
        cbar.outline.set_linewidth(0.5)

    note = f"shared colour scale, 0 – {vmax/1000:.0f} kW/m²; panel labels are peak kW/m²"
    fig.text(0.5, 0.005, note, ha="center", fontsize=7.5, color="#777777")

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight", facecolor="white")
    return fig


def single_vs_field(store, cfg, key, heliostat_row=0, heliostat_id=None,
                    crop_mm=None, save_path=None):
    """One heliostat beside the whole field: 'the focus is not just a point'."""
    import matplotlib.pyplot as plt

    single = store.heliostat_flux(key, heliostat_row, cfg=cfg)
    combined = store.field_flux(key, cfg=cfg)

    window = cfg.receiver.window_mm
    extent = [-window, window, -window, window]
    cmap = flux_colormap()

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.0))
    for ax, data, label in (
        (axes[0], single, "single heliostat"
                          + (f" (id {heliostat_id})" if heliostat_id is not None else "")),
        (axes[1], combined, "whole field"),
    ):
        im = ax.imshow(data / 1000.0, origin="lower", cmap=cmap, extent=extent,
                       aspect="equal", vmin=0.0, interpolation="nearest")
        ax.set_title(f"{label}\npeak {data.max()/1000:.1f} kW/m²", fontsize=10)
        ax.set_xlabel("x (mm)", fontsize=9)
        if crop_mm:
            ax.set_xlim(-crop_mm, crop_mm); ax.set_ylim(-crop_mm, crop_mm)
        _style_axis(ax)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("kW/m²", fontsize=8)
        cb.ax.tick_params(labelsize=7)
    axes[0].set_ylabel("y (mm)", fontsize=9)
    fig.suptitle(f"Receiver flux — {key}", fontsize=11)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight", facecolor="white")
    return fig


def encircled_energy_axes(ax, radii_mm, power_w, aperture_mm=None,
                          colour="#2166ac", label=None, fractions=(0.5, 0.9)):
    """Draw one cumulative-power-vs-aperture-radius curve onto ``ax``.

    Plotted in absolute kW rather than as a normalised fraction, so a change that
    scales the whole spot -- shading, blocking, DNI -- visibly moves the curve.
    A fraction-only plot would be invariant to exactly the weights the explorer
    exists to let you toggle.

    Radii are measured about the receiver axis, not the spot centroid, because
    that is what a physical aperture intercepts; the caller is expected to say so
    in the axis label.
    """
    radii_mm = np.asarray(radii_mm, float)
    power_w = np.asarray(power_w, float)
    total = float(power_w[-1]) if power_w.size else 0.0

    ax.plot(radii_mm, power_w / 1000.0, "-", color=colour, linewidth=2.0,
            label=label, zorder=3)

    marks = []
    for frac in fractions:
        if total <= 0:
            continue
        r = float(np.interp(frac * total, power_w, radii_mm))
        ax.plot([r], [frac * total / 1000.0], "o", color=colour, markersize=5,
                zorder=4)
        ax.annotate(f"r{int(frac*100)} {r:.0f} mm", (r, frac * total / 1000.0),
                    textcoords="offset points", xytext=(7, -3),
                    fontsize=8, color=colour)
        marks.append(r)

    if aperture_mm:
        captured = float(np.interp(aperture_mm, radii_mm, power_w))
        note = f"Ø{aperture_mm:.0f}"
        if total > 0:
            note += f"\n{captured/1000:.1f} kW\nspill {1 - captured/total:.1%}"
        ax.axvline(aperture_mm, color="#d6604d", linewidth=1.4, linestyle="--",
                   zorder=2)
        ax.annotate(note, (aperture_mm, captured / 1000.0),
                    textcoords="offset points", xytext=(8, 8), fontsize=8,
                    color="#d6604d", va="bottom")

    _style_axis(ax)
    ax.grid(True, color="#eeeeee", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    return marks


def annual_energy_axes(ax, daily, sine_fit=None, traced_dates=(), hourly=None,
                       title="Daily collected energy over the year"):
    """Draw one year of daily collected energy onto ``ax``, sine fit overlaid.

    ``daily`` is ``energy.annual_energy(...)["daily"]``: one row per day of
    the year, almost all of them interpolated from a handful of traced dates
    rather than ray-traced themselves. ``traced_dates`` (an iterable of
    ``datetime.date``) marks the days that carry real ray-trace data --
    everything else rides on the (declination, hour angle) interpolation
    described at the top of ``energy.py``, so the distinction matters.

    ``sine_fit`` is ``energy.fit_annual_sine(daily)``; passing it overlays the
    fitted curve and states its R² so the fit is never presented as more than
    it is.

    ``hourly`` is ``energy.annual_energy(...)["hourly"]``, used only to shade
    days where more than half the daylight hours needed extrapolation beyond
    the traced declination hull -- the least trustworthy part of the curve.
    Pass ``None`` to skip the shading if it clutters the figure.
    """
    import matplotlib.dates as mdates

    daily = daily.sort_values("date")
    x = pd.to_datetime(daily["date"])
    y_mwh = daily["energy_kwh"].to_numpy(float) / 1000.0
    traced = set(traced_dates)
    is_traced = daily["date"].isin(traced).to_numpy()

    if hourly is not None and "extrapolated" in hourly.columns:
        frac = hourly.groupby("date")["extrapolated"].mean()
        for d, f in frac.items():
            if f > 0.5:
                t = pd.Timestamp(d)
                ax.axvspan(t - pd.Timedelta(hours=12), t + pd.Timedelta(hours=12),
                          color="#f4a261", alpha=0.12, zorder=0, linewidth=0)

    ax.plot(x, y_mwh, "-", color="#9ecae1", linewidth=1.0, zorder=1)
    ax.scatter(x[~is_traced], y_mwh[~is_traced], s=8, color="#4393c3", zorder=2,
              label="interpolated day")
    if is_traced.any():
        ax.scatter(x[is_traced], y_mwh[is_traced], s=46, color="#d6604d",
                  edgecolors="white", linewidths=0.7, zorder=4, label="traced day")

    fit_note = ""
    if sine_fit is not None:
        doy_grid = np.linspace(1.0, 366.0, 400)
        y_fit = sine_fit["predict"](doy_grid) / 1000.0
        year = int(x.dt.year.iloc[0]) if len(x) else _dt.date.today().year
        x_fit = pd.Timestamp(f"{year}-01-01") + pd.to_timedelta(doy_grid - 1.0, unit="D")
        ax.plot(x_fit, y_fit, "--", color="#4d4d4d", linewidth=1.6, zorder=3,
               label="fitted sinusoid")
        peak = pd.Timestamp(f"{year}-01-01") + pd.Timedelta(days=sine_fit["phase_day_of_year"] - 1.0)
        fit_note = (f"  ·  fitted sinusoid: mean {sine_fit['mean']/1000:.2f} MWh, "
                    f"amplitude {sine_fit['amplitude']/1000:.2f} MWh, "
                    f"peak ~{peak:%b %d}, R²={sine_fit['r_squared']:.2f}")

    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.set_ylabel("Energy collected (MWh/day)", fontsize=9)
    ax.set_title(title + fit_note, fontsize=9.5)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.set_ylim(bottom=0.0)
    _style_axis(ax)
    ax.grid(True, color="#eeeeee", linewidth=0.5)
    ax.set_axisbelow(True)
    return ax


def annual_energy_figure(daily, sine_fit=None, traced_dates=(), hourly=None,
                         title="Daily collected energy over the year", save_path=None):
    """Standalone figure wrapping :func:`annual_energy_axes`, for direct use outside the GUI."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    annual_energy_axes(ax, daily, sine_fit=sine_fit, traced_dates=traced_dates,
                       hourly=hourly, title=title)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight", facecolor="white")
    return fig


def field_scatter(summary, column, timestep=None, cfg=None, title=None,
                  cmap="viridis", save_path=None):
    """All heliostats, positioned as in the field, coloured by any metric."""
    import matplotlib.pyplot as plt

    data = summary[summary.timestep == timestep] if timestep else (
        summary.groupby("heliostat_id", as_index=False).agg(
            x_m=("x_m", "first"), y_m=("y_m", "first"), **{column: (column, "mean")}
        )
    )

    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    sc = ax.scatter(data.x_m, data.y_m, c=data[column], s=26, cmap=cmap,
                    edgecolors="none")
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)", fontsize=9)
    ax.set_ylabel("y (m)", fontsize=9)
    ax.set_title(title or f"{column}" + (f" — {timestep}" if timestep else " (mean over time)"),
                 fontsize=11)
    _style_axis(ax)
    ax.grid(True, color="#eeeeee", linewidth=0.5)
    ax.set_axisbelow(True)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(column, fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight", facecolor="white")
    return fig


def power_through_day(summary, value="power_w", ylabel="Receiver power (kW)",
                      scale=1e-3, save_path=None, title="Delivered power through the day"):
    """One line per date. Direct-labelled, so no legend box is needed."""
    import matplotlib.pyplot as plt

    grouped = (summary.groupby(["date", "hour"], as_index=False)[value].sum())
    dates = sorted(grouped["date"].unique())

    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    colors = ["#2166ac", "#d6604d", "#4d9221", "#8073ac", "#b8860b",
              "#0f7b7b", "#a6522c", "#7b3294"]
    for i, date in enumerate(dates):
        sub = grouped[grouped["date"] == date].sort_values("hour")
        color = colors[i % len(colors)]
        ax.plot(sub["hour"], sub[value] * scale, "-o", color=color,
                linewidth=2.0, markersize=4.5)
        last = sub.iloc[-1]
        ax.annotate(f" {date:%b %d}" if hasattr(date, "strftime") else f" {date}",
                    (last["hour"], last[value] * scale),
                    color=color, fontsize=8.5, va="center")

    ax.set_xlabel("Local hour", fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=11)
    ax.set_xlim(right=ax.get_xlim()[1] + 1.4)
    _style_axis(ax)
    ax.grid(True, axis="y", color="#eeeeee", linewidth=0.5)
    ax.set_axisbelow(True)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight", facecolor="white")
    return fig


def efficiency_breakdown(summary, save_path=None):
    """Cosine, shading and blocking through the day -- where the losses are."""
    import matplotlib.pyplot as plt

    grouped = summary.groupby(["date", "hour"], as_index=False).agg(
        cosine=("cosine_efficiency", "mean"),
        shade=("eta_shade", "mean"),
        block=("eta_block", "mean"),
        transmission=("transmission", "mean"),
    )
    dates = sorted(grouped["date"].unique())

    fig, axes = plt.subplots(1, len(dates), figsize=(3.4 * len(dates), 3.6),
                             sharey=True, squeeze=False)
    series = [("cosine", "#2166ac", "cosine"), ("shade", "#d6604d", "shading"),
              ("block", "#4d9221", "blocking"), ("transmission", "#8073ac", "optical")]

    for i, date in enumerate(dates):
        ax = axes[0][i]
        sub = grouped[grouped["date"] == date].sort_values("hour")
        for col, color, label in series:
            ax.plot(sub["hour"], sub[col], "-", color=color, linewidth=2.0,
                    label=label if i == 0 else None)
        ax.set_title(f"{date:%b %d}" if hasattr(date, "strftime") else str(date), fontsize=10)
        ax.set_xlabel("Local hour", fontsize=9)
        ax.set_ylim(0, 1.02)
        _style_axis(ax)
        ax.grid(True, axis="y", color="#eeeeee", linewidth=0.5)
        ax.set_axisbelow(True)
    axes[0][0].set_ylabel("Efficiency", fontsize=9)
    axes[0][0].legend(frameon=False, fontsize=8.5, loc="lower center")
    fig.suptitle("Loss mechanisms through the day", fontsize=11)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight", facecolor="white")
    return fig
