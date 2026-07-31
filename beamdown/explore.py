"""Interactive heliostat explorer.

    python -m beamdown.explore --output analysis_output/demo25

Click any heliostat in the field map to see the receiver spot it produces at the
selected time, and its delivered power across the whole day. The slider moves
through timesteps; the dropdown recolours the field by any summary metric.

This is the tool that makes a long sweep usable: the summary table answers
"which heliostat is worst", and this answers "why".

Needs an interactive matplotlib backend (the default on Windows is fine). For
static output use ``python -m beamdown figures`` instead.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

COLOUR_OPTIONS = [
    "power_w", "eta_shade", "eta_block", "transmission",
    "cosine_efficiency", "r90_mm", "peak_flux_w_m2", "spillage",
]


class Explorer:
    def __init__(self, store, cfg, crop_mm=None):
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        self.store = store
        self.cfg = cfg
        self.crop_mm = crop_mm
        self.summary = store.summary()
        self.keys = store.timestep_keys()
        self.ids = sorted(self.summary.heliostat_id.unique())
        self.step_i = len(self.keys) // 2
        self.colour_by = "power_w"
        self.selected = self.ids[0]

        # heliostat_id -> row index within a timestep's counts array
        self.row_of = {hid: i for i, hid in enumerate(self.ids)}

        self.fig = plt.figure(figsize=(14.5, 7.6))
        gs = GridSpec(2, 2, figure=self.fig, width_ratios=[1.15, 1.0],
                      height_ratios=[1.0, 0.85], hspace=0.30, wspace=0.22,
                      left=0.06, right=0.97, top=0.93, bottom=0.16)
        self.ax_field = self.fig.add_subplot(gs[:, 0])
        self.ax_spot = self.fig.add_subplot(gs[0, 1])
        self.ax_curve = self.fig.add_subplot(gs[1, 1])

        self._build_widgets()
        self._draw_all()

        self.fig.canvas.mpl_connect("button_press_event", self._on_click)
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)

    # -- widgets ----------------------------------------------------------
    def _build_widgets(self):
        from matplotlib.widgets import RadioButtons, Slider

        ax_slider = self.fig.add_axes([0.06, 0.06, 0.45, 0.03])
        self.slider = Slider(
            ax_slider, "timestep", 0, max(1, len(self.keys) - 1),
            valinit=self.step_i, valstep=1,
        )
        self.slider.on_changed(self._on_slider)

        ax_radio = self.fig.add_axes([0.60, 0.015, 0.16, 0.115])
        ax_radio.set_facecolor("#fafafa")
        available = [c for c in COLOUR_OPTIONS if c in self.summary.columns]
        self.radio = RadioButtons(ax_radio, available[:5],
                                  active=available.index(self.colour_by))
        for label in self.radio.labels:
            label.set_fontsize(8)
        self.radio.on_clicked(self._on_colour)

        self.fig.text(0.79, 0.10, "click a heliostat", fontsize=8.5, color="#555555")
        self.fig.text(0.79, 0.075, "left / right arrows: step time", fontsize=8.5,
                      color="#555555")

    # -- data helpers -----------------------------------------------------
    @property
    def key(self) -> str:
        return self.keys[self.step_i]

    def _step_rows(self):
        return self.summary[self.summary.timestep == self.key].sort_values("heliostat_id")

    # -- drawing ----------------------------------------------------------
    def _draw_all(self):
        self._draw_field()
        self._draw_spot()
        self._draw_curve()

    def _draw_field(self):
        from .plots import _style_axis

        ax = self.ax_field
        ax.clear()
        rows = self._step_rows()
        vals = rows[self.colour_by].to_numpy(float)

        sc = ax.scatter(rows.x_m, rows.y_m, c=vals, s=48, cmap="viridis",
                        edgecolors="none", picker=True)
        sel = rows[rows.heliostat_id == self.selected]
        if len(sel):
            ax.scatter(sel.x_m, sel.y_m, s=190, facecolors="none",
                       edgecolors="#d6604d", linewidths=2.0, zorder=5)

        ax.set_aspect("equal")
        ax.set_xlabel("x (m)", fontsize=9)
        ax.set_ylabel("y (m)", fontsize=9)
        row = self.summary[self.summary.timestep == self.key].iloc[0]
        ax.set_title(
            f"{row['date']}  {row['hour']:g}:00   sun az {row['solar_az_deg']:.0f}°"
            f"  el {row['solar_el_deg']:.0f}°\ncoloured by {self.colour_by}",
            fontsize=10,
        )
        _style_axis(ax)
        ax.grid(True, color="#eeeeee", linewidth=0.5)
        ax.set_axisbelow(True)

        if getattr(self, "_field_cbar", None) is not None:
            try:
                self._field_cbar.remove()
            except Exception:
                pass
        self._field_cbar = self.fig.colorbar(sc, ax=ax, fraction=0.043, pad=0.02)
        self._field_cbar.set_label(self.colour_by, fontsize=8.5)
        self._field_cbar.ax.tick_params(labelsize=7.5)
        self._field_xy = rows[["x_m", "y_m"]].to_numpy()
        self._field_ids = rows.heliostat_id.to_numpy()

    def _draw_spot(self):
        from .plots import flux_colormap, _style_axis

        ax = self.ax_spot
        ax.clear()
        row_idx = self.row_of.get(self.selected, 0)
        flux = self.store.heliostat_flux(self.key, row_idx, cfg=self.cfg) / 1000.0

        w = self.cfg.receiver.window_mm
        im = ax.imshow(flux, origin="lower", cmap=flux_colormap(),
                       extent=[-w, w, -w, w], aspect="equal", vmin=0.0,
                       interpolation="nearest")
        if self.crop_mm:
            ax.set_xlim(-self.crop_mm, self.crop_mm)
            ax.set_ylim(-self.crop_mm, self.crop_mm)

        rows = self._step_rows()
        r = rows[rows.heliostat_id == self.selected]
        subtitle = ""
        if len(r):
            r = r.iloc[0]
            subtitle = (f"shade {r['eta_shade']:.2f}  block {r['eta_block']:.2f}  "
                        f"r90 {r['r90_mm']:.0f} mm")
        ax.set_title(f"heliostat {self.selected} — peak {flux.max():.1f} kW/m²\n{subtitle}",
                     fontsize=9.5)
        ax.set_xlabel("x (mm)", fontsize=8.5)
        ax.set_ylabel("y (mm)", fontsize=8.5)
        _style_axis(ax)

        if getattr(self, "_spot_cbar", None) is not None:
            try:
                self._spot_cbar.remove()
            except Exception:
                pass
        self._spot_cbar = self.fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        self._spot_cbar.set_label("kW/m²", fontsize=8)
        self._spot_cbar.ax.tick_params(labelsize=7.5)

    def _draw_curve(self):
        from .plots import _style_axis

        ax = self.ax_curve
        ax.clear()
        mine = self.summary[self.summary.heliostat_id == self.selected]
        current_date = self.summary[self.summary.timestep == self.key].iloc[0]["date"]

        for date, grp in mine.groupby("date"):
            grp = grp.sort_values("hour")
            is_current = date == current_date
            ax.plot(grp.hour, grp.power_w / 1000.0, "-o",
                    linewidth=2.2 if is_current else 1.0,
                    markersize=4.5 if is_current else 2.5,
                    color="#d6604d" if is_current else "#bbbbbb",
                    zorder=3 if is_current else 1,
                    label=str(date) if is_current else None)

        row = self.summary[self.summary.timestep == self.key].iloc[0]
        ax.axvline(row["hour"], color="#888888", linewidth=0.8, linestyle="--", zorder=0)
        ax.set_xlabel("Local hour", fontsize=8.5)
        ax.set_ylabel("Power (kW)", fontsize=8.5)
        ax.set_title(f"heliostat {self.selected} through the day "
                     f"(highlighted: {current_date})", fontsize=9.5)
        _style_axis(ax)
        ax.grid(True, axis="y", color="#eeeeee", linewidth=0.5)
        ax.set_axisbelow(True)

    # -- events -----------------------------------------------------------
    def _on_slider(self, val):
        self.step_i = int(val)
        self._draw_all()
        self.fig.canvas.draw_idle()

    def _on_colour(self, label):
        self.colour_by = label
        self._draw_field()
        self.fig.canvas.draw_idle()

    def _on_click(self, event):
        if event.inaxes is not self.ax_field or event.xdata is None:
            return
        d = np.hypot(self._field_xy[:, 0] - event.xdata,
                     self._field_xy[:, 1] - event.ydata)
        self.selected = int(self._field_ids[int(np.argmin(d))])
        self._draw_all()
        self.fig.canvas.draw_idle()

    def _on_key(self, event):
        if event.key in ("right", "left"):
            delta = 1 if event.key == "right" else -1
            self.step_i = int(np.clip(self.step_i + delta, 0, len(self.keys) - 1))
            self.slider.set_val(self.step_i)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="beamdown.explore")
    ap.add_argument("--config", default=None)
    ap.add_argument("--output", default=None, help="sweep output directory")
    ap.add_argument("--crop", type=float, default=None, help="crop the spot to +/- mm")
    args = ap.parse_args(argv)

    from .config import load_config
    from .store import RunStore

    cfg = load_config(args.config)
    if args.output:
        object.__setattr__(cfg.storage, "root", args.output)

    store = RunStore(cfg.output_root, cfg=cfg, mode="r")
    import matplotlib.pyplot as plt

    Explorer(store, cfg, crop_mm=args.crop)
    plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
