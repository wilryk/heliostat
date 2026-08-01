"""One matplotlib style for every figure this project produces, and one way to save it.

The brief, in the owner's words: white background, visible line weights, no
wasted white space, simple designs that communicate effectively, and figures
that survive being dropped into a single column of a paper at ~3.5 in wide.

What that means concretely
-------------------------
**White, everywhere.** Figure, axes and saved-file facecolors are all white and
opaque. A transparent PNG on a white page looks fine and on a dark viewer looks
broken; a grey axes patch prints as grey. Neither is what a journal wants.

**Type sized for 3.5 in, not for a 27 in monitor.** Body text is 9 pt with 8 pt
ticks. In a 3.5 in single-column figure that is the conventional size and stays
legible at print scale; in the GUI's ~9 in figures it reads small but correct.
Sizing for the screen instead would produce figures whose labels are 20 pt once
scaled down, which is the usual reason a paper figure looks amateur.

**Lines you can see.** ``lines.linewidth`` is 2.0. A 1.0 pt curve reduced to
3.5 in and printed is a hairline. Axis spines go the other way -- 0.8 pt and
grey -- so the data dominates the frame rather than competing with it.

**Constrained layout, with the pads turned down.** ``tight_layout`` is a
one-shot fit that has to be re-run after anything changes; constrained layout is
a solver that runs at draw time, handles colorbars and suptitles properly, and
here has its pads cut to roughly a third of the default. That is the single
biggest cure for "more white space than needed". Use :func:`finish` rather than
calling ``fig.tight_layout()``: on a figure that already has a layout engine,
``tight_layout`` *replaces* it and warns.

**No decoration.** Top and right spines off for line plots (:func:`style_axes`),
grids off by default, legends without frames. Nothing is drawn that does not
carry information.

Usage::

    from beamdown import plot_style
    plot_style.apply()                      # once, before any Figure is made
    fig, ax = plot_style.paper_figure()     # 3.5 in single column
    ...
    plot_style.save_figure(fig, "out/annual_energy")   # -> .png (600 dpi) + .pdf

Every figure is written twice on purpose: the PDF is vector, which is what a
journal actually wants and what survives being scaled; the PNG is the one you
can paste into a slide or an email. Rasterised content inside the PDF (an
``imshow`` flux map, for instance) is written at the same 600 dpi, so the PDF is
not quietly the lower-resolution copy.
"""

from __future__ import annotations

from pathlib import Path

# Column widths in inches. 3.5 in is the near-universal single-column figure
# width (IEEE 3.5, Elsevier 3.54, most others within a few hundredths); 7.16 in
# is the matching full-width. Quoted here so a script never has to guess.
SINGLE_COLUMN_IN = 3.5
DOUBLE_COLUMN_IN = 7.16

# 600 dpi for line art and text-bearing raster output. 300 is the usual floor
# for photographs; text and thin lines alias visibly at it.
EXPORT_DPI = 600

# Margin left around the figure when saving with a tight bounding box. Small but
# not zero -- at exactly zero, descenders and minus signs get clipped by the
# rounding in the bbox calculation.
SAVE_PAD_IN = 0.02

# Greys used for anything that is not data.
SPINE_COLOUR = "#8a8a8a"
TICK_COLOUR = "#4d4d4d"
GRID_COLOUR = "#e6e6e6"

_RCPARAMS = {
    # -- white, opaque, everywhere -------------------------------------------
    "figure.facecolor": "white",
    "figure.edgecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white",
    "savefig.transparent": False,

    # -- type: 9 pt body, legible at 3.5 in ----------------------------------
    "font.size": 9.0,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9.0,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "legend.fontsize": 8.0,
    "figure.titlesize": 9.5,
    "axes.titlepad": 4.0,
    "axes.labelpad": 2.5,

    # -- data heavy, frame light ---------------------------------------------
    "lines.linewidth": 2.0,
    "lines.markersize": 4.0,
    "lines.markeredgewidth": 0.0,
    "patch.linewidth": 0.8,
    "axes.linewidth": 0.8,
    "axes.edgecolor": SPINE_COLOUR,
    "axes.labelcolor": "#222222",
    "text.color": "#222222",
    "xtick.color": TICK_COLOUR,
    "ytick.color": TICK_COLOUR,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.minor.size": 1.6,
    "ytick.minor.size": 1.6,

    # -- decoration off by default -------------------------------------------
    "axes.grid": False,
    "grid.color": GRID_COLOUR,
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "legend.handlelength": 1.8,
    "legend.borderpad": 0.3,
    "legend.labelspacing": 0.3,

    # -- kill the white space -------------------------------------------------
    #
    # Constrained layout rather than tight_layout, with h/w pads cut from the
    # default 0.04167 in (3 pt) to 0.015 in and the inter-axes space from 0.02
    # to 0.01 of the figure. See finish().
    "figure.constrained_layout.use": True,
    "figure.constrained_layout.h_pad": 0.015,
    "figure.constrained_layout.w_pad": 0.015,
    "figure.constrained_layout.hspace": 0.01,
    "figure.constrained_layout.wspace": 0.01,

    # -- saving ----------------------------------------------------------------
    "savefig.dpi": EXPORT_DPI,
    "savefig.bbox": "tight",
    "savefig.pad_inches": SAVE_PAD_IN,
    # Real text in the PDF, not outlines: it stays searchable and selectable,
    # and a copy-editor can find a typo in an axis label.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def apply(**overrides) -> dict:
    """Install the paper style into the global rcParams. Returns what it set.

    Call once, **before** any ``Figure`` is constructed: several of these --
    ``figure.constrained_layout.use`` above all -- are read when the figure is
    made, not when it is drawn, so a figure created earlier keeps the old
    behaviour. Cheap and idempotent, so calling it again is harmless.
    """
    import matplotlib

    params = dict(_RCPARAMS)
    params.update(overrides)
    matplotlib.rcParams.update(params)
    return params


def style_axes(ax, spines=("left", "bottom"), grid: str | None = None):
    """Minimal frame: keep only the named spines, mute the ticks.

    ``grid`` is ``None`` (no grid), ``"y"``, ``"x"`` or ``"both"``. A grid is a
    reading aid for a value plot and clutter on an image or a map, so it is
    opt-in rather than a default.
    """
    for side in ("top", "right", "bottom", "left"):
        keep = side in spines
        ax.spines[side].set_visible(keep)
        if keep:
            ax.spines[side].set_linewidth(0.8)
            ax.spines[side].set_color(SPINE_COLOUR)
    ax.tick_params(colors=TICK_COLOUR, labelsize=8.0, width=0.8, length=3.0)
    if grid:
        ax.grid(True, axis=("both" if grid == "both" else grid),
                color=GRID_COLOUR, linewidth=0.6)
        ax.set_axisbelow(True)
    return ax


def finish(fig) -> None:
    """Lay the figure out, whichever engine it has.

    A figure built after :func:`apply` already carries the constrained-layout
    engine and needs nothing -- calling ``fig.tight_layout()`` on it would
    *replace* that engine and emit "The figure layout has changed to tight". A
    figure built before ``apply`` (or by third-party code) has no engine, and
    still wants its margins pulled in. This tells them apart instead of
    guessing.
    """
    try:
        if fig.get_layout_engine() is not None:
            return
    except AttributeError:            # matplotlib < 3.6
        pass
    fig.tight_layout(pad=0.4)


def paper_figure(width_in: float = SINGLE_COLUMN_IN, aspect: float = 0.72,
                 nrows: int = 1, ncols: int = 1, **kwargs):
    """``(fig, ax_or_axes)`` sized for a paper column. ``aspect`` is height/width."""
    import matplotlib.pyplot as plt

    apply()
    return plt.subplots(nrows, ncols,
                        figsize=(width_in, width_in * aspect), **kwargs)


def save_figure(fig, basepath, dpi: int = EXPORT_DPI,
                formats=("png", "pdf")) -> list[Path]:
    """Write ``fig`` once per format beside a common stem. Returns the paths.

    ``basepath`` may carry an image extension or not -- ``plots/annual`` and
    ``plots/annual.png`` both produce ``plots/annual.png`` and
    ``plots/annual.pdf``, so a caller that came from a "Save as..." dialog does
    not end up with ``annual.png.pdf``.

    The bounding box is tight and the facecolor is forced white and opaque here
    as well as in the rcParams, because a caller may have styled the figure by
    hand and this is the last place to get it right.
    """
    base = Path(basepath)
    if base.suffix.lower() in (".png", ".pdf", ".svg", ".eps", ".jpg", ".jpeg",
                               ".tif", ".tiff"):
        base = base.with_suffix("")
    if base.parent and str(base.parent) not in ("", "."):
        base.parent.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    for ext in formats:
        path = base.with_suffix("." + ext.lstrip("."))
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=SAVE_PAD_IN,
                    facecolor="white", edgecolor="white", transparent=False)
        written.append(path)
    return written


def describe() -> str:
    """One line naming what a saved figure will be, for a status bar or a log."""
    return (f"PNG at {EXPORT_DPI} dpi + vector PDF, white background, tight bbox "
            f"({SAVE_PAD_IN:g} in pad)")
