"""Shared matplotlib style presets used by both CLI plots and ``draw_fit``.

* ``_ERRBAR_STYLE`` / ``_MARKER_STYLE`` — kwargs dicts so all panels render
  data points and crossing markers with identical visual weight.
* ``_setup_plot_style`` — rcParams preset for publication-quality output.
* ``_journal_style_ticks`` — inward major+minor ticks on every side, matched
  to journal figure conventions.
* ``_smart_ylim`` — log-flux y-limit picker capped at 12 decades.
"""

import numpy as np

# Errorbar style for plotted data points (Fitter.draw_fit + CLI plots).
_ERRBAR_STYLE = dict(
    linestyle="none",
    linewidth=0.5,
    markeredgecolor="black",
    fmt="o",
    markersize=4,
    markeredgewidth=0.2,
    capsize=3,
    capthick=3,
)

# Bare marker style (no errorbars) for the bottom panel's crossing points.
# Smaller than the upper panel's data markers so they don't compete visually.
_MARKER_STYLE = dict(
    linestyle="none",
    markersize=3,
    markeredgecolor=_ERRBAR_STYLE["markeredgecolor"],
    markeredgewidth=_ERRBAR_STYLE["markeredgewidth"],
)


def _setup_plot_style(font=None):
    """Configure matplotlib for publication-quality output.

    Uses mathtext with STIX fonts (no LaTeX subprocess) for fast rendering.
    """
    import matplotlib as mpl

    _SANS_SERIF = {
        "helvetica",
        "arial",
        "liberation sans",
        "dejavu sans",
        "gill sans",
        "futura",
        "optima",
        "verdana",
        "tahoma",
    }

    if font is None:
        family = "sans-serif"
        font_list = ["Helvetica"]
        math_font = "Helvetica"
    elif font.lower() in _SANS_SERIF:
        family = "sans-serif"
        font_list = [font]
        math_font = font
    else:
        family = "serif"
        font_list = [font]
        math_font = font

    style = {
        "text.usetex": False,
        "mathtext.fontset": "custom",
        "mathtext.rm": math_font,
        "mathtext.it": f"{math_font}:italic",
        "mathtext.bf": f"{math_font}:bold",
        "font.family": family,
        f"font.{family}": font_list,
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "legend.fontsize": 5.5,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "axes.linewidth": 0.5,
        "lines.linewidth": 1.0,
        "axes.grid": True,
        "axes.grid.which": "both",
        "grid.color": "0.7",
        "grid.linestyle": ":",
        "grid.linewidth": 0.3,
        "grid.alpha": 0.4,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
    }

    mpl.rcParams.update(style)


def _journal_style_ticks(ax_top, ax_top_R, ax_bot):
    """Apply journal-paper tick styling: inward major + minor ticks on all sides.

    When ``ax_top_R`` is present it owns the right ticks; otherwise ``ax_top``
    shows them on its own right edge.
    """
    common = dict(direction="in", which="both", labelsize=10, length=4, width=0.8)
    minor_kw = dict(direction="in", which="minor", length=2.5, width=0.6)
    right_on_top = ax_top_R is None

    ax_top.minorticks_on()
    ax_top.tick_params(**common, top=True, bottom=True, left=True, right=right_on_top)
    ax_top.tick_params(**minor_kw, top=True, bottom=True, left=True, right=right_on_top)

    if ax_top_R is not None:
        ax_top_R.minorticks_on()
        ax_top_R.tick_params(**common, right=True, left=False, top=False, bottom=False)
        ax_top_R.tick_params(
            **minor_kw, right=True, left=False, top=False, bottom=False
        )

    if ax_bot is not None:
        ax_bot.minorticks_on()
        ax_bot.tick_params(**common, top=True, bottom=True, left=True, right=True)
        ax_bot.tick_params(**minor_kw, top=True, bottom=True, left=True, right=True)


def _smart_ylim(ax, flux_arrays, scale=1.0):
    """Set smart log y-axis limits from flux arrays, capped at 12 decades."""
    all_pos = []
    for arr in flux_arrays:
        scaled = arr / scale
        pos = scaled[scaled > 0]
        if pos.size > 0:
            all_pos.append(pos)
    if all_pos:
        combined = np.concatenate(all_pos)
        f_max, f_min = np.max(combined), np.min(combined)
        y_top = f_max * 10
        y_bot = f_min / 10
        if y_top / y_bot > 1e12:
            y_bot = y_top * 1e-12
        ax.set_ylim(bottom=y_bot, top=y_top)
