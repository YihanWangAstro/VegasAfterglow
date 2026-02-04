"""Regression test visualization functions for power-law scaling verification."""

import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from validation.visualization.common import (BAND_COLORS, COLORS, PAGE_PORTRAIT, PHASE_COLORS, PHASE_NAMES, QTY_SYMBOLS, format_slope, to_float)
from validation.regression.run_regression import SPECTRAL_REGIMES


def _smooth_slopes(log_x, log_y, window=7):
    """Compute smoothed local power-law slopes from log-spaced data."""
    d_log_x = np.diff(log_x)
    d_log_y = np.diff(log_y)
    slopes = d_log_y / d_log_x
    x_mid = 10 ** (0.5 * (log_x[:-1] + log_x[1:]))
    if len(slopes) >= window:
        kernel = np.ones(window) / window
        trim = (window - 1) // 2
        slopes = np.convolve(slopes, kernel, mode="valid")
        x_mid = x_mid[trim:trim + len(slopes)]
    return x_mid, slopes


def _get_cell_style(result: Dict) -> Tuple[str, str, str, str]:
    if not result:
        return "#E8E8E8", "#888888", "\u2014", ""
    if result.get("measured") is None:
        return "#CCCCCC", "black", "N/A", ""
    exp_str = result.get("expected_expr") or (format_slope(result.get("expected")) if result.get("expected") is not None else "")
    color, tc = ("#90EE90", "darkgreen") if result.get("passed", False) else ("#FFB6C1", "darkred")
    return color, tc, f"{result['measured']:.2f}", exp_str


def _draw_cell(fig, x, y, width, height, color, text_color, line1, line2, fontsize=8,
               measured=None, expected=None, tolerance=None, passed=None):
    """Draw a cell with optional deviation indicator bar."""
    rect = plt.Rectangle((x + 0.003, y + 0.003), width * 0.94, height * 0.9, facecolor=color, edgecolor="white",
                          linewidth=1.5, transform=fig.transFigure, clip_on=False)
    fig.add_artist(rect)

    has_bar = (measured is not None and expected is not None and tolerance is not None and tolerance > 0)

    if has_bar:
        # Deviation bar centered vertically in cell
        bar_y = y + height * 0.38
        bar_h = height * 0.14
        bar_left = x + width * 0.10
        bar_width = width * 0.74
        bar_range = 3.0 * tolerance  # full range: expected ± 3*tol

        # Tolerance zone (middle third)
        tol_frac = tolerance / bar_range
        tol_rect = plt.Rectangle(
            (bar_left + bar_width * (0.5 - tol_frac), bar_y),
            bar_width * 2 * tol_frac, bar_h,
            facecolor="#B0B0B0", edgecolor="none", alpha=0.35,
            transform=fig.transFigure, clip_on=False)
        fig.add_artist(tol_rect)

        # Center tick at expected value
        cx = bar_left + bar_width * 0.5
        center_line = plt.Line2D([cx, cx], [bar_y, bar_y + bar_h], color="#999999",
                                  linewidth=0.7, transform=fig.transFigure, clip_on=False)
        fig.add_artist(center_line)

        # Measured marker
        deviation = measured - expected
        frac = 0.5 + (deviation / bar_range) * 0.5
        frac = max(0.03, min(0.97, frac))
        mx = bar_left + bar_width * frac
        marker_color = "#2E7D32" if passed else "#C62828"
        marker = plt.Line2D([mx], [bar_y + bar_h / 2], marker='o', color=marker_color,
                             markersize=3.5, linestyle='None', transform=fig.transFigure,
                             clip_on=False, zorder=5)
        fig.add_artist(marker)

        # Measured value text: above bar, anchored at dot x-position, same color as dot
        text_x_min = x + width * 0.08
        text_x_max = x + width * 0.86
        mx_clamped = max(text_x_min, min(text_x_max, mx))
        fig.text(mx_clamped, y + height * 0.78, line1, ha="center", va="center", fontsize=fontsize,
                 fontweight="bold", color=marker_color)

        # Expected value text: below bar, anchored at center line x-position, same color as center line
        if line2:
            fig.text(cx, y + height * 0.18, line2, ha="center", va="center",
                     fontsize=fontsize - 1, color="#999999")
    else:
        # Original layout without deviation bar
        fig.text(x + width / 2, y + height * 0.6, line1, ha="center", va="center", fontsize=fontsize,
                 fontweight="bold", color=text_color)
        if line2:
            fig.text(x + width / 2, y + height * 0.28, line2, ha="center", va="center",
                     fontsize=fontsize - 1, color=text_color)


def _count_pass_fail(grid_data: Dict) -> Tuple[int, int]:
    results = [r for md in grid_data.values() for pd in md.values() for r in pd.values() if r.get("measured") is not None]
    return sum(1 for r in results if r.get("passed")), len(results)


def _count_spectrum_pass_fail(spectrum_grid: Dict) -> Tuple[int, int]:
    """Count pass/fail for spectrum shape tests."""
    results = [s for md in spectrum_grid.values() for rd in md.values() for s in rd.get("segments", []) if s.get("measured") is not None]
    return sum(1 for r in results if r.get("passed")), len(results)


def _draw_grid_section(fig, grid_data, quantities, title, top, bottom, grid_left, cell_width, phases, media):
    """Draw a grid section and return (n_pass, n_total, medium_rects).

    medium_rects is a dict mapping medium name to (x, y, w, h) in figure coordinates.
    """
    n_pass, n_total = _count_pass_fail(grid_data)
    rate = n_pass / n_total * 100 if n_total > 0 else 0
    title_color = COLORS["pass"] if rate >= 80 else COLORS["fail"]
    cell_height = (top - bottom) / len(quantities)

    # Create separate clickable rectangles for each medium (ISM, Wind)
    medium_rects = {}
    n_phase_cols = len(phases)
    for k, medium in enumerate(media):
        x_start = grid_left + k * n_phase_cols * cell_width
        rect_width = n_phase_cols * cell_width
        medium_rects[medium] = (x_start - 0.01, bottom, rect_width + 0.02, top - bottom + 0.05)

    fig.text(0.5, top + 0.03, f"{title} ({n_pass}/{n_total})", ha="center", fontsize=11,
             fontweight="bold", color=title_color)

    for k, medium in enumerate(media):
        for j, phase in enumerate(phases):
            col = k * len(phases) + j
            x = grid_left + col * cell_width + cell_width / 2
            fig.text(x, top + 0.01, f"{medium}|{PHASE_NAMES[phase][:5]}", ha="center", va="bottom",
                     fontsize=7, color="#333333")

    for i, qty in enumerate(quantities):
        y_center = top - (i + 0.5) * cell_height
        fig.text(grid_left - 0.015, y_center, QTY_SYMBOLS.get(qty, qty), ha="right", va="center",
                 fontsize=9, fontweight="bold")

        for k, medium in enumerate(media):
            for j, phase in enumerate(phases):
                col = k * len(phases) + j
                x = grid_left + col * cell_width
                y = top - (i + 1) * cell_height

                phase_data = grid_data.get(medium, {}).get(phase, {})
                qty_result = phase_data.get(qty, {})

                color, text_color, line1, line2 = _get_cell_style(qty_result)
                _draw_cell(fig, x, y, cell_width, cell_height, color, text_color, line1, line2, fontsize=8,
                           measured=to_float(qty_result.get("measured")),
                           expected=to_float(qty_result.get("expected")),
                           tolerance=to_float(qty_result.get("tolerance")),
                           passed=qty_result.get("passed"))

    return n_pass, n_total, medium_rects


def _format_segment_label(seg_name: str) -> str:
    """Convert segment name to readable label like 'ν < ν_a' or 'ν_a < ν < ν_m'."""
    _fs = {"nu_a": r"$\nu_a$", "nu_m": r"$\nu_m$", "nu_c": r"$\nu_c$"}
    for prefix, fmt in [("below_", r"$\nu<${}"), ("above_", r"$\nu>${}")]:
        if seg_name.startswith(prefix):
            freq = seg_name[len(prefix):]
            return fmt.format(_fs.get(freq, freq))
    if "_to_" in seg_name and len(parts := seg_name.split("_to_")) == 2:
        return _fs.get(parts[0], parts[0]) + r"$<\nu<$" + _fs.get(parts[1], parts[1])
    return seg_name


def _draw_spectrum_grid_section(fig, spectrum_grid, top, bottom, grid_left, grid_right,
                                section_title="Spectrum Shapes"):
    """Draw spectrum shapes grid with regimes as rows and segments as columns."""
    regimes = ["I", "II", "III", "IV", "V"]
    n_cols = 4

    n_pass, n_total = _count_spectrum_pass_fail(spectrum_grid)
    rate = n_pass / n_total * 100 if n_total > 0 else 0
    title_color = COLORS["pass"] if rate >= 80 else COLORS["fail"]

    # Layout: reserve space for column headers above each row
    header_height = 0.018  # Height for column header row per regime
    row_height = (top - bottom) / len(regimes)
    cell_height = row_height - header_height
    cell_width = (grid_right - grid_left) / n_cols

    fig.text(0.5, top + 0.02, f"{section_title} ({n_pass}/{n_total})", ha="center", fontsize=11,
             fontweight="bold", color=title_color)

    # Draw each regime row with its own column headers
    for i, regime in enumerate(regimes):
        row_top = top - i * row_height
        row_bottom = row_top - row_height

        # Row label on the left
        y_center = row_bottom + row_height / 2
        fig.text(grid_left - 0.015, y_center, f"Regime {regime}", ha="right", va="center",
                 fontsize=8, fontweight="bold")

        regime_data = next((md[regime] for md in spectrum_grid.values() if regime in md), None)
        segments = regime_data.get("segments", []) if regime_data else []

        # Get expected segment names from SPECTRAL_REGIMES (not from actual data)
        expected_segments = SPECTRAL_REGIMES.get(regime, {}).get("segments", [])

        # Draw column headers for this row (using expected segment names)
        for j in range(n_cols):
            x = grid_left + j * cell_width
            if j < len(expected_segments):
                seg_name = expected_segments[j][0]  # First element is segment name
                label = _format_segment_label(seg_name)
            else:
                label = ""
            fig.text(x + cell_width / 2, row_top - 0.003, label,
                     ha="center", va="top", fontsize=9, color="#333333", fontweight="bold")

        # Build lookup from segment name to result
        seg_lookup = {s.get("segment", ""): s for s in segments}

        # Draw cells below headers (matching by expected segment name)
        for j in range(n_cols):
            x, y = grid_left + j * cell_width, row_bottom
            seg_result = seg_lookup.get(expected_segments[j][0]) if j < len(expected_segments) else None
            color, text_color, line1, line2 = _get_cell_style(seg_result) if seg_result else ("#E8E8E8", "#888888", "\u2014", "")
            _draw_cell(fig, x, y, cell_width, cell_height, color, text_color, line1, line2, fontsize=7,
                       measured=to_float(seg_result.get("measured")) if seg_result else None,
                       expected=to_float(seg_result.get("expected")) if seg_result else None,
                       tolerance=to_float(seg_result.get("tolerance")) if seg_result else None,
                       passed=seg_result.get("passed") if seg_result else None)

    return n_pass, n_total


_FWD_PHASES = ["coasting", "BM", "deep_newtonian"]
_RVS_PHASES = ["crossing", "BM", "deep_newtonian"]
_MEDIA = ["ISM", "Wind"]
_GRID_LEFT = 0.18


def _plot_multi_grid_summary(results, page_title, quantities, grids, suptitle_y=0.96):
    """Shared layout for dynamics and frequencies summary grids.

    Args:
        results: Full results dict.
        page_title: Figure suptitle prefix.
        quantities: List of quantity names for the grid rows.
        grids: List of (grid_key, section_title, phases, top, bottom) tuples.
        suptitle_y: Vertical position of the suptitle.

    Returns:
        (fig, section_rects) where section_rects is a dict mapping grid_key to
        a dict of {medium: (x, y, w, h)} for clickable regions per medium.
    """
    fig = plt.figure(figsize=PAGE_PORTRAIT)
    cell_width = (0.92 - _GRID_LEFT) / (len(_FWD_PHASES) * len(_MEDIA))

    total_pass, total_tests = 0, 0
    section_rects = {}
    for grid_key, title, phases, top, bottom in grids:
        grid_data = results.get(grid_key, {})
        if grid_data:
            n_p, n_t, medium_rects = _draw_grid_section(
                fig, grid_data, quantities, title, top=top, bottom=bottom,
                grid_left=_GRID_LEFT, cell_width=cell_width, phases=phases, media=_MEDIA)
            total_pass += n_p
            total_tests += n_t
            section_rects[grid_key] = medium_rects

    total_rate = total_pass / total_tests * 100 if total_tests > 0 else 0
    status_color = COLORS["pass"] if total_rate >= 80 else COLORS["fail"]
    fig.suptitle(f"{page_title}: {total_pass}/{total_tests} passed ({total_rate:.0f}%)",
                 fontsize=13, fontweight="bold", y=suptitle_y, color=status_color)
    return fig, section_rects


def plot_dynamics_summary_grid(results: Dict) -> Tuple[plt.Figure, Dict[str, Tuple]]:
    """Plot dynamics summary: FWD + RVS thin + RVS thick shock dynamics grids.

    Returns:
        (fig, section_rects) where section_rects maps grid_key to (x, y, w, h).
    """
    return _plot_multi_grid_summary(results, "Dynamics Summary", ["u", "r", "B", "N_p"], [
        ("shock_grid",           "Forward Shock",                _FWD_PHASES, 0.86, 0.66),
        ("rvs_shock_grid_thin",  "Reverse Shock \u2014 Thin Shell",  _RVS_PHASES, 0.56, 0.36),
        ("rvs_shock_grid_thick", "Reverse Shock \u2014 Thick Shell", _RVS_PHASES, 0.26, 0.06),
    ])


def plot_spectrum_summary_grid(results: Dict) -> Tuple[plt.Figure, Dict[str, Dict]]:
    """Plot standalone spectrum shapes pass/fail grid.

    Returns:
        (fig, section_rects) where section_rects["spectrum_grid"] is a single rect
        (spectrum has no ISM/Wind split).
    """
    fig = plt.figure(figsize=PAGE_PORTRAIT)

    spectrum_grid = results.get("spectrum_grid", {})
    grid_right = 0.92
    top, bottom = 0.78, 0.40

    n_pass, n_total = 0, 0
    section_rects = {}
    if spectrum_grid:
        n_pass, n_total = _draw_spectrum_grid_section(
            fig, spectrum_grid, top=top, bottom=bottom, grid_left=_GRID_LEFT, grid_right=grid_right,
            section_title="Synchrotron Spectrum Shapes",
        )
        # Clickable region covers the entire grid (no ISM/Wind split for spectrum)
        full_rect = (_GRID_LEFT - 0.02, bottom, grid_right - _GRID_LEFT + 0.04, top - bottom + 0.05)
        section_rects["spectrum_grid"] = {"_all": full_rect}

    total_rate = n_pass / n_total * 100 if n_total > 0 else 0
    status_color = COLORS["pass"] if total_rate >= 80 else COLORS["fail"]
    fig.suptitle(f"Synchrotron Spectrum Shapes: {n_pass}/{n_total} passed ({total_rate:.0f}%)",
                 fontsize=13, fontweight="bold", y=0.92, color=status_color)
    return fig, section_rects


def plot_frequencies_summary_grid(results: Dict) -> Tuple[plt.Figure, Dict[str, Tuple]]:
    """Plot frequencies summary: FWD + RVS thin + RVS thick frequency grids.

    Returns:
        (fig, section_rects) where section_rects maps grid_key to (x, y, w, h).
    """
    return _plot_multi_grid_summary(results, "Frequencies Summary", ["nu_m", "nu_c", "nu_M"], [
        ("freq_grid",           "Forward Shock",                _FWD_PHASES, 0.86, 0.70),
        ("rvs_freq_grid_thin",  "Reverse Shock \u2014 Thin Shell",  _RVS_PHASES, 0.56, 0.40),
        ("rvs_freq_grid_thick", "Reverse Shock \u2014 Thick Shell", _RVS_PHASES, 0.26, 0.10),
    ])


def _smooth_1d(arr, size):
    """Uniform 1D smoothing with nearest-edge padding."""
    return np.convolve(np.pad(arr, size // 2, mode="edge"), np.ones(size) / size, mode="valid")


def _plot_quantity_panel(ax_qty, ax_deriv, t, y, qty, phases):

    if len(y) == 0:
        ax_qty.text(0.5, 0.5, f"No {qty} data", ha="center", va="center", transform=ax_qty.transAxes)
        ax_qty.set_ylabel(QTY_SYMBOLS.get(qty, qty), fontsize=9)
        ax_deriv.set_xlabel("Observer Time [s]", fontsize=7)
        return

    valid = (y > 0) & np.isfinite(y)
    if np.sum(valid) < 3:
        ax_qty.text(0.5, 0.5, f"No valid {qty} data", ha="center", va="center", transform=ax_qty.transAxes)
        return

    ax_qty.loglog(t[valid], y[valid], "k-", alpha=0.4, linewidth=1)
    t_min, t_max = t[valid].min(), t[valid].max()

    for phase, phase_data in phases.items():
        if qty not in phase_data.get("fits", {}):
            continue
        t_range, fit_info = phase_data["t_range"], phase_data["fits"][qty]
        alpha_exp = fit_info.get("expected")
        mask = (t >= t_range[0]) & (t <= t_range[1]) & valid
        if np.sum(mask) < 2:
            continue
        color = PHASE_COLORS.get(phase, "gray")
        t_center = np.sqrt(t_range[0] * t_range[1])
        if alpha_exp is not None:
            t_fit = np.logspace(np.log10(max(t_min, t_range[0] / 10)), np.log10(min(t_max, t_range[1] * 10)), 80)
            y_exp = np.interp(t_center, t[mask], y[mask]) * (t_fit / t_center) ** to_float(alpha_exp)
            ax_qty.loglog(t_fit, y_exp, "--", color=color, linewidth=2, label=PHASE_NAMES[phase])
        ax_qty.axvspan(t_range[0], t_range[1], alpha=0.06, color=color)

    ax_qty.set_ylabel(QTY_SYMBOLS.get(qty, qty), fontsize=9)
    ax_qty.legend(loc="best", fontsize=7)
    ax_qty.grid(True, alpha=0.3)
    plt.setp(ax_qty.get_xticklabels(), visible=False)

    t_mid, local_slopes = _smooth_slopes(np.log10(t[valid]), np.log10(y[valid]))
    ax_deriv.semilogx(t_mid, local_slopes, "k-", alpha=0.5, linewidth=1)

    all_vals = []
    for phase, phase_data in phases.items():
        if qty not in phase_data.get("fits", {}):
            continue
        t_range, fit_info = phase_data["t_range"], phase_data["fits"][qty]
        alpha_exp, alpha_fit = fit_info.get("expected"), fit_info["measured"]
        color = PHASE_COLORS.get(phase, "gray")
        t_ext = (max(t_min, t_range[0] / 10), min(t_max, t_range[1] * 10))
        t_text = np.sqrt(t_range[0] * t_range[1])
        ax_deriv.axvspan(t_range[0], t_range[1], alpha=0.08, color=color)
        if alpha_exp is not None:
            alpha_exp_float = to_float(alpha_exp)
            ax_deriv.hlines(alpha_exp_float, *t_ext, colors=color, linestyles="--", linewidth=2)
            ax_deriv.text(t_text, alpha_exp_float - 0.12, format_slope(alpha_exp), color=color,
                          fontsize=7, fontweight="bold", va="top", ha="center")
            all_vals.append(alpha_exp_float)
        if not np.isnan(alpha_fit):
            ax_deriv.plot(t_text, alpha_fit, "o", color=color, markersize=5, zorder=5)
        all_vals.append(alpha_fit)

    ax_deriv.set_ylabel("slope", fontsize=7)
    ax_deriv.set_xlabel("Observer Time [s]", fontsize=7)
    ax_deriv.grid(True, alpha=0.3)
    if all_vals:
        margin = max(0.3, (max(all_vals) - min(all_vals)) * 0.3)
        ax_deriv.set_ylim(min(all_vals) - margin, max(all_vals) + margin)


def _plot_quantities_combined(viz_data: Dict, medium: str, quantities: List[str], title: str) -> plt.Figure:
    fig = plt.figure(figsize=PAGE_PORTRAIT)

    data = viz_data.get(medium, {})
    if not data:
        fig.text(0.5, 0.5, f"No data for {medium}", ha="center", va="center")
        return fig

    t = np.array(data.get("t", []))
    phases = data.get("phases", {})

    positions = [
        (0.12, 0.54, 0.35, 0.32),  # top-left
        (0.54, 0.54, 0.35, 0.32),  # top-right
        (0.12, 0.12, 0.35, 0.32),  # bottom-left
        (0.54, 0.12, 0.35, 0.32),  # bottom-right
    ]

    for idx, qty in enumerate(quantities):
        left, bottom, width, height = positions[idx]
        y = np.array(data.get(qty, []))

        ax_qty = fig.add_axes([left, bottom + height * 0.35, width, height * 0.65])
        ax_deriv = fig.add_axes([left, bottom, width, height * 0.30], sharex=ax_qty)

        _plot_quantity_panel(ax_qty, ax_deriv, t, y, qty, phases)

    fig.suptitle(f"{title}: {medium}", fontsize=12, fontweight="bold", y=0.96)
    return fig


def plot_shock_quantities_combined(viz_data: Dict, medium: str = "ISM") -> plt.Figure:
    return _plot_quantities_combined(viz_data, medium, ["u", "r", "B", "N_p"], "Forward Shock Dynamics")


def plot_frequency_quantities_combined(viz_data: Dict, medium: str = "ISM") -> plt.Figure:
    return _plot_quantities_combined(viz_data, medium, ["nu_m", "nu_c", "nu_a", "nu_M"], "Forward Shock Frequencies")


REGIME_TITLES = {
    "I": r"$\nu_a < \nu_m < \nu_c$",
    "II": r"$\nu_m < \nu_a < \nu_c$",
    "III": r"$\nu_a < \nu_c < \nu_m$",
    "IV": r"$\nu_c < \nu_a < \nu_m$",
    "V": r"$\nu_c < \nu_m < \nu_a$",
}


def _format_beta_expr(expr: str) -> str:
    """Format spectral index expression for display (LaTeX-style)."""
    if not expr:
        return ""
    expr = str(expr)
    if "p" in expr:
        return f"${expr.replace('(p-1)', '(p{-}1)')}$"
    if "/" in expr and len(parts := expr.split("/")) == 2:
        return f"${parts[0]}/{parts[1]}$"
    return expr

SEGMENT_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]


def _plot_spectrum_panel(ax_spec, ax_slope, nu, flux, nu_a, nu_m, nu_c, segments, regime,
                         nu_range=(1e6, 1e20)):
    """Plot spectrum panel with slope derivative, similar to shock dynamics panels."""

    valid = (flux > 0) & np.isfinite(flux) & (nu > 0)
    if np.sum(valid) < 3:
        ax_spec.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax_spec.transAxes)
        ax_spec.set_xlim(*nu_range)
        ax_slope.set_xlim(*nu_range)
        return

    nu_valid, flux_valid = nu[valid], flux[valid]
    sort_idx = np.argsort(nu_valid)
    nu_valid, flux_valid = nu_valid[sort_idx], flux_valid[sort_idx]

    # Upper panel: spectrum
    ax_spec.loglog(nu_valid, flux_valid, "k-", alpha=0.6, linewidth=1.5)

    # Mark characteristic frequencies
    freq_info = [("nu_a", nu_a, "firebrick"), ("nu_m", nu_m, "yellowgreen"), ("nu_c", nu_c, "royalblue")]
    for fname, freq, color in freq_info:
        if freq > 0 and nu_valid.min() < freq < nu_valid.max():
            ax_spec.axvline(freq, color=color, linestyle="--", alpha=0.7, linewidth=1)

    ax_spec.set_ylabel(r"$F_\nu$ [mJy]", fontsize=8)
    ax_spec.grid(True, alpha=0.3)
    ax_spec.tick_params(labelsize=7)
    plt.setp(ax_spec.get_xticklabels(), visible=False)

    # Lower panel: slope (d log F / d log ν)
    log_nu = np.log10(nu_valid)
    log_flux = np.log10(flux_valid)
    nu_mid, local_slopes = _smooth_slopes(log_nu, log_flux)

    ax_slope.semilogx(nu_mid, local_slopes, "k-", alpha=0.5, linewidth=1)

    # Mark expected slopes for each segment
    all_vals = []
    for i, seg in enumerate(segments):
        nu_low, nu_high = seg.get("nu_low"), seg.get("nu_high")
        if nu_low is None or nu_high is None:
            continue
        expected, measured = seg.get("expected"), seg.get("measured")
        color = SEGMENT_COLORS[i % len(SEGMENT_COLORS)]
        nu_center = np.sqrt(nu_low * nu_high)
        for ax in (ax_spec, ax_slope):
            ax.axvspan(nu_low, nu_high, alpha=0.08, color=color)
        if expected is not None:
            ax_slope.hlines(expected, nu_low, nu_high, colors=color, linestyles="--", linewidth=2)
            label = _format_beta_expr(seg.get("expected_expr")) if seg.get("expected_expr") else f"{expected:.2f}"
            ax_slope.text(nu_center, expected + 0.15, label, color=color, fontsize=6,
                          fontweight="bold", va="bottom", ha="center")
            all_vals.append(expected)
        if measured is not None and not np.isnan(measured):
            ax_slope.plot(nu_center, measured, "o", color=color, markersize=4, zorder=5)
            all_vals.append(measured)

    # Mark frequencies on slope panel
    for fname, freq, color in freq_info:
        if freq > 0 and nu_valid.min() < freq < nu_valid.max():
            ax_slope.axvline(freq, color=color, linestyle="--", alpha=0.5, linewidth=1)

    ax_slope.set_ylabel(r"$\beta$", fontsize=8)
    ax_slope.set_xlabel(r"$\nu$ [Hz]", fontsize=8)
    ax_slope.grid(True, alpha=0.3)
    ax_slope.tick_params(labelsize=7)

    # Fix x-axis range
    ax_spec.set_xlim(*nu_range)
    ax_slope.set_xlim(*nu_range)

    if all_vals:
        margin = max(0.5, (max(all_vals) - min(all_vals)) * 0.4)
        ax_slope.set_ylim(min(all_vals) - margin, max(all_vals) + margin)


def plot_spectrum_shapes(viz_data: Dict) -> plt.Figure:
    """Plot spectrum shapes with slope derivative panels for all tested regimes."""
    fig = plt.figure(figsize=PAGE_PORTRAIT)

    regime_names = ["I", "II", "III", "IV", "V"]

    # Layout: 3 rows x 2 columns, each with spectrum + slope subpanels
    positions = [
        (0.08, 0.68, 0.38, 0.22),   # Regime I (top-left)
        (0.54, 0.68, 0.38, 0.22),   # Regime II (top-right)
        (0.08, 0.38, 0.38, 0.22),   # Regime III (middle-left)
        (0.54, 0.38, 0.38, 0.22),   # Regime IV (middle-right)
        (0.08, 0.08, 0.38, 0.22),   # Regime V (bottom-left)
    ]

    for i, regime in enumerate(regime_names):
        key = f"regime_{regime}"
        data = viz_data.get(key, {})

        left, bottom, width, height = positions[i]

        # Create spectrum (upper) and slope (lower) axes
        ax_spec = fig.add_axes([left, bottom + height * 0.4, width, height * 0.6])
        ax_slope = fig.add_axes([left, bottom, width, height * 0.35], sharex=ax_spec)

        if not data:
            ax_spec.text(0.5, 0.5, f"No data", ha="center", va="center", transform=ax_spec.transAxes)
            ax_spec.set_title(f"Regime {regime}", fontsize=9)
            continue

        nu = np.array(data.get("nu", []))
        flux = np.array(data.get("flux", []))
        nu_a = data.get("nu_a", 0)
        nu_m = data.get("nu_m", 0)
        nu_c = data.get("nu_c", 0)
        actual_regime = data.get("actual_regime", regime)
        medium = data.get("medium", "ISM")
        t = data.get("t", 0)
        segments = data.get("expected_segments", [])

        _plot_spectrum_panel(ax_spec, ax_slope, nu, flux, nu_a, nu_m, nu_c, segments, regime)

        title = f"Regime {regime}: {REGIME_TITLES.get(regime, '')}"
        if actual_regime != regime:
            title += f" (actual: {actual_regime})"
        ax_spec.set_title(title, fontsize=9, fontweight="bold")
        ax_spec.text(0.02, 0.95, f"{medium}, t={t:.0e}s", transform=ax_spec.transAxes, fontsize=7, va="top")

    fig.suptitle("Synchrotron Spectrum Shapes: Spectral Index Verification", fontsize=12, fontweight="bold", y=0.97)
    return fig


# --- Reverse Shock Visualization ---


def _rvs_suffix(regime): return f" — {regime.title()} Shell" if regime else ""

def plot_rvs_shock_quantities_combined(viz_data: Dict, medium: str = "ISM", regime: str = "") -> plt.Figure:
    return _plot_quantities_combined(viz_data, medium, ["u", "r", "B", "N_p"], f"Reverse Shock Dynamics{_rvs_suffix(regime)}")

def plot_rvs_frequency_quantities_combined(viz_data: Dict, medium: str = "ISM", regime: str = "") -> plt.Figure:
    return _plot_quantities_combined(viz_data, medium, ["nu_m", "nu_c", "nu_a", "nu_M"], f"Reverse Shock Frequencies{_rvs_suffix(regime)}")
