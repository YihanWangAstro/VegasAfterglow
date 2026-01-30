"""Regression test visualization functions for power-law scaling verification."""

import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from validation.visualization.common import (BAND_COLORS, COLORS, PAGE_PORTRAIT, PHASE_COLORS, PHASE_NAMES, QTY_SYMBOLS, format_slope, to_float)
from validation.regression.run_regression import SPECTRAL_REGIMES


def _get_cell_style(result: Dict) -> Tuple[str, str, str, str]:
    if result:
        measured = result.get("measured")
        expected = result.get("expected")
        expected_expr = result.get("expected_expr")
        passed = result.get("passed", False)

        if measured is not None:
            # Use expected_expr if available (e.g., "-(p-1)/2"), else fall back to format_slope
            exp_str = expected_expr if expected_expr is not None else (format_slope(expected) if expected is not None else "")
            if passed:
                return "#90EE90", "darkgreen", f"{measured:.2f}", exp_str
            else:
                return "#FFB6C1", "darkred", f"{measured:.2f}", exp_str
        else:
            return "#CCCCCC", "black", "N/A", ""
    return "#E8E8E8", "#888888", "\u2014", ""


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
    n_pass = 0
    n_total = 0

    for medium_data in grid_data.values():
        for phase_data in medium_data.values():
            for result in phase_data.values():
                if result.get("measured") is not None:
                    n_total += 1
                    if result.get("passed", False):
                        n_pass += 1

    return n_pass, n_total


def _count_spectrum_pass_fail(spectrum_grid: Dict) -> Tuple[int, int]:
    """Count pass/fail for spectrum shape tests."""
    n_pass = 0
    n_total = 0

    for medium_data in spectrum_grid.values():
        for regime_data in medium_data.values():
            for seg_result in regime_data.get("segments", []):
                if seg_result.get("measured") is not None:
                    n_total += 1
                    if seg_result.get("passed", False):
                        n_pass += 1

    return n_pass, n_total


def _draw_grid_section(fig, grid_data, quantities, title, top, bottom, grid_left, cell_width, phases, media):
    """Draw a grid section and return pass/fail counts."""
    n_pass, n_total = _count_pass_fail(grid_data)
    rate = n_pass / n_total * 100 if n_total > 0 else 0
    title_color = COLORS["pass"] if rate >= 80 else COLORS["fail"]
    cell_height = (top - bottom) / len(quantities)

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

    return n_pass, n_total


def _format_segment_label(seg_name: str) -> str:
    """Convert segment name to readable label like 'ν < ν_a' or 'ν_a < ν < ν_m'."""
    freq_symbols = {"nu_a": r"$\nu_a$", "nu_m": r"$\nu_m$", "nu_c": r"$\nu_c$"}

    if seg_name.startswith("below_"):
        freq = seg_name.replace("below_", "")
        return r"$\nu<$" + freq_symbols.get(freq, freq)
    elif seg_name.startswith("above_"):
        freq = seg_name.replace("above_", "")
        return r"$\nu>$" + freq_symbols.get(freq, freq)
    elif "_to_" in seg_name:
        parts = seg_name.split("_to_")
        if len(parts) == 2:
            f1, f2 = parts
            return freq_symbols.get(f1, f1) + r"$<\nu<$" + freq_symbols.get(f2, f2)
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

        # Find regime data (could be in ISM or Wind)
        regime_data = None
        for medium_data in spectrum_grid.values():
            if regime in medium_data:
                regime_data = medium_data[regime]
                break

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
            x = grid_left + j * cell_width
            y = row_bottom

            seg_result = None
            if j < len(expected_segments):
                expected_seg_name = expected_segments[j][0]
                seg_result = seg_lookup.get(expected_seg_name)
                if seg_result:
                    color, text_color, line1, line2 = _get_cell_style(seg_result)
                else:
                    color, text_color, line1, line2 = "#E8E8E8", "#888888", "\u2014", ""
            else:
                color, text_color, line1, line2 = "#E8E8E8", "#888888", "\u2014", ""

            _draw_cell(fig, x, y, cell_width, cell_height, color, text_color, line1, line2, fontsize=7,
                       measured=to_float(seg_result.get("measured")) if seg_result else None,
                       expected=to_float(seg_result.get("expected")) if seg_result else None,
                       tolerance=to_float(seg_result.get("tolerance")) if seg_result else None,
                       passed=seg_result.get("passed") if seg_result else None)

    return n_pass, n_total


def plot_dynamics_summary_grid(results: Dict) -> plt.Figure:
    """Plot dynamics summary: FWD + RVS thin + RVS thick shock dynamics grids."""
    fig = plt.figure(figsize=PAGE_PORTRAIT)

    fwd_phases = ["coasting", "BM", "deep_newtonian"]
    rvs_phases = ["crossing", "BM", "deep_newtonian"]
    media = ["ISM", "Wind"]
    quantities = ["u", "r", "B", "N_p"]

    grid_left = 0.18
    cell_width = (0.92 - grid_left) / (len(fwd_phases) * len(media))

    total_pass, total_tests = 0, 0

    grids = [
        (results.get("shock_grid", {}), "Forward Shock", fwd_phases, 0.86, 0.66),
        (results.get("rvs_shock_grid_thin", {}), "Reverse Shock — Thin Shell", rvs_phases, 0.56, 0.36),
        (results.get("rvs_shock_grid_thick", {}), "Reverse Shock — Thick Shell", rvs_phases, 0.26, 0.06),
    ]

    for grid_data, title, phases, top, bottom in grids:
        if grid_data:
            n_p, n_t = _draw_grid_section(
                fig, grid_data, quantities, title, top=top, bottom=bottom,
                grid_left=grid_left, cell_width=cell_width, phases=phases, media=media)
            total_pass += n_p
            total_tests += n_t

    total_rate = total_pass / total_tests * 100 if total_tests > 0 else 0
    status_color = COLORS["pass"] if total_rate >= 80 else COLORS["fail"]
    fig.suptitle(f"Dynamics Summary: {total_pass}/{total_tests} passed ({total_rate:.0f}%)",
                 fontsize=13, fontweight="bold", y=0.96, color=status_color)
    return fig


def plot_spectrum_summary_grid(results: Dict) -> plt.Figure:
    """Plot standalone spectrum shapes pass/fail grid."""
    fig = plt.figure(figsize=PAGE_PORTRAIT)

    spectrum_grid = results.get("spectrum_grid", {})
    grid_left, grid_right = 0.18, 0.92

    n_pass, n_total = 0, 0
    if spectrum_grid:
        n_pass, n_total = _draw_spectrum_grid_section(
            fig, spectrum_grid, top=0.78, bottom=0.40, grid_left=grid_left, grid_right=grid_right,
            section_title="Synchrotron Spectrum Shapes",
        )

    total_rate = n_pass / n_total * 100 if n_total > 0 else 0
    status_color = COLORS["pass"] if total_rate >= 80 else COLORS["fail"]
    fig.suptitle(f"Synchrotron Spectrum Shapes: {n_pass}/{n_total} passed ({total_rate:.0f}%)",
                 fontsize=13, fontweight="bold", y=0.92, color=status_color)
    return fig


def plot_frequencies_summary_grid(results: Dict) -> plt.Figure:
    """Plot frequencies summary: FWD + RVS thin + RVS thick frequency grids."""
    fig = plt.figure(figsize=PAGE_PORTRAIT)

    fwd_phases = ["coasting", "BM", "deep_newtonian"]
    rvs_phases = ["crossing", "BM", "deep_newtonian"]
    media = ["ISM", "Wind"]
    quantities = ["nu_m", "nu_c", "nu_M"]

    grid_left = 0.18
    cell_width = (0.92 - grid_left) / (len(fwd_phases) * len(media))

    total_pass, total_tests = 0, 0

    grids = [
        (results.get("freq_grid", {}), "Forward Shock", fwd_phases, 0.86, 0.70),
        (results.get("rvs_freq_grid_thin", {}), "Reverse Shock — Thin Shell", rvs_phases, 0.56, 0.40),
        (results.get("rvs_freq_grid_thick", {}), "Reverse Shock — Thick Shell", rvs_phases, 0.26, 0.10),
    ]

    for grid_data, title, phases, top, bottom in grids:
        if grid_data:
            n_p, n_t = _draw_grid_section(
                fig, grid_data, quantities, title, top=top, bottom=bottom,
                grid_left=grid_left, cell_width=cell_width, phases=phases, media=media)
            total_pass += n_p
            total_tests += n_t

    total_rate = total_pass / total_tests * 100 if total_tests > 0 else 0
    status_color = COLORS["pass"] if total_rate >= 80 else COLORS["fail"]
    fig.suptitle(f"Frequencies Summary: {total_pass}/{total_tests} passed ({total_rate:.0f}%)",
                 fontsize=13, fontweight="bold", y=0.96, color=status_color)
    return fig


def _smooth_1d(arr, size):
    """Uniform 1D smoothing with nearest-edge padding (replaces scipy.ndimage.uniform_filter1d)."""
    pad = size // 2
    padded = np.pad(arr, pad, mode="edge")
    kernel = np.ones(size) / size
    return np.convolve(padded, kernel, mode="valid")


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

        t_range = phase_data["t_range"]
        fit_info = phase_data["fits"][qty]
        alpha_fit = fit_info["measured"]
        alpha_exp = fit_info.get("expected")

        mask = (t >= t_range[0]) & (t <= t_range[1]) & valid
        if np.sum(mask) < 2:
            continue

        color = PHASE_COLORS.get(phase, "gray")
        t_center = np.sqrt(t_range[0] * t_range[1])
        y_center = np.interp(t_center, t[mask], y[mask])

        t_ext_min = max(t_min, t_range[0] / 10)
        t_ext_max = min(t_max, t_range[1] * 10)
        t_fit = np.logspace(np.log10(t_ext_min), np.log10(t_ext_max), 80)

        if alpha_exp is not None:
            y_exp = y_center * (t_fit / t_center) ** to_float(alpha_exp)
            label = f"{PHASE_NAMES[phase]}"
            ax_qty.loglog(t_fit, y_exp, "--", color=color, linewidth=2, label=label)

        ax_qty.axvspan(t_range[0], t_range[1], alpha=0.06, color=color)

    ax_qty.set_ylabel(QTY_SYMBOLS.get(qty, qty), fontsize=9)
    ax_qty.legend(loc="best", fontsize=7)
    ax_qty.grid(True, alpha=0.3)
    plt.setp(ax_qty.get_xticklabels(), visible=False)

    log_t = np.log10(t[valid])
    log_y = np.log10(y[valid])
    d_log_t = np.diff(log_t)
    d_log_y = np.diff(log_y)
    local_slopes = d_log_y / d_log_t
    t_mid = 10 ** (0.5 * (log_t[:-1] + log_t[1:]))

    ax_deriv.semilogx(t_mid, local_slopes, "k-", alpha=0.5, linewidth=1)

    all_vals = []
    for phase, phase_data in phases.items():
        if qty not in phase_data.get("fits", {}):
            continue

        t_range = phase_data["t_range"]
        fit_info = phase_data["fits"][qty]
        alpha_exp = fit_info.get("expected")
        alpha_fit = fit_info["measured"]
        color = PHASE_COLORS.get(phase, "gray")

        t_ext_min = max(t_min, t_range[0] / 10)
        t_ext_max = min(t_max, t_range[1] * 10)

        ax_deriv.axvspan(t_range[0], t_range[1], alpha=0.08, color=color)
        t_text = np.sqrt(t_range[0] * t_range[1])

        if alpha_exp is not None:
            alpha_exp_float = to_float(alpha_exp)
            ax_deriv.hlines(alpha_exp_float, t_ext_min, t_ext_max, colors=color, linestyles="--", linewidth=2)
            slope_str = format_slope(alpha_exp)
            ax_deriv.text(t_text, alpha_exp_float - 0.12, slope_str, color=color, fontsize=7,
                          fontweight="bold", va="top", ha="center")
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
    if expr is None:
        return ""
    expr = str(expr)
    # Handle p-related expressions
    if "p" in expr:
        # Convert -(p-1)/2 -> $-(p-1)/2$
        expr = expr.replace("(p-1)", "(p{-}1)")
        return f"${expr}$"
    # Handle fractions like "1/3", "-1/2", "5/2"
    if "/" in expr:
        parts = expr.split("/")
        if len(parts) == 2:
            num, den = parts
            return f"${num}/{den}$"
    # Plain number
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
    d_log_nu = np.diff(log_nu)
    d_log_flux = np.diff(log_flux)
    local_slopes = d_log_flux / d_log_nu
    nu_mid = 10 ** (0.5 * (log_nu[:-1] + log_nu[1:]))

    ax_slope.semilogx(nu_mid, local_slopes, "k-", alpha=0.5, linewidth=1)

    # Mark expected slopes for each segment
    all_vals = []
    for i, seg in enumerate(segments):
        expected = seg.get("expected")
        expected_expr = seg.get("expected_expr")  # Original expression (fraction or p-related)
        measured = seg.get("measured")
        color = SEGMENT_COLORS[i % len(SEGMENT_COLORS)]

        # Use stored frequency ranges (preferred) or fall back to parsing segment name
        nu_low = seg.get("nu_low")
        nu_high = seg.get("nu_high")
        if nu_low is None or nu_high is None:
            continue

        # Shade segment region
        ax_spec.axvspan(nu_low, nu_high, alpha=0.08, color=color)
        ax_slope.axvspan(nu_low, nu_high, alpha=0.08, color=color)

        nu_center = np.sqrt(nu_low * nu_high)

        if expected is not None:
            ax_slope.hlines(expected, nu_low, nu_high, colors=color, linestyles="--", linewidth=2)
            # Display expression if available, otherwise decimal
            label = _format_beta_expr(expected_expr) if expected_expr else f"{expected:.2f}"
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


def plot_rvs_shock_quantities_combined(viz_data: Dict, medium: str = "ISM", regime: str = "") -> plt.Figure:
    """Plot reverse shock dynamics (mirrors plot_shock_quantities_combined)."""
    suffix = f" — {regime.title()} Shell" if regime else ""
    return _plot_quantities_combined(viz_data, medium, ["u", "r", "B", "N_p"], f"Reverse Shock Dynamics{suffix}")


def plot_rvs_frequency_quantities_combined(viz_data: Dict, medium: str = "ISM", regime: str = "") -> plt.Figure:
    """Plot reverse shock frequencies (mirrors plot_frequency_quantities_combined)."""
    suffix = f" — {regime.title()} Shell" if regime else ""
    return _plot_quantities_combined(viz_data, medium, ["nu_m", "nu_c", "nu_a", "nu_M"], f"Reverse Shock Frequencies{suffix}")
