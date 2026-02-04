"""Benchmark visualization functions for convergence analysis and performance timing."""

import multiprocessing as mp
import os
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from validation.visualization.common import (BAND_COLORS, COLORS, FIDUCIAL_VALUES, MAX_ERROR_THRESHOLD,
                                             MEAN_ERROR_THRESHOLD, PAGE_PORTRAIT, find_fiducial_index,
                                             get_flux_time, setup_plot_style)


def _get_fiducial_idx(dim_key, values):
    """Look up the fiducial index for a convergence dimension."""
    fid = FIDUCIAL_VALUES.get(dim_key.replace("_convergence", ""), 0)
    return next(((i, fid) for i, v in enumerate(values) if abs(v - fid) < 1e-9), (-1, fid))


def _classify_convergence(worst_mean_err, worst_max_err):
    """Classify convergence status based on error thresholds."""
    if worst_mean_err < MEAN_ERROR_THRESHOLD and worst_max_err < MAX_ERROR_THRESHOLD:
        return "PASS", "green"
    if worst_mean_err < MEAN_ERROR_THRESHOLD:
        return "ACCEPTABLE", "blue"
    return "FAIL", "red"


def _get_combinations(configs, jets, media):
    """Get jet/medium combinations for a set of configs."""
    combinations = []
    combo_labels = []
    for jet in jets:
        for medium in media:
            matching = [c for c in configs if c.get("jet_type") == jet and c.get("medium") == medium]
            if matching:
                combinations.append((jet, medium, matching))
                combo_labels.append(f"{jet}/{medium}")
    return combinations, combo_labels


def _plot_stage_breakdown_panel(ax, configs, jets, media):
    """Plot stage breakdown stacked bar chart for a set of configs. Bars = jet × medium."""
    has_stage_data = any(c.get("timing", {}).get("stage_breakdown") for c in configs)

    combinations, combo_labels = _get_combinations(configs, jets, media)

    if not combinations:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes, fontsize=10, color="gray")
        return

    if has_stage_data:
        # Collect all stage names across configs, ordered by total time
        stage_totals = {}
        for _, _, matching in combinations:
            for c in matching:
                for stage, ms in c.get("timing", {}).get("stage_breakdown", {}).items():
                    if stage != "total":
                        stage_totals[stage] = stage_totals.get(stage, 0) + ms
        stage_names = sorted(stage_totals, key=lambda s: stage_totals[s], reverse=True)

        stage_labels = {
            "mesh": "grid generation", "EAT_grid": "EAT grid",
            "shock_dynamics": "shock dynamics", "observer": "observer",
            "syn_electrons": "syn electron gen", "syn_photons": "syn photon gen",
            "cooling": "cooling", "sync_flux": "flux integration",
            "ic_photons": "IC photon gen", "ssc_flux": "SSC flux",
        }
        stage_colors = {
            "mesh": "#2ECC71", "EAT_grid": "#27AE60",
            "shock_dynamics": "#3498DB", "observer": "#1ABC9C",
            "syn_electrons": "#9B59B6", "syn_photons": "#E74C3C", "cooling": "#F39C12",
            "sync_flux": "#E67E22", "ic_photons": "#8E44AD", "ssc_flux": "#C0392B",
        }
        fallback_colors = plt.cm.Set2(np.linspace(0, 1, max(len(stage_names), 1)))

        x = np.arange(len(combinations))
        width = 0.7
        bottom = np.zeros(len(combinations))

        for i, stage in enumerate(stage_names):
            values = []
            for _, _, matching in combinations:
                # Average across bands (per-config stage breakdown), not across configs
                band_times = []
                for c in matching:
                    t = c.get("timing", {}).get("stage_breakdown", {}).get(stage, 0)
                    if t > 0:
                        band_times.append(t)
                values.append(np.mean(band_times) if band_times else 0)
            color = stage_colors.get(stage, fallback_colors[i % len(fallback_colors)])
            label = stage_labels.get(stage, stage)
            ax.bar(x, values, width, bottom=bottom, label=label, color=color, alpha=0.8)
            bottom += np.array(values)

        ax.set_xticks(x)
        ax.set_xticklabels(combo_labels, rotation=45, ha="right")
        ax.set_ylabel("Time [ms]")
        ax.legend(fontsize=6, loc="upper right", ncol=2)
        ax.grid(axis="y", alpha=0.3)
        ax.minorticks_on()
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.grid(axis="y", which="minor", alpha=0.15, linestyle=":")
    else:
        # Fallback: total flux time bar chart
        x = np.arange(len(combinations))
        width = 0.7
        values = []
        for _, _, matching in combinations:
            band_times = [get_flux_time(c) for c in matching if get_flux_time(c) > 0]
            values.append(np.mean(band_times) if band_times else 0)
        ax.bar(x, values, width, color="#3498DB", alpha=0.8)
        for xi, v in zip(x, values):
            if v > 0:
                ax.text(xi, v + 0.5, f"{v:.0f}", ha="center", va="bottom", fontsize=7)
        ax.set_xticks(x)
        ax.set_xticklabels(combo_labels, rotation=45, ha="right")
        ax.set_ylabel("Time [ms]")
        ax.grid(axis="y", alpha=0.3)
        ax.minorticks_on()
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.grid(axis="y", which="minor", alpha=0.15, linestyle=":")


def plot_benchmark_overview(session: Dict, angle_filter: str = "all") -> plt.Figure:
    """Create benchmark overview page with 4 panels, one per radiation config.

    Each panel shows a stage breakdown stacked bar chart with bars = jet type × medium.
    Currently only synchrotron is populated; the other 3 panels are placeholders.
    angle_filter: 'all', 'on-axis', or 'off-axis'.
    """
    fig = plt.figure(figsize=PAGE_PORTRAIT)
    gs = GridSpec(2, 2, figure=fig, hspace=0.40, wspace=0.35, top=0.92, bottom=0.08, left=0.10, right=0.95)

    all_configs = session.get("configs", session.get("results", []))

    if angle_filter == "on-axis":
        configs = [c for c in all_configs if c.get("theta_obs_ratio", 0) <= 1]
        title_suffix = " (On-axis, θ_v=0)"
    elif angle_filter == "off-axis":
        configs = [c for c in all_configs if c.get("theta_obs_ratio", 0) > 1]
        title_suffix = " (Off-axis, θ_v/θ_c>1)"
    else:
        configs = all_configs
        title_suffix = ""

    if not configs:
        fig.text(0.5, 0.5, f"No {angle_filter} data available", ha="center", va="center", fontsize=16, color="gray")
        fig.suptitle(f"Benchmark Overview{title_suffix}", fontsize=14, fontweight="bold", y=0.98)
        return fig

    jets = sorted(set(c.get("jet_type", "unknown") for c in configs))
    media = sorted(set(c.get("medium", "unknown") for c in configs))
    # Define the 4 radiation config panels
    rad_panels = [
        ("synchrotron", "Fwd Synchrotron"),
        ("rvs_sync_thin", "Rvs Thin Shell"),
        ("rvs_sync_thick", "Rvs Thick Shell"),
        ("fwd_ic", "Fwd IC"),
    ]

    for idx, (rad_key, rad_label) in enumerate(rad_panels):
        row, col = divmod(idx, 2)
        ax = fig.add_subplot(gs[row, col])
        ax.set_title(rad_label, fontsize=10, fontweight="bold")

        rad_configs = [c for c in configs if c.get("radiation") == rad_key]
        if rad_configs:
            _plot_stage_breakdown_panel(ax, rad_configs, jets, media)
        else:
            ax.text(0.5, 0.5, "TBD", ha="center", va="center", transform=ax.transAxes,
                    fontsize=14, color="#CCCCCC", fontweight="bold")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_color("#EEEEEE")

    fig.suptitle(f"Benchmark Overview{title_suffix}", fontsize=12, fontweight="bold", y=0.96)
    return fig


def plot_single_model_convergence_page(config: Dict, model_id: int) -> plt.Figure:
    """Create convergence page for a single model (4 rows x 3 cols: errors, times, light curves)."""
    fig = plt.figure(figsize=PAGE_PORTRAIT)
    gs = GridSpec(4, 3, figure=fig, hspace=0.40, wspace=0.35, top=0.91, bottom=0.05, left=0.10, right=0.95)

    jet_type = config.get("jet_type", "?")
    medium = config.get("medium", "?")
    radiation = config.get("radiation", "?")
    theta_obs = config.get("theta_obs", 0)
    theta_obs_ratio = config.get("theta_obs_ratio", theta_obs / 0.1 if theta_obs else 0)

    dimensions = [
        ("phi_convergence", "\u03c6 Resolution", "\u03c6 [ppd]"),
        ("theta_convergence", "\u03b8 Resolution", "\u03b8 [ppd]"),
        ("t_convergence", "t Resolution", "t [ppd]"),
    ]

    all_fiducial_max_errors = []
    all_fiducial_mean_errors = []

    def plot_error_row(row_idx, error_key, ylabel, threshold, collected_errors):
        for col, (dim_key, title, xlabel) in enumerate(dimensions):
            ax = fig.add_subplot(gs[row_idx, col])
            has_any_data = False
            all_valid_errors = []
            fiducial_errors = []

            conv = config.get(dim_key) or {}
            vals = conv.get("values", [])
            errors_by_band = conv.get(error_key, {}) or {}

            fiducial_idx, fiducial_val = _get_fiducial_idx(dim_key, vals)

            if vals:
                for band_name, errors in errors_by_band.items():
                    if not errors:
                        continue
                    color = BAND_COLORS.get(band_name, "gray")
                    plot_errs = []
                    for e in errors:
                        if e is not None and np.isfinite(e):
                            abs_e = abs(e)
                            plot_errs.append(max(abs_e, 1e-10))
                            if abs_e > 0:
                                all_valid_errors.append(abs_e)
                        else:
                            plot_errs.append(np.nan)

                    if fiducial_idx >= 0 and fiducial_idx < len(errors):
                        fid_err = errors[fiducial_idx]
                        if fid_err is not None and np.isfinite(fid_err):
                            fiducial_errors.append(abs(fid_err))

                    if any(np.isfinite(e) for e in plot_errs):
                        has_any_data = True
                        ax.semilogy(vals, plot_errs, marker="o", linestyle="-", color=color, markersize=3,
                                    linewidth=1.0, alpha=0.8, label=band_name)
                        if fiducial_idx >= 0 and fiducial_idx < len(plot_errs):
                            fid_err = plot_errs[fiducial_idx]
                            if np.isfinite(fid_err):
                                ax.semilogy([fiducial_val], [fid_err], marker="*", color=color, markersize=8, zorder=10)

            ax.set_xlabel(xlabel, fontsize=8)
            ax.set_ylabel(ylabel, fontsize=8)
            ax.set_title(title, fontsize=9)
            ax.tick_params(labelsize=7)
            ax.grid(True, alpha=0.3)
            ax.axhline(threshold, color="black", linestyle="--", linewidth=1, alpha=0.7)

            if has_any_data:
                ax.legend(fontsize=6, loc="upper right")
                if all_valid_errors:
                    y_min = min(min(all_valid_errors), threshold * 0.1)
                    y_max = max(max(all_valid_errors), threshold * 10)
                    if y_max / max(y_min, 1e-10) < 100:
                        y_center = np.sqrt(y_min * y_max)
                        y_min = y_center / 10
                        y_max = y_center * 10
                    ax.set_ylim(y_min, y_max)

                if fiducial_errors:
                    max_fid_err = max(fiducial_errors)
                    fid_color = "red" if max_fid_err > threshold else "green"
                    ax.axhline(max_fid_err, color=fid_color, linestyle="-", linewidth=1.0, alpha=0.8)
            else:
                ax.set_ylim(threshold * 0.01, threshold * 100)
                ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes, fontsize=9, color="gray")

            collected_errors.extend(fiducial_errors)

    plot_error_row(0, "errors_by_band", "Max Rel. Error", MAX_ERROR_THRESHOLD, all_fiducial_max_errors)
    plot_error_row(1, "mean_errors_by_band", "Mean Rel. Error", MEAN_ERROR_THRESHOLD, all_fiducial_mean_errors)

    # Row 3: CPU time
    for col, (dim_key, title, xlabel) in enumerate(dimensions):
        ax = fig.add_subplot(gs[2, col])
        has_any_data = False

        conv = config.get(dim_key) or {}
        vals = conv.get("values", [])
        times_by_band = conv.get("times_by_band", {}) or {}

        fiducial_idx, fiducial_val = _get_fiducial_idx(dim_key, vals)

        if not times_by_band:
            times = conv.get("times_ms", [])
            if vals and times:
                times_by_band = {"Optical": times}

        if vals:
            for band_name, times in times_by_band.items():
                if not times:
                    continue
                color = BAND_COLORS.get(band_name, "gray")
                valid_times = [t if (t is not None and np.isfinite(t)) else np.nan for t in times]
                if any(np.isfinite(t) for t in valid_times):
                    has_any_data = True
                    ax.plot(vals, valid_times, linestyle="-", marker="o", color=color, alpha=0.8, linewidth=1.0,
                            markersize=4, label=band_name)
                    if fiducial_idx >= 0 and fiducial_idx < len(valid_times):
                        fid_time = valid_times[fiducial_idx]
                        if np.isfinite(fid_time):
                            ax.plot([fiducial_val], [fid_time], marker="*", color=color, markersize=10, zorder=10)

        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel("CPU Time [ms]", fontsize=8)
        ax.set_title(title, fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)

        if has_any_data:
            ax.legend(fontsize=6, loc="upper left")
        else:
            ax.text(0.5, 0.5, "No timing data", ha="center", va="center", transform=ax.transAxes, fontsize=9, color="gray")

    # Row 4: Light curve panels
    for col, (dim_key, title, xlabel) in enumerate(dimensions):
        ax = fig.add_subplot(gs[3, col])
        has_any_data = False

        conv = config.get(dim_key) or {}
        t_array = conv.get("t_array", [])
        flux_by_band = conv.get("flux_by_band", {}) or {}
        res_values = conv.get("values", [])
        dim_name = dim_key.replace("_convergence", "")
        fiducial_val = FIDUCIAL_VALUES.get(dim_name, 0)

        is_rvs = config.get("radiation", "").startswith("rvs_")
        ref_flux_by_band = conv.get("ref_flux_by_band", {}) or {}
        flux_peak = 0
        flux_floor = np.inf
        if t_array and flux_by_band:
            for band_name, fluxes_list in flux_by_band.items():
                if not fluxes_list:
                    continue
                # For rvs radiation configs, only plot the reverse shock component
                if is_rvs and "(rvs)" not in band_name:
                    continue
                color = BAND_COLORS.get(band_name, "gray")

                # Plot reference resolution curve first (dashed, behind sweep curves)
                ref_flux = ref_flux_by_band.get(band_name)
                if ref_flux and len(ref_flux) == len(t_array):
                    ref_arr = np.asarray(ref_flux)
                    pos = ref_arr[ref_arr > 0]
                    if len(pos) > 0:
                        flux_peak = max(flux_peak, np.max(pos))
                        flux_floor = min(flux_floor, np.min(pos))
                    ax.loglog(t_array, ref_flux, linestyle="--", color=color, alpha=0.4, linewidth=0.8,
                              zorder=1)

                n_res = len(fluxes_list)
                for i, flux in enumerate(fluxes_list):
                    if not flux or len(flux) != len(t_array):
                        continue
                    has_any_data = True
                    flux_arr = np.asarray(flux)
                    pos = flux_arr[flux_arr > 0]
                    if len(pos) > 0:
                        flux_peak = max(flux_peak, np.max(pos))
                        flux_floor = min(flux_floor, np.min(pos))
                    res_val = res_values[i] if i < len(res_values) else fiducial_val
                    ls = "-"
                    label = band_name if i == n_res - 1 else None
                    ax.loglog(t_array, flux, linestyle=ls, color=color, alpha=1, linewidth=0.5, label=label)

        ax.set_xlabel("Time [s]", fontsize=8)
        ax.set_ylabel("Flux [erg/cm\u00b2/s/Hz]", fontsize=8)
        ax.set_title(f"Light Curves ({title.split()[0]})", fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)

        if has_any_data:
            # Clamp y-axis: use data range if tighter than 1e-15 * peak
            if flux_peak > 0:
                y_lo = flux_floor * 0.1 if np.isfinite(flux_floor) else flux_peak * 1e-15
                ax.set_ylim(max(y_lo, flux_peak * 1e-15), flux_peak * 10)

            ax.legend(fontsize=6, loc=0)
        else:
            ax.text(0.5, 0.5, "No light curve data", ha="center", va="center", transform=ax.transAxes,
                    fontsize=9, color="gray")

    # Determine pass/fail status
    worst_max_err = max(all_fiducial_max_errors) if all_fiducial_max_errors else 0
    worst_mean_err = max(all_fiducial_mean_errors) if all_fiducial_mean_errors else 0
    status_text, status_color = _classify_convergence(worst_mean_err, worst_max_err)

    model_info = f"#{model_id}: {jet_type} / {medium} / {radiation} / \u03b8_v/\u03b8_c={theta_obs_ratio:.1f}"
    fig.text(0.08, 0.97, f"[{status_text}]", ha="left", fontsize=11, fontweight="bold", color=status_color)
    fig.text(0.5, 0.97, model_info, ha="center", fontsize=10, fontweight="bold")

    return fig


def _collect_fiducial_errors(configs, error_key, base_band, shock_filter):
    """Collect fiducial errors for a given band and shock type across all configs/dimensions.

    shock_filter: "fwd" collects plain + (fwd) bands, "rvs" collects (rvs) bands.
    Returns list of absolute error values.
    """
    errors = []
    for config in configs:
        for dim_key in ["phi_convergence", "theta_convergence", "t_convergence"]:
            conv = config.get(dim_key) or {}
            errors_by_band = conv.get(error_key, {}) or {}
            tested_values = conv.get("values", [])

            dim_name = dim_key.replace("_convergence", "")
            fiducial_value = FIDUCIAL_VALUES.get(dim_name, 0)
            fiducial_idx = find_fiducial_index(tested_values, fiducial_value)

            for band_name, errs in errors_by_band.items():
                if not (band_name == base_band or band_name.startswith(base_band + " (")):
                    continue
                is_rvs = "(rvs)" in band_name
                is_fwd = "(fwd)" in band_name
                is_plain = not is_rvs and not is_fwd
                if shock_filter == "fwd" and not (is_fwd or is_plain):
                    continue
                if shock_filter == "rvs" and not is_rvs:
                    continue
                if errs and 0 <= fiducial_idx < len(errs):
                    err = errs[fiducial_idx]
                    if err is not None and np.isfinite(err) and abs(err) > 0:
                        errors.append(abs(err))
    return errors


def _plot_error_hist_panel(ax, errors, threshold, color, label):
    """Plot a single error histogram panel."""
    if not errors:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="gray")
        return

    e_min = min(errors) * 0.5
    e_max = max(errors) * 2
    bins = np.logspace(np.log10(max(e_min, 1e-10)), np.log10(e_max), 20)

    ax.hist(errors, bins=bins, alpha=0.6, color=color, edgecolor="white", linewidth=0.5,
            label=f"{label} (n={len(errors)})")
    ax.axvline(threshold, color="black", linestyle="--", linewidth=1.2, alpha=0.8)
    ax.set_xscale("log")

    n_fail = sum(1 for e in errors if e >= threshold)
    if n_fail > 0:
        ax.text(0.97, 0.97, f"{n_fail}/{len(errors)} exceed", ha="right", va="top",
                transform=ax.transAxes, fontsize=7, color="red", fontweight="bold")
    else:
        ax.text(0.97, 0.97, f"all {len(errors)} pass", ha="right", va="top",
                transform=ax.transAxes, fontsize=7, color="green", fontweight="bold")

    ax.legend(fontsize=6, loc="upper left")


def plot_error_distribution(session: Dict) -> plt.Figure:
    """Create error distribution page: histograms of max/mean relative errors by band.

    Layout: 4 rows x 3 columns (Radio, Optical, X-ray).
    Rows: Fwd Max Error, Fwd Mean Error, Rvs Max Error, Rvs Mean Error.
    """
    fig = plt.figure(figsize=PAGE_PORTRAIT)
    gs = GridSpec(4, 3, figure=fig, hspace=0.40, wspace=0.30, top=0.90, bottom=0.05, left=0.12, right=0.95)

    configs = session.get("configs", session.get("results", []))
    if not configs:
        fig.text(0.5, 0.5, "No convergence data", ha="center", va="center", fontsize=16, color="gray")
        fig.suptitle("Error Distribution", fontsize=12, fontweight="bold", y=0.96)
        return fig

    base_bands = ["Radio", "Optical", "X-ray"]

    # 4 rows: (error_key, row_label, shock_filter, threshold)
    rows = [
        ("errors_by_band", "Fwd Max Error", "fwd", MAX_ERROR_THRESHOLD),
        ("mean_errors_by_band", "Fwd Mean Error", "fwd", MEAN_ERROR_THRESHOLD),
        ("errors_by_band", "Rvs Max Error", "rvs", MAX_ERROR_THRESHOLD),
        ("mean_errors_by_band", "Rvs Mean Error", "rvs", MEAN_ERROR_THRESHOLD),
    ]

    for row, (error_key, row_label, shock_filter, threshold) in enumerate(rows):
        for col, base_band in enumerate(base_bands):
            ax = fig.add_subplot(gs[row, col])

            errors = _collect_fiducial_errors(configs, error_key, base_band, shock_filter)
            color_key = f"{base_band} ({shock_filter})" if shock_filter == "rvs" else base_band
            color = BAND_COLORS.get(color_key, BAND_COLORS.get(base_band, "gray"))
            _plot_error_hist_panel(ax, errors, threshold, color, shock_filter)

            ax.set_xlabel("Relative Error" if row == 3 else "", fontsize=8)
            ax.set_ylabel("Count", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.grid(True, alpha=0.3)

            if row == 0:
                ax.set_title(base_band, fontsize=10, fontweight="bold", color=BAND_COLORS.get(base_band, "gray"))

        # Row label on the left
        row_y = 0.90 - (row + 0.5) * (0.85 / 4)
        fig.text(0.02, row_y, row_label, fontsize=8, fontweight="bold", rotation=90, va="center", ha="center")

    fig.suptitle("Error Distribution Across All Configurations", fontsize=12, fontweight="bold", y=0.96)
    return fig


def plot_convergence_summary(session: Dict) -> Tuple[plt.Figure, Dict[str, int], List[Dict]]:
    """Create convergence pass/fail grid. Returns (fig, model_id_map, models_list)."""
    configs = session.get("configs", session.get("results", []))

    if not configs:
        fig = plt.figure(figsize=PAGE_PORTRAIT)
        fig.text(0.5, 0.5, "No convergence data", ha="center", va="center", fontsize=16, color="gray")
        return fig, {}, []

    models = []
    for config in configs:
        jet = config.get("jet_type", "?")
        medium = config.get("medium", "?")
        radiation = config.get("radiation", "?")
        theta_obs = config.get("theta_obs", 0)
        theta_obs_ratio = config.get("theta_obs_ratio", theta_obs / 0.1 if theta_obs else 0)

        dim_results = {}
        max_errors = {}
        mean_errors = {}

        def _worst_at_fid(errors_by_band, fid_idx):
            worst = 0.0
            for errs in errors_by_band.values():
                if errs and 0 <= fid_idx < len(errs) and errs[fid_idx] is not None and np.isfinite(errs[fid_idx]):
                    worst = max(worst, abs(errs[fid_idx]))
            return worst

        for dim_key in ["phi_convergence", "theta_convergence", "t_convergence"]:
            conv = config.get(dim_key) or {}
            dim_name = dim_key.replace("_convergence", "")
            fid_idx = find_fiducial_index(conv.get("values", []), FIDUCIAL_VALUES.get(dim_name, 0))
            max_errors[dim_name] = _worst_at_fid(conv.get("errors_by_band", {}), fid_idx)
            mean_errors[dim_name] = _worst_at_fid(conv.get("mean_errors_by_band", {}), fid_idx)

        worst_max_err = max(max_errors.values()) if max_errors else 0
        worst_mean_err = max(mean_errors.values()) if mean_errors else 0
        status, _ = _classify_convergence(worst_mean_err, worst_max_err)

        has_data = any(config.get(k) for k in ["phi_convergence", "theta_convergence", "t_convergence"])

        models.append({
            "jet": jet, "medium": medium, "radiation": radiation, "theta_obs": theta_obs,
            "theta_obs_ratio": theta_obs_ratio, "status": status, "has_data": has_data,
            "dim_results": dim_results, "max_errors": max_errors, "mean_errors": mean_errors, "config": config,
        })

    model_id_map = {}
    for i, model in enumerate(models):
        model["id"] = i + 1
        key = f"{model['jet']}_{model['medium']}_{model['radiation']}_{model['theta_obs_ratio']:.1f}"
        model_id_map[key] = model["id"]

    fig = plt.figure(figsize=PAGE_PORTRAIT)
    n_models = len(models)

    if n_models <= 12:
        n_cols = 4
    elif n_models <= 24:
        n_cols = 6
    elif n_models <= 48:
        n_cols = 8
    else:
        n_cols = 10

    n_rows = (n_models + n_cols - 1) // n_cols

    grid_top = 0.88
    grid_bottom = 0.08
    grid_height = grid_top - grid_bottom
    cell_height = grid_height / max(n_rows, 1)
    cell_width = 0.90 / n_cols
    x_start = 0.05

    # cell_rects[i] = (x, y, w, h) in figure coordinates for each model cell
    cell_rects = {}

    for i, model in enumerate(models):
        row = i // n_cols
        col = i % n_cols
        x = x_start + col * cell_width
        y = grid_top - (row + 1) * cell_height

        cw, ch = cell_width * 0.95, cell_height * 0.9
        cell_rects[i] = (x, y, cw, ch)

        _status_colors = {"PASS": ("#90EE90", "darkgreen"), "ACCEPTABLE": ("#ADD8E6", "darkblue"), "FAIL": ("#FFB6C1", "darkred")}
        color, text_color = ("#CCCCCC", "black") if not model["has_data"] else _status_colors.get(model["status"], ("#FFB6C1", "darkred"))

        rect = plt.Rectangle((x, y), cw, ch, facecolor=color, edgecolor="white",
                              linewidth=1.5, transform=fig.transFigure, clip_on=False)
        fig.add_artist(rect)

        fig.text(x + cell_width * 0.475, y + cell_height * 0.70, f"#{model['id']}", ha="center", va="center",
                 fontsize=9, fontweight="bold", color=text_color)

        theta_ratio = model.get("theta_obs_ratio", 0)
        fig.text(x + cell_width * 0.475, y + cell_height * 0.52,
                 f"{model['jet'][:5]}/{model['medium'][:3]}/\u03b8{theta_ratio:.0f}",
                 ha="center", va="center", fontsize=6, color=text_color)
        fig.text(x + cell_width * 0.475, y + cell_height * 0.38,
                 f"{model['radiation']}", ha="center", va="center", fontsize=5, color=text_color)

        max_err = max(model["max_errors"].values()) if model["max_errors"] else 0
        err_str = f"\u03b5={max_err:.1e}" if max_err > 0 else "no data"
        fig.text(x + cell_width * 0.475, y + cell_height * 0.22, err_str, ha="center", va="center",
                 fontsize=6, color=text_color)

    n_pass = sum(1 for m in models if m["status"] == "PASS" and m["has_data"])
    n_acceptable = sum(1 for m in models if m["status"] == "ACCEPTABLE" and m["has_data"])
    n_fail = sum(1 for m in models if m["status"] == "FAIL" and m["has_data"])
    n_no_data = sum(1 for m in models if not m["has_data"])

    summary = f"Pass: {n_pass}  |  Acceptable: {n_acceptable}  |  Fail: {n_fail}  |  No Data: {n_no_data}"
    fig.text(0.5, 0.93, summary, ha="center", fontsize=10, color="#333333")
    fig.text(0.5, 0.97, "Resolution Convergence Summary", ha="center", fontsize=12, fontweight="bold")

    return fig, model_id_map, models, cell_rects


def _init_plot_worker():
    matplotlib.use("Agg")
    setup_plot_style()


def _generate_convergence_page_worker(args: Tuple) -> Tuple[int, Optional[str], bool, str]:
    model_idx, model_info, temp_dir = args

    try:
        config = model_info["config"]
        model_id = model_info["id"]
        jet = config.get("jet_type", "?")
        medium = config.get("medium", "?")
        desc = f"#{model_id} {jet}/{medium}"

        has_conv = any(config.get(k) for k in ["phi_convergence", "theta_convergence", "t_convergence"])
        if not has_conv:
            return (model_idx, None, True, desc)

        fig = plot_single_model_convergence_page(config, model_id)
        temp_file = os.path.join(temp_dir, f"conv_page_{model_idx:04d}.pdf")
        fig.savefig(temp_file, format="pdf")
        plt.close(fig)

        return (model_idx, temp_file, True, desc)

    except Exception as e:
        jet = model_info.get("config", {}).get("jet_type", "?")
        medium = model_info.get("config", {}).get("medium", "?")
        return (model_idx, None, False, f"#{model_info.get('id', '?')} {jet}/{medium} - Error: {e}")


def generate_convergence_pages(models_list: List[Dict], n_workers: int = 0) -> List[str]:
    """Generate per-model convergence pages as individual PDF files.

    Returns list of PDF file paths in model order. Caller is responsible
    for embedding these into the final report and cleaning up the temp directory.
    """
    import shutil

    n_models = len(models_list)
    persistent_temp_dir = tempfile.mkdtemp(prefix="convergence_pages_")
    pdf_files = []

    if n_workers > 0 and n_models > 1:
        print(f"  - Per-model convergence pages ({n_models} models, {n_workers} workers)")

        with tempfile.TemporaryDirectory() as temp_dir:
            worker_args = [(i, model_info, temp_dir) for i, model_info in enumerate(models_list)]

            temp_files = {}
            completed = 0
            with mp.Pool(n_workers, initializer=_init_plot_worker) as pool:
                for result in pool.imap_unordered(_generate_convergence_page_worker, worker_args):
                    completed += 1
                    model_idx, temp_file, success, desc = result
                    if success and temp_file:
                        temp_files[model_idx] = temp_file
                    print(f"    [{completed}/{n_models}] {desc}" + (" FAILED" if not success else ""))

            for i in sorted(temp_files.keys()):
                dst = os.path.join(persistent_temp_dir, f"page_{i:04d}.pdf")
                shutil.copy2(temp_files[i], dst)
                pdf_files.append(dst)
    else:
        print(f"  - Per-model convergence pages ({n_models} models)")
        for i, model_info in enumerate(models_list):
            config = model_info["config"]
            model_id = model_info["id"]
            has_conv = any(config.get(k) for k in ["phi_convergence", "theta_convergence", "t_convergence"])
            if has_conv:
                jet = config.get("jet_type", "?")
                medium = config.get("medium", "?")
                print(f"    [{i+1}/{n_models}] #{model_id} {jet}/{medium}")
                fig = plot_single_model_convergence_page(config, model_id)
                dst = os.path.join(persistent_temp_dir, f"page_{i:04d}.pdf")
                fig.savefig(dst, format="pdf")
                plt.close(fig)
                pdf_files.append(dst)

    return pdf_files
