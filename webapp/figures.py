"""Plotly figure builders and trace helpers."""

from collections import OrderedDict

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from VegasAfterglow.cli import _freq_colors

from .constants import (
    AXIS_COMMON,
    COMP_LABELS,
    FBAND_TITLE,
    FLUX_SCALES,
    FREQ_SCALES,
    INSTRUMENTS,
    LAYOUT_COMMON,
    LEGEND_COMMON,
    OBS_COLORS,
    PLOTLY_DASH,
    PLOTLY_FLUX_LABELS,
    TIME_LABELS,
    TIME_SCALES,
)
from .helpers import (
    band_label,
    cgs_to_ab_mag,
    format_time_label,
    freq_label,
    mag_to_cgs,
    smart_ylim,
    time_colors,
)


_SUPER = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def _pow10(v):
    """Format a number as '1.23×10⁴' using Unicode superscripts."""
    if v == 0 or not np.isfinite(v):
        return "0"
    exp = int(np.floor(np.log10(np.abs(v))))
    coeff = v / 10**exp
    sup = str(exp).translate(_SUPER)
    if abs(coeff - 1.0) < 0.005:
        return f"10{sup}"
    return f"{coeff:.2f}×10{sup}"


def _add_traces(fig, x, groups, x_name, x_unit, secondary_y=None):
    """Add trace groups to a figure.

    groups: list of (label, color, components, hover_y_unit)
    components: list of (comp_name, y_array)
    """
    for label, color, comps, hover_y_unit in groups:
        for comp_name, y in comps:
            is_total = comp_name == "total"
            name = label if is_total else f"{label} ({COMP_LABELS.get(comp_name, comp_name)})"
            mask = np.isfinite(y)
            x_m, y_m = x[mask], y[mask]
            hover = [
                f"{_pow10(xv)} {x_unit}<br>{_pow10(yv)} {hover_y_unit}"
                for xv, yv in zip(x_m, y_m)
            ]
            fig.add_trace(
                go.Scatter(
                    x=x_m, y=y_m, mode="lines", name=name,
                    line=dict(
                        color=color,
                        width=1.2 if is_total else 0.9,
                        dash=PLOTLY_DASH.get(comp_name, "solid"),
                    ),
                    opacity=1.0 if is_total else 0.75,
                    legendgroup=name,
                    hoverlabel=dict(bgcolor="white", font_color="black",
                                    bordercolor="#ccc"),
                    hovertemplate="%{text}<extra>" + name + "</extra>",
                    text=hover,
                ),
                secondary_y=secondary_y,
            )


def _add_sensitivity_traces(fig, instruments, mode, freq_scale=1.0,
                            flux_scale=1.0, x_range=None, nufnu=False,
                            has_secondary=False):
    """Overlay instrument sensitivity limits on a figure."""
    fnu_ys, fband_ys = [], []

    for name in instruments:
        nu_min, nu_max, sens, kind, color = INSTRUMENTS[name]
        is_fnu = (kind == "Fnu")

        if mode == "lightcurve":
            if not x_range:
                continue
            xs = x_range
        else:
            xs = [nu_min / freq_scale, nu_max / freq_scale]

        if is_fnu:
            y_val = sens / flux_scale
            if mode == "spectrum" and nufnu:
                y_val = np.sqrt(nu_min * nu_max) * sens
            fnu_ys.append(y_val)
        else:
            y_val = sens
            fband_ys.append(y_val)

        fig.add_trace(
            go.Scatter(
                x=xs, y=[y_val, y_val], mode="lines",
                line=dict(color=color, width=1, dash="dash" if is_fnu else "solid"),
                name=name, legendgroup="inst",
                legendgrouptitle_text="Instruments",
                showlegend=True,
                hovertemplate=(
                    f"{name}<br>F<sub>\u03bd</sub>=%{{y:.1e}}<extra></extra>" if is_fnu
                    else f"{name}<br>F=%{{y:.1e}} erg/cm\u00b2/s<extra></extra>"
                ),
            ),
            secondary_y=(False if has_secondary else None) if is_fnu else True,
        )

    if fnu_ys:
        min_s = min(fnu_ys)
        cur = fig.layout.yaxis.range
        if cur and min_s > 0 and np.log10(min_s) < cur[0]:
            fig.update_layout(yaxis_range=[np.log10(min_s * 0.3), cur[1]])
    if fband_ys and has_secondary:
        min_s = min(fband_ys)
        cur2 = fig.layout.yaxis2.range if fig.layout.yaxis2 else None
        if cur2 and min_s > 0 and np.log10(min_s) < cur2[0]:
            fig.update_layout(yaxis2_range=[np.log10(min_s * 0.3), cur2[1]])


def _add_obs_traces(fig, obs_data, flux_unit, x_unit, has_secondary,
                    mode="lightcurve", nufnu=False):
    """Overlay observational data points with error bars.

    obs_data: tuple of (label, x_val, x_unit_row, y_val, err_val, y_unit) rows.
    mode: "lightcurve" (x = time) or "spectrum" (x = frequency).
    """
    is_lc = (mode == "lightcurve")
    x_scales = TIME_SCALES if is_lc else FREQ_SCALES

    fnu_groups = OrderedDict()
    fband_groups = OrderedDict()

    for label, x_val, x_unit_row, y_val, err_val, y_unit in obs_data:
        if not (np.isfinite(x_val) and np.isfinite(y_val) and x_val > 0):
            continue
        x_phys = x_val * x_scales.get(x_unit_row, 1.0)
        err_val = abs(err_val) if np.isfinite(err_val) else 0.0
        if y_unit == "erg/cm\u00b2/s":
            fband_groups.setdefault(label, []).append((x_phys, y_val, err_val))
        elif y_unit == "AB mag":
            f_cgs, e_cgs = mag_to_cgs(y_val, err_val)
            fnu_groups.setdefault(label, []).append((x_phys, f_cgs, e_cgs))
        else:
            scale = FLUX_SCALES.get(y_unit, 1.0)
            fnu_groups.setdefault(label, []).append((x_phys, y_val * scale, err_val * scale))

    x_disp_scale = (TIME_SCALES if is_lc else FREQ_SCALES)[x_unit]
    is_mag = (flux_unit == "AB mag")
    f_scale = 1.0 if is_mag else FLUX_SCALES[flux_unit]
    if is_lc:
        x_hover = f"t=%{{x:.2e}} {TIME_LABELS[x_unit]}"
    else:
        x_hover = f"\u03bd=%{{x:.2e}} {x_unit}"

    all_labels = list(dict.fromkeys(list(fnu_groups) + list(fband_groups)))
    all_fnu_ys = []

    for idx, label in enumerate(all_labels):
        color = OBS_COLORS[idx % len(OBS_COLORS)]

        if label in fnu_groups:
            pts = fnu_groups[label]
            xs = [r[0] / x_disp_scale for r in pts]
            if is_mag:
                ys = [float(cgs_to_ab_mag(np.array([r[1]]))) for r in pts]
                errs = [r[2] / r[1] * 2.5 / np.log(10) if r[1] > 0 else 0.0 for r in pts]
                hover_y = "mag=%{y:.2f}"
            elif not is_lc and nufnu:
                ys = [r[0] * r[1] for r in pts]
                errs = [r[0] * r[2] for r in pts]
                hover_y = "\u03bdF\u03bd=%{y:.2e} erg/cm\u00b2/s"
            else:
                ys = [r[1] / f_scale for r in pts]
                errs = [r[2] / f_scale for r in pts]
                hover_y = f"F\u03bd=%{{y:.2e}} {PLOTLY_FLUX_LABELS[flux_unit]}"
            all_fnu_ys.extend(ys)
            fig.add_trace(
                go.Scatter(
                    x=xs, y=ys, mode="markers",
                    name=label, legendgroup=f"obs_{label}",
                    marker=dict(color=color, size=6, symbol="circle"),
                    error_y=dict(type="data", array=errs, visible=True,
                                 color=color, thickness=1.0, width=3),
                    hovertemplate=f"{x_hover}<br>{hover_y}<extra>{label}</extra>",
                ),
                secondary_y=False if has_secondary else None,
            )

        if label in fband_groups and has_secondary:
            pts = fband_groups[label]
            xs = [r[0] / x_disp_scale for r in pts]
            ys = [r[1] for r in pts]
            errs = [r[2] for r in pts]
            fig.add_trace(
                go.Scatter(
                    x=xs, y=ys, mode="markers",
                    name=label, legendgroup=f"obs_{label}",
                    marker=dict(color=color, size=6, symbol="diamond"),
                    error_y=dict(type="data", array=errs, visible=True,
                                 color=color, thickness=1.0, width=3),
                    hovertemplate=f"{x_hover}<br>F=%{{y:.2e}} erg/cm\u00b2/s<extra>{label}</extra>",
                ),
                secondary_y=True,
            )

    if is_mag:
        finite = [v for v in all_fnu_ys if np.isfinite(v)]
        if finite:
            cur = fig.layout.yaxis.range
            if cur:
                new_hi = max(cur[0], max(finite) + 1)  # faint end (big number)
                new_lo = min(cur[1], min(finite) - 1)  # bright end (small number)
                fig.update_layout(yaxis_range=[new_hi, new_lo])
    else:
        pos = [v for v in all_fnu_ys if v > 0]
        if pos:
            cur = fig.layout.yaxis.range
            if cur:
                new_lo = min(cur[0], np.log10(min(pos) * 0.3))
                new_hi = max(cur[1], np.log10(max(pos) * 3.0))
                fig.update_layout(yaxis_range=[new_lo, new_hi])

    all_fband_ys = [r[1] for pts in fband_groups.values() for r in pts if r[1] > 0]
    if all_fband_ys and has_secondary:
        cur2 = fig.layout.yaxis2.range if fig.layout.yaxis2 else None
        if cur2:
            new_lo = min(cur2[0], np.log10(min(all_fband_ys) * 0.3))
            new_hi = max(cur2[1], np.log10(max(all_fband_ys) * 3.0))
            fig.update_layout(yaxis2_range=[new_lo, new_hi])


def make_figure(data, flux_unit, time_unit, t_min, t_max,
                need_secondary_y=False):
    """Build an interactive Plotly figure from computation results."""
    times = data["times"]
    freqs = data["frequencies"]
    pt_components = data["pt_components"]
    band_data = data["band_data"]

    is_mag = (flux_unit == "AB mag")
    f_scale = 1.0 if is_mag else FLUX_SCALES[flux_unit]
    t_scale = TIME_SCALES[time_unit]
    t_disp = times / t_scale

    has_points = len(pt_components) > 0 and freqs.size > 0
    has_bands = len(band_data) > 0
    dual = has_points and has_bands
    use_secondary = dual or need_secondary_y

    all_freqs = list(freqs) if has_points else []
    all_freqs += [np.sqrt(b[0] * b[1]) for b in band_data]
    all_colors = _freq_colors(np.array(all_freqs)) if all_freqs else []
    pt_colors = all_colors[:len(freqs)] if has_points else []
    bd_colors = all_colors[len(freqs):] if has_bands else []

    if use_secondary:
        fig = make_subplots(specs=[[{"secondary_y": True}]])
    else:
        fig = go.Figure()

    t_unit_label = TIME_LABELS[time_unit]

    if has_points:
        if is_mag:
            pt_groups = [
                (freq_label(nu), pt_colors[i],
                 [(cn, cgs_to_ab_mag(cf[i])) for cn, cf in pt_components],
                 "mag")
                for i, nu in enumerate(freqs)
            ]
        else:
            pt_groups = [
                (freq_label(nu), pt_colors[i],
                 [(cn, cf[i] / f_scale) for cn, cf in pt_components],
                 PLOTLY_FLUX_LABELS[flux_unit])
                for i, nu in enumerate(freqs)
            ]
        _add_traces(fig, t_disp, pt_groups, "t", t_unit_label,
                     secondary_y=False if use_secondary else None)

    if has_bands:
        bd_groups = [
            (band_label(nu_min, nu_max, blabel), bd_colors[b_idx],
             comps, "erg/cm\u00b2/s")
            for b_idx, (nu_min, nu_max, blabel, comps) in enumerate(band_data)
        ]
        if dual:
            _add_traces(fig, t_disp, bd_groups, "t", t_unit_label,
                         secondary_y=True)
        else:
            _add_traces(fig, t_disp, bd_groups, "t", t_unit_label,
                         secondary_y=False if use_secondary else None)

    x_lo = np.log10(t_min / t_scale)
    x_hi = np.log10(t_max / t_scale)

    if has_points:
        if is_mag:
            all_mag = np.concatenate([cgs_to_ab_mag(cf[i]) for _, cf in pt_components
                                      for i in range(len(freqs))])
            finite = all_mag[np.isfinite(all_mag)]
            if finite.size > 0:
                y_range_pt = [np.max(finite) + 1, np.min(finite) - 1]  # reversed
            else:
                y_range_pt = None
        else:
            pt_flux = [cf / f_scale for _, cf in pt_components]
            y_bot, y_top = smart_ylim(pt_flux)
            y_range_pt = [np.log10(y_bot), np.log10(y_top)] if y_bot else None
    else:
        y_range_pt = None

    if has_bands:
        bd_flux = []
        for _, _, _, comps in band_data:
            bd_flux.extend([cf for _, cf in comps])
        bd_bot, bd_top = smart_ylim(bd_flux)
    else:
        bd_bot, bd_top = None, None

    x_range = [x_lo, x_hi]

    x_axis = dict(
        type="log",
        title=f"t<sub>obs</sub> ({TIME_LABELS[time_unit]})",
        range=x_range,
        **AXIS_COMMON,
    )

    if is_mag:
        y_axis_pt = dict(
            type="linear",
            title="AB mag",
            range=y_range_pt,
            autorange="reversed" if y_range_pt is None else None,
            **AXIS_COMMON,
        )
    else:
        y_axis_pt = dict(
            type="log",
            title=f"F<sub>\u03bd</sub> ({PLOTLY_FLUX_LABELS[flux_unit]})",
            range=y_range_pt,
            **AXIS_COMMON,
        )

    if dual:
        y_range_bd = [np.log10(bd_bot), np.log10(bd_top)] if bd_bot else None
        y_axis_bd = {
            "type": "log", "title": FBAND_TITLE, "range": y_range_bd,
            **AXIS_COMMON, "showgrid": False,
        }
        fig.update_layout(xaxis=x_axis, yaxis=y_axis_pt, yaxis2=y_axis_bd)
    elif has_bands and not has_points:
        y_range_bd = [np.log10(bd_bot), np.log10(bd_top)] if bd_bot else None
        y_axis_bd = dict(
            type="log", title=FBAND_TITLE, range=y_range_bd, **AXIS_COMMON,
        )
        if use_secondary:
            y_axis_inst = {
                "type": "log", "title": FBAND_TITLE,
                **AXIS_COMMON, "showgrid": False,
            }
            fig.update_layout(xaxis=x_axis, yaxis=y_axis_bd, yaxis2=y_axis_inst)
        else:
            fig.update_layout(xaxis=x_axis, yaxis=y_axis_bd)
    elif use_secondary:
        y_axis_bd = {
            "type": "log", "title": FBAND_TITLE,
            **AXIS_COMMON, "showgrid": False,
        }
        fig.update_layout(xaxis=x_axis, yaxis=y_axis_pt, yaxis2=y_axis_bd)
    else:
        fig.update_layout(xaxis=x_axis, yaxis=y_axis_pt)

    fig.update_layout(legend=LEGEND_COMMON, **LAYOUT_COMMON)

    return fig


def make_sed_figure(data, flux_unit, freq_unit, nufnu=False,
                    need_secondary_y=False):
    """Build an interactive Plotly spectrum/SED figure from computation results."""
    t_snapshots = data["t_snapshots"]
    freqs = data["frequencies"]
    components = data["components"]

    is_mag = (flux_unit == "AB mag")
    f_scale = 1.0 if is_mag else FLUX_SCALES[flux_unit]
    nu_scale = FREQ_SCALES[freq_unit]
    nu_disp = freqs / nu_scale

    colors = time_colors(t_snapshots)

    if need_secondary_y:
        fig = make_subplots(specs=[[{"secondary_y": True}]])
    else:
        fig = go.Figure()

    if is_mag:
        y_label = "AB mag"
        hover_unit = "mag"
    elif nufnu:
        y_label = f"\u03bdF<sub>\u03bd</sub> (erg cm<sup>\u22122</sup> s<sup>\u22121</sup>)"
        hover_unit = "erg/cm\u00b2/s"
    else:
        y_label = f"F<sub>\u03bd</sub> ({PLOTLY_FLUX_LABELS[flux_unit]})"
        hover_unit = PLOTLY_FLUX_LABELS[flux_unit]

    groups = []
    for j, t_snap in enumerate(t_snapshots):
        if is_mag:
            comps = [(cn, cgs_to_ab_mag(cf[:, j])) for cn, cf in components]
        elif nufnu:
            comps = [(cn, freqs * cf[:, j]) for cn, cf in components]
        else:
            comps = [(cn, cf[:, j] / f_scale) for cn, cf in components]
        groups.append((format_time_label(t_snap), colors[j], comps, hover_unit))
    _add_traces(fig, nu_disp, groups, "\u03bd", freq_unit,
                secondary_y=False if need_secondary_y else None)

    nu_min, nu_max = freqs[0], freqs[-1]
    x_lo = np.log10(nu_min / nu_scale)
    x_hi = np.log10(nu_max / nu_scale)

    if is_mag:
        all_mag = np.concatenate([cgs_to_ab_mag(cf[:, j])
                                  for _, cf in components
                                  for j in range(len(t_snapshots))])
        finite = all_mag[np.isfinite(all_mag)]
        if finite.size > 0:
            y_range = [np.max(finite) + 1, np.min(finite) - 1]  # reversed
        else:
            y_range = None
    else:
        if nufnu:
            all_flux = [freqs[:, None] * cf for _, cf in components]
        else:
            all_flux = [cf / f_scale for _, cf in components]
        y_bot, y_top = smart_ylim(all_flux)
        y_range = [np.log10(y_bot), np.log10(y_top)] if y_bot else None

    x_axis = dict(
        type="log",
        title=f"\u03bd ({freq_unit})",
        range=[x_lo, x_hi],
        **AXIS_COMMON,
    )
    if is_mag:
        y_axis = dict(
            type="linear",
            title=y_label,
            range=y_range,
            autorange="reversed" if y_range is None else None,
            **AXIS_COMMON,
        )
    else:
        y_axis = dict(
            type="log",
            title=y_label,
            range=y_range,
            **AXIS_COMMON,
        )
    fig.update_layout(xaxis=x_axis, yaxis=y_axis)
    if need_secondary_y:
        y_axis_bd = {
            "type": "log",
            "title": FBAND_TITLE,
            **AXIS_COMMON,
            "showgrid": False,
        }
        fig.update_layout(yaxis2=y_axis_bd)
    fig.update_layout(legend=LEGEND_COMMON, **LAYOUT_COMMON)

    return fig


def _skymap_heatmap(image, x_centers, y_centers, zmin=None, zmax=None):
    """Create a Heatmap trace for a single sky image frame."""
    img_plot = np.where(image > 0, image, np.nan)
    z = np.log10(img_plot.T)
    kw = {}
    if zmin is not None:
        kw["zmin"] = zmin
        kw["zmax"] = zmax
    return go.Heatmap(
        z=z, x=x_centers, y=y_centers,
        colorscale="Inferno",
        colorbar=dict(
            title=dict(text="log<sub>10</sub>(I)", side="right",
                       font=dict(size=11, color="#000")),
            tickfont=dict(size=9, color="#000"),
        ),
        hovertemplate=(
            "\u0394x=%{x:.1f} \u03bcas<br>"
            "\u0394y=%{y:.1f} \u03bcas<br>"
            "log\u2081\u2080(I)=%{z:.2f}<extra></extra>"
        ),
        **kw,
    )


def _skymap_layout(extent, t_label, nu_label):
    """Common layout dict for sky map figures."""
    return dict(
        xaxis=dict(
            title="\u0394x (\u03bcas)", scaleanchor="y",
            constrain="domain", **AXIS_COMMON,
        ),
        yaxis=dict(title="\u0394y (\u03bcas)", **AXIS_COMMON),
        title=dict(
            text=f"t = {t_label}, \u03bd = {nu_label}",
            x=0.5, font=dict(size=12, color="#333"),
        ),
        height=550, width=600,
        margin=dict(l=60, r=80, t=40, b=60),
        plot_bgcolor="#ffffff", paper_bgcolor="#ffffff",
    )


def make_skymap_figure(data):
    """Build an interactive Plotly sky map figure from computation results."""
    image = data["images"][0]
    extent = data["extent_uas"]
    t_obs = data["t_obs_array"][0]
    nu_obs = data["nu_obs"]

    x_min, x_max, y_min, y_max = extent
    npixel = image.shape[0]
    dx = (x_max - x_min) / npixel
    dy = (y_max - y_min) / npixel
    x_centers = np.linspace(x_min + dx / 2, x_max - dx / 2, npixel)
    y_centers = np.linspace(y_min + dy / 2, y_max - dy / 2, npixel)

    fig = go.Figure(data=_skymap_heatmap(image, x_centers, y_centers))
    fig.update_layout(**_skymap_layout(extent, format_time_label(t_obs), freq_label(nu_obs)))
    return fig
