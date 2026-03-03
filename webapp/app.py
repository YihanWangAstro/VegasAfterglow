"""VegasAfterglow Interactive Light Curve Generator.

Launch locally:
    pip install -e ".[webapp]"
    streamlit run webapp/app.py
"""

import base64
import pathlib
import re

from PIL import Image as _PILImage

import numpy as np
import streamlit as st
import streamlit.components.v1 as _stc

from webapp.constants import (
    FLUX_SCALES,
    FREQ_SCALES,
    INSTRUMENTS,
    OBS_FLUX_UNITS,
    TIME_SCALES,
)
from webapp.compute import compute_model, compute_sed, compute_skymap
from webapp.exports import (
    export_csv,
    export_json,
    export_sed_csv,
    export_sed_json,
    export_skymap_json,
)
from webapp.figures import (
    _add_obs_traces,
    _add_sensitivity_traces,
    make_figure,
    make_sed_figure,
    make_skymap_figure,
)
from webapp.helpers import log_slider, parse_entry, z_from_lumi_dist_mpc
from webapp.style import CLIPBOARD_JS, SIDEBAR_CSS

_ASSETS = pathlib.Path(__file__).resolve().parent.parent / "assets"
_LOGO_B64 = base64.b64encode((_ASSETS / "logo-horizontal.svg").read_bytes()).decode()
_FAVICON = _ASSETS / "favicon.png"

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="VegasAfterglow",
    page_icon=_PILImage.open(_FAVICON),
    layout="wide",
)

# ---------------------------------------------------------------------------
# Sidebar CSS
# ---------------------------------------------------------------------------
st.markdown(SIDEBAR_CSS, unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Sidebar: parameters (no form — every change triggers recompute)
# ---------------------------------------------------------------------------

with st.sidebar:
    st.markdown(
        f'<a href="https://github.com/YihanWangAstro/VegasAfterglow" target="_blank">'
        f'<img src="data:image/svg+xml;base64,{_LOGO_B64}" style="width:100%;" />'
        f'</a>',
        unsafe_allow_html=True,
    )

    st.markdown(
        '<p style="text-align:right;font-size:0.8em;color:#888;margin:0;">'
        'You can use \u2190/\u2192 to scroll the sidebar.</p>',
        unsafe_allow_html=True,
    )

    # -- Mode toggle -------------------------------------------------
    plot_mode = st.radio("Mode", ["Light Curve", "Spectrum", "Sky Image"], horizontal=True, key="plot_mode")

    # -- Mode-dependent inputs ---------------------------------------
    if plot_mode == "Light Curve":
        nu_str = st.text_input(
            "Frequencies",
            value="1e9, R, 1keV",
            key="nu_str",
            placeholder="e.g. 1e9, R, 1keV, XRT, [0.3keV,10keV]",
        )
        c1, c2 = st.columns(2)
        with c1:
            t_min = log_slider(r"$t_{\rm min}$ (s)", -2.0, 6.0, 0.0, key="t_min")
        with c2:
            t_max = log_slider(r"$t_{\rm max}$ (s)", 3.0, 10.0, 8.0, key="t_max")
    elif plot_mode == "Spectrum":
        t_snap_str = st.text_input(
            "Observer times (s)",
            value="1e3, 1e4, 1e5, 1e6",
            key="t_snap_str",
            placeholder="e.g. 1e3, 1e4, 1e5",
        )
        c1, c2 = st.columns(2)
        with c1:
            sed_nu_min = log_slider(r"$\nu_{\rm min}$ (Hz)", 6.0, 20.0, 8.0, key="sed_nu_min")
        with c2:
            sed_nu_max = log_slider(r"$\nu_{\rm max}$ (Hz)", 10.0, 35.0, 20.0, key="sed_nu_max")
        c1, c2, c3 = st.columns(3)
        with c1:
            num_nu = st.slider("Freq points", 50, 500, 200, 10, key="num_nu")
        with c2:
            freq_unit = st.selectbox("Freq unit", list(FREQ_SCALES.keys()), key="freq_unit")
        with c3:
            show_nufnu = st.checkbox(r"$\nu F_\nu$", key="show_nufnu")
    else:  # Sky Image
        sky_t_obs = log_slider(r"$t_{\rm obs}$ (s)", 2.0, 9.0, 6.0, key="sky_t_obs")
        sky_nu_str = st.text_input(
            "Frequency",
            value="1e9",
            key="sky_nu_str",
            placeholder="e.g. 1e9, 1GHz, 1keV",
        )
        c1, c2 = st.columns(2)
        with c1:
            sky_fov = log_slider("FOV (\u03bcas)", 1.0, 5.0, np.log10(500), key="sky_fov")
        with c2:
            sky_npixel = st.select_slider("Pixels", [64, 128, 256, 512, 1024, 2048], value=256, key="sky_npixel")

    # -- Shared observer + units ------------------------------------
    c1, c2 = st.columns(2)
    with c1:
        d_L_mpc = st.number_input(r"$d_L$ (Mpc)", min_value=0.1, max_value=1e6, value=100.0, step=10.0, key="d_L_mpc")
    with c2:
        theta_obs = st.slider(r"$\theta_{\rm obs}$ (rad)", 0.0, 1.57, 0.0, 0.01, key="theta_obs")
    z = z_from_lumi_dist_mpc(d_L_mpc)
    d_L_cm = d_L_mpc * 3.0856775814913673e24  # Mpc -> cm
    if plot_mode != "Sky Image":
        c1, c2 = st.columns(2)
        with c1:
            flux_unit = st.selectbox("Flux", list(FLUX_SCALES.keys()), key="flux_unit")
        with c2:
            if plot_mode == "Light Curve":
                time_unit = st.selectbox("Time", list(TIME_SCALES.keys()), key="time_unit")
            else:
                time_unit = "s"
    else:
        flux_unit = "cgs"
        time_unit = "s"

    # -- Jet --------------------------------------------------------
    with st.container(border=True):
        c1, c2 = st.columns([3, 2])
        with c1:
            jet_type = st.selectbox(
                "Jet", ["Top-hat", "Gaussian", "Power-law", "Two-component"],
                key="jet_type",
            )
        with c2:
            theta_c = st.slider(r"$\theta_c$", 0.01, 1.0, 0.1, 0.01, key="theta_c")
        c1, c2 = st.columns(2)
        with c1:
            E_iso = log_slider(r"$E_{\rm iso}$ (erg)", 48.0, 57.0, 52.0, key="E_iso")
        with c2:
            Gamma0 = log_slider(r"$\Gamma_0$", 1.0, 3.5, np.log10(300), key="Gamma0")
        spreading = st.checkbox("Spreading", key="spreading")

        if jet_type == "Power-law":
            c1, c2 = st.columns(2)
            with c1:
                k_e = st.slider(r"$k_e$", 0.5, 10.0, 2.0, 0.1, key="k_e")
            with c2:
                k_g = st.slider(r"$k_g$", 0.5, 10.0, 2.0, 0.1, key="k_g")
        else:
            k_e, k_g = 2.0, 2.0

        if jet_type == "Two-component":
            c1, c2 = st.columns(2)
            with c1:
                theta_w = st.slider(r"$\theta_w$", 0.05, 1.5, 0.3, 0.01, key="theta_w")
            with c2:
                E_iso_w = log_slider(r"$E_{\rm iso,w}$", 48.0, 55.0, 51.0, key="E_iso_w")
            Gamma0_w = log_slider(r"$\Gamma_{0,w}$", 1.0, 3.0, 2.0, key="Gamma0_w")
        else:
            theta_w, E_iso_w, Gamma0_w = 0.3, 1e51, 100

    # -- Medium -----------------------------------------------------
    with st.container(border=True):
        c1, c2 = st.columns([2, 3])
        with c1:
            medium_type = st.selectbox("Medium", ["ISM", "Wind"], key="medium_type")
        with c2:
            if medium_type == "ISM":
                n_ism = log_slider(r"$n_{\rm ism}$ (cm⁻³)", -5.0, 5.0, 0.0, key="n_ism")
                A_star, k_m = 0.1, 2.0
            else:
                A_star = log_slider(r"$A_*$", -3.0, 2.0, -1.0, key="A_star")
                n_ism = 1.0
        if medium_type == "Wind":
            k_m = st.slider(r"$k_m$", 0.0, 4.0, 2.0, 0.1, key="k_m")
        else:
            k_m = 2.0

    # -- Radiation --------------------------------------------------
    with st.container(border=True):
        c1, c2 = st.columns(2)
        with c1:
            eps_e = log_slider(r"$\varepsilon_e$", -5.0, 0.0, -1.0, key="eps_e")
        with c2:
            eps_B = log_slider(r"$\varepsilon_B$", -6.0, 0.0, -3.0, key="eps_B")
        c1, c2 = st.columns(2)
        with c1:
            p_val = st.slider("p", 2.01, 3.0, 2.3, 0.01, key="p_val")
        with c2:
            xi_e = log_slider(r"$\xi_e$", -3.0, 0.0, 0.0, key="xi_e")
        c1, c2 = st.columns(2)
        with c1:
            ssc = st.checkbox("SSC", key="ssc")
        with c2:
            kn = st.checkbox("KN", key="kn")

    # -- Reverse shock (collapsed expander) -------------------------
    with st.expander("Reverse Shock", expanded=False):
        enable_rvs = st.checkbox("Enable", key="enable_rvs")
        if enable_rvs:
            duration = log_slider("Duration (s)", -1.0, 4.0, 0.0, key="duration")
            c1, c2 = st.columns(2)
            with c1:
                eps_e_r = log_slider(r"$\varepsilon_e$", -5.0, 0.0, -1.0, key="eps_e_r")
            with c2:
                eps_B_r = log_slider(r"$\varepsilon_B$", -6.0, 0.0, -3.0, key="eps_B_r")
            c1, c2 = st.columns(2)
            with c1:
                p_r = st.slider("p", 2.01, 3.0, 2.3, 0.01, key="p_r")
            with c2:
                xi_e_r = log_slider(r"$\xi_e$", -3.0, 0.0, 0.0, key="xi_e_r")
            c1, c2 = st.columns(2)
            with c1:
                rvs_ssc = st.checkbox("SSC", key="rvs_ssc")
            with c2:
                rvs_kn = st.checkbox("KN", key="rvs_kn")
        else:
            duration = 1.0
            eps_e_r, eps_B_r, p_r, xi_e_r = eps_e, eps_B, p_val, xi_e
            rvs_ssc, rvs_kn = False, False

    # -- Instruments (collapsed) ----------------------------------------
    if plot_mode != "Sky Image":
        with st.expander("Instruments", expanded=False):
            selected_instruments = st.multiselect(
                "Show sensitivity",
                list(INSTRUMENTS.keys()),
                default=[],
                key="instruments",
            )
    else:
        selected_instruments = []

    # -- Observations (Light Curve & Spectrum) ---------------------------
    _is_lc = (plot_mode == "Light Curve")
    _obs_state_key = "obs_groups_lc" if _is_lc else "obs_groups_sed"
    if _obs_state_key not in st.session_state:
        st.session_state[_obs_state_key] = []
    obs_data_tuple = ()
    if plot_mode in ("Light Curve", "Spectrum"):
        # Handle pending group removal before rendering widgets
        _obs_rm_key = f"{_obs_state_key}_rm"
        if _obs_rm_key in st.session_state:
            idx = st.session_state.pop(_obs_rm_key)
            if 0 <= idx < len(st.session_state[_obs_state_key]):
                st.session_state[_obs_state_key].pop(idx)

        with st.expander("Observation Data", expanded=False):
            groups = st.session_state[_obs_state_key]
            if st.button("\u2795 Add group", key="obs_add_group"):
                groups.append({
                    "legend": f"data {len(groups) + 1}",
                    "x_unit": "day" if _is_lc else "Hz",
                    "y_unit": "mJy",
                    "text": "",
                    "visible": True,
                })

            # Render each group
            if groups:
                x_name = "t" if _is_lc else "\u03bd"
                tab_labels = [g.get("legend", f"Group {i+1}") or f"Group {i+1}"
                              for i, g in enumerate(groups)]
                tabs = st.tabs(tab_labels)
                for i, tab in enumerate(tabs):
                    g = groups[i]
                    with tab:
                        c_vis, c_del = st.columns([5, 1])
                        with c_vis:
                            g["visible"] = st.checkbox(
                                "Show", value=g.get("visible", True),
                                key=f"obs_vis_{i}")
                        with c_del:
                            if st.button("\u2716", key=f"obs_del_{i}",
                                         help="Remove this group",
                                         type="tertiary"):
                                st.session_state[_obs_rm_key] = i
                                st.rerun()
                        g["legend"] = st.text_input(
                            "Legend", value=g.get("legend", ""),
                            key=f"obs_legend_{i}")
                        c1, c2 = st.columns(2)
                        if _is_lc:
                            x_units = list(TIME_SCALES.keys())
                            x_label = "t unit"
                            x_default = "day"
                        else:
                            x_units = list(FREQ_SCALES.keys())
                            x_label = "\u03bd unit"
                            x_default = "Hz"
                        with c1:
                            g["x_unit"] = st.selectbox(
                                x_label, x_units,
                                index=x_units.index(
                                    g.get("x_unit", x_default)),
                                key=f"obs_xunit_{i}")
                        with c2:
                            g["y_unit"] = st.selectbox(
                                "Flux unit", OBS_FLUX_UNITS,
                                index=OBS_FLUX_UNITS.index(
                                    g.get("y_unit", "mJy")),
                                key=f"obs_yunit_{i}")
                        g["text"] = st.text_area(
                            f"{x_name}, flux, err",
                            value=g.get("text", ""),
                            height=120,
                            placeholder=f"e.g.\n1e4  0.5  0.1\n1e5  0.3  0.05",
                            help="One row per line. Columns: "
                                 f"{x_name}, flux, err (optional). "
                                 "Separated by space, tab, or comma.",
                            key=f"obs_text_{i}",
                        )

            st.session_state[_obs_state_key] = groups

        # Parse all groups into flat tuple for caching
        # Format: (label, x_val, x_unit, y_val, err_val, y_unit)
        x_col = "t" if _is_lc else "nu"
        rows = []
        for g in groups:
            if not g.get("visible", True):
                continue
            label = g.get("legend", "") or "data"
            x_unit = g.get("x_unit", "day" if _is_lc else "Hz")
            y_unit = g.get("y_unit", "mJy")
            for line in g.get("text", "").strip().splitlines():
                parts = re.split(r"[,\s\t]+", line.strip())
                if len(parts) < 2:
                    continue
                try:
                    x_val = float(parts[0])
                    y_val = float(parts[1])
                except (ValueError, IndexError):
                    continue
                try:
                    err_val = float(parts[2]) if len(parts) > 2 else 0.0
                except ValueError:
                    err_val = 0.0
                rows.append((label, x_val, x_unit, y_val, err_val, y_unit))
        obs_data_tuple = tuple(rows)

    # -- Config (collapsed) -----------------------------------------
    with st.expander("Config", expanded=False):
        num_t = st.slider("Time points", 50, 300, 100, 10, key="num_t")
        c1, c2, c3 = st.columns(3)
        with c1:
            res_phi = st.slider(r"$\phi$ ppd", 0.05, 1.0, 0.1, 0.05, key="res_phi")
        with c2:
            res_theta = st.slider(r"$\theta$ ppd", 0.1, 2.0, 0.5, 0.1, key="res_theta")
        with c3:
            res_t = st.slider("t ppd", 1.0, 20.0, 5.0, 0.5, key="res_t")

# ---------------------------------------------------------------------------
# Shared physics params
# ---------------------------------------------------------------------------

_physics_params = dict(
    jet_type=jet_type,
    theta_c=theta_c,
    E_iso=E_iso,
    Gamma0=Gamma0,
    spreading=spreading,
    duration=duration,
    k_e=k_e,
    k_g=k_g,
    theta_w=theta_w,
    E_iso_w=E_iso_w,
    Gamma0_w=Gamma0_w,
    medium_type=medium_type,
    n_ism=n_ism,
    A_star=A_star,
    k_m=k_m,
    d_L_cm=d_L_cm,
    z=z,
    theta_obs=theta_obs,
    eps_e=eps_e,
    eps_B=eps_B,
    p=p_val,
    xi_e=xi_e,
    ssc=ssc,
    kn=kn,
    enable_rvs=enable_rvs,
    eps_e_r=eps_e_r,
    eps_B_r=eps_B_r,
    p_r=p_r,
    xi_e_r=xi_e_r,
    rvs_ssc=rvs_ssc,
    rvs_kn=rvs_kn,
    res_phi=res_phi,
    res_theta=res_theta,
    res_t=res_t,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _show_plot(fig):
    """Render via fixed-height iframe — no layout reflow on rerun."""
    html = fig.to_html(
        include_plotlyjs="cdn", full_html=False,
        config={"toImageButtonOptions": {"format": "png", "scale": 3}},
    )
    _stc.html(f'<div style="width:100%;background:#fff;">{html}</div>',
              height=545, scrolling=False)


def _download_row(fig, downloads, prefix):
    """Show download buttons + Save Figure + Cite in a row."""
    cols = st.columns(len(downloads) + 2)
    for i, (label, content, ext, mime) in enumerate(downloads):
        with cols[i]:
            st.download_button(label, content, file_name=f"{prefix}.{ext}", mime=mime)
    with cols[-2]:
        png_key = f"_png_{prefix}"
        if st.button("Save Figure", key=f"save_{prefix}"):
            st.session_state[png_key] = fig.to_image(format="png", scale=3)
        if png_key in st.session_state:
            st.download_button("Download PNG", st.session_state[png_key],
                               file_name=f"{prefix}.png", mime="image/png",
                               key=f"dl_png_{prefix}")
    with cols[-1]:
        if st.button("Cite"):
            _stc.html(CLIPBOARD_JS, height=0)
            st.toast("\u2714 BibTeX copied to clipboard!")


# ---------------------------------------------------------------------------
# Main area: compute + display
# ---------------------------------------------------------------------------

if plot_mode == "Light Curve":
    frequencies, bands = [], []
    tokens = [t.strip() for t in re.split(r",(?![^\[]*\])", nu_str) if t.strip()]
    for part in tokens:
        try:
            entry = parse_entry(part)
            (bands if isinstance(entry, tuple) else frequencies).append(entry)
        except Exception:
            st.sidebar.warning(f"Unknown frequency or filter: '{part}'")
    if not frequencies and not bands:
        frequencies = [1e9]

    params = dict(**_physics_params, frequencies=frequencies, bands=bands,
                  t_min=t_min, t_max=t_max, num_t=num_t)
    try:
        data = compute_model(params)
    except Exception as e:
        st.error(f"Computation failed: {e}")
        st.stop()

    has_fband_inst = any(INSTRUMENTS[n][3] == "Fband" for n in selected_instruments)
    has_fband_obs = any(len(r) >= 6 and r[5] == "erg/cm\u00b2/s" for r in obs_data_tuple)
    need_sec = has_fband_inst or has_fband_obs
    fig = make_figure(data, flux_unit, time_unit, t_min, t_max, need_secondary_y=need_sec)
    use_sec = need_sec or len(data["band_data"]) > 0
    if selected_instruments:
        t_scale = TIME_SCALES[time_unit]
        _add_sensitivity_traces(fig, selected_instruments, "lightcurve",
                                flux_scale=FLUX_SCALES[flux_unit],
                                x_range=[t_min / t_scale, t_max / t_scale],
                                has_secondary=use_sec)
    if obs_data_tuple:
        _add_obs_traces(fig, obs_data_tuple, flux_unit, time_unit, has_secondary=use_sec)

    _show_plot(fig)
    _download_row(fig, [
        ("Download CSV", export_csv(data, flux_unit, time_unit), "csv", "text/csv"),
        ("Download JSON", export_json(data, flux_unit, time_unit), "json", "application/json"),
    ], "afterglow_lightcurve")

elif plot_mode == "Spectrum":
    t_snapshots = []
    for tok in t_snap_str.split(","):
        tok = tok.strip()
        if tok:
            try:
                t_snapshots.append(float(tok))
            except ValueError:
                st.sidebar.warning(f"Invalid time value: '{tok}'")
    if not t_snapshots:
        t_snapshots = [1e4]
    if sed_nu_min >= sed_nu_max:
        st.sidebar.warning("nu_min must be less than nu_max")
        st.stop()

    params = dict(**_physics_params, t_snapshots=sorted(t_snapshots),
                  nu_min=sed_nu_min, nu_max=sed_nu_max, num_nu=num_nu)
    try:
        data = compute_sed(params)
    except Exception as e:
        st.error(f"Computation failed: {e}")
        st.stop()

    has_fband_inst = any(INSTRUMENTS[n][3] == "Fband" for n in selected_instruments)
    has_fband_obs = any(len(r) >= 6 and r[5] == "erg/cm\u00b2/s" for r in obs_data_tuple)
    need_sec = has_fband_inst or has_fband_obs
    fig = make_sed_figure(data, flux_unit, freq_unit, nufnu=show_nufnu,
                          need_secondary_y=need_sec)
    if selected_instruments:
        _add_sensitivity_traces(fig, selected_instruments, "spectrum",
                                freq_scale=FREQ_SCALES[freq_unit],
                                flux_scale=FLUX_SCALES[flux_unit],
                                nufnu=show_nufnu, has_secondary=need_sec)
    if obs_data_tuple:
        _add_obs_traces(fig, obs_data_tuple, flux_unit, freq_unit,
                        has_secondary=need_sec, mode="spectrum", nufnu=show_nufnu)

    _show_plot(fig)
    _download_row(fig, [
        ("Download CSV", export_sed_csv(data, flux_unit, freq_unit), "csv", "text/csv"),
        ("Download JSON", export_sed_json(data, flux_unit, freq_unit), "json", "application/json"),
    ], "afterglow_sed")

else:  # Sky Image
    sky_nu_str_val = sky_nu_str.strip()
    try:
        sky_nu = parse_entry(sky_nu_str_val) if sky_nu_str_val else 1e9
        if isinstance(sky_nu, tuple):
            st.sidebar.warning("Sky image requires a single frequency, not a band")
            st.stop()
    except Exception:
        st.sidebar.warning(f"Invalid frequency: '{sky_nu_str_val}'")
        sky_nu = 1e9

    params = dict(**_physics_params, t_obs=sky_t_obs, nu_obs=float(sky_nu),
                  fov=sky_fov, npixel=sky_npixel)
    try:
        data = compute_skymap(params)
    except Exception as e:
        st.error(f"Computation failed: {e}")
        st.stop()

    fig = make_skymap_figure(data)
    _show_plot(fig)
    _download_row(fig, [
        ("Download JSON", export_skymap_json(data), "json", "application/json"),
    ], "afterglow_skymap")
