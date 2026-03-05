"""VegasAfterglow Interactive Light Curve Generator.

Launch locally:
    pip install -e ".[webapp]"
    streamlit run webapp/app.py
"""

import base64
import pathlib
import re
import sys

# Ensure repo root is on sys.path so `webapp.*` imports resolve
# when the file is run outside of `streamlit run`.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

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
)
from webapp.helpers import (
    apply_query_params,
    build_share_url,
    log_slider,
    parse_entry,
    z_from_lumi_dist_mpc,
)
from webapp.style import CLIPBOARD_JS, SIDEBAR_CSS

import VegasAfterglow as _va


def _skymap_mpl_style(plt):
    """Apply matplotlib rcParams matching the Plotly skymap style."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 8,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.minor.size": 2.5,
        "ytick.minor.size": 2.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "axes.grid": True,
        "axes.axisbelow": False,
        "grid.linestyle": ":",
        "grid.linewidth": 0.3,
        "grid.alpha": 0.3,
    })


def _parse_obs_file(uploaded_file):
    """Parse uploaded CSV/TXT/DAT into text_area format (x  y  err)."""
    import csv, io as _io
    content = uploaded_file.getvalue().decode("utf-8", errors="replace")
    for delim in [",", "\t", " "]:
        try:
            rows = []
            for row in csv.reader(_io.StringIO(content), delimiter=delim):
                cleaned = [c.strip() for c in row if c.strip()]
                if len(cleaned) < 2:
                    continue
                try:
                    float(cleaned[0]); float(cleaned[1])
                except ValueError:
                    continue
                rows.append("  ".join(cleaned[:3]))
            if rows:
                return "\n".join(rows)
        except Exception:
            continue
    return ""


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
# Restore params from URL on first load
# ---------------------------------------------------------------------------
apply_query_params()

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
    _FREQ_HELP = (
        "Accepts Hz values (1e9), unit suffixes (1GHz, 1keV), "
        "filter names, instrument bands, or custom ranges ([0.3keV,10keV]).\n\n"
        "**Filters:** U B V R I J H Ks | u b v uvw1 uvm2 uvw2 | "
        "g r i z | w VT_B VT_R | "
        "F225W F275W F336W F438W F475W F555W F606W F625W "
        "F775W F814W F850LP F105W F110W F125W F140W F160W\n\n"
        "**Bands:** XRT BAT FXT WXT MXT ECLAIRs LAT GBM\n\n"
        "**Units:** Hz kHz MHz GHz | eV keV MeV GeV"
    )
    if plot_mode == "Light Curve":
        nu_str = st.text_input(
            "Frequencies",
            value="1e9, R, 1keV",
            key="nu_str",
            placeholder="e.g. 1e9, R, 1keV, XRT, [0.3keV,10keV]",
            help=_FREQ_HELP,
        )
        c1, c2 = st.columns(2)
        with c1:
            t_min = log_slider(r"$t_{\rm min}$ (s)", -1.0, 6.0, 0.0, key="t_min")
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
        sky_animate = st.checkbox("Animate", key="sky_animate")
        if sky_animate:
            c1, c2 = st.columns(2)
            with c1:
                sky_t_min = log_slider(r"$t_{\rm min}$ (s)", 2.0, 9.0, 4.0, key="sky_t_min")
            with c2:
                sky_t_max = log_slider(r"$t_{\rm max}$ (s)", 2.0, 9.0, 7.0, key="sky_t_max")
            sky_nframes = st.slider("Frames", 3, 60, 10, 1, key="sky_nframes")
            sky_t_obs = sky_t_min  # for single-frame fallback
        else:
            sky_t_obs = log_slider(r"$t_{\rm obs}$ (s)", 2.0, 9.0, 6.0, key="sky_t_obs")
        sky_nu_str = st.text_input(
            "Frequency",
            value="1e9",
            key="sky_nu_str",
            placeholder="e.g. 1e9, 1GHz, 1keV",
            help=_FREQ_HELP,
        )
        c1, c2 = st.columns(2)
        with c1:
            sky_fov = log_slider("FOV (\u03bcas)", 1.0, 5.0, np.log10(500), key="sky_fov")
        with c2:
            sky_npixel = st.select_slider("Pixels", [64, 128, 256, 512, 1024], value=256, key="sky_npixel")

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
            _FLUX_CHOICES = list(FLUX_SCALES.keys()) + ["AB mag"]
            flux_unit = st.selectbox("Flux", _FLUX_CHOICES, key="flux_unit")
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
                        if _is_lc:
                            x_units = list(TIME_SCALES.keys())
                            x_default = "day"
                        else:
                            x_units = list(FREQ_SCALES.keys())
                            x_default = "Hz"
                        c1, c2 = st.columns(2)
                        with c1:
                            g["x_unit"] = st.selectbox(
                                f"{x_name} unit", x_units,
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
                            height=100,
                            placeholder=f"{x_name}  flux  err\n1e4  0.5  0.1\n1e5  0.3  0.05",
                            help="One row per line. Columns: "
                                 f"{x_name}, flux, err (optional). "
                                 "Separated by space, tab, or comma.",
                            key=f"obs_text_{i}",
                        )
                        with st.popover("Upload file"):
                            uploaded = st.file_uploader(
                                "Upload file",
                                type=["csv", "txt", "dat"],
                                key=f"obs_file_{i}",
                            )
                            if uploaded is not None:
                                parsed = _parse_obs_file(uploaded)
                                if parsed:
                                    g["text"] = parsed
                                    st.session_state[f"obs_text_{i}"] = parsed

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

    # -- Resolutions (collapsed) ------------------------------------
    with st.expander("Resolutions", expanded=False):
        num_t = st.slider("Time points", 50, 300, 100, 10, key="num_t")
        c1, c2, c3 = st.columns(3)
        with c1:
            res_phi = st.slider(r"$\phi$ ppd", 0.05, 1.0, 0.1, 0.05, key="res_phi")
        with c2:
            res_theta = st.slider(r"$\theta$ ppd", 0.1, 2.0, 0.25, 0.1, key="res_theta")
        with c3:
            res_t = st.slider("t ppd", 1.0, 20.0, 10.0, 0.5, key="res_t")

    if st.button("\U0001F517 Share URL", key="share_btn", help="Update URL with current parameters"):
        build_share_url()
        st.toast("\u2714 URL updated — copy from browser address bar")

    st.markdown(
        f'<p style="text-align:center;font-size:0.7em;color:#888;margin-top:0.5rem;">'
        f'VegasAfterglow v{_va.__version__.split(".dev")[0]}<br>'
        f'This interactive tool provides a subset of VegasAfterglow features. '
        f'For MCMC fitting, custom jet/medium models, and the full API, see the '
        f'<a href="https://github.com/YihanWangAstro/VegasAfterglow" target="_blank" '
        f'style="color:#4ecdc4;">Python package</a>.</p>',
        unsafe_allow_html=True,
    )

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
    """Render Plotly chart natively with auto-sized width."""
    fig.update_layout(width=None)
    st.plotly_chart(fig, use_container_width=True,
                    config={"toImageButtonOptions": {"format": "png", "scale": 3}})


def _download_row(fig, downloads, prefix):
    """Show download buttons + Save PNG + Cite in a centered row."""
    n_btns = len(downloads) + 2
    pad, *btn_cols, _ = st.columns([0.5] + [1] * n_btns + [0.5])
    for i, (label, content, ext, mime) in enumerate(downloads):
        with btn_cols[i]:
            st.download_button(label, content, file_name=f"{prefix}.{ext}", mime=mime,
                               use_container_width=True)
    with btn_cols[-2]:
        png_key = f"_png_{prefix}"
        if png_key in st.session_state:
            st.download_button("Download PNG", st.session_state[png_key],
                               file_name=f"{prefix}.png", mime="image/png",
                               key=f"dl_png_{prefix}", use_container_width=True)
        else:
            if st.button("Save PNG", key=f"save_{prefix}", use_container_width=True):
                st.session_state[png_key] = fig.to_image(format="png", scale=3)
                st.rerun()
    with btn_cols[-1]:
        if st.button("Cite", use_container_width=True):
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
    fig = make_figure(data, flux_unit, time_unit, t_min, t_max,
                      need_secondary_y=need_sec)
    use_sec = need_sec or len(data["band_data"]) > 0
    is_mag = (flux_unit == "AB mag")
    if selected_instruments and not is_mag:
        t_scale = TIME_SCALES[time_unit]
        _add_sensitivity_traces(fig, selected_instruments, "lightcurve",
                                flux_scale=FLUX_SCALES[flux_unit],
                                x_range=[t_min / t_scale, t_max / t_scale],
                                has_secondary=use_sec)
    if obs_data_tuple:
        _add_obs_traces(fig, obs_data_tuple, flux_unit, time_unit,
                        has_secondary=use_sec)

    _show_plot(fig)
    export_unit = "cgs" if is_mag else flux_unit
    _download_row(fig, [
        ("CSV", export_csv(data, export_unit, time_unit), "csv", "text/csv"),
        ("JSON", export_json(data, export_unit, time_unit), "json", "application/json"),
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

    is_mag = (flux_unit == "AB mag")
    if is_mag and show_nufnu:
        st.warning("AB mag is not compatible with \u03bdF\u03bd. Showing F\u03bd in AB mag.")
        show_nufnu = False

    has_fband_inst = any(INSTRUMENTS[n][3] == "Fband" for n in selected_instruments)
    has_fband_obs = any(len(r) >= 6 and r[5] == "erg/cm\u00b2/s" for r in obs_data_tuple)
    need_sec = has_fband_inst or has_fband_obs
    fig = make_sed_figure(data, flux_unit, freq_unit, nufnu=show_nufnu,
                          need_secondary_y=need_sec)
    if selected_instruments and not is_mag:
        _add_sensitivity_traces(fig, selected_instruments, "spectrum",
                                freq_scale=FREQ_SCALES[freq_unit],
                                flux_scale=FLUX_SCALES[flux_unit],
                                nufnu=show_nufnu, has_secondary=need_sec)
    if obs_data_tuple:
        _add_obs_traces(fig, obs_data_tuple, flux_unit, freq_unit,
                        has_secondary=need_sec, mode="spectrum",
                        nufnu=show_nufnu)

    _show_plot(fig)
    export_unit = "cgs" if is_mag else flux_unit
    _download_row(fig, [
        ("CSV", export_sed_csv(data, export_unit, freq_unit), "csv", "text/csv"),
        ("JSON", export_sed_json(data, export_unit, freq_unit), "json", "application/json"),
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
    if sky_animate:
        t_arr = np.logspace(np.log10(sky_t_min), np.log10(sky_t_max), sky_nframes).tolist()
        params["t_obs_array"] = t_arr
    try:
        data = compute_skymap(params)
    except Exception as e:
        st.error(f"Computation failed: {e}")
        st.stop()

    if sky_animate and len(data["images"]) > 1:
        import io as _io
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from webapp.helpers import format_time_label, freq_label

        n_frames = len(data["images"])
        frame_ms = max(50, 2000 // n_frames)
        extent = data["extent_uas"]
        ext_arr = [extent[0], extent[1], extent[2], extent[3]]
        all_pos = np.concatenate([img[img > 0] for img in data["images"]
                                  if np.any(img > 0)])
        vmin = float(np.log10(all_pos.min())) if all_pos.size > 0 else 0
        vmax = float(np.log10(all_pos.max())) if all_pos.size > 0 else 1

        _skymap_mpl_style(plt)

        # Create figure once, reuse for all frames
        fig_mpl, ax = plt.subplots(figsize=(4.0, 3.3), dpi=200, facecolor="white")
        log_img0 = np.where(data["images"][0] > 0,
                            np.log10(data["images"][0]), np.nan)
        import matplotlib.cm as _mcm
        cmap = _mcm.get_cmap("inferno").copy()
        cmap.set_bad(color="white")
        im = ax.imshow(log_img0.T, origin="lower", extent=ext_arr,
                       cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
        ax.set_xlabel("\u0394x (\u03bcas)")
        ax.set_ylabel("\u0394y (\u03bcas)")
        title = ax.set_title("")
        cb = fig_mpl.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(r"$\log_{10}\,I$ (erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ sr$^{-1}$)")
        fig_mpl.subplots_adjust(left=0.15, right=0.88, bottom=0.13, top=0.92)

        # Full render once (no title) to get static overlay
        nu_label = freq_label(data["nu_obs"])
        title.set_visible(False)
        canvas = fig_mpl.canvas
        canvas.draw()
        w, h = canvas.get_width_height()

        # Axes pixel bbox (canvas y=0 at bottom, image y=0 at top)
        bbox = ax.get_window_extent(canvas.get_renderer())
        x0, y0 = int(round(bbox.x0)), int(round(bbox.y0))
        x1, y1 = int(round(bbox.x1)), int(round(bbox.y1))
        img_w, img_h = x1 - x0, y1 - y0

        # Capture overlay (axes, ticks, colorbar — no image data, no title)
        # Make axes interior transparent so we can composite image underneath
        im.set_visible(False)
        ax.set_facecolor("none")
        fig_mpl.patch.set_alpha(0)
        canvas.draw()
        overlay_pil = _PILImage.frombuffer("RGBA", (w, h), canvas.buffer_rgba()).copy()
        im.set_visible(True)
        plt.close(fig_mpl)

        # Pre-compute colormapped RGBA for each frame (pure numpy)
        from matplotlib.colors import Normalize
        norm = Normalize(vmin=vmin, vmax=vmax)
        mapped_frames = []
        for i in range(n_frames):
            img_data = data["images"][i]
            log_img = np.where(img_data > 0, np.log10(img_data), np.nan)
            rgba = (cmap(norm(log_img.T)) * 255).astype(np.uint8)
            rgba = rgba[::-1]  # flip for canvas y-down
            pil_img = _PILImage.fromarray(rgba)
            if pil_img.size != (img_w, img_h):
                pil_img = pil_img.resize((img_w, img_h), _PILImage.NEAREST)
            mapped_frames.append(pil_img)

        # Title via PIL ImageDraw (no matplotlib per frame)
        from PIL import ImageDraw, ImageFont
        try:
            title_font = ImageFont.truetype("Helvetica", int(9 * 200 / 72))
        except OSError:
            try:
                title_font = ImageFont.truetype("Arial", int(9 * 200 / 72))
            except OSError:
                title_font = ImageFont.load_default()
        title_y = h - y1 - int(9 * 200 / 72) - 4  # above axes top
        title_color = (51, 51, 51, 255)

        # Composite: white bg + image + overlay (ticks on top) + title
        frames_pil = []
        for i in range(n_frames):
            frame = _PILImage.new("RGBA", (w, h), (255, 255, 255, 255))
            frame.paste(mapped_frames[i], (x0, h - y1))
            frame = _PILImage.alpha_composite(frame, overlay_pil)
            # Draw title centered
            t_label = format_time_label(data["t_obs_array"][i])
            txt = f"t = {t_label},  \u03bd = {nu_label}"
            draw = ImageDraw.Draw(frame)
            tw = draw.textlength(txt, font=title_font)
            draw.text(((w - tw) / 2, title_y), txt, fill=title_color, font=title_font)
            frames_pil.append(frame)
        buf = _io.BytesIO()
        frames_pil[0].save(buf, format="GIF", save_all=True,
                           append_images=frames_pil[1:],
                           duration=frame_ms, loop=0)
        gif_bytes = buf.getvalue()
        gif_b64 = base64.b64encode(gif_bytes).decode()
        st.markdown(
            f'<div style="text-align:center;">'
            f'<img src="data:image/gif;base64,{gif_b64}" '
            f'style="max-width:70%;"></div>',
            unsafe_allow_html=True,
        )
        _, c1, c2, c3, _ = st.columns([1, 1, 1, 1, 1])
        with c1:
            st.download_button("JSON", export_skymap_json(data),
                               file_name="afterglow_skymap.json", mime="application/json",
                               use_container_width=True)
        with c2:
            st.download_button("GIF", gif_bytes,
                               file_name="afterglow_skymap.gif", mime="image/gif",
                               key="dl_gif_skymap", use_container_width=True)
        with c3:
            if st.button("Cite", key="cite_skymap_anim", use_container_width=True):
                _stc.html(CLIPBOARD_JS, height=0)
                st.toast("\u2714 BibTeX copied to clipboard!")
    else:
        import io as _io
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from webapp.helpers import format_time_label, freq_label

        extent = data["extent_uas"]
        ext_arr = [extent[0], extent[1], extent[2], extent[3]]
        image = data["images"][0]
        log_img = np.where(image > 0, np.log10(image), np.nan)
        pos = image[image > 0]
        vmin = float(np.log10(pos.min())) if pos.size > 0 else 0
        vmax = float(np.log10(pos.max())) if pos.size > 0 else 1

        _skymap_mpl_style(plt)
        import matplotlib.cm as _mcm
        cmap = _mcm.get_cmap("inferno").copy()
        cmap.set_bad(color="white")
        fig_mpl, ax = plt.subplots(figsize=(4.0, 3.3), dpi=200, facecolor="white")
        im = ax.imshow(log_img.T, origin="lower", extent=ext_arr,
                       cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
        ax.set_xlabel("\u0394x (\u03bcas)")
        ax.set_ylabel("\u0394y (\u03bcas)")
        t_label = format_time_label(data["t_obs_array"][0])
        nu_label = freq_label(data["nu_obs"])
        ax.set_title(f"t = {t_label},  \u03bd = {nu_label}")
        cb = fig_mpl.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(r"$\log_{10}\,I$ (erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ sr$^{-1}$)")
        fig_mpl.subplots_adjust(left=0.15, right=0.88, bottom=0.13, top=0.92)

        buf = _io.BytesIO()
        fig_mpl.savefig(buf, format="png", facecolor="white",
                        bbox_inches="tight", pad_inches=0.03)
        plt.close(fig_mpl)
        buf.seek(0)
        png_bytes = buf.getvalue()
        png_b64 = base64.b64encode(png_bytes).decode()
        st.markdown(
            f'<div style="text-align:center;">'
            f'<img src="data:image/png;base64,{png_b64}" '
            f'style="max-width:70%;"></div>',
            unsafe_allow_html=True,
        )
        _, c1, c2, c3, _ = st.columns([1, 1, 1, 1, 1])
        with c1:
            st.download_button("JSON", export_skymap_json(data),
                               file_name="afterglow_skymap.json", mime="application/json",
                               use_container_width=True)
        with c2:
            st.download_button("PNG", png_bytes,
                               file_name="afterglow_skymap.png", mime="image/png",
                               key="dl_png_skymap_static", use_container_width=True)
        with c3:
            if st.button("Cite", key="cite_skymap_static", use_container_width=True):
                _stc.html(CLIPBOARD_JS, height=0)
                st.toast("\u2714 BibTeX copied to clipboard!")
