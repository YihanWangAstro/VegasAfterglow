"""VegasAfterglow Interactive Light Curve Generator.

Launch locally:
    pip install -e ".[webapp]"
    streamlit run webapp/app.py
"""

import base64
from collections import OrderedDict
import hashlib
import json
import pathlib
import re
import sys
import time as time_mod

# On Streamlit Cloud the local VegasAfterglow/ source tree (no compiled
# C++ extension) shadows the pip-installed wheel.  Fix: strip the repo
# root from sys.path, import VegasAfterglow from site-packages (caches
# in sys.modules), then add the repo root back for webapp.* imports.
# Subsequent VegasAfterglow submodule imports use the cached package's
# __path__ (site-packages), not sys.path.
_repo_root = str(pathlib.Path(__file__).resolve().parent.parent)
sys.path = [p for p in sys.path if p not in ("", ".", _repo_root)]
import VegasAfterglow  # noqa: E402, F401  — from site-packages
sys.path.insert(0, _repo_root)

from PIL import Image as _PILImage

import numpy as np
import plotly.graph_objects as go
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


_ASSETS = pathlib.Path(__file__).resolve().parent.parent / "assets"
_LOGO_B64 = base64.b64encode((_ASSETS / "logo-horizontal.svg").read_bytes()).decode()
_FAVICON = _ASSETS / "favicon.png"
_PLOTLY_COMPONENT_DIR = pathlib.Path(__file__).resolve().parent / "components" / "plotly_live"
_PLOTLY_LIVE = _stc.declare_component("plotly_live_v2", path=str(_PLOTLY_COMPONENT_DIR))

_PLOT_HEIGHT = 560
_FIG_CACHE_MAX = 8
_COMPUTE_SEEN_MAX = 256


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------


def _perf_debug_enabled():
    raw = st.query_params.get("debug_perf", "0")
    if isinstance(raw, list):
        raw = raw[0] if raw else "0"
    return str(raw).strip().lower() in {"1", "true", "yes", "on"}


def _to_jsonable(value):
    if isinstance(value, dict):
        return {str(k): _to_jsonable(v) for k, v in sorted(value.items())}
    if isinstance(value, (list, tuple)):
        return [_to_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def _decode_plotly_binary_arrays(value):
    if isinstance(value, dict):
        if (
            "dtype" in value
            and "bdata" in value
            and isinstance(value["dtype"], str)
            and isinstance(value["bdata"], str)
        ):
            try:
                arr = np.frombuffer(
                    base64.b64decode(value["bdata"]),
                    dtype=np.dtype(value["dtype"]),
                )
                shape = value.get("shape")
                if isinstance(shape, (list, tuple)) and len(shape) > 0:
                    arr = arr.reshape(tuple(int(s) for s in shape))
                return arr.tolist()
            except Exception:
                pass
        return {k: _decode_plotly_binary_arrays(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_decode_plotly_binary_arrays(v) for v in value]
    return value


def _plotly_component_figure(fig):
    # Plotly.py v6 serializes arrays as {"dtype","bdata"} blobs.
    # The custom JS path needs plain arrays for reliable updates.
    return _decode_plotly_binary_arrays(json.loads(fig.to_json()))


def _canonical_json(value):
    return json.dumps(_to_jsonable(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _make_key(prefix, payload):
    digest = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
    return f"{prefix}:{digest}"


def _get_lru_cache(name):
    key = f"_lru_{name}"
    cache = st.session_state.get(key)
    if not isinstance(cache, OrderedDict):
        cache = OrderedDict()
        st.session_state[key] = cache
    return cache


def _figure_cache_get(mode, view_key):
    cache = _get_lru_cache(f"fig_{mode}")
    fig_json = cache.get(view_key)
    if fig_json is None:
        return None, False
    cache.move_to_end(view_key)
    return go.Figure(json.loads(fig_json)), True


def _figure_cache_put(mode, view_key, fig):
    cache = _get_lru_cache(f"fig_{mode}")
    cache[view_key] = fig.to_json()
    cache.move_to_end(view_key)
    while len(cache) > _FIG_CACHE_MAX:
        cache.popitem(last=False)


def _seen_compute_key(mode, compute_key):
    cache = _get_lru_cache(f"seen_compute_{mode}")
    hit = compute_key in cache
    cache[compute_key] = None
    cache.move_to_end(compute_key)
    while len(cache) > _COMPUTE_SEEN_MAX:
        cache.popitem(last=False)
    return hit


def _lock_y_limits(fig, mode, axis_key):
    """Persist and reapply y-axis ranges to prevent auto-rescaling flicker."""
    state_key = f"_ylims_{mode}"
    state = st.session_state.get(state_key)
    if not isinstance(state, dict) or state.get("axis_key") != axis_key:
        state = {"axis_key": axis_key}

    for axis_name in ("yaxis", "yaxis2"):
        axis_obj = getattr(fig.layout, axis_name, None)
        range_key = f"{axis_name}_range"
        if axis_obj is None:
            state.pop(range_key, None)
            continue

        current_range = list(axis_obj.range) if axis_obj.range is not None else None
        saved_range = state.get(range_key)

        if saved_range is None:
            if current_range is not None:
                state[range_key] = current_range
        else:
            fig.update_layout(
                **{
                    f"{axis_name}_range": saved_range,
                    f"{axis_name}_autorange": False,
                }
            )

    st.session_state[state_key] = state


def _timed_export_payload(metric_key, fn):
    t0 = time_mod.perf_counter()
    payload = fn()
    dt = time_mod.perf_counter() - t0
    d = st.session_state.setdefault("_perf_export_sec", {})
    d[metric_key] = dt
    return payload


def _show_perf_debug(mode_prefix, perf):
    if not _perf_debug_enabled():
        return
    export_sec = st.session_state.get("_perf_export_sec", {})
    mode_exports = {
        k: v for k, v in export_sec.items() if k.startswith(f"{mode_prefix}_")
    }
    with st.expander("Perf Debug", expanded=False):
        st.caption(
            f"key_build={perf.get('key_build_s', 0.0):.4f}s | "
            f"compute={perf.get('compute_s', 0.0):.4f}s "
            f"({'hit' if perf.get('compute_hit') else 'miss'}) | "
            f"figure={perf.get('figure_s', 0.0):.4f}s "
            f"({'hit' if perf.get('figure_hit') else 'miss'})"
        )
        if "obs_parse_s" in perf:
            st.caption(
                f"obs_parse={perf.get('obs_parse_s', 0.0):.4f}s "
                f"({'hit' if perf.get('obs_parse_hit') else 'miss'})"
            )
        if mode_exports:
            txt = " | ".join(f"{k}={v:.4f}s" for k, v in sorted(mode_exports.items()))
            st.caption(f"last_exports: {txt}")


# ---------------------------------------------------------------------------
# Cached compute wrappers
# ---------------------------------------------------------------------------


@st.cache_data(show_spinner=False, max_entries=64)
def _cached_compute_model(params_json):
    return compute_model(json.loads(params_json))


@st.cache_data(show_spinner=False, max_entries=64)
def _cached_compute_sed(params_json):
    return compute_sed(json.loads(params_json))


# ---------------------------------------------------------------------------
# UI helpers
# ---------------------------------------------------------------------------


def _render_logo_block(show_mode=True):
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
    if show_mode:
        return st.radio(
            "Mode",
            ["Light Curve", "Spectrum", "Sky Image"],
            horizontal=True,
            key="plot_mode",
        )
    return None


def _render_footer(mode):
    mode_key = mode.lower().replace(" ", "_")
    if st.button("\U0001F517 Share URL", key=f"share_btn_{mode_key}", help="Update URL with current parameters"):
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


def _parse_obs_rows(groups, is_lc):
    rows = []
    x_default = "day" if is_lc else "Hz"
    for g in groups:
        if not g.get("visible", True):
            continue
        label = g.get("legend", "") or "data"
        x_unit = g.get("x_unit", x_default)
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
    return tuple(rows)


def _obs_signature(groups):
    sig = []
    for g in groups:
        sig.append(
            (
                bool(g.get("visible", True)),
                str(g.get("legend", "")),
                str(g.get("x_unit", "")),
                str(g.get("y_unit", "")),
                str(g.get("text", "")),
            )
        )
    return tuple(sig)


def _parse_obs_with_cache(groups, is_lc, perf):
    cache_key = "_obs_parse_lc" if is_lc else "_obs_parse_sed"
    sig = _obs_signature(groups)

    t0 = time_mod.perf_counter()
    cached = st.session_state.get(cache_key)
    if cached and cached.get("sig") == sig:
        perf["obs_parse_s"] = time_mod.perf_counter() - t0
        perf["obs_parse_hit"] = True
        return cached.get("rows", ())

    rows = _parse_obs_rows(groups, is_lc)
    st.session_state[cache_key] = {"sig": sig, "rows": rows}
    perf["obs_parse_s"] = time_mod.perf_counter() - t0
    perf["obs_parse_hit"] = False
    return rows


def _parse_obs_file(uploaded_file):
    """Parse uploaded CSV/TXT/DAT into text_area format (x  y  err)."""
    import csv
    import io as _io

    content = uploaded_file.getvalue().decode("utf-8", errors="replace")
    for delim in [",", "\t", " "]:
        try:
            rows = []
            for row in csv.reader(_io.StringIO(content), delimiter=delim):
                cleaned = [c.strip() for c in row if c.strip()]
                if len(cleaned) < 2:
                    continue
                try:
                    float(cleaned[0])
                    float(cleaned[1])
                except ValueError:
                    continue
                rows.append("  ".join(cleaned[:3]))
            if rows:
                return "\n".join(rows)
        except Exception:
            continue
    return ""


def _render_observation_section(is_lc, perf):
    _obs_state_key = "obs_groups_lc" if is_lc else "obs_groups_sed"
    prefix = "lc" if is_lc else "sed"
    if _obs_state_key not in st.session_state:
        st.session_state[_obs_state_key] = []

    _obs_rm_key = f"{_obs_state_key}_rm"
    if _obs_rm_key in st.session_state:
        idx = st.session_state.pop(_obs_rm_key)
        if 0 <= idx < len(st.session_state[_obs_state_key]):
            st.session_state[_obs_state_key].pop(idx)

    with st.expander("Observation Data", expanded=False):
        groups = st.session_state[_obs_state_key]
        if st.button("\u2795 Add group", key=f"obs_add_group_{prefix}"):
            groups.append(
                {
                    "legend": f"data {len(groups) + 1}",
                    "x_unit": "day" if is_lc else "Hz",
                    "y_unit": "mJy",
                    "text": "",
                    "visible": True,
                }
            )

        if groups:
            x_name = "t" if is_lc else "\u03bd"
            tab_labels = [
                g.get("legend", f"Group {i+1}") or f"Group {i+1}"
                for i, g in enumerate(groups)
            ]
            tabs = st.tabs(tab_labels)
            for i, tab in enumerate(tabs):
                g = groups[i]
                with tab:
                    c_vis, c_del = st.columns([5, 1])
                    with c_vis:
                        g["visible"] = st.checkbox(
                            "Show",
                            value=g.get("visible", True),
                            key=f"obs_vis_{prefix}_{i}",
                        )
                    with c_del:
                        if st.button(
                            "\u2716",
                            key=f"obs_del_{prefix}_{i}",
                            help="Remove this group",
                            type="tertiary",
                        ):
                            st.session_state[_obs_rm_key] = i
                            st.rerun()
                    g["legend"] = st.text_input(
                        "Legend",
                        value=g.get("legend", ""),
                        key=f"obs_legend_{prefix}_{i}",
                    )
                    if is_lc:
                        x_units = list(TIME_SCALES.keys())
                        x_default = "day"
                    else:
                        x_units = list(FREQ_SCALES.keys())
                        x_default = "Hz"

                    c1, c2 = st.columns(2)
                    with c1:
                        current_x = g.get("x_unit", x_default)
                        x_index = x_units.index(current_x) if current_x in x_units else 0
                        g["x_unit"] = st.selectbox(
                            f"{x_name} unit",
                            x_units,
                            index=x_index,
                            key=f"obs_xunit_{prefix}_{i}",
                        )
                    with c2:
                        current_y = g.get("y_unit", "mJy")
                        y_index = OBS_FLUX_UNITS.index(current_y) if current_y in OBS_FLUX_UNITS else 0
                        g["y_unit"] = st.selectbox(
                            "Flux unit",
                            OBS_FLUX_UNITS,
                            index=y_index,
                            key=f"obs_yunit_{prefix}_{i}",
                        )
                    g["text"] = st.text_area(
                        f"{x_name}, flux, err",
                        value=g.get("text", ""),
                        height=100,
                        placeholder=f"{x_name}  flux  err\n1e4  0.5  0.1\n1e5  0.3  0.05",
                        help="One row per line. Columns: "
                        f"{x_name}, flux, err (optional). "
                        "Separated by space, tab, or comma.",
                        key=f"obs_text_{prefix}_{i}",
                    )
                    with st.popover("Upload file"):
                        uploaded = st.file_uploader(
                            "Upload file",
                            type=["csv", "txt", "dat"],
                            key=f"obs_file_{prefix}_{i}",
                        )
                        if uploaded is not None:
                            parsed = _parse_obs_file(uploaded)
                            if parsed:
                                g["text"] = parsed
                                st.session_state[f"obs_text_{prefix}_{i}"] = parsed

        st.session_state[_obs_state_key] = groups

    return _parse_obs_with_cache(st.session_state[_obs_state_key], is_lc, perf)


def _render_shared_controls(plot_mode, perf):
    is_lc = plot_mode == "Light Curve"
    is_sky = plot_mode == "Sky Image"

    c1, c2 = st.columns(2)
    with c1:
        d_L_mpc = st.number_input(
            r"$d_L$ (Mpc)",
            min_value=0.1,
            max_value=1e6,
            value=100.0,
            step=10.0,
            key="d_L_mpc",
        )
    with c2:
        theta_obs = st.slider(
            r"$\theta_{\rm obs}$ (rad)",
            0.0,
            1.57,
            0.0,
            0.01,
            key="theta_obs",
        )
    z = z_from_lumi_dist_mpc(d_L_mpc)
    d_L_cm = d_L_mpc * 3.0856775814913673e24

    if not is_sky:
        c1, c2 = st.columns(2)
        with c1:
            _flux_choices = list(FLUX_SCALES.keys()) + ["AB mag"]
            flux_unit = st.selectbox("Flux", _flux_choices, key="flux_unit")
        with c2:
            if is_lc:
                time_unit = st.selectbox("Time", list(TIME_SCALES.keys()), key="time_unit")
            else:
                time_unit = "s"
    else:
        flux_unit = "cgs"
        time_unit = "s"

    with st.container(border=True):
        c1, c2 = st.columns([3, 2])
        with c1:
            jet_type = st.selectbox(
                "Jet",
                ["Top-hat", "Gaussian", "Power-law", "Two-component"],
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

    with st.container(border=True):
        c1, c2 = st.columns([2, 3])
        with c1:
            medium_type = st.selectbox("Medium", ["ISM", "Wind"], key="medium_type")
        with c2:
            if medium_type == "ISM":
                n_ism = log_slider(r"$n_{\rm ism}$ (cm⁻³)", -5.0, 5.0, 0.0, key="n_ism")
                A_star = 0.1
            else:
                A_star = log_slider(r"$A_*$", -3.0, 2.0, -1.0, key="A_star")
                n_ism = 1.0
        if medium_type == "Wind":
            k_m = st.slider(r"$k_m$", 0.0, 4.0, 2.0, 0.1, key="k_m")
        else:
            k_m = 2.0

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

    if not is_sky:
        with st.expander("Instruments", expanded=False):
            selected_instruments = st.multiselect(
                "Show sensitivity",
                list(INSTRUMENTS.keys()),
                default=[],
                key="instruments",
            )
    else:
        selected_instruments = []

    obs_data_tuple = ()
    if plot_mode in ("Light Curve", "Spectrum"):
        obs_data_tuple = _render_observation_section(is_lc, perf)

    with st.expander("Resolutions", expanded=False):
        num_t = st.slider("Time points", 50, 300, 100, 10, key="num_t")
        c1, c2, c3 = st.columns(3)
        with c1:
            res_phi = st.slider(r"$\phi$ ppd", 0.05, 1.0, 0.1, 0.05, key="res_phi")
        with c2:
            res_theta = st.slider(r"$\theta$ ppd", 0.1, 2.0, 0.25, 0.1, key="res_theta")
        with c3:
            res_t = st.slider("t ppd", 1.0, 20.0, 10.0, 0.5, key="res_t")

    _render_footer(plot_mode)

    return {
        "d_L_cm": d_L_cm,
        "z": z,
        "theta_obs": theta_obs,
        "flux_unit": flux_unit,
        "time_unit": time_unit,
        "jet_type": jet_type,
        "theta_c": theta_c,
        "E_iso": E_iso,
        "Gamma0": Gamma0,
        "spreading": spreading,
        "duration": duration,
        "k_e": k_e,
        "k_g": k_g,
        "theta_w": theta_w,
        "E_iso_w": E_iso_w,
        "Gamma0_w": Gamma0_w,
        "medium_type": medium_type,
        "n_ism": n_ism,
        "A_star": A_star,
        "k_m": k_m,
        "eps_e": eps_e,
        "eps_B": eps_B,
        "p": p_val,
        "xi_e": xi_e,
        "ssc": ssc,
        "kn": kn,
        "enable_rvs": enable_rvs,
        "eps_e_r": eps_e_r,
        "eps_B_r": eps_B_r,
        "p_r": p_r,
        "xi_e_r": xi_e_r,
        "rvs_ssc": rvs_ssc,
        "rvs_kn": rvs_kn,
        "selected_instruments": selected_instruments,
        "obs_data_tuple": obs_data_tuple,
        "num_t": num_t,
        "res_phi": res_phi,
        "res_theta": res_theta,
        "res_t": res_t,
    }


def _physics_params_from_shared(shared):
    return {
        "jet_type": shared["jet_type"],
        "theta_c": shared["theta_c"],
        "E_iso": shared["E_iso"],
        "Gamma0": shared["Gamma0"],
        "spreading": shared["spreading"],
        "duration": shared["duration"],
        "k_e": shared["k_e"],
        "k_g": shared["k_g"],
        "theta_w": shared["theta_w"],
        "E_iso_w": shared["E_iso_w"],
        "Gamma0_w": shared["Gamma0_w"],
        "medium_type": shared["medium_type"],
        "n_ism": shared["n_ism"],
        "A_star": shared["A_star"],
        "k_m": shared["k_m"],
        "d_L_cm": shared["d_L_cm"],
        "z": shared["z"],
        "theta_obs": shared["theta_obs"],
        "eps_e": shared["eps_e"],
        "eps_B": shared["eps_B"],
        "p": shared["p"],
        "xi_e": shared["xi_e"],
        "ssc": shared["ssc"],
        "kn": shared["kn"],
        "enable_rvs": shared["enable_rvs"],
        "eps_e_r": shared["eps_e_r"],
        "eps_B_r": shared["eps_B_r"],
        "p_r": shared["p_r"],
        "xi_e_r": shared["xi_e_r"],
        "rvs_ssc": shared["rvs_ssc"],
        "rvs_kn": shared["rvs_kn"],
        "res_phi": shared["res_phi"],
        "res_theta": shared["res_theta"],
        "res_t": shared["res_t"],
    }


def _render_plot_actions(fig, downloads, prefix, cite_key, main_slot, render_token):
    fig.update_layout(height=_PLOT_HEIGHT)
    # Keep the plotting area size stable across rerenders.
    fig.update_layout(margin=dict(autoexpand=False))
    fig.for_each_xaxis(lambda ax: ax.update(automargin=False))
    fig.for_each_yaxis(lambda ax: ax.update(automargin=False))
    with main_slot:
        _PLOTLY_LIVE(
            figure=_plotly_component_figure(fig),
            height=_PLOT_HEIGHT,
            revision=render_token,
            config={"toImageButtonOptions": {"format": "png", "scale": 3}},
            key=f"plot_{prefix}",
            default=None,
        )

    cols = st.columns(len(downloads) + 1)
    for i, spec in enumerate(downloads):
        label, metric_key, producer, ext, mime = spec
        with cols[i]:
            st.download_button(
                label,
                data=lambda fn=producer, mk=metric_key: _timed_export_payload(mk, fn),
                file_name=f"{prefix}.{ext}",
                mime=mime,
                key=f"dl_{metric_key}_{prefix}",
                use_container_width=True,
                on_click="ignore",
            )

    with cols[-1]:
        if st.button("Cite", key=f"cite_{cite_key}", use_container_width=True):
            _stc.html(CLIPBOARD_JS, height=0)
            st.toast("\u2714 BibTeX copied to clipboard!")


# ---------------------------------------------------------------------------
# Page config and static style
# ---------------------------------------------------------------------------


st.set_page_config(
    page_title="VegasAfterglow",
    page_icon=_PILImage.open(_FAVICON),
    layout="wide",
)

apply_query_params()
st.markdown(SIDEBAR_CSS, unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# Light Curve fragment
# ---------------------------------------------------------------------------


@st.fragment
def render_lightcurve_fragment(main_slot):
    perf = {}

    if True:
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

        shared = _render_shared_controls("Light Curve", perf)

    physics_params = _physics_params_from_shared(shared)

    frequencies, bands = [], []
    tokens = [t.strip() for t in re.split(r",(?![^\[]*\])", nu_str) if t.strip()]
    for part in tokens:
        try:
            entry = parse_entry(part)
            (bands if isinstance(entry, tuple) else frequencies).append(entry)
        except Exception:
            st.warning(f"Unknown frequency or filter: '{part}'")
    if not frequencies and not bands:
        frequencies = [1e9]

    tk0 = time_mod.perf_counter()
    params = {
        **physics_params,
        "frequencies": sorted(float(nu) for nu in frequencies),
        "bands": [list(b) for b in bands],
        "t_min": float(t_min),
        "t_max": float(t_max),
        "num_t": int(shared["num_t"]),
    }
    params_json = _canonical_json(params)
    compute_key = _make_key("lc_compute", params)
    perf["key_build_s"] = time_mod.perf_counter() - tk0

    try:
        perf["compute_hit"] = _seen_compute_key("lc", compute_key)
        tc0 = time_mod.perf_counter()
        data = _cached_compute_model(params_json)
        perf["compute_s"] = time_mod.perf_counter() - tc0
    except Exception as e:
        st.error(f"Computation failed: {e}")
        st.stop()

    selected_instruments = shared["selected_instruments"]
    obs_data_tuple = shared["obs_data_tuple"]
    flux_unit = shared["flux_unit"]
    time_unit = shared["time_unit"]

    has_fband_inst = any(INSTRUMENTS[n][3] == "Fband" for n in selected_instruments)
    has_fband_obs = any(len(r) >= 6 and r[5] == "erg/cm²/s" for r in obs_data_tuple)
    need_sec = has_fband_inst or has_fband_obs
    use_sec = need_sec or len(data["band_data"]) > 0
    is_mag = flux_unit == "AB mag"

    view_payload = {
        "compute_key": compute_key,
        "flux_unit": flux_unit,
        "time_unit": time_unit,
        "t_min": float(t_min),
        "t_max": float(t_max),
        "selected_instruments": list(selected_instruments),
        "obs_data_tuple": [list(r) for r in obs_data_tuple],
    }
    view_key = _make_key("lc_view", view_payload)

    tf0 = time_mod.perf_counter()
    fig, fig_hit = _figure_cache_get("lc", view_key)
    if fig is None:
        fig = make_figure(
            data,
            flux_unit,
            time_unit,
            t_min,
            t_max,
            need_secondary_y=need_sec,
        )
        if selected_instruments and not is_mag:
            t_scale = TIME_SCALES[time_unit]
            _add_sensitivity_traces(
                fig,
                selected_instruments,
                "lightcurve",
                flux_scale=FLUX_SCALES[flux_unit],
                x_range=[t_min / t_scale, t_max / t_scale],
                has_secondary=use_sec,
            )
        if obs_data_tuple:
            _add_obs_traces(
                fig,
                obs_data_tuple,
                flux_unit,
                time_unit,
                has_secondary=use_sec,
            )
        _figure_cache_put("lc", view_key, fig)
    perf["figure_s"] = time_mod.perf_counter() - tf0
    perf["figure_hit"] = fig_hit

    axis_key = _make_key(
        "lc_axis",
        {
            "flux_unit": flux_unit,
            "time_unit": time_unit,
            "need_secondary": bool(need_sec),
            "is_mag": bool(is_mag),
            "t_min": float(t_min),
            "t_max": float(t_max),
        },
    )
    fig.update_layout(uirevision=axis_key)
    _lock_y_limits(fig, "lc", axis_key)

    export_unit = "cgs" if is_mag else flux_unit
    fig_json = fig.to_json()
    downloads = [
        (
            "CSV",
            "lc_csv",
            lambda d=data, u=export_unit, tu=time_unit: export_csv(d, u, tu),
            "csv",
            "text/csv",
        ),
        (
            "JSON",
            "lc_json",
            lambda d=data, u=export_unit, tu=time_unit: export_json(d, u, tu),
            "json",
            "application/json",
        ),
        (
            "PNG",
            "lc_png",
            lambda fj=fig_json: go.Figure(json.loads(fj)).to_image(format="png", scale=3),
            "png",
            "image/png",
        ),
    ]
    _render_plot_actions(fig, downloads, "afterglow_lightcurve", "lightcurve", main_slot, view_key)
    _show_perf_debug("lc", perf)


# ---------------------------------------------------------------------------
# Spectrum fragment
# ---------------------------------------------------------------------------


@st.fragment
def render_spectrum_fragment(main_slot):
    perf = {}

    if True:
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

        shared = _render_shared_controls("Spectrum", perf)

    t_snapshots = []
    for tok in t_snap_str.split(","):
        tok = tok.strip()
        if tok:
            try:
                t_snapshots.append(float(tok))
            except ValueError:
                st.warning(f"Invalid time value: '{tok}'")

    if not t_snapshots:
        t_snapshots = [1e4]

    if sed_nu_min >= sed_nu_max:
        st.warning("nu_min must be less than nu_max")
        st.stop()

    physics_params = _physics_params_from_shared(shared)

    tk0 = time_mod.perf_counter()
    params = {
        **physics_params,
        "t_snapshots": sorted(t_snapshots),
        "nu_min": float(sed_nu_min),
        "nu_max": float(sed_nu_max),
        "num_nu": int(num_nu),
    }
    params_json = _canonical_json(params)
    compute_key = _make_key("sed_compute", params)
    perf["key_build_s"] = time_mod.perf_counter() - tk0

    try:
        perf["compute_hit"] = _seen_compute_key("sed", compute_key)
        tc0 = time_mod.perf_counter()
        data = _cached_compute_sed(params_json)
        perf["compute_s"] = time_mod.perf_counter() - tc0
    except Exception as e:
        st.error(f"Computation failed: {e}")
        st.stop()

    flux_unit = shared["flux_unit"]
    selected_instruments = shared["selected_instruments"]
    obs_data_tuple = shared["obs_data_tuple"]

    is_mag = flux_unit == "AB mag"
    if is_mag and show_nufnu:
        st.warning("AB mag is not compatible with νFν. Showing Fν in AB mag.")
        show_nufnu = False

    has_fband_inst = any(INSTRUMENTS[n][3] == "Fband" for n in selected_instruments)
    has_fband_obs = any(len(r) >= 6 and r[5] == "erg/cm²/s" for r in obs_data_tuple)
    need_sec = has_fband_inst or has_fband_obs

    view_payload = {
        "compute_key": compute_key,
        "flux_unit": flux_unit,
        "freq_unit": freq_unit,
        "show_nufnu": bool(show_nufnu),
        "selected_instruments": list(selected_instruments),
        "obs_data_tuple": [list(r) for r in obs_data_tuple],
    }
    view_key = _make_key("sed_view", view_payload)

    tf0 = time_mod.perf_counter()
    fig, fig_hit = _figure_cache_get("sed", view_key)
    if fig is None:
        fig = make_sed_figure(
            data,
            flux_unit,
            freq_unit,
            nufnu=show_nufnu,
            need_secondary_y=need_sec,
        )
        if selected_instruments and not is_mag:
            _add_sensitivity_traces(
                fig,
                selected_instruments,
                "spectrum",
                freq_scale=FREQ_SCALES[freq_unit],
                flux_scale=FLUX_SCALES[flux_unit],
                nufnu=show_nufnu,
                has_secondary=need_sec,
            )
        if obs_data_tuple:
            _add_obs_traces(
                fig,
                obs_data_tuple,
                flux_unit,
                freq_unit,
                has_secondary=need_sec,
                mode="spectrum",
                nufnu=show_nufnu,
            )
        _figure_cache_put("sed", view_key, fig)
    perf["figure_s"] = time_mod.perf_counter() - tf0
    perf["figure_hit"] = fig_hit

    axis_key = _make_key(
        "sed_axis",
        {
            "flux_unit": flux_unit,
            "freq_unit": freq_unit,
            "show_nufnu": bool(show_nufnu),
            "need_secondary": bool(need_sec),
            "nu_min": float(sed_nu_min),
            "nu_max": float(sed_nu_max),
        },
    )
    fig.update_layout(uirevision=axis_key)
    _lock_y_limits(fig, "sed", axis_key)

    export_unit = "cgs" if is_mag else flux_unit
    fig_json = fig.to_json()
    downloads = [
        (
            "CSV",
            "sed_csv",
            lambda d=data, u=export_unit, fu=freq_unit: export_sed_csv(d, u, fu),
            "csv",
            "text/csv",
        ),
        (
            "JSON",
            "sed_json",
            lambda d=data, u=export_unit, fu=freq_unit: export_sed_json(d, u, fu),
            "json",
            "application/json",
        ),
        (
            "PNG",
            "sed_png",
            lambda fj=fig_json: go.Figure(json.loads(fj)).to_image(format="png", scale=3),
            "png",
            "image/png",
        ),
    ]
    _render_plot_actions(fig, downloads, "afterglow_sed", "sed", main_slot, view_key)
    _show_perf_debug("sed", perf)


# ---------------------------------------------------------------------------
# Sky image mode (non-fragment, unchanged behavior)
# ---------------------------------------------------------------------------


def _skymap_mpl_style(plt):
    """Apply matplotlib rcParams matching the Plotly skymap style."""
    plt.rcParams.update(
        {
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
        }
    )


def render_skymap_mode():
    with st.sidebar:
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

        sky_animate = st.checkbox("Animate", key="sky_animate")
        if sky_animate:
            c1, c2 = st.columns(2)
            with c1:
                sky_t_min = log_slider(r"$t_{\rm min}$ (s)", 2.0, 9.0, 4.0, key="sky_t_min")
            with c2:
                sky_t_max = log_slider(r"$t_{\rm max}$ (s)", 2.0, 9.0, 7.0, key="sky_t_max")
            sky_nframes = st.slider("Frames", 3, 60, 10, 1, key="sky_nframes")
            sky_t_obs = sky_t_min
        else:
            sky_t_obs = log_slider(r"$t_{\rm obs}$ (s)", 2.0, 9.0, 6.0, key="sky_t_obs")
            sky_t_min, sky_t_max, sky_nframes = sky_t_obs, sky_t_obs, 1

        sky_nu_str = st.text_input(
            "Frequency",
            value="1e9",
            key="sky_nu_str",
            placeholder="e.g. 1e9, 1GHz, 1keV",
            help=_FREQ_HELP,
        )

        c1, c2 = st.columns(2)
        with c1:
            sky_fov = log_slider("FOV (μas)", 1.0, 5.0, np.log10(500), key="sky_fov")
        with c2:
            sky_npixel = st.select_slider("Pixels", [64, 128, 256, 512, 1024], value=256, key="sky_npixel")

        shared = _render_shared_controls("Sky Image", perf={})

    physics_params = _physics_params_from_shared(shared)

    sky_nu_str_val = sky_nu_str.strip()
    try:
        sky_nu = parse_entry(sky_nu_str_val) if sky_nu_str_val else 1e9
        if isinstance(sky_nu, tuple):
            st.sidebar.warning("Sky image requires a single frequency, not a band")
            st.stop()
    except Exception:
        st.sidebar.warning(f"Invalid frequency: '{sky_nu_str_val}'")
        sky_nu = 1e9

    params = {
        **physics_params,
        "t_obs": sky_t_obs,
        "nu_obs": float(sky_nu),
        "fov": sky_fov,
        "npixel": sky_npixel,
    }
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
        all_pos = np.concatenate([img[img > 0] for img in data["images"] if np.any(img > 0)])
        vmin = float(np.log10(all_pos.min())) if all_pos.size > 0 else 0
        vmax = float(np.log10(all_pos.max())) if all_pos.size > 0 else 1

        _skymap_mpl_style(plt)

        fig_mpl, ax = plt.subplots(figsize=(4.0, 3.3), dpi=200, facecolor="white")
        log_img0 = np.where(data["images"][0] > 0, np.log10(data["images"][0]), np.nan)
        import matplotlib.cm as _mcm

        cmap = _mcm.get_cmap("inferno").copy()
        cmap.set_bad(color="white")
        im = ax.imshow(log_img0.T, origin="lower", extent=ext_arr, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
        ax.set_xlabel("Δx (μas)")
        ax.set_ylabel("Δy (μas)")
        title = ax.set_title("")
        cb = fig_mpl.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(r"$\log_{10}\,I$ (erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ sr$^{-1}$)")
        fig_mpl.subplots_adjust(left=0.15, right=0.88, bottom=0.13, top=0.92)

        nu_label = freq_label(data["nu_obs"])
        title.set_visible(False)
        canvas = fig_mpl.canvas
        canvas.draw()
        w, h = canvas.get_width_height()

        bbox = ax.get_window_extent(canvas.get_renderer())
        x0, y0 = int(round(bbox.x0)), int(round(bbox.y0))
        x1, y1 = int(round(bbox.x1)), int(round(bbox.y1))
        img_w, img_h = x1 - x0, y1 - y0

        im.set_visible(False)
        ax.set_facecolor("none")
        fig_mpl.patch.set_alpha(0)
        canvas.draw()
        overlay_pil = _PILImage.frombuffer("RGBA", (w, h), canvas.buffer_rgba()).copy()
        im.set_visible(True)
        plt.close(fig_mpl)

        from matplotlib.colors import Normalize

        norm = Normalize(vmin=vmin, vmax=vmax)
        mapped_frames = []
        for i in range(n_frames):
            img_data = data["images"][i]
            log_img = np.where(img_data > 0, np.log10(img_data), np.nan)
            rgba = (cmap(norm(log_img.T)) * 255).astype(np.uint8)
            rgba = rgba[::-1]
            pil_img = _PILImage.fromarray(rgba)
            if pil_img.size != (img_w, img_h):
                pil_img = pil_img.resize((img_w, img_h), _PILImage.NEAREST)
            mapped_frames.append(pil_img)

        from PIL import ImageDraw, ImageFont

        _title_size = int(9 * 200 / 72)
        _font_candidates = [
            "Helvetica",
            "Arial",
            "DejaVuSans",
            "DejaVu Sans",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/TTF/DejaVuSans.ttf",
            "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        ]
        for _fname in _font_candidates:
            try:
                title_font = ImageFont.truetype(_fname, _title_size)
                break
            except OSError:
                continue
        else:
            title_font = ImageFont.load_default(size=_title_size)
        title_y = h - y1 - int(9 * 200 / 72) - 4
        title_color = (51, 51, 51, 255)

        frames_pil = []
        for i in range(n_frames):
            frame = _PILImage.new("RGBA", (w, h), (255, 255, 255, 255))
            frame.paste(mapped_frames[i], (x0, h - y1))
            frame = _PILImage.alpha_composite(frame, overlay_pil)
            t_label = format_time_label(data["t_obs_array"][i])
            txt = f"t = {t_label},  freq = {nu_label}"
            draw = ImageDraw.Draw(frame)
            tw = draw.textlength(txt, font=title_font)
            draw.text(((w - tw) / 2, title_y), txt, fill=title_color, font=title_font)
            frames_pil.append(frame)

        buf = _io.BytesIO()
        frames_pil[0].save(
            buf,
            format="GIF",
            save_all=True,
            append_images=frames_pil[1:],
            duration=frame_ms,
            loop=0,
        )
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
            st.download_button(
                "JSON",
                export_skymap_json(data),
                file_name="afterglow_skymap.json",
                mime="application/json",
                width="stretch",
            )
        with c2:
            st.download_button(
                "GIF",
                gif_bytes,
                file_name="afterglow_skymap.gif",
                mime="image/gif",
                key="dl_gif_skymap",
                width="stretch",
            )
        with c3:
            if st.button("Cite", key="cite_skymap_anim", width="stretch"):
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
        im = ax.imshow(log_img.T, origin="lower", extent=ext_arr, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
        ax.set_xlabel("Δx (μas)")
        ax.set_ylabel("Δy (μas)")
        t_label = format_time_label(data["t_obs_array"][0])
        nu_label = freq_label(data["nu_obs"])
        ax.set_title(f"t = {t_label},  ν = {nu_label}")
        cb = fig_mpl.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(r"$\log_{10}\,I$ (erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ sr$^{-1}$)")
        fig_mpl.subplots_adjust(left=0.15, right=0.88, bottom=0.13, top=0.92)

        buf = _io.BytesIO()
        fig_mpl.savefig(buf, format="png", facecolor="white", bbox_inches="tight", pad_inches=0.03)
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
            st.download_button(
                "JSON",
                export_skymap_json(data),
                file_name="afterglow_skymap.json",
                mime="application/json",
                width="stretch",
            )
        with c2:
            st.download_button(
                "PNG",
                png_bytes,
                file_name="afterglow_skymap.png",
                mime="image/png",
                key="dl_png_skymap_static",
                width="stretch",
            )
        with c3:
            if st.button("Cite", key="cite_skymap_static", width="stretch"):
                _stc.html(CLIPBOARD_JS, height=0)
                st.toast("\u2714 BibTeX copied to clipboard!")


# ---------------------------------------------------------------------------
# App shell routing
# ---------------------------------------------------------------------------


with st.sidebar:
    plot_mode = _render_logo_block(show_mode=True)

main_slot = st.empty()

if plot_mode == "Light Curve":
    with st.sidebar:
        render_lightcurve_fragment(main_slot)
elif plot_mode == "Spectrum":
    with st.sidebar:
        render_spectrum_fragment(main_slot)
else:
    render_skymap_mode()
