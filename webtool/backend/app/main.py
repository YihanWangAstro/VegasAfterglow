from __future__ import annotations

import json
import os
import re
import sys
import tempfile
import time as time_mod
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import plotly.graph_objects as go
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from pydantic import BaseModel, Field

# Avoid matplotlib cache permission issues during frequency color generation.
_mpl_config_dir = os.path.join(tempfile.gettempdir(), "matplotlib")
os.makedirs(_mpl_config_dir, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", _mpl_config_dir)

# Ensure repo root is importable for local dev runs where VegasAfterglow source
# lives in this repo. In container deployments the package is installed via pip.
_HERE = Path(__file__).resolve()
_repo_root = _HERE.parents[3] if len(_HERE.parents) > 3 else None
if _repo_root and (_repo_root / "VegasAfterglow").exists():
    if str(_repo_root) not in sys.path:
        sys.path.insert(0, str(_repo_root))

from .compute import compute_model, compute_sed, compute_skymap
from .constants import FLUX_SCALES, FREQ_SCALES, INSTRUMENTS, OBS_FLUX_UNITS, TIME_SCALES
from .exports import (
    export_csv,
    export_json,
    export_sed_csv,
    export_sed_json,
    export_skymap_gif_base64,
    export_skymap_json,
)
from .figures import _add_obs_traces, _add_sensitivity_traces, make_figure, make_sed_figure
from .helpers import format_time_label, freq_label, parse_entry, z_from_lumi_dist_mpc

try:
    import VegasAfterglow as _va

    APP_VERSION = _va.__version__.split(".dev")[0]
except Exception:
    APP_VERSION = "unknown"


class ObservationGroup(BaseModel):
    legend: str = "data"
    x_unit: str | None = None
    y_unit: str = "mJy"
    text: str = ""
    visible: bool = True


class SharedParams(BaseModel):
    d_L_mpc: float = 100.0
    theta_obs: float = 0.0
    flux_unit: str = "mJy"
    time_unit: str = "s"
    jet_type: str = "Top-hat"
    theta_c: float = 0.1
    E_iso: float = 1e52
    Gamma0: float = 300.0
    spreading: bool = False
    duration: float = 1.0
    k_e: float = 2.0
    k_g: float = 2.0
    theta_w: float = 0.3
    E_iso_w: float = 1e51
    Gamma0_w: float = 100.0
    medium_type: str = "ISM"
    n_ism: float = 1.0
    A_star: float = 0.1
    k_m: float = 2.0
    eps_e: float = 0.1
    eps_B: float = 1e-3
    p: float = 2.3
    xi_e: float = 1.0
    ssc: bool = False
    kn: bool = False
    enable_rvs: bool = False
    eps_e_r: float = 0.1
    eps_B_r: float = 1e-3
    p_r: float = 2.3
    xi_e_r: float = 1.0
    rvs_ssc: bool = False
    rvs_kn: bool = False
    num_t: int = 100
    res_phi: float = 0.1
    res_theta: float = 0.25
    res_t: float = 10.0


class LightCurveRequest(BaseModel):
    shared: SharedParams = Field(default_factory=SharedParams)
    frequencies_input: str = "1e9, R, 1keV"
    t_min: float = 1.0
    t_max: float = 1e8
    selected_instruments: list[str] = Field(default_factory=list)
    observation_groups: list[ObservationGroup] = Field(default_factory=list)
    include_figure: bool = True
    include_exports: bool = True
    export_kinds: list[str] = Field(default_factory=list)


class SpectrumRequest(BaseModel):
    shared: SharedParams = Field(default_factory=SharedParams)
    t_snapshots_input: str = "1e3, 1e4, 1e5, 1e6"
    nu_min: float = 1e8
    nu_max: float = 1e20
    num_nu: int = 200
    freq_unit: str = "Hz"
    show_nufnu: bool = False
    selected_instruments: list[str] = Field(default_factory=list)
    observation_groups: list[ObservationGroup] = Field(default_factory=list)
    include_figure: bool = True
    include_exports: bool = True
    export_kinds: list[str] = Field(default_factory=list)


class SkyMapRequest(BaseModel):
    shared: SharedParams = Field(default_factory=SharedParams)
    animate: bool = False
    t_obs: float = 1e6
    t_min: float = 1e4
    t_max: float = 1e7
    n_frames: int = 10
    nu_input: str = "1e9"
    fov: float = 500.0
    npixel: int = 256
    include_figure: bool = True
    include_exports: bool = True
    export_kinds: list[str] = Field(default_factory=list)


app = FastAPI(title="webtool-api", version="0.2.0")

allowed_origins = [
    origin.strip()
    for origin in os.getenv(
        "ALLOWED_ORIGINS",
        "http://localhost:3000,http://127.0.0.1:3000",
    ).split(",")
    if origin.strip()
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
app.add_middleware(GZipMiddleware, minimum_size=1000, compresslevel=5)


def _default_shared() -> dict[str, Any]:
    return SharedParams().model_dump()


def _default_lightcurve() -> dict[str, Any]:
    return {
        "frequencies_input": "1e9, R, 1keV",
        "t_min": 1.0,
        "t_max": 1e8,
        "selected_instruments": [],
        "observation_groups": [],
    }


def _default_spectrum() -> dict[str, Any]:
    return {
        "t_snapshots_input": "1e3, 1e4, 1e5, 1e6",
        "nu_min": 1e8,
        "nu_max": 1e20,
        "num_nu": 200,
        "freq_unit": "Hz",
        "show_nufnu": False,
        "selected_instruments": [],
        "observation_groups": [],
    }


def _default_skymap() -> dict[str, Any]:
    return {
        "animate": False,
        "t_obs": 1e6,
        "t_min": 1e4,
        "t_max": 1e7,
        "n_frames": 10,
        "nu_input": "1e9",
        "fov": 500.0,
        "npixel": 256,
    }


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


@lru_cache(maxsize=64)
def _cached_compute_model(params_json: str):
    return compute_model(json.loads(params_json))


@lru_cache(maxsize=64)
def _cached_compute_sed(params_json: str):
    return compute_sed(json.loads(params_json))


@lru_cache(maxsize=32)
def _cached_compute_skymap(params_json: str):
    return compute_skymap(json.loads(params_json))


def _shared_to_physics(shared: SharedParams) -> dict[str, Any]:
    d_L_cm = shared.d_L_mpc * 3.0856775814913673e24
    z_val = z_from_lumi_dist_mpc(shared.d_L_mpc)
    return {
        "jet_type": shared.jet_type,
        "theta_c": shared.theta_c,
        "E_iso": shared.E_iso,
        "Gamma0": shared.Gamma0,
        "spreading": shared.spreading,
        "duration": shared.duration,
        "k_e": shared.k_e,
        "k_g": shared.k_g,
        "theta_w": shared.theta_w,
        "E_iso_w": shared.E_iso_w,
        "Gamma0_w": shared.Gamma0_w,
        "medium_type": shared.medium_type,
        "n_ism": shared.n_ism,
        "A_star": shared.A_star,
        "k_m": shared.k_m,
        "d_L_cm": d_L_cm,
        "z": z_val,
        "theta_obs": shared.theta_obs,
        "eps_e": shared.eps_e,
        "eps_B": shared.eps_B,
        "p": shared.p,
        "xi_e": shared.xi_e,
        "ssc": shared.ssc,
        "kn": shared.kn,
        "enable_rvs": shared.enable_rvs,
        "eps_e_r": shared.eps_e_r,
        "eps_B_r": shared.eps_B_r,
        "p_r": shared.p_r,
        "xi_e_r": shared.xi_e_r,
        "rvs_ssc": shared.rvs_ssc,
        "rvs_kn": shared.rvs_kn,
        "res_phi": shared.res_phi,
        "res_theta": shared.res_theta,
        "res_t": shared.res_t,
    }


def _parse_frequency_input(frequencies_input: str) -> tuple[list[float], list[list[Any]], list[str]]:
    frequencies: list[float] = []
    bands: list[list[Any]] = []
    warnings: list[str] = []
    tokens = [tok.strip() for tok in re.split(r",(?![^\[]*\])", frequencies_input) if tok.strip()]

    for token in tokens:
        try:
            entry = parse_entry(token)
            if isinstance(entry, tuple):
                bands.append([float(entry[0]), float(entry[1]), str(entry[2])])
            else:
                frequencies.append(float(entry))
        except Exception:
            warnings.append(f"Unknown frequency or filter: '{token}'")

    if not frequencies and not bands:
        frequencies = [1e9]

    return sorted(frequencies), bands, warnings


def _parse_observation_rows(groups: list[ObservationGroup], is_lc: bool) -> tuple:
    rows = []
    x_default = "day" if is_lc else "Hz"
    for group in groups:
        g = group.model_dump()
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


def _to_fig_json(fig: go.Figure) -> dict[str, Any]:
    # Avoid JSON encode/decode roundtrip on every request.
    return fig.to_plotly_json()


def _validate_shared(shared: SharedParams) -> None:
    if shared.flux_unit not in list(FLUX_SCALES.keys()) + ["AB mag"]:
        raise HTTPException(status_code=400, detail=f"Unsupported flux_unit: {shared.flux_unit}")
    if shared.time_unit not in TIME_SCALES:
        raise HTTPException(status_code=400, detail=f"Unsupported time_unit: {shared.time_unit}")


def _validate_instruments(selected_instruments: list[str]) -> list[str]:
    unknown = [name for name in selected_instruments if name not in INSTRUMENTS]
    if unknown:
        raise HTTPException(status_code=400, detail=f"Unknown instruments: {', '.join(unknown)}")
    return selected_instruments


def _resolve_export_kinds(requested: list[str], defaults: set[str], allowed: set[str]) -> set[str]:
    if not requested:
        return defaults
    normalized = {kind.strip().lower() for kind in requested if kind.strip()}
    unknown = normalized - allowed
    if unknown:
        raise HTTPException(status_code=400, detail=f"Unknown export kind(s): {', '.join(sorted(unknown))}")
    return normalized


def _make_skymap_figure(data: dict[str, Any]) -> go.Figure:
    images = data["images"]
    extent = data["extent_uas"]
    t_arr = data["t_obs_array"]
    nu_obs = data["nu_obs"]
    nu_label = freq_label(float(nu_obs))

    if not images:
        raise HTTPException(status_code=500, detail="Sky map returned no images")

    # Plotly heatmap x/y are cell centers. Build centers from extent edges to avoid
    # half-cell clipping, which is most visible at small FOV values.
    nx = int(images[0].shape[0])
    ny = int(images[0].shape[1])
    x_edges = np.linspace(float(extent[0]), float(extent[1]), nx + 1)
    y_edges = np.linspace(float(extent[2]), float(extent[3]), ny + 1)
    x = 0.5 * (x_edges[:-1] + x_edges[1:])
    y = 0.5 * (y_edges[:-1] + y_edges[1:])

    log_frames = [np.where(img > 0, np.log10(img), np.nan).T for img in images]
    finite_chunks = [frame[np.isfinite(frame)] for frame in log_frames if np.any(np.isfinite(frame))]
    finite_vals = np.concatenate(finite_chunks) if finite_chunks else np.array([])
    if finite_vals.size > 0:
        z_min = float(np.min(finite_vals))
        z_max = float(np.max(finite_vals))
    else:
        z_min, z_max = 0.0, 1.0

    colorbar_title = "log<sub>10</sub> I (erg cm<sup>-2</sup> s<sup>-1</sup> Hz<sup>-1</sup> sr<sup>-1</sup>)"

    fig = go.Figure(
        data=[
            go.Heatmap(
                z=log_frames[0],
                x=x,
                y=y,
                colorscale="Inferno",
                zmin=z_min,
                zmax=z_max,
                colorbar={"title": {"text": colorbar_title, "side": "right"}},
                hovertemplate="Δx=%{x:.2f} μas<br>Δy=%{y:.2f} μas<br>log10 I=%{z:.3f}<extra></extra>",
            )
        ]
    )

    fig.update_layout(
        title=f"t = {format_time_label(float(t_arr[0]))}, ν = {nu_label}",
        xaxis_title="Δx (μas)",
        yaxis_title="Δy (μas)",
        xaxis={"range": [float(extent[0]), float(extent[1])]},
        yaxis={"range": [float(extent[2]), float(extent[3])]},
        template="none",
        plot_bgcolor="#ffffff",
        paper_bgcolor="#ffffff",
        margin={"l": 60, "r": 20, "t": 50, "b": 60},
    )
    axis_style = dict(
        showline=True,
        linewidth=0.8,
        linecolor="#000",
        mirror=True,
        ticks="inside",
        ticklen=5,
        tickwidth=0.8,
        tickcolor="#000",
        showgrid=True,
        gridcolor="rgba(0,0,0,0.10)",
        griddash="dot",
        gridwidth=0.3,
    )
    fig.update_xaxes(constrain="domain", **axis_style)
    fig.update_yaxes(scaleanchor="x", scaleratio=1, constrain="domain", **axis_style)

    if len(log_frames) > 1:
        frame_objs = []
        slider_steps = []
        for i, frame in enumerate(log_frames):
            frame_name = f"frame_{i}"
            frame_objs.append(
                go.Frame(
                    name=frame_name,
                    data=[
                        go.Heatmap(
                            z=frame,
                            x=x,
                            y=y,
                            colorscale="Inferno",
                            zmin=z_min,
                            zmax=z_max,
                            colorbar={"title": {"text": colorbar_title, "side": "right"}},
                        )
                    ],
                    layout={"title": f"t = {format_time_label(float(t_arr[i]))}, ν = {nu_label}"},
                )
            )
            slider_steps.append(
                {
                    "label": format_time_label(float(t_arr[i])),
                    "method": "animate",
                    "args": [[frame_name], {"mode": "immediate", "frame": {"duration": 0, "redraw": True}, "transition": {"duration": 0}}],
                }
            )

        fig.frames = frame_objs
        fig.update_layout(
            margin={"b": 230},
            updatemenus=[
                {
                    "type": "buttons",
                    "showactive": True,
                    "active": 0,
                    "x": 0.53,
                    "xanchor": "center",
                    "y": -0.44,
                    "yanchor": "top",
                    "buttons": [
                        {
                            "label": "Pause/Play",
                            "method": "animate",
                            "args": [
                                None,
                                {
                                    "frame": {"duration": 300, "redraw": True},
                                    "fromcurrent": True,
                                    "transition": {"duration": 0},
                                    "mode": "immediate",
                                },
                            ],
                            "args2": [[None], {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"}],
                        },
                    ],
                }
            ],
            sliders=[
                {
                    "active": 0,
                    "x": 0.13,
                    "y": -0.08,
                    "len": 0.8,
                    "pad": {"t": 8, "b": 0},
                    "currentvalue": {"prefix": "t: "},
                    "steps": slider_steps,
                }
            ],
        )

    return fig


@app.get("/api/health")
def health() -> dict[str, str]:
    return {"status": "ok", "service": "fastapi"}


@app.get("/api/options")
def options() -> dict[str, Any]:
    return {
        "version": APP_VERSION,
        "flux_units": list(FLUX_SCALES.keys()) + ["AB mag"],
        "time_units": list(TIME_SCALES.keys()),
        "freq_units": list(FREQ_SCALES.keys()),
        "obs_flux_units": OBS_FLUX_UNITS,
        "instruments": sorted(list(INSTRUMENTS.keys())),
    }


@app.get("/api/defaults")
def defaults() -> dict[str, Any]:
    return {
        "shared": _default_shared(),
        "lightcurve": _default_lightcurve(),
        "spectrum": _default_spectrum(),
        "skymap": _default_skymap(),
    }


@app.post("/api/lightcurve")
def lightcurve(req: LightCurveRequest) -> dict[str, Any]:
    _validate_shared(req.shared)
    selected_instruments = _validate_instruments(req.selected_instruments)

    if req.t_min <= 0 or req.t_max <= 0:
        raise HTTPException(status_code=400, detail="t_min and t_max must be positive")
    if req.t_min >= req.t_max:
        raise HTTPException(status_code=400, detail="t_min must be less than t_max")

    frequencies, bands, parse_warnings = _parse_frequency_input(req.frequencies_input)

    params = {
        **_shared_to_physics(req.shared),
        "frequencies": frequencies,
        "bands": bands,
        "t_min": float(req.t_min),
        "t_max": float(req.t_max),
        "num_t": int(req.shared.num_t),
    }

    params_json = _canonical_json(params)

    t0 = time_mod.perf_counter()
    try:
        data = _cached_compute_model(params_json)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Light curve computation failed: {exc}") from exc
    compute_s = time_mod.perf_counter() - t0

    obs_rows = _parse_observation_rows(req.observation_groups, is_lc=True)

    has_fband_inst = any(INSTRUMENTS[name][3] == "Fband" for name in selected_instruments)
    has_fband_obs = any(len(row) >= 6 and row[5] == "erg/cm²/s" for row in obs_rows)
    need_secondary = has_fband_inst or has_fband_obs
    use_secondary = need_secondary or len(data["band_data"]) > 0
    is_mag = req.shared.flux_unit == "AB mag"

    export_unit = "cgs" if is_mag else req.shared.flux_unit
    response: dict[str, Any] = {
        "meta": {
            "compute_seconds": compute_s,
            "frequencies_hz": frequencies,
            "bands": bands,
            "warnings": parse_warnings,
        },
    }
    if req.include_figure:
        fig = make_figure(
            data,
            req.shared.flux_unit,
            req.shared.time_unit,
            req.t_min,
            req.t_max,
            need_secondary_y=need_secondary,
        )

        if selected_instruments and not is_mag:
            t_scale = TIME_SCALES[req.shared.time_unit]
            _add_sensitivity_traces(
                fig,
                selected_instruments,
                "lightcurve",
                flux_scale=FLUX_SCALES[req.shared.flux_unit],
                x_range=[req.t_min / t_scale, req.t_max / t_scale],
                has_secondary=use_secondary,
            )

        if obs_rows:
            _add_obs_traces(fig, obs_rows, req.shared.flux_unit, req.shared.time_unit, has_secondary=use_secondary)
        response["figure"] = _to_fig_json(fig)

    if req.include_exports:
        kinds = _resolve_export_kinds(req.export_kinds, defaults={"csv", "json"}, allowed={"csv", "json"})
        exports: dict[str, str] = {}
        if "csv" in kinds:
            exports["csv"] = export_csv(data, export_unit, req.shared.time_unit)
        if "json" in kinds:
            exports["json"] = export_json(data, export_unit, req.shared.time_unit)
        response["exports"] = exports

    return response


@app.post("/api/spectrum")
def spectrum(req: SpectrumRequest) -> dict[str, Any]:
    _validate_shared(req.shared)
    selected_instruments = _validate_instruments(req.selected_instruments)

    if req.freq_unit not in FREQ_SCALES:
        raise HTTPException(status_code=400, detail=f"Unsupported freq_unit: {req.freq_unit}")

    if req.nu_min <= 0 or req.nu_max <= 0:
        raise HTTPException(status_code=400, detail="nu_min and nu_max must be positive")
    if req.nu_min >= req.nu_max:
        raise HTTPException(status_code=400, detail="nu_min must be less than nu_max")

    t_snapshots = []
    warnings = []
    for tok in req.t_snapshots_input.split(","):
        tok = tok.strip()
        if tok:
            try:
                t_snapshots.append(float(tok))
            except ValueError:
                warnings.append(f"Invalid time value: '{tok}'")
    if not t_snapshots:
        t_snapshots = [1e4]

    params = {
        **_shared_to_physics(req.shared),
        "t_snapshots": sorted(t_snapshots),
        "nu_min": float(req.nu_min),
        "nu_max": float(req.nu_max),
        "num_nu": int(req.num_nu),
    }
    params_json = _canonical_json(params)

    t0 = time_mod.perf_counter()
    try:
        data = _cached_compute_sed(params_json)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Spectrum computation failed: {exc}") from exc
    compute_s = time_mod.perf_counter() - t0

    obs_rows = _parse_observation_rows(req.observation_groups, is_lc=False)

    is_mag = req.shared.flux_unit == "AB mag"
    show_nufnu = req.show_nufnu
    if is_mag and show_nufnu:
        show_nufnu = False
        warnings.append("AB mag is not compatible with νFν. Showing Fν in AB mag.")

    has_fband_inst = any(INSTRUMENTS[name][3] == "Fband" for name in selected_instruments)
    has_fband_obs = any(len(row) >= 6 and row[5] == "erg/cm²/s" for row in obs_rows)
    need_secondary = has_fband_inst or has_fband_obs

    export_unit = "cgs" if is_mag else req.shared.flux_unit
    response: dict[str, Any] = {
        "meta": {
            "compute_seconds": compute_s,
            "t_snapshots_s": sorted(t_snapshots),
            "warnings": warnings,
        },
    }
    if req.include_figure:
        fig = make_sed_figure(
            data,
            req.shared.flux_unit,
            req.freq_unit,
            nufnu=show_nufnu,
            need_secondary_y=need_secondary,
        )

        if selected_instruments and not is_mag:
            _add_sensitivity_traces(
                fig,
                selected_instruments,
                "spectrum",
                freq_scale=FREQ_SCALES[req.freq_unit],
                flux_scale=FLUX_SCALES[req.shared.flux_unit],
                nufnu=show_nufnu,
                has_secondary=need_secondary,
            )

        if obs_rows:
            _add_obs_traces(
                fig,
                obs_rows,
                req.shared.flux_unit,
                req.freq_unit,
                has_secondary=need_secondary,
                mode="spectrum",
                nufnu=show_nufnu,
            )
        response["figure"] = _to_fig_json(fig)

    if req.include_exports:
        kinds = _resolve_export_kinds(req.export_kinds, defaults={"csv", "json"}, allowed={"csv", "json"})
        exports: dict[str, str] = {}
        if "csv" in kinds:
            exports["csv"] = export_sed_csv(data, export_unit, req.freq_unit)
        if "json" in kinds:
            exports["json"] = export_sed_json(data, export_unit, req.freq_unit)
        response["exports"] = exports

    return response


@app.post("/api/skymap")
def skymap(req: SkyMapRequest) -> dict[str, Any]:
    _validate_shared(req.shared)

    if req.fov <= 0:
        raise HTTPException(status_code=400, detail="fov must be positive")
    if req.npixel <= 1:
        raise HTTPException(status_code=400, detail="npixel must be greater than 1")

    try:
        nu_entry = parse_entry(req.nu_input.strip()) if req.nu_input.strip() else 1e9
        if isinstance(nu_entry, tuple):
            raise HTTPException(status_code=400, detail="Sky image requires a single frequency, not a band")
        nu_obs = float(nu_entry)
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Invalid frequency: '{req.nu_input}'") from exc

    params = {
        **_shared_to_physics(req.shared),
        "t_obs": float(req.t_obs),
        "nu_obs": nu_obs,
        "fov": float(req.fov),
        "npixel": int(req.npixel),
    }

    if req.animate:
        if req.t_min <= 0 or req.t_max <= 0:
            raise HTTPException(status_code=400, detail="t_min and t_max must be positive")
        if req.t_min >= req.t_max:
            raise HTTPException(status_code=400, detail="t_min must be less than t_max")
        if req.n_frames < 2:
            raise HTTPException(status_code=400, detail="n_frames must be >= 2 when animate=true")
        params["t_obs_array"] = np.logspace(np.log10(req.t_min), np.log10(req.t_max), int(req.n_frames)).tolist()

    params_json = _canonical_json(params)

    t0 = time_mod.perf_counter()
    try:
        data = _cached_compute_skymap(params_json)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Sky image computation failed: {exc}") from exc
    compute_s = time_mod.perf_counter() - t0

    response: dict[str, Any] = {
        "meta": {
            "compute_seconds": compute_s,
            "n_frames": len(data.get("images", [])),
        },
    }
    if req.include_figure:
        fig = _make_skymap_figure(data)
        response["figure"] = _to_fig_json(fig)
    if req.include_exports:
        kinds = _resolve_export_kinds(req.export_kinds, defaults={"json"}, allowed={"json", "gif"})
        exports: dict[str, str] = {}
        if "json" in kinds:
            exports["json"] = export_skymap_json(data)
        if "gif" in kinds:
            exports["gif"] = export_skymap_gif_base64(data)
        response["exports"] = exports
    return response
