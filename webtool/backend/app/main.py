from __future__ import annotations

import os
import socket
import sys
import tempfile
import time as time_mod
from contextlib import asynccontextmanager
from functools import lru_cache
from pathlib import Path
from typing import Any
from urllib.error import URLError
from urllib.request import Request, urlopen

import numpy as np
import orjson
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import ORJSONResponse

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

from .constants import FLUX_SCALES, FREQ_SCALES, INSTRUMENTS, OBS_FLUX_UNITS, TIME_SCALES
from .exports import (
    export_csv,
    export_json,
    export_sed_csv,
    export_sed_json,
    export_skymap_gif_base64,
    export_skymap_json,
)
from .figures import _add_obs_traces, _add_sensitivity_traces, make_figure, make_sed_figure, make_skymap_figure
from .helpers import parse_entry
from .schemas import LightCurveRequest, SharedParams, SkyMapRequest, SpectrumRequest
from .services.cache import cached_compute_model, cached_compute_sed, cached_compute_skymap, canonical_json
from .services.parsing import parse_frequency_input, parse_observation_rows, shared_to_physics
from .services.validation import resolve_export_kinds, validate_instruments, validate_shared

try:
    import VegasAfterglow as _va

    APP_VERSION = _va.__version__.split(".dev")[0]
except Exception:
    APP_VERSION = "unknown"

_FLUX_UNITS = list(FLUX_SCALES.keys()) + ["AB mag"]
_TIME_UNITS = list(TIME_SCALES.keys())
_FREQ_UNITS = list(FREQ_SCALES.keys())


@asynccontextmanager
async def _lifespan(app: FastAPI):
    # Warm up: run a minimal computation so cold-start latency is paid at boot,
    # not on the first user request.
    try:
        _warmup_params = {
            **shared_to_physics(SharedParams()),
            "frequencies": [1e9],
            "bands": [],
            "t_min": 1.0,
            "t_max": 1e6,
            "num_t": 2,
        }
        cached_compute_model(canonical_json(_warmup_params))
    except Exception:
        pass
    yield


app = FastAPI(title="webtool-api", version="0.2.0", default_response_class=ORJSONResponse, lifespan=_lifespan)

_default_origins = {
    "http://localhost:3000",
    "http://127.0.0.1:3000",
    "https://vegasafterglow.com",
    "https://www.vegasafterglow.com",
    "https://vegasafterglow.vercel.app",
}
_env_origins = {
    origin.strip()
    for origin in os.getenv(
        "ALLOWED_ORIGINS",
        "",
    ).split(",")
    if origin.strip()
}
allowed_origins = sorted(_default_origins | _env_origins)

app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
app.add_middleware(GZipMiddleware, minimum_size=200, compresslevel=6)


@lru_cache(maxsize=1)
def _detect_region() -> str:
    for key in ("SERVER_REGION", "K_REGION", "CLOUD_RUN_REGION", "GOOGLE_CLOUD_REGION"):
        value = os.getenv(key, "").strip()
        if value:
            return value

    # Cloud Run may not always expose region in env vars; use metadata server.
    metadata_url = "http://metadata.google.internal/computeMetadata/v1/instance/region"
    request = Request(metadata_url, headers={"Metadata-Flavor": "Google"})
    try:
        with urlopen(request, timeout=0.5) as response:
            raw = response.read().decode("utf-8").strip()
            if raw:
                return raw.rsplit("/", 1)[-1]
    except (URLError, TimeoutError, OSError, ValueError):
        pass
    return ""



@app.get("/api/health")
def health() -> dict[str, str]:
    region = _detect_region()
    service_name = os.getenv("K_SERVICE", "")
    revision = os.getenv("K_REVISION", "")
    instance = socket.gethostname()

    payload: dict[str, str] = {"status": "ok", "service": "fastapi"}
    if region:
        payload["region"] = region
        region_to_location = {
            "us-west2": "Los Angeles",
            "us-east4": "Northern Virginia",
            "europe-west4": "Netherlands",
            "asia-east2": "Hong Kong",
            "asia-northeast1": "Tokyo",
        }
        payload["location_label"] = region_to_location.get(region, region)
    if service_name:
        payload["service_name"] = service_name
    if revision:
        payload["revision"] = revision
    if instance and instance not in {"localhost", "127.0.0.1"}:
        payload["instance"] = instance
    return payload


@app.get("/api/options")
def options() -> dict[str, Any]:
    return {
        "version": APP_VERSION,
        "flux_units": _FLUX_UNITS,
        "time_units": _TIME_UNITS,
        "freq_units": _FREQ_UNITS,
        "obs_flux_units": OBS_FLUX_UNITS,
        "instruments": sorted(INSTRUMENTS),
    }


@app.get("/api/defaults")
def defaults() -> dict[str, Any]:
    return {
        "shared": SharedParams().model_dump(),
        "lightcurve": {
            "frequencies_input": "1e9, R, 1keV",
            "t_min": 1.0,
            "t_max": 1e8,
            "selected_instruments": [],
            "observation_groups": [],
        },
        "spectrum": {
            "t_snapshots_input": "1e3, 1e4, 1e5, 1e6",
            "nu_min": 1e8,
            "nu_max": 1e20,
            "num_nu": 200,
            "freq_unit": "Hz",
            "show_nufnu": False,
            "selected_instruments": [],
            "observation_groups": [],
        },
        "skymap": {
            "animate": False,
            "t_obs": 1e6,
            "t_min": 1e4,
            "t_max": 1e7,
            "n_frames": 15,
            "nu_input": "1e9",
            "fov": 500.0,
            "npixel": 256,
        },
    }


@app.post("/api/lightcurve")
def lightcurve(req: LightCurveRequest) -> dict[str, Any]:
    validate_shared(req.shared)
    selected_instruments = validate_instruments(req.selected_instruments)

    if req.t_min <= 0 or req.t_max <= 0:
        raise HTTPException(status_code=400, detail="t_min and t_max must be positive")
    if req.t_min >= req.t_max:
        raise HTTPException(status_code=400, detail="t_min must be less than t_max")

    frequencies, bands, parse_warnings = parse_frequency_input(req.frequencies_input)

    params = {
        **shared_to_physics(req.shared),
        "frequencies": frequencies,
        "bands": bands,
        "t_min": req.t_min,
        "t_max": req.t_max,
        "num_t": req.shared.num_t,
    }

    params_json = canonical_json(params)

    t0 = time_mod.perf_counter()
    try:
        data = cached_compute_model(params_json)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Light curve computation failed: {exc}") from exc
    compute_s = time_mod.perf_counter() - t0

    obs_rows = parse_observation_rows(req.observation_groups, is_lc=True)

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
        response["figure"] = orjson.loads(fig.to_json(engine="orjson"))

    if req.include_exports:
        kinds = resolve_export_kinds(req.export_kinds, defaults={"csv", "json"}, allowed={"csv", "json"})
        exports: dict[str, str] = {}
        if "csv" in kinds:
            exports["csv"] = export_csv(data, export_unit, req.shared.time_unit)
        if "json" in kinds:
            exports["json"] = export_json(data, export_unit, req.shared.time_unit)
        response["exports"] = exports

    return response


@app.post("/api/spectrum")
def spectrum(req: SpectrumRequest) -> dict[str, Any]:
    validate_shared(req.shared)
    selected_instruments = validate_instruments(req.selected_instruments)

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
        raise HTTPException(status_code=400, detail="t_snapshots_input must include at least one time")

    params = {
        **shared_to_physics(req.shared),
        "t_snapshots": sorted(t_snapshots),
        "nu_min": req.nu_min,
        "nu_max": req.nu_max,
        "num_nu": req.num_nu,
    }
    params_json = canonical_json(params)

    t0 = time_mod.perf_counter()
    try:
        data = cached_compute_sed(params_json)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Spectrum computation failed: {exc}") from exc
    compute_s = time_mod.perf_counter() - t0

    obs_rows = parse_observation_rows(req.observation_groups, is_lc=False)

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
        response["figure"] = orjson.loads(fig.to_json(engine="orjson"))

    if req.include_exports:
        kinds = resolve_export_kinds(req.export_kinds, defaults={"csv", "json"}, allowed={"csv", "json"})
        exports: dict[str, str] = {}
        if "csv" in kinds:
            exports["csv"] = export_sed_csv(data, export_unit, req.freq_unit)
        if "json" in kinds:
            exports["json"] = export_sed_json(data, export_unit, req.freq_unit)
        response["exports"] = exports

    return response


@app.post("/api/skymap")
def skymap(req: SkyMapRequest) -> dict[str, Any]:
    validate_shared(req.shared)

    if req.fov <= 0:
        raise HTTPException(status_code=400, detail="fov must be positive")
    if req.npixel <= 1:
        raise HTTPException(status_code=400, detail="npixel must be greater than 1")
    if req.npixel > 1024:
        raise HTTPException(status_code=400, detail="npixel must be <= 1024")

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
        **shared_to_physics(req.shared),
        "t_obs": req.t_obs,
        "nu_obs": nu_obs,
        "fov": req.fov,
        "npixel": req.npixel,
    }

    if req.animate:
        if req.t_min <= 0 or req.t_max <= 0:
            raise HTTPException(status_code=400, detail="t_min and t_max must be positive")
        if req.t_min >= req.t_max:
            raise HTTPException(status_code=400, detail="t_min must be less than t_max")
        if req.n_frames < 2:
            raise HTTPException(status_code=400, detail="n_frames must be >= 2 when animate=true")
        if req.n_frames > 30:
            raise HTTPException(status_code=400, detail="n_frames must be <= 30 when animate=true")
        if req.npixel > 512:
            raise HTTPException(status_code=400, detail="npixel must be <= 512 when animate=true")
        params["t_obs_array"] = np.logspace(np.log10(req.t_min), np.log10(req.t_max), req.n_frames).tolist()

    params_json = canonical_json(params)

    t0 = time_mod.perf_counter()
    try:
        data = cached_compute_skymap(params_json)
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
        fig = make_skymap_figure(data)
        response["figure"] = orjson.loads(fig.to_json(engine="orjson"))
    if req.include_exports:
        kinds = resolve_export_kinds(req.export_kinds, defaults={"json"}, allowed={"json", "gif"})
        exports: dict[str, str] = {}
        if "json" in kinds:
            exports["json"] = export_skymap_json(data)
        if "gif" in kinds:
            exports["gif"] = export_skymap_gif_base64(data)
        response["exports"] = exports
    return response
