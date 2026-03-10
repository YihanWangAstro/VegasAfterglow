from __future__ import annotations

import os
import socket
import sys
import tempfile
import threading
import time as time_mod
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeoutError
from contextlib import asynccontextmanager
from functools import lru_cache
from pathlib import Path
from typing import Any
from urllib.error import URLError
from urllib.request import Request, urlopen

import numpy as np
from fastapi import FastAPI, HTTPException, Request as FastAPIRequest, Response
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

from .compute import compute_model, compute_sed, compute_skymap
from .helpers import OBS_FLUX_UNITS, parse_entry, parse_frequency_input, shared_to_physics
from .schemas import MAX_FREQUENCIES, MAX_SNAPSHOTS, LightCurveRequest, SharedParams, SkyMapRequest, SpectrumRequest
from .services.plot_data import build_lc_plot_data, build_sed_plot_data, build_skymap_plot_data

try:
    import VegasAfterglow as _va

    APP_VERSION = _va.__version__.split(".dev")[0]
except Exception:
    APP_VERSION = "unknown"



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
        compute_model(_warmup_params)
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
    allow_credentials=False,
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["Content-Type", "Accept"],
)
app.add_middleware(GZipMiddleware, minimum_size=1024, compresslevel=6)

# ---------------------------------------------------------------------------
# Rate limiting (per-IP, in-memory)
# ---------------------------------------------------------------------------

RATE_LIMIT_WINDOW = 10  # seconds
RATE_LIMIT_MAX = int(os.getenv("RATE_LIMIT_MAX", "300"))  # requests per window
_RATE_EVICT_INTERVAL = 60  # seconds between stale-IP eviction sweeps

_rate_lock = threading.Lock()
_rate_buckets: dict[str, list[float]] = defaultdict(list)
_rate_last_evict = 0.0


def _get_client_ip(request: FastAPIRequest) -> str:
    forwarded = request.headers.get("x-forwarded-for")
    if forwarded:
        return forwarded.split(",")[0].strip()
    return request.client.host if request.client else "unknown"


def _check_rate_limit(request: FastAPIRequest) -> None:
    ip = _get_client_ip(request)
    now = time_mod.monotonic()
    with _rate_lock:
        global _rate_last_evict
        if now - _rate_last_evict > _RATE_EVICT_INTERVAL:
            cutoff = now - RATE_LIMIT_WINDOW
            stale = [k for k, v in _rate_buckets.items() if not v or v[-1] <= cutoff]
            for k in stale:
                del _rate_buckets[k]
            _rate_last_evict = now

        bucket = _rate_buckets[ip]
        cutoff = now - RATE_LIMIT_WINDOW
        _rate_buckets[ip] = bucket = [t for t in bucket if t > cutoff]
        if len(bucket) >= RATE_LIMIT_MAX:
            raise HTTPException(status_code=429, detail="Too many requests. Please wait before trying again.")
        bucket.append(now)


# ---------------------------------------------------------------------------
# Computation timeout
# ---------------------------------------------------------------------------

COMPUTE_TIMEOUT = int(os.getenv("COMPUTE_TIMEOUT", "120"))  # seconds
_compute_pool = ThreadPoolExecutor(max_workers=4)


def _run_with_timeout(fn, params: dict, label: str):
    """Run *fn(params)* in the compute pool with a timeout. Raises HTTPException on failure."""
    future = _compute_pool.submit(fn, params)
    try:
        return future.result(timeout=COMPUTE_TIMEOUT)
    except FuturesTimeoutError:
        raise HTTPException(status_code=504, detail=f"Computation timed out (limit: {COMPUTE_TIMEOUT}s)")
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"{label} failed: {exc}") from exc


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
def options(response: Response) -> dict[str, Any]:
    response.headers["Cache-Control"] = "public, max-age=300"
    return {
        "version": APP_VERSION,
        "obs_flux_units": OBS_FLUX_UNITS,
    }


@app.get("/api/defaults")
def defaults(response: Response) -> dict[str, Any]:
    response.headers["Cache-Control"] = "public, max-age=300"
    return {
        "shared": SharedParams().model_dump(),
        "lightcurve": {
            "frequencies_input": "1e9, R, 1keV",
            "t_min": 1.0,
            "t_max": 1e8,
            "observation_groups": [],
        },
        "spectrum": {
            "t_snapshots_input": "1e3, 1e4, 1e5, 1e6",
            "nu_min": 1e8,
            "nu_max": 1e20,
            "num_nu": 200,
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
def lightcurve(req: LightCurveRequest, request: FastAPIRequest) -> dict[str, Any]:
    _check_rate_limit(request)

    if req.t_min <= 0 or req.t_max <= 0:
        raise HTTPException(status_code=400, detail="t_min and t_max must be positive")
    if req.t_min >= req.t_max:
        raise HTTPException(status_code=400, detail="t_min must be less than t_max")

    frequencies, bands, parse_warnings = parse_frequency_input(req.frequencies_input)
    if len(frequencies) + len(bands) > MAX_FREQUENCIES:
        raise HTTPException(status_code=400, detail=f"Too many frequencies (max {MAX_FREQUENCIES})")

    params = {
        **shared_to_physics(req.shared),
        "frequencies": frequencies,
        "bands": bands,
        "t_min": req.t_min,
        "t_max": req.t_max,
        "num_t": req.shared.num_t,
    }

    t0 = time_mod.perf_counter()
    data = _run_with_timeout(compute_model, params, "Light curve computation")
    compute_s = time_mod.perf_counter() - t0

    response: dict[str, Any] = {
        "meta": {
            "compute_seconds": compute_s,
            "frequencies_hz": frequencies,
            "bands": bands,
            "warnings": parse_warnings,
        },
    }

    plot_data = build_lc_plot_data(data, req.t_min, req.t_max)
    response["plot_data"] = plot_data

    return ORJSONResponse(response)


@app.post("/api/spectrum")
def spectrum(req: SpectrumRequest, request: FastAPIRequest) -> dict[str, Any]:
    _check_rate_limit(request)

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
    if len(t_snapshots) > MAX_SNAPSHOTS:
        raise HTTPException(status_code=400, detail=f"Too many time snapshots (max {MAX_SNAPSHOTS})")

    params = {
        **shared_to_physics(req.shared),
        "t_snapshots": sorted(t_snapshots),
        "nu_min": req.nu_min,
        "nu_max": req.nu_max,
        "num_nu": req.num_nu,
    }
    t0 = time_mod.perf_counter()
    data = _run_with_timeout(compute_sed, params, "Spectrum computation")
    compute_s = time_mod.perf_counter() - t0

    response: dict[str, Any] = {
        "meta": {
            "compute_seconds": compute_s,
            "t_snapshots_s": sorted(t_snapshots),
            "warnings": warnings,
        },
    }

    plot_data = build_sed_plot_data(data)
    response["plot_data"] = plot_data

    return ORJSONResponse(response)


@app.post("/api/skymap")
def skymap(req: SkyMapRequest, request: FastAPIRequest) -> dict[str, Any]:
    _check_rate_limit(request)

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

    t0 = time_mod.perf_counter()
    data = _run_with_timeout(compute_skymap, params, "Sky image computation")
    compute_s = time_mod.perf_counter() - t0

    response: dict[str, Any] = {
        "meta": {
            "compute_seconds": compute_s,
            "n_frames": len(data.get("images", [])),
        },
    }
    response["plot_data"] = build_skymap_plot_data(data)
    return ORJSONResponse(response)
