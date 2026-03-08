"""Raw plot data builders for client-side rendering."""

import base64

import numpy as np
from fastapi import HTTPException

from ..constants import (
    INSTRUMENTS,
)


def _build_instruments_list(selected_instruments, t_lo_s=None, t_hi_s=None):
    """Build the instruments list for plot_data payloads.

    Optionally includes t_lo_s/t_hi_s for lightcurve mode.
    """
    if not selected_instruments:
        return []
    result = []
    for name in selected_instruments:
        nu_min_i, nu_max_i, sens, kind = INSTRUMENTS[name]
        entry = {
            "name": name,
            "nu_min": float(nu_min_i),
            "nu_max": float(nu_max_i),
            "sensitivity": float(sens),
            "kind": kind,
        }
        if t_lo_s is not None:
            entry["t_lo_s"] = float(t_lo_s)
            entry["t_hi_s"] = float(t_hi_s)
        result.append(entry)
    return result


def build_lc_plot_data(data, t_min, t_max, flux_unit, time_unit,
                       selected_instruments=None):
    """Return raw LC data for client-side rendering.

    Returns dict with times_s, pt (per-frequency components in CGS), bands,
    and instruments — all in CGS/seconds so the frontend can convert.
    """
    times = data["times"]
    freqs = data["frequencies"]
    pt_components = data["pt_components"]
    band_data = data["band_data"]

    has_points = len(pt_components) > 0 and freqs.size > 0

    # Build pt block: {name: [[n_t floats] per freq]}
    pt_block = None
    if has_points:
        components_dict = {}
        for comp_name, cf in pt_components:
            components_dict[comp_name] = [cf[i].tolist() for i in range(len(freqs))]
        pt_block = {
            "freq_hz": [float(nu) for nu in freqs],
            "components": components_dict,
        }

    # Build bands block
    bands_list = []
    for nu_min, nu_max, blabel, comps in band_data:
        band_components = {}
        for comp_name, flux_1d in comps:
            band_components[comp_name] = flux_1d.tolist()
        bands_list.append({
            "nu_min": float(nu_min),
            "nu_max": float(nu_max),
            "name": blabel,
            "nu_cen": float(np.sqrt(nu_min * nu_max)),
            "components": band_components,
        })

    instruments_list = _build_instruments_list(selected_instruments, t_lo_s=t_min, t_hi_s=t_max)

    return {
        "flux_unit": flux_unit,
        "time_unit": time_unit,
        "t_min_s": float(t_min),
        "t_max_s": float(t_max),
        "times_s": times.tolist(),
        "pt": pt_block,
        "bands": bands_list,
        "instruments": instruments_list,
    }


def build_sed_plot_data(data, flux_unit, freq_unit, nufnu,
                        selected_instruments=None):
    """Return raw spectrum/SED data for client-side rendering.

    Returns dict with freq_hz, t_snapshots_s, components
    (all in CGS erg/cm²/s/Hz), and instruments.
    """
    t_snapshots = data["t_snapshots"]
    freqs = data["frequencies"]
    components = data["components"]

    # components: list of (name, cf) where cf.shape = (n_nu, n_t)
    # Store as {name: [cf[:, j].tolist() for j in range(n_t)]} — [t_idx][nu_idx]
    components_dict = {}
    for comp_name, cf in components:
        components_dict[comp_name] = [cf[:, j].tolist() for j in range(len(t_snapshots))]

    instruments_list = _build_instruments_list(selected_instruments)

    return {
        "flux_unit": flux_unit,
        "freq_unit": freq_unit,
        "nufnu": nufnu,
        "freq_hz": freqs.tolist(),
        "t_snapshots_s": [float(t) for t in t_snapshots],
        "components": components_dict,
        "instruments": instruments_list,
    }


def build_skymap_plot_data(data: dict) -> dict:
    """Return raw skymap data for client-side rendering.

    Frames are encoded as base64 float32 (row-major, NaN for non-positive pixels).
    Shape is [ny, nx] after the .T transpose applied during log10 conversion.
    """
    images = data["images"]
    extent = data["extent_uas"]
    t_arr = data["t_obs_array"]
    nu_obs = data["nu_obs"]

    if not images:
        raise HTTPException(status_code=500, detail="Sky map returned no images")

    nx, ny = int(images[0].shape[0]), int(images[0].shape[1])
    x_min, x_max, y_min, y_max = float(extent[0]), float(extent[1]), float(extent[2]), float(extent[3])
    dx = (x_max - x_min) / max(nx, 1)
    dy = (y_max - y_min) / max(ny, 1)
    x0 = x_min + 0.5 * dx
    y0 = y_min + 0.5 * dy

    def _log10_frame(img):
        frame = np.asarray(img, dtype=np.float32).T  # shape (ny, nx)
        out = np.full(frame.shape, np.nan, dtype=np.float32)
        np.log10(frame, out=out, where=frame > 0)
        return out

    log_frames = [_log10_frame(img) for img in images]
    stacked = np.stack(log_frames)
    has_finite = bool(np.any(np.isfinite(stacked)))
    z_min = float(np.nanmin(stacked)) if has_finite else 0.0
    z_max = float(np.nanmax(stacked)) if has_finite else 1.0

    frames_b64 = [base64.b64encode(f.tobytes()).decode("ascii") for f in log_frames]

    return {
        "nx": nx,
        "ny": ny,
        "frames_b64f32": frames_b64,
        "extent_uas": [x_min, x_max, y_min, y_max],
        "t_obs_s": [float(t) for t in t_arr],
        "nu_obs_hz": float(nu_obs),
        "dx": dx,
        "dy": dy,
        "x0": x0,
        "y0": y0,
        "z_min": z_min,
        "z_max": z_max,
    }
