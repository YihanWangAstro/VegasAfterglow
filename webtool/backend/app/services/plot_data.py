"""Raw plot data builders for client-side rendering."""

import base64

import numpy as np
from fastapi import HTTPException


def build_lc_plot_data(data, t_min, t_max):
    """Return raw LC data for client-side rendering.

    Returns dict with times_s, pt (per-frequency components in CGS), and bands
    — all in CGS/seconds so the frontend can convert.
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

    return {
        "t_min_s": float(t_min),
        "t_max_s": float(t_max),
        "times_s": times.tolist(),
        "pt": pt_block,
        "bands": bands_list,
    }


def build_sed_plot_data(data):
    """Return raw spectrum/SED data for client-side rendering.

    Returns dict with freq_hz, t_snapshots_s, components
    — all in CGS erg/cm²/s/Hz.
    """
    t_snapshots = data["t_snapshots"]
    freqs = data["frequencies"]
    components = data["components"]

    # components: list of (name, cf) where cf.shape = (n_nu, n_t)
    # Store as {name: [cf[:, j].tolist() for j in range(n_t)]} — [t_idx][nu_idx]
    components_dict = {}
    for comp_name, cf in components:
        components_dict[comp_name] = [cf[:, j].tolist() for j in range(len(t_snapshots))]

    return {
        "freq_hz": freqs.tolist(),
        "t_snapshots_s": [float(t) for t in t_snapshots],
        "components": components_dict,
    }


def build_skymap_plot_data(data: dict) -> dict:
    """Return raw skymap data for client-side rendering.

    Zero-border is trimmed before log10 to reduce computation and transfer size.
    Frames are encoded as base64 float32 (row-major, NaN for non-positive pixels).
    x0/y0/nx/ny are adjusted to reflect the cropped region.
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

    # Stack transposed frames into single (n_frames, ny, nx) array
    stacked = np.stack([np.asarray(img, dtype=np.float32).T for img in images])

    # Crop to bounding box of positive pixels (union across all frames).
    # np.max avoids materializing a full boolean temporary like `stacked > 0`.
    max_frame = np.max(stacked, axis=0)  # (ny, nx)
    rows_any = np.any(max_frame > 0, axis=1)
    cols_any = np.any(max_frame > 0, axis=0)
    if rows_any.any() and cols_any.any():
        row_idx = np.where(rows_any)[0]
        col_idx = np.where(cols_any)[0]
        r0, r1 = int(row_idx[0]), int(row_idx[-1]) + 1
        c0, c1 = int(col_idx[0]), int(col_idx[-1]) + 1
        stacked = stacked[:, r0:r1, c0:c1]
        x0 += c0 * dx
        y0 += r0 * dy
        nx = c1 - c0
        ny = r1 - r0

    # log10 in-place on the cropped array, NaN for non-positive pixels
    stacked[stacked <= 0] = np.nan
    np.log10(stacked, out=stacked, where=np.isfinite(stacked))

    has_finite = bool(np.any(np.isfinite(stacked)))
    z_min = float(np.nanmin(stacked)) if has_finite else 0.0
    z_max = float(np.nanmax(stacked)) if has_finite else 1.0

    frames_b64 = [base64.b64encode(stacked[i].tobytes()).decode("ascii") for i in range(stacked.shape[0])]

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
