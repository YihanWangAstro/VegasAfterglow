"""CSV and JSON export functions."""

import base64
import io

import orjson

import numpy as np
from VegasAfterglow.cli import _format_energy
from matplotlib import cm, colors
from PIL import Image

from .constants import FLUX_SCALES, FREQ_SCALES, TIME_SCALES
from .helpers import band_label, format_time_label

_ORJSON_DUMP_OPTIONS = orjson.OPT_INDENT_2 | getattr(orjson, "OPT_SERIALIZE_NUMPY", 0)


def _dump_json(obj):
    return orjson.dumps(obj, option=_ORJSON_DUMP_OPTIONS).decode()


def export_csv(data, flux_unit, time_unit):
    """Generate CSV string from results."""
    times = data["times"]
    freqs = data["frequencies"]
    pt_components = data["pt_components"]
    band_data = data["band_data"]
    f_scale = FLUX_SCALES[flux_unit]
    t_scale = TIME_SCALES[time_unit]

    buf = io.StringIO()
    cols = [f"t({time_unit})"]
    for nu in freqs:
        nu_label = _format_energy(nu)
        for comp_name, _ in pt_components:
            cols.append(f"F_{comp_name}({nu_label})[{flux_unit}]")
    for nu_min, nu_max, blabel, comps in band_data:
        bl = band_label(nu_min, nu_max, blabel)
        for comp_name, _ in comps:
            cols.append(f"F_{comp_name}({bl})[erg/cm2/s]")
    buf.write(",".join(cols) + "\n")
    for j in range(len(times)):
        row = [f"{times[j] / t_scale:.6e}"]
        for i in range(len(freqs)):
            for _, comp_flux in pt_components:
                row.append(f"{comp_flux[i, j] / f_scale:.6e}")
        for _, _, _, comps in band_data:
            for _, comp_flux in comps:
                row.append(f"{comp_flux[j]:.6e}")
        buf.write(",".join(row) + "\n")
    return buf.getvalue()


def export_json(data, flux_unit, time_unit):
    """Generate JSON string from results."""
    times = data["times"]
    freqs = data["frequencies"]
    pt_components = data["pt_components"]
    band_data = data["band_data"]
    f_scale = FLUX_SCALES[flux_unit]
    t_scale = TIME_SCALES[time_unit]

    obj = {
        "units": {"time": time_unit, "flux_density": flux_unit},
        "times": times / t_scale,
        "frequencies_Hz": freqs,
        "flux_density": {},
    }
    for nu_idx, nu in enumerate(freqs):
        nu_label = _format_energy(nu)
        obj["flux_density"][nu_label] = {
            comp_name: comp_flux[nu_idx] / f_scale
            for comp_name, comp_flux in pt_components
        }
    if band_data:
        obj["units"]["band_flux"] = "erg/cm2/s"
        obj["bands"] = {}
        for nu_min, nu_max, blabel, comps in band_data:
            bl = band_label(nu_min, nu_max, blabel)
            obj["bands"][bl] = {
                "nu_min_Hz": nu_min,
                "nu_max_Hz": nu_max,
                "flux": {comp_name: comp_flux for comp_name, comp_flux in comps},
            }
    return _dump_json(obj)


def export_sed_csv(data, flux_unit, freq_unit):
    """Generate CSV string from SED results."""
    t_snapshots = data["t_snapshots"]
    freqs = data["frequencies"]
    components = data["components"]
    f_scale = FLUX_SCALES[flux_unit]
    nu_scale = FREQ_SCALES[freq_unit]

    buf = io.StringIO()
    cols = [f"nu({freq_unit})"]
    for t in t_snapshots:
        t_label = format_time_label(t)
        for comp_name, _ in components:
            cols.append(f"F_{comp_name}(t={t_label})[{flux_unit}]")
    buf.write(",".join(cols) + "\n")
    for i in range(len(freqs)):
        row = [f"{freqs[i] / nu_scale:.6e}"]
        for j in range(len(t_snapshots)):
            for _, comp_flux in components:
                row.append(f"{comp_flux[i, j] / f_scale:.6e}")
        buf.write(",".join(row) + "\n")
    return buf.getvalue()


def export_sed_json(data, flux_unit, freq_unit):
    """Generate JSON string from SED results."""
    t_snapshots = data["t_snapshots"]
    freqs = data["frequencies"]
    components = data["components"]
    f_scale = FLUX_SCALES[flux_unit]
    nu_scale = FREQ_SCALES[freq_unit]

    obj = {
        "units": {"frequency": freq_unit, "flux_density": flux_unit},
        "frequencies": freqs / nu_scale,
        "frequencies_Hz": freqs,
        "t_snapshots_s": t_snapshots,
        "flux_density": {},
    }
    for j, t in enumerate(t_snapshots):
        t_label = format_time_label(t)
        obj["flux_density"][t_label] = {
            comp_name: comp_flux[:, j] / f_scale
            for comp_name, comp_flux in components
        }
    return _dump_json(obj)


def export_skymap_json(data):
    """Generate JSON string from sky map results."""
    obj = {
        "t_obs_s": data["t_obs_array"],
        "nu_obs_Hz": data["nu_obs"],
        "fov_uas": data["fov_uas"],
        "extent_uas": data["extent_uas"],
        "units": {"image": "erg/cm2/s/Hz/sr", "extent": "uas"},
        "images": data["images"],
    }
    return _dump_json(obj)


def export_skymap_gif_base64(data):
    """Generate GIF (base64-encoded) from sky map frames."""
    images = data.get("images", [])
    if not images:
        raise ValueError("Sky map returned no frames")

    all_pos = np.concatenate([img[img > 0] for img in images if np.any(img > 0)])
    vmin = float(np.log10(all_pos.min())) if all_pos.size > 0 else 0.0
    vmax = float(np.log10(all_pos.max())) if all_pos.size > 0 else 1.0

    cmap = cm.get_cmap("inferno").copy()
    cmap.set_bad(color="white")
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    pil_frames = []
    for img in images:
        log_img = np.where(img > 0, np.log10(img), np.nan)
        rgba = (cmap(norm(log_img.T)) * 255).astype(np.uint8)
        rgba = rgba[::-1]
        pil_frames.append(Image.fromarray(rgba, mode="RGBA"))

    frame_ms = max(50, 2000 // max(1, len(pil_frames)))
    buf = io.BytesIO()
    pil_frames[0].save(
        buf,
        format="GIF",
        save_all=True,
        append_images=pil_frames[1:],
        duration=frame_ms,
        loop=0,
        disposal=2,
    )
    return base64.b64encode(buf.getvalue()).decode("ascii")
