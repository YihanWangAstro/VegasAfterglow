"""CSV and JSON export functions."""

import io
import json

from VegasAfterglow.cli import _format_energy

from .constants import FLUX_SCALES, FREQ_SCALES, TIME_SCALES
from .helpers import band_label, format_time_label


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
        "times": (times / t_scale).tolist(),
        "frequencies_Hz": freqs.tolist(),
        "flux_density": {},
    }
    for nu_idx, nu in enumerate(freqs):
        nu_label = _format_energy(nu)
        obj["flux_density"][nu_label] = {
            comp_name: (comp_flux[nu_idx] / f_scale).tolist()
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
                "flux": {
                    comp_name: comp_flux.tolist() for comp_name, comp_flux in comps
                },
            }
    return json.dumps(obj, indent=2)


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
        "frequencies": (freqs / nu_scale).tolist(),
        "frequencies_Hz": freqs.tolist(),
        "t_snapshots_s": t_snapshots.tolist(),
        "flux_density": {},
    }
    for j, t in enumerate(t_snapshots):
        t_label = format_time_label(t)
        obj["flux_density"][t_label] = {
            comp_name: (comp_flux[:, j] / f_scale).tolist()
            for comp_name, comp_flux in components
        }
    return json.dumps(obj, indent=2)


def export_skymap_json(data):
    """Generate JSON string from sky map results."""
    obj = {
        "t_obs_s": data["t_obs_array"].tolist(),
        "nu_obs_Hz": data["nu_obs"],
        "fov_uas": data["fov_uas"],
        "extent_uas": data["extent_uas"].tolist(),
        "units": {"image": "erg/cm2/s/Hz/sr", "extent": "uas"},
        "images": [img.tolist() for img in data["images"]],
    }
    return json.dumps(obj, indent=2)
