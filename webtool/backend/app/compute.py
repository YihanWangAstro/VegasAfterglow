"""Model building and computation functions."""

import time as time_mod

import numpy as np

from VegasAfterglow import (
    GaussianJet,
    ISM,
    Model,
    Observer,
    PowerLawJet,
    Radiation,
    TophatJet,
    TwoComponentJet,
    Wind,
)
from VegasAfterglow.cli import _get_components
from VegasAfterglow.units import uas


def build_jet(p):
    common = dict(
        theta_c=p["theta_c"],
        E_iso=p["E_iso"],
        Gamma0=p["Gamma0"],
        spreading=p["spreading"],
        duration=p["duration"],
    )
    jt = p["jet_type"]
    if jt == "Top-hat":
        return TophatJet(**common)
    if jt == "Gaussian":
        return GaussianJet(**common)
    if jt == "Power-law":
        return PowerLawJet(**common, k_e=p["k_e"], k_g=p["k_g"])
    if jt == "Two-component":
        return TwoComponentJet(
            **common,
            theta_w=p["theta_w"],
            E_iso_w=p["E_iso_w"],
            Gamma0_w=p["Gamma0_w"],
        )


def build_medium(p):
    if p["medium_type"] == "ISM":
        return ISM(n_ism=p["n_ism"])
    if p["medium_type"] == "Wind bubble":
        return Wind(A_star=p["A_star"], n_ism=p["n_ism"], k_m=p["k_m"])
    return Wind(A_star=p["A_star"], k_m=p["k_m"])


def build_observer(p):
    return Observer(lumi_dist=p["d_L_cm"], z=p["z"], theta_obs=p["theta_obs"])


def build_radiation(p):
    fwd = Radiation(
        eps_e=p["eps_e"],
        eps_B=p["eps_B"],
        p=p["p"],
        xi_e=p["xi_e"],
        ssc=p["ssc"],
        kn=p["kn"],
    )
    rvs = None
    if p["enable_rvs"]:
        rvs = Radiation(
            eps_e=p["eps_e_r"],
            eps_B=p["eps_B_r"],
            p=p["p_r"],
            xi_e=p["xi_e_r"],
            ssc=p.get("rvs_ssc", False),
            kn=p.get("rvs_kn", False),
        )
    return fwd, rvs


def _build_model(p):
    """Construct a Model from a params dict."""
    jet = build_jet(p)
    medium = build_medium(p)
    observer = build_observer(p)
    fwd_rad, rvs_rad = build_radiation(p)
    return Model(
        jet=jet,
        medium=medium,
        observer=observer,
        fwd_rad=fwd_rad,
        rvs_rad=rvs_rad,
        resolutions=(p["res_phi"], p["res_theta"], p["res_t"]),
    )


def compute_model(p):
    """Build model and compute light curve."""
    model = _build_model(p)
    bands = p.get("bands", [])

    times = np.logspace(np.log10(p["t_min"]), np.log10(p["t_max"]), p["num_t"])

    t0 = time_mod.time()

    pt_components = []
    freqs = np.array(sorted(p["frequencies"])) if p["frequencies"] else np.array([])
    if freqs.size > 0:
        result = model.flux_density_grid(times, freqs)
        pt_components = [(name, np.array(flux)) for name, flux in _get_components(result)]

    band_data = []
    for nu_min, nu_max, label in bands:
        num_nu_band = max(5, int(np.log10(nu_max / nu_min) * 5))
        br = model.flux(times, nu_min, nu_max, num_nu_band)
        comps = [(name, np.array(flux)) for name, flux in _get_components(br)]
        band_data.append((nu_min, nu_max, label, comps))

    elapsed = time_mod.time() - t0

    return {
        "times": times,
        "frequencies": freqs,
        "pt_components": pt_components,
        "band_data": band_data,
        "elapsed": elapsed,
    }


def compute_sed(p):
    """Build model and compute SED at snapshot times."""
    t_snapshots = np.array(sorted(p["t_snapshots"]))
    nu_array = np.logspace(np.log10(p["nu_min"]), np.log10(p["nu_max"]), p["num_nu"])

    model = _build_model(p)

    t0 = time_mod.time()
    result = model.flux_density_grid(t_snapshots, nu_array)
    components = [(name, np.array(flux)) for name, flux in _get_components(result)]
    elapsed = time_mod.time() - t0

    return {
        "t_snapshots": t_snapshots,
        "frequencies": nu_array,
        "components": components,
        "elapsed": elapsed,
    }


def compute_skymap(p):
    """Build model and compute a single sky image."""
    model = _build_model(p)
    t_obs = p["t_obs"]
    t0 = time_mod.time()
    result = model.sky_image(
        [t_obs],
        nu_obs=p["nu_obs"],
        fov=p["fov"] * uas,
        npixel=p["npixel"],
    )
    elapsed = time_mod.time() - t0
    return {
        "image": np.array(result.image[0]),
        "extent_uas": np.array(result.extent) / uas,
        "t_obs": t_obs,
        "nu_obs": p["nu_obs"],
        "fov_uas": p["fov"],
        "elapsed": elapsed,
    }
