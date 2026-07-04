"""Regenerate the golden-baseline .npz files used by test_golden.py.

Run from the repo root:

    python tests/python/golden/regenerate.py

WARNING: only run this when a physics change is INTENDED. The baselines
define the expected output of the code; regenerating them silently blesses
whatever the current build produces. Inspect the diff against the old
baselines first, then commit the new .npz files together with the physics
change that motivated them.
"""

import json
import os

import numpy as np

import VegasAfterglow as va

GOLDEN_DIR = os.path.dirname(os.path.abspath(__file__))

T = np.logspace(2, 8, 40)
NUS = np.array([1e9, 1e14, 1e17, 1e22])

CONFIGS = {
    "tophat_ism": {
        "jet": {"type": "TophatJet", "theta_c": 0.1, "E_iso": 1e53, "Gamma0": 300},
        "medium": {"type": "ISM", "n_ism": 1.0},
        "observer": {"lumi_dist": 3e28, "z": 0.5, "theta_obs": 0.0},
        "fwd_rad": {"eps_e": 0.1, "eps_B": 1e-3, "p": 2.3},
    },
    "rs_thick": {
        "jet": {"type": "TophatJet", "theta_c": 0.1, "E_iso": 1e53, "Gamma0": 100, "duration": 1000},
        "medium": {"type": "ISM", "n_ism": 1.0},
        "observer": {"lumi_dist": 3e28, "z": 0.5, "theta_obs": 0.0},
        "fwd_rad": {"eps_e": 0.1, "eps_B": 1e-3, "p": 2.3},
        "rvs_rad": {"eps_e": 0.1, "eps_B": 1e-2, "p": 2.5},
    },
    "gauss_wind_ssc": {
        "jet": {"type": "GaussianJet", "theta_c": 0.1, "E_iso": 1e53, "Gamma0": 300},
        "medium": {"type": "Wind", "A_star": 0.1},
        "observer": {"lumi_dist": 3e28, "z": 1.0, "theta_obs": 0.2},
        "fwd_rad": {"eps_e": 0.1, "eps_B": 1e-4, "p": 2.3, "ssc": True, "kn": True},
    },
    # Structured-jet reverse shock, viewed off-axis at a fine theta grid: the jet
    # wings put rows at Gamma0(theta) ~ 20-60, where the region-3 seed state is
    # most fragile at early times (see the physical-domain projection in
    # FRShockEqn::operator()). No other golden reaches this regime, and the full
    # validation benchmark that does only runs on the weekly deploy workflow.
    "gauss_ism_rs": {
        "jet": {"type": "GaussianJet", "theta_c": 0.1, "E_iso": 1e52, "Gamma0": 300, "duration": 1.0},
        "medium": {"type": "ISM", "n_ism": 1.0},
        "observer": {"lumi_dist": 1e28, "z": 1.0, "theta_obs": 0.4},
        "fwd_rad": {"eps_e": 0.1, "eps_B": 0.01, "p": 2.3},
        "rvs_rad": {"eps_e": 0.1, "eps_B": 0.01, "p": 2.3},
        "resolutions": (0.1, 1.2, 10),
    },
    # Magnetized ejecta reverse shock: the only golden exercising the sigma > 0
    # jump conditions (cubic root of the downstream four-velocity), the shell
    # magnetization evolution, and the magnetosonic crossing-rate cap.
    "tophat_sigma_rs": {
        "jet": {"type": "MagnetizedTophat", "theta_c": 0.1, "E_iso": 1e52, "Gamma0": 300,
                "sigma0": 0.1, "duration": 1.0},
        "medium": {"type": "ISM", "n_ism": 1.0},
        "observer": {"lumi_dist": 3e28, "z": 0.5, "theta_obs": 0.0},
        "fwd_rad": {"eps_e": 0.1, "eps_B": 1e-3, "p": 2.3},
        "rvs_rad": {"eps_e": 0.1, "eps_B": 1e-2, "p": 2.5},
    },
    # Power-law jet in a wind: covers the power-law structure profile and the
    # only reverse shock evolved in a wind (k = 2) density gradient.
    "powerlaw_wind_rs": {
        "jet": {"type": "PowerLawJet", "theta_c": 0.1, "E_iso": 1e52, "Gamma0": 300,
                "k_e": 2.0, "k_g": 2.0, "duration": 1.0},
        "medium": {"type": "Wind", "A_star": 0.1},
        "observer": {"lumi_dist": 1e28, "z": 1.0, "theta_obs": 0.3},
        "fwd_rad": {"eps_e": 0.1, "eps_B": 0.01, "p": 2.3},
        "rvs_rad": {"eps_e": 0.1, "eps_B": 0.01, "p": 2.3},
    },
    # Two-component jet viewed between the core and the wing: covers the
    # two-component structure profile (narrow fast core + wide slow sheath).
    "two_component_ism": {
        "jet": {"type": "TwoComponentJet", "theta_c": 0.05, "E_iso": 1e52, "Gamma0": 300,
                "theta_w": 0.3, "E_iso_w": 1e50, "Gamma0_w": 50},
        "medium": {"type": "ISM", "n_ism": 1.0},
        "observer": {"lumi_dist": 1e28, "z": 1.0, "theta_obs": 0.15},
        "fwd_rad": {"eps_e": 0.1, "eps_B": 0.01, "p": 2.3},
    },
}


def _magnetized_tophat(theta_c, E_iso, Gamma0, sigma0, duration=1.0):
    # Named jet factories take no magnetization; build the tophat profile on the
    # generic Ejecta with a constant sigma0 (scalars stay JSON-serializable).
    return va.Ejecta(
        E_iso=lambda phi, theta: E_iso if theta <= theta_c else 0.0,
        Gamma0=lambda phi, theta: Gamma0 if theta <= theta_c else 1.0,
        sigma0=lambda phi, theta: sigma0,
        duration=duration,
    )


def build_model(config):
    jet_kwargs = dict(config["jet"])
    jet_type = jet_kwargs.pop("type")
    if jet_type == "MagnetizedTophat":
        jet = _magnetized_tophat(**jet_kwargs)
    else:
        jet = getattr(va, jet_type)(**jet_kwargs)
    medium_kwargs = dict(config["medium"])
    medium = getattr(va, medium_kwargs.pop("type"))(**medium_kwargs)
    observer = va.Observer(**config["observer"])
    fwd_rad = va.Radiation(**config["fwd_rad"])
    rvs_rad = va.Radiation(**config["rvs_rad"]) if "rvs_rad" in config else None
    kwargs = {}
    if "resolutions" in config:
        kwargs["resolutions"] = tuple(config["resolutions"])
    return va.Model(jet=jet, medium=medium, observer=observer, fwd_rad=fwd_rad, rvs_rad=rvs_rad, **kwargs)


def compute_components(model):
    flux = model.flux_density_grid(T, NUS)
    return {
        "total": np.asarray(flux.total),
        "fwd_sync": np.asarray(flux.fwd.sync),
        "fwd_ssc": np.asarray(flux.fwd.ssc),
        "rvs_sync": np.asarray(flux.rvs.sync),
        "rvs_ssc": np.asarray(flux.rvs.ssc),
    }


def main():
    version = getattr(va, "__version__", "unknown")
    print(f"VegasAfterglow version: {version}")
    print("WARNING: regenerating golden baselines. Only do this when a physics")
    print("change is INTENDED; the new .npz files must be committed with it.")
    for name, config in CONFIGS.items():
        model = build_model(config)
        components = compute_components(model)
        path = os.path.join(GOLDEN_DIR, f"{name}.npz")
        np.savez(path, t=T, nus=NUS, config=json.dumps(config), **components)
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
