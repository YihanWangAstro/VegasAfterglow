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
}


def build_model(config):
    jet_kwargs = dict(config["jet"])
    jet = getattr(va, jet_kwargs.pop("type"))(**jet_kwargs)
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
