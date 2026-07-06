"""Absolute accuracy of the shipped default resolutions, across all grid axes.

The validation suite pins its own fiducial resolution and the golden baselines
bless whatever state they are regenerated from, so neither gates the accuracy
of the *defaults* users actually get. This test does: it compares the default
grid against a reference converged in phi, theta, AND time simultaneously, on
the families that pin the calibrated defaults:

- tophat/wind reverse shock viewed at theta_obs = theta_c: the wind thin-shell
  family whose rvs component pins the time default, with the jet edge on the
  line of sight (the phi-worst geometry) and the bright-limb theta gradient.
- powerlaw + reverse shock on-axis: the theta-worst family (wing rows drive
  both shocks' theta convergence).

Calibrated 2026-07-05 at the (0.06, 0.15, 6)/(0.06, 0.2, 10) per-mode
defaults; the asserted gates and masking are the validation suite's own
(imported from validation.common).
"""
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from validation.common import (  # noqa: E402
    MAX_ERROR_THRESHOLD,
    MEAN_ERROR_THRESHOLD,
    masked_relative_errors,
)


def _load_regenerate_module():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden", "regenerate.py")
    spec = importlib.util.spec_from_file_location("golden_regenerate", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


regen = _load_regenerate_module()

pytestmark = pytest.mark.physics

T = np.logspace(1, 7, 100)
LOG_T = np.log10(T)
BANDS = np.array([1e9, 4.84e14, 1e18])
CONVERGED = (0.3, 0.8, 40)

RAD = {"eps_e": 0.1, "eps_B": 0.01, "p": 2.3, "xi_e": 1.0}

# The forward-shock-only and reverse-shock modes have separate default
# resolutions, so BOTH modes need their pinning families gated here. Configs
# use the golden-baseline schema and builder (regenerate.build_model).
FAMILIES = {
    # --- reverse-shock mode: pins (0.06, 0.2, 10) ---
    "rs_tophat_wind_edge": {
        "jet": {"type": "TophatJet", "theta_c": 0.1, "E_iso": 1e52, "Gamma0": 300, "duration": 1.0},
        "medium": {"type": "Wind", "A_star": 0.1},
        "observer": {"lumi_dist": 1e28, "z": 1.0, "theta_obs": 0.1},
        "fwd_rad": RAD,
        "rvs_rad": RAD,
    },
    "rs_powerlaw_onaxis": {
        "jet": {"type": "PowerLawJet", "theta_c": 0.1, "E_iso": 1e52, "Gamma0": 300,
                "k_e": 2.0, "k_g": 2.0, "duration": 1.0},
        "medium": {"type": "ISM", "n_ism": 1.0},
        "observer": {"lumi_dist": 1e28, "z": 1.0, "theta_obs": 0.0},
        "fwd_rad": RAD,
        "rvs_rad": RAD,
    },
    # --- forward-shock-only mode: pins (0.06, 0.15, 6) ---
    "fs_two_component": {
        "jet": {"type": "TwoComponentJet", "theta_c": 0.05, "E_iso": 1e52, "Gamma0": 300,
                "theta_w": 0.3, "E_iso_w": 1e50, "Gamma0_w": 50},
        "medium": {"type": "ISM", "n_ism": 1.0},
        "observer": {"lumi_dist": 1e28, "z": 1.0, "theta_obs": 0.075},
        "fwd_rad": RAD,
    },
    "fs_tophat_wind_edge": {
        "jet": {"type": "TophatJet", "theta_c": 0.1, "E_iso": 1e52, "Gamma0": 300},
        "medium": {"type": "Wind", "A_star": 0.1},
        "observer": {"lumi_dist": 1e28, "z": 1.0, "theta_obs": 0.1},
        "fwd_rad": RAD,
    },
}


@pytest.mark.parametrize("family", sorted(FAMILIES))
def test_default_resolution_holds_convergence_gates(family):
    """Every emission component at the (mode-specific) default resolution stays within both validation gates (mean error < 5%, max error < 15%) against a phi+theta+time-converged reference, on the families that pin the calibrated defaults of each mode."""
    cfg = FAMILIES[family]
    ref = regen.build_model({**cfg, "resolutions": CONVERGED}).flux_density_grid(T, BANDS)
    got = regen.build_model(cfg).flux_density_grid(T, BANDS)

    # Every emission component is gated individually: total flux is dominated
    # by the forward shock and hides multi-percent reverse-shock errors (this
    # exact blind spot let an under-resolved reverse-shock default ship in
    # 2026-07).
    modes = ("fwd", "rvs") if "rvs_rad" in cfg else ("fwd",)
    for mode in modes:
        got_c = np.asarray(getattr(got, mode).sync)
        ref_c = np.asarray(getattr(ref, mode).sync)
        mean_err, max_err = 0.0, 0.0
        for i in range(len(BANDS)):
            rel = masked_relative_errors(got_c[i], ref_c[i], LOG_T)
            if rel is None:
                continue
            mean_err = max(mean_err, float(rel.mean()))
            max_err = max(max_err, float(rel.max()))
        assert mean_err < MEAN_ERROR_THRESHOLD, (
            f"{family}/{mode}.sync: default-resolution mean error {mean_err:.2%} "
            f">= {MEAN_ERROR_THRESHOLD:.0%} gate")
        assert max_err < MAX_ERROR_THRESHOLD, (
            f"{family}/{mode}.sync: default-resolution max error {max_err:.2%} "
            f">= {MAX_ERROR_THRESHOLD:.0%} gate")
