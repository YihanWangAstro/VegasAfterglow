"""Utility functions for regression tests."""

import numpy as np

STANDARD_PARAMS = {
    "E_iso": 1e52, "Gamma0": 300, "theta_c": 0.1, "n_ism": 1.0, "A_star": 0.1,
    "eps_e": 0.1, "eps_B": 0.01, "p": 2.2, "xi_e": 1.0, "lumi_dist": 1e28, "z": 1.0,
}


def fit_powerlaw(t, f):
    """Compute average power-law slope from local derivatives in log-log space."""
    valid = (f > 0) & np.isfinite(f) & (t > 0) & np.isfinite(t)
    if np.sum(valid) < 3:
        return np.nan
    t_valid, f_valid = t[valid], f[valid]
    sort_idx = np.argsort(t_valid)
    log_t, log_f = np.log10(t_valid[sort_idx]), np.log10(f_valid[sort_idx])
    d_log_t, d_log_f = np.diff(log_t), np.diff(log_f)
    valid_diff = np.abs(d_log_t) > 1e-10
    if np.sum(valid_diff) < 2:
        return np.nan
    local_slopes = d_log_f[valid_diff] / d_log_t[valid_diff]
    return np.mean(local_slopes[1:]) if len(local_slopes) > 1 else np.mean(local_slopes)
