//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "smooth-power-law-syn.h"

#include "../util/utilities.h"
//========================================================================================================
//                                  SmoothPowerLawSyn Class Methods
//========================================================================================================
inline Real log2_broken_power_ratio(Real log2_nu, Real log2_nu_b, Real beta_lo, Real beta_hi, Real s) {
    return -log2_softplus(s * (beta_lo - beta_hi) * (log2_nu - log2_nu_b)) / s;
}

inline Real log2_broken_power(Real log2_nu, Real log2_nu_b, Real beta_lo, Real beta_hi, Real s) {
    return beta_lo * (log2_nu - log2_nu_b) + log2_broken_power_ratio(log2_nu, log2_nu_b, beta_lo, beta_hi, s);
}

inline Real log2_smooth_one(Real log2_a, Real log2_b) {
    return log2_a - log2_softplus(log2_a - log2_b);
}

Real SmoothPowerLawSyn::log2_optical_thin(Real log2_nu) const {
    if (log2_nu_m < log2_nu_c) {
        // Use precomputed p-dependent slopes
        return log2_broken_power(log2_nu, log2_nu_m, 1. / 3, half_one_minus_p_, smooth_m_slow_) +
               log2_broken_power_ratio(log2_nu, log2_nu_c, half_one_minus_p_, minus_half_p_, smooth_c_slow_);
    } else {
        return log2_broken_power(log2_nu, log2_nu_c, 1. / 3, -0.5, smooth_m_fast_) +
               log2_broken_power_ratio(log2_nu, log2_nu_m, -0.5, minus_half_p_, smooth_c_fast_);
    }
}

Real SmoothPowerLawSyn::log2_optical_thick(Real log2_nu) const {
    const Real log2_x = log2_nu - log2_nu_m;
    const Real s = -smooth_a_ * std::exp2(2. / 3 * log2_x);

    return 2.5 * log2_x + log2_softplus(-0.5 * log2_x + s);
}

Real SmoothPowerLawSyn::log2_optical_thick_sharp(Real log2_nu) const {
    if (log2_nu < log2_nu_m) {
        return 2. * (log2_nu - log2_nu_m);
    } else {
        return 2.5 * (log2_nu - log2_nu_m);
    }
}

Real SmoothPowerLawSyn::log2_optical_thin_sharp(Real log2_nu) const {
    if (log2_nu_m < log2_nu_c) {
        if (log2_nu < log2_nu_m) {
            return (log2_nu - log2_nu_m) / 3.0;
        } else if (log2_nu < log2_nu_c) {
            return half_one_minus_p_ * (log2_nu - log2_nu_m);
        } else {
            return half_one_minus_p_ * (log2_nu_c - log2_nu_m) + minus_half_p_ * (log2_nu - log2_nu_c);
        }
    } else {
        if (log2_nu < log2_nu_c) {
            return (log2_nu - log2_nu_c) / 3.0;
        } else if (log2_nu < log2_nu_m) {
            return -0.5 * (log2_nu - log2_nu_c);
        } else {
            return -0.5 * (log2_nu_m - log2_nu_c) + minus_half_p_ * (log2_nu - log2_nu_m);
        }
    }
}

Real SmoothPowerLawSyn::compute_spectrum(Real nu) const {
    return fast_exp2(compute_log2_spectrum(fast_log2(nu)));
}

Real SmoothPowerLawSyn::compute_log2_spectrum(Real log2_nu) const {
    Real log2_f_thin = log2_optical_thin(log2_nu);
    Real log2_f_thick = log2_optical_thick(log2_nu);
    return log2_norm_ + log2_smooth_one(log2_f_thin, log2_f_thick + log2_thick_norm_);
}

void SmoothPowerLawSyn::update_constant() {
    // Precompute smoothing parameters
    smooth_m_slow_ = 1.8 - 0.4 * p;
    smooth_c_slow_ = 1 - 0.04 * p;

    smooth_m_fast_ = 3.5 - 0.85 * p;
    smooth_c_fast_ = 0.6;

    smooth_a_ = 3.5 * p - 1.5;

    // Precompute power-law slopes (used repeatedly in spectrum calculations)
    half_one_minus_p_ = 0.5 * (1.0 - p);
    minus_half_p_ = -0.5 * p;

    // Precompute inverse for division optimization
    inv_nu_M_ = 1.0 / nu_M;
    ln2_div_nu_M_ = 1.442695040888963407359924681001892137 * inv_nu_M_;

    // Compute normalization
    if (nu_m < nu_c) {
        log2_norm_ = 1.0 / smooth_m_slow_;
    } else {
        log2_norm_ = 1.0 / smooth_c_fast_;
    }
    log2_thick_norm_ = log2_optical_thin_sharp(log2_nu_a) - log2_optical_thick_sharp(log2_nu_a);
}

Real SmoothPowerLawSyn::compute_I_nu(Real nu) const {
    if (nu <= nu_c) { // Below cooling frequency, simple scaling
        return fast_exp(-nu / nu_M) * I_nu_max * compute_spectrum(nu);
    } else {
        return fast_exp(-nu / nu_M) * I_nu_max * compute_spectrum(nu) * inverse_compton_correction(*this, nu);
    }
}

Real SmoothPowerLawSyn::compute_log2_I_nu(Real log2_nu) const {
    if (log2_nu <= log2_nu_c) { // Below cooling frequency, simple scaling
        return log2_I_nu_max + compute_log2_spectrum(log2_nu) - ln2_div_nu_M_ * fast_exp2(log2_nu);
    } else {
        const Real nu = fast_exp2(log2_nu);
        return log2_I_nu_max + compute_log2_spectrum(log2_nu) - ln2_div_nu_M_ * nu +
               fast_log2(inverse_compton_correction(*this, nu));
    }
}
