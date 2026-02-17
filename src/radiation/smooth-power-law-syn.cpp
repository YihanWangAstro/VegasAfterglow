//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "smooth-power-law-syn.h"

//========================================================================================================
//                                  SmoothPowerLawSyn Class Methods
//========================================================================================================
inline Real log2_smooth_one(Real log2_a, Real log2_b) noexcept {
    return log2_a - log2_softplus(log2_a - log2_b);
}

Real SmoothPowerLawSyn::log2_optical_thin(Real log2_nu) const noexcept {
    if (log2_nu_m < log2_nu_c) {
        return (log2_nu - log2_nu_m) / 3.0 +
               log2_broken_power_ratio(log2_nu, log2_nu_m, diff_slope_m_slow_, smooth_m_slow_) +
               log2_broken_power_ratio(log2_nu, log2_nu_c, diff_slope_c_slow_, smooth_c_slow_);
    } else {
        return (log2_nu - log2_nu_c) / 3.0 +
               log2_broken_power_ratio(log2_nu, log2_nu_c, diff_slope_c_fast_, smooth_c_fast_) +
               log2_broken_power_ratio(log2_nu, log2_nu_m, diff_slope_m_fast_, smooth_m_fast_);
    }
}

Real SmoothPowerLawSyn::log2_optical_thick(Real log2_nu) const noexcept {
    const Real log2_x = log2_nu - log2_nu_m;
    const Real s = -smooth_a_ * fast_exp2(2. / 3 * log2_x);

    return 2.5 * log2_x + log2_softplus(-0.5 * log2_x + s);
}

Real SmoothPowerLawSyn::log2_optical_thick_sharp(Real log2_nu) const noexcept {
    if (log2_nu < log2_nu_m) {
        return 2. * (log2_nu - log2_nu_m);
    } else {
        return 2.5 * (log2_nu - log2_nu_m);
    }
}

Real SmoothPowerLawSyn::log2_optical_thin_sharp(Real log2_nu) const noexcept {
    if (log2_nu_m < log2_nu_c) {
        if (log2_nu < log2_nu_m) {
            return (log2_nu - log2_nu_m) / 3.0;
        } else if (log2_nu < log2_nu_c) {
            return 0.5 * (1.0 - p) * (log2_nu - log2_nu_m);
        } else {
            return 0.5 * (1.0 - p) * (log2_nu_c - log2_nu_m) - 0.5 * p * (log2_nu - log2_nu_c);
        }
    } else {
        if (log2_nu < log2_nu_c) {
            return (log2_nu - log2_nu_c) / 3.0;
        } else if (log2_nu < log2_nu_m) {
            return -0.5 * (log2_nu - log2_nu_c);
        } else {
            return -0.5 * (log2_nu_m - log2_nu_c) - 0.5 * p * (log2_nu - log2_nu_m);
        }
    }
}

Real SmoothPowerLawSyn::compute_spectrum(Real nu) const noexcept {
    return fast_exp2(compute_log2_spectrum(fast_log2(nu)));
}

Real SmoothPowerLawSyn::compute_log2_spectrum(Real log2_nu) const noexcept {
    Real log2_f_thin = log2_optical_thin(log2_nu);
    Real log2_f_thick = log2_optical_thick(log2_nu);
    return log2_norm_ + log2_smooth_one(log2_f_thin, log2_f_thick + log2_thick_norm_);
}

void SmoothPowerLawSyn::build() noexcept {
    log2_I_nu_max = fast_log2(I_nu_max);
    log2_nu_m = fast_log2(nu_m);
    log2_nu_c = fast_log2(nu_c);
    log2_nu_a = fast_log2(nu_a);
    log2_nu_M = fast_log2(nu_M);

    // Uniform smoothing (spectral-index-dependent tuning disabled)
    smooth_m_slow_ = 1;
    smooth_c_slow_ = 1;
    smooth_m_fast_ = 1;
    smooth_c_fast_ = 1;
    constexpr Real ln2 = std::numbers::ln2;
    smooth_a_ = (3.5 * p - 1.5) / ln2;

    // Precompute s*(beta_lo - beta_hi) for broken power law transitions
    diff_slope_m_slow_ = smooth_m_slow_ * (0.5 * p - 1.0 / 6.0);
    diff_slope_c_slow_ = smooth_c_slow_ * 0.5;
    diff_slope_m_fast_ = smooth_m_fast_ * 0.5 * (p - 1.0);
    diff_slope_c_fast_ = smooth_c_fast_ * 5.0 / 6.0;

    inv_nu_M_ = 1.0 / nu_M;

    // Compute normalization
    if (nu_m < nu_c) {
        log2_norm_ = 1.0 / smooth_m_slow_;
    } else {
        log2_norm_ = 1.0 / smooth_c_fast_;
    }
    log2_thick_norm_ = log2_optical_thin_sharp(log2_nu_a) - log2_optical_thick_sharp(log2_nu_a);
}

Real SmoothPowerLawSyn::compute_I_nu(Real nu) const noexcept {
    if (nu <= nu_c) { // Below cooling frequency, simple scaling
        return fast_exp(-nu * inv_nu_M_) * I_nu_max * compute_spectrum(nu);
    } else {
        return fast_exp(-nu * inv_nu_M_) * I_nu_max * compute_spectrum(nu) * inverse_compton_correction(*this, nu);
    }
}

Real SmoothPowerLawSyn::compute_log2_I_nu(Real log2_nu) const noexcept {
    constexpr Real log2e = std::numbers::log2e;
    if (log2_nu <= log2_nu_c) { // Below cooling frequency, simple scaling
        return log2_I_nu_max + compute_log2_spectrum(log2_nu) - log2e * inv_nu_M_ * fast_exp2(log2_nu);
    } else {
        const Real nu = fast_exp2(log2_nu);
        return log2_I_nu_max + compute_log2_spectrum(log2_nu) - log2e * inv_nu_M_ * nu +
               fast_log2(inverse_compton_correction(*this, nu));
    }
}
