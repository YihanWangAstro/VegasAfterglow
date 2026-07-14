//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "smooth-power-law-syn.h"

#include <algorithm>

//========================================================================================================
//                                  SmoothPowerLawSyn Class Methods
//========================================================================================================
inline Real log2_smooth_one(Real log2_a, Real log2_b, Real s) noexcept {
    return log2_a - log2_softplus(s * (log2_a - log2_b)) / s;
}

Real SmoothPowerLawSyn::log2_optical_thin(Real log2_nu) const noexcept {
    return (log2_nu - log2_nu_lo_) / 3.0 + log2_broken_power_ratio(log2_nu, log2_nu_lo_, diff_lo_, smooth_lo_) +
           log2_broken_power_ratio(log2_nu, log2_nu_hi_, diff_hi_, smooth_hi_);
}

Real SmoothPowerLawSyn::log2_optical_thick(Real log2_nu) const noexcept {
    const Real log2_x = log2_nu - log2_nu_m;

    // Beyond log2_x_far_ the softplus argument -0.5*log2_x + s is below the -20
    // shortcut of log2_softplus, which then returns exactly 0 — skip its exp2.
    if (log2_x > log2_x_far_) {
        return 2.5 * log2_x;
    }

    const Real s = -smooth_thick_ * fast_exp2(2. / 3 * log2_x);
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
    // IC steepens the spectrum above nu_c; the SSA thick branch is set by the source
    // function, so the correction applies to the thin component only, before the blend.
    if (log2_nu > log2_nu_c && has_IC_correction(*this)) {
        const Real nu = fast_exp2(log2_nu);
        log2_f_thin += fast_log2(inverse_compton_correction(*this, nu));
    }
    return log2_norm_ + log2_smooth_one(log2_f_thin, log2_f_thick + log2_thick_norm_, s_a_blend_);
}

void SmoothPowerLawSyn::build() noexcept {
    log2_I_nu_max = fast_log2(I_nu_max);
    log2_nu_m = fast_log2(nu_m);
    log2_nu_c = fast_log2(nu_c);
    log2_nu_a = fast_log2(nu_a);
    log2_nu_M = fast_log2(nu_M);
    inv_nu_M_ = 1.0 / nu_M;

    constexpr Real ln2 = std::numbers::ln2;
    // 2 -> 5/2 SSA turnover sharpness (pre-G&S convention; p is fixed per model, so
    // it cannot vary along a light curve), /ln2: the thick form uses 2^ not exp.
    smooth_thick_ = (3.5 * p - 1.5) / ln2;
    // Where the optical-thick softplus provably hits its -20 shortcut:
    // smooth_thick_ * 2^(2*log2_x/3) >= 20 alone pushes the argument below -20
    // for log2_x > 0 (the -0.5*log2_x term only helps).
    log2_x_far_ = 1.5 * fast_log2(20.0 / smooth_thick_);

    // ----- Cooling regime by break order (hard branch, original convention).
    // Slow cooling: slopes 1/3 -> -(p-1)/2 -> -p/2, breaks at nu_m then nu_c.
    // Fast cooling: slopes 1/3 -> -1/2     -> -p/2, breaks at nu_c then nu_m.
    log2_nu_lo_ = std::min(log2_nu_m, log2_nu_c);
    log2_nu_hi_ = std::max(log2_nu_m, log2_nu_c);

    // Uniform break smoothing (s = 1): break-dependent sharpness would tie the
    // spectral shape to the evolving breaks and imprint jumps on light curves.
    smooth_lo_ = 1;
    smooth_hi_ = 1;
    const Real alpha_mid = (log2_nu_m <= log2_nu_c) ? -0.5 * (p - 1.0) : -0.5;
    diff_lo_ = smooth_lo_ * (1.0 / 3.0 - alpha_mid);
    diff_hi_ = smooth_hi_ * (alpha_mid + 0.5 * p);

    // Thin/thick join at nu_a: uniform sharpness as well.
    s_a_blend_ = 1;

    // Analytic 1/smooth_lo_ (G&S F_b convention): sub-break flux stays consistent
    // across nu_c positions; the peak drops naturally when breaks are close.
    log2_norm_ = 1.0 / smooth_lo_;
    log2_thick_norm_ = log2_optical_thin_sharp(log2_nu_a) - log2_optical_thick_sharp(log2_nu_a);
}

Real SmoothPowerLawSyn::compute_I_nu(Real nu) const noexcept {
    return fast_exp(-nu * inv_nu_M_) * I_nu_max * compute_spectrum(nu);
}

Real SmoothPowerLawSyn::compute_log2_I_nu(Real log2_nu) const noexcept {
    const Real spec = log2_I_nu_max + compute_log2_spectrum(log2_nu);
    // Skip the exponential-cutoff term when nu/nu_M < 2^-20 (negligible).
    if (log2_nu - log2_nu_M < -20)
        return spec;
    constexpr Real log2e = std::numbers::log2e;
    return spec - log2e * inv_nu_M_ * fast_exp2(log2_nu);
}
