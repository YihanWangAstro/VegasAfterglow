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

inline Real sigmoid2(Real x) noexcept {
    return 1.0 / (1.0 + fast_exp2(-x));
}
inline Real blend(Real w, Real a, Real b) noexcept {
    return w * a + (1.0 - w) * b;
}

Real SmoothPowerLawSyn::log2_optical_thin(Real log2_nu) const noexcept {
    // Single unified double-smoothed expression, no regime branching. The
    // cached log2_nu_lo_/log2_nu_hi_ and the smooth_lo_/smooth_hi_ blend in
    // build() smoothly interpolate between slow-cooling and fast-cooling
    // Granot & Sari (2002) Table 2 limits as nu_m and nu_c approach each other.
    return (log2_nu - log2_nu_lo_) / 3.0 + log2_broken_power_ratio(log2_nu, log2_nu_lo_, diff_lo_, smooth_lo_) +
           log2_broken_power_ratio(log2_nu, log2_nu_hi_, diff_hi_, smooth_hi_);
}

Real SmoothPowerLawSyn::log2_optical_thick(Real log2_nu) const noexcept {
    const Real log2_x = log2_nu - log2_nu_m;
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
    if (log2_nu > log2_nu_c) {
        // IC steepens the emission spectrum above nu_c. The SSA-dominated thick
        // branch is set by the source function and is not directly modified, so
        // apply the IC correction to the thin component only, before the blend.
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
    // G&S Table 2 b=4 (ISM, k=0): 2 -> 5/2 transition at nu_m, s = 3.44p - 1.41.
    // Divided by ln2 because the optical-thick form uses 2^(...) where G&S uses exp(...).
    smooth_thick_ = (3.44 * p - 1.41) / ln2;

    constexpr Real s_swap = 4.0;  // sigmoid sharpness for all regime blends (~1 octave window)
    constexpr Real s_floor = 0.1; // protect against atypical p flipping G&S linear fits negative

    // ----- Cooling-regime blend: deep slow (nu_m << nu_c) vs deep fast (nu_c << nu_m).
    // Slow cooling: slopes 1/3 -> -(p-1)/2 -> -p/2, breaks at nu_m (b=2) then nu_c (b=3).
    // Fast cooling: slopes 1/3 -> -1/2     -> -p/2, breaks at nu_c (b=11) then nu_m (b=9).
    const Real w_slow = sigmoid2(s_swap * (log2_nu_c - log2_nu_m));

    // Soft min/max of (log2_nu_m, log2_nu_c); at equality, separated by 2/s_swap.
    const Real soft_offset = log2_softplus(-s_swap * std::abs(log2_nu_c - log2_nu_m)) / s_swap;
    log2_nu_lo_ = std::min(log2_nu_m, log2_nu_c) - soft_offset;
    log2_nu_hi_ = std::max(log2_nu_m, log2_nu_c) + soft_offset;

    // G&S Table 2 ISM (k=0) smoothing parameters for the four optically-thin breaks.
    const Real s_m_slow = std::max(1.84 - 0.40 * p, s_floor); // b=2
    const Real s_c_slow = std::max(1.15 - 0.06 * p, s_floor); // b=3
    constexpr Real s_c_fast = 0.597;                          // b=11
    const Real s_m_fast = std::max(3.34 - 0.82 * p, s_floor); // b=9

    smooth_lo_ = blend(w_slow, s_m_slow, s_c_fast); // lo break: slow nu_m <-> fast nu_c
    smooth_hi_ = blend(w_slow, s_c_slow, s_m_fast); // hi break: slow nu_c <-> fast nu_m
    const Real alpha_mid = blend(w_slow, -0.5 * (p - 1.0), -0.5);
    diff_lo_ = smooth_lo_ * (1.0 / 3.0 - alpha_mid);
    diff_hi_ = smooth_hi_ * (alpha_mid + 0.5 * p);

    // ----- nu_a position blend: below/between/above (nu_m, nu_c).
    // G&S Table 2: b=1 (below), b=5 (between, slow), b=6 (above).
    // Regime 4 (nu_c < nu_a < nu_m, fast) reuses b=5 as approximation.
    const Real u = sigmoid2(s_swap * (log2_nu_a - log2_nu_m)); // -> 1 when nu_a > nu_m
    const Real v = sigmoid2(s_swap * (log2_nu_a - log2_nu_c)); // -> 1 when nu_a > nu_c
    const Real w_below = (1.0 - u) * (1.0 - v);
    const Real w_above = u * v;
    constexpr Real s_a_below = 1.64;                           // b=1
    const Real s_a_mid = std::max(1.47 - 0.21 * p, s_floor);   // b=5
    const Real s_a_above = std::max(0.94 - 0.14 * p, s_floor); // b=6
    s_a_blend_ = w_below * s_a_below + w_above * s_a_above + (1.0 - w_below - w_above) * s_a_mid;

    // Self-normalise: peak value of optical-thin at log2_nu_lo_ becomes 1.0 linear.
    log2_norm_ = -log2_optical_thin(log2_nu_lo_);
    log2_thick_norm_ = log2_optical_thin_sharp(log2_nu_a) - log2_optical_thick_sharp(log2_nu_a);
}

Real SmoothPowerLawSyn::compute_I_nu(Real nu) const noexcept {
    // IC correction folded into compute_log2_spectrum (thin branch only).
    return fast_exp(-nu * inv_nu_M_) * I_nu_max * compute_spectrum(nu);
}

Real SmoothPowerLawSyn::compute_log2_I_nu(Real log2_nu) const noexcept {
    constexpr Real log2e = std::numbers::log2e;
    return log2_I_nu_max + compute_log2_spectrum(log2_nu) - log2e * inv_nu_M_ * fast_exp2(log2_nu);
}
