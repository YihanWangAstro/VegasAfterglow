//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <bit>
#include <cmath>
#include <cstdint>
#include <numbers>

#include "macros.h"

/**
 * <!-- ************************************************************************************** -->
 * @defgroup FastMath Fast Math Functions
 * @brief Optimized versions of common mathematical functions.
 * @details These functions provide fast approximations of exponential and logarithm functions using
 *          polynomial approximations when AFTERGLOW_FAST_MATH is defined.
 * <!-- ************************************************************************************** -->
 */

constexpr uint64_t FRAC_MASK = (1ull << 52) - 1;

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast approximation of the base-2 logarithm
 * @param x The input value
 * @return log2(x)
 * @details Bit manipulation + degree-3 polynomial via atanh substitution. Branchless.
 *          Raw rel error ~4e-6; tuned for the GRB afterglow pipeline accuracy budget.
 * <!-- ************************************************************************************** -->
 */
inline double fast_log2(double x) noexcept {
#ifdef AFTERGLOW_FAST_MATH
    uint64_t bits = std::bit_cast<uint64_t>(x);

    // Extract exponent
    int64_t e = static_cast<int64_t>((bits >> 52) & 0x7FF) - 1023;

    // Normalize mantissa to [1, 2)
    bits = (bits & FRAC_MASK) | (1023ull << 52);
    double m = std::bit_cast<double>(bits);

    // Branchless centering around sqrt(2) so |u| stays small for fast polynomial convergence.
    // If m > sqrt(2), use m/2 and shift exponent up by 1.
    const bool large = m > std::numbers::sqrt2;
    m = large ? m * 0.5 : m;
    e += large;

    // m is now in [sqrt(2)/2, sqrt(2)] ≈ [0.707, 1.414]; u in [-0.172, 0.172]
    double u = (m - 1.0) / (m + 1.0);
    double u2 = u * u;

    // log2(m) = 2 * INV_LN2 * (u + u^3/3); raw error ~3.6e-6.
    constexpr double c1 = 2.8853900817779268; // 2/ln(2)
    constexpr double c3 = 0.9617966939259756; // 2/(3*ln(2))

    double poly = c1 + u2 * c3;

    return static_cast<double>(e) + u * poly;
#else
    return std::log2(x);
#endif
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast approximation of 2 raised to a power
 * @param x The exponent
 * @return 2^x
 * @details Branchless clamp + degree-3 Taylor on fractional part + direct exponent injection.
 *          Raw rel error ~8e-3; tuned for the GRB afterglow pipeline accuracy budget.
 * <!-- ************************************************************************************** -->
 */
inline double fast_exp2(double x) noexcept {
#ifdef AFTERGLOW_FAST_MATH
    // Branchless clamp. Bottom is -1023 (not -1074): below that the result rounds to 0 anyway
    // via exp_bits = 0, and clamping here keeps i + 1023 non-negative for the exponent injection.
    x = x < -1023.0 ? -1023.0 : (x > 1023.0 ? 1023.0 : x);

    // Split into integer and fractional parts (round to nearest keeps |f| small).
    double i_part = std::floor(x + 0.5);
    double f = x - i_part; // f in [-0.5, 0.5]
    int64_t i = static_cast<int64_t>(i_part);

    // Polynomial for 2^f, f in [-0.5, 0.5]. Order 3 Taylor; raw error ~8e-3.
    constexpr double c1 = std::numbers::ln2;
    constexpr double c2 = 0.2402265069591007; // ln(2)^2 / 2!
    constexpr double c3 = 0.0555041086648216; // ln(2)^3 / 3!

    double poly = 1.0 + f * (c1 + f * (c2 + f * c3));

    // Multiply by 2^i via direct exponent injection. i in [-1023, 1023] -> i + 1023 in [0, 2046],
    // a valid biased exponent (0 encodes 0.0 / denormals, 2046 encodes the largest normal).
    const uint64_t exp_bits = static_cast<uint64_t>(i + 1023) << 52;
    const double pow2_i = std::bit_cast<double>(exp_bits);
    return poly * pow2_i;
#else
    return std::exp2(x);
#endif
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast approximation of a raised to the power of b
 * @param a The base
 * @param b The exponent
 * @return a^b
 * <!-- ************************************************************************************** -->
 */
inline Real fast_pow(Real a, Real b) noexcept {
    return fast_exp2(b * fast_log2(a));
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast approximation of the exponential function
 * @param x The exponent
 * @return e^x
 * <!-- ************************************************************************************** -->
 */
inline Real fast_exp(Real x) noexcept {
#ifdef AFTERGLOW_FAST_MATH
    return fast_exp2(x * std::numbers::log2e);
#else
    return std::exp(x);
#endif
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast approximation of the natural logarithm
 * @param x The input value
 * @return ln(x)
 * <!-- ************************************************************************************** -->
 */
inline double fast_log(double x) noexcept {
#ifdef AFTERGLOW_FAST_MATH
    return fast_log2(x) * std::numbers::ln2;
#else
    return std::log(x);
#endif
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes log2(1+2^x)
 * @param x The input value
 * @return log2(1+2^x)
 * <!-- ************************************************************************************** -->
 */
inline Real log2_softplus(Real x) noexcept {
    if (x > 20.0)
        return x;
    if (x < -20.0)
        return 0.0;
    return fast_log2(1.0 + fast_exp2(x));
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Smooth broken power-law transition in log2 space.
 * @details Returns a correction term that smoothly changes the slope by delta_beta at the break.
 *          Far below break: correction -> 0 (original slope preserved).
 *          Far above break: correction -> -delta_beta * (log2_x - log2_x_break) / s (slope shifted).
 * @param log2_x        Log2 of the evaluation point
 * @param log2_x_break  Log2 of the break point
 * @param s_delta_beta  s * (slope_below - slope_above), where s controls transition sharpness
 * @param s             Smoothing sharpness parameter (larger = sharper transition)
 * @return Log2-space correction to add to the base power law
 * <!-- ************************************************************************************** -->
 */
inline Real log2_broken_power_ratio(Real log2_x, Real log2_x_break, Real s_delta_beta, Real s) noexcept {
    return -log2_softplus(s_delta_beta * (log2_x - log2_x_break)) / s;
}
