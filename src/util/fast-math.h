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
 * @param x The input value (positive normal double)
 * @return log2(x)
 * @details Subtract-trick range reduction (mantissa centered on sqrt(2) with no compare/select)
 *          + degree-3 minimax polynomial in u = (z-1)/(z+1). Branchless, NEON/SSE
 *          autovectorizable (the division pipelines; it is not the throughput limiter).
 *          Max abs error 5.6e-6 (measured exhaustively over the mantissa range).
 * <!-- ************************************************************************************** -->
 */
inline double fast_log2(double x) noexcept {
#ifdef AFTERGLOW_FAST_MATH
    // Subtracting the bit pattern of sqrt(2)/2 makes the exponent field of `tmp` carry
    // round-to-nearest(log2(x)), so z = x / 2^k lands in [sqrt(2)/2, sqrt(2)) directly —
    // replaces the mask/or + compare + two selects of the classic reduction.
    constexpr uint64_t off = 0x3FE6A09E00000000ull; // bits of ~sqrt(2)/2
    constexpr uint64_t exp_mask_hi = 0xFFF0000000000000ull;

    uint64_t ix = std::bit_cast<uint64_t>(x);
    uint64_t tmp = ix - off;
    int64_t k = static_cast<int64_t>(tmp) >> 52;
    double z = std::bit_cast<double>(ix - (tmp & exp_mask_hi));
    double e = static_cast<double>(k);

    // z in [sqrt(2)/2, sqrt(2)); u in (-0.172, 0.172)
    double u = (z - 1.0) / (z + 1.0);
    double u2 = u * u;

    // Degree-3 minimax (Remez) refit of the atanh form; 16x more accurate than the
    // Taylor coefficients at identical cost. Max abs err 5.568e-6.
    constexpr double c1 = 2.885228583309733;
    constexpr double c3 = 0.9835335869591781;

    double poly = c1 + u2 * c3;

    return e + u * poly;
#else
    return std::log2(x);
#endif
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast approximation of 2 raised to a power
 * @param x The exponent
 * @return 2^x
 * @details Branchless clamp + round-to-nearest split + degree-3 minimax polynomial for 2^f
 *          constrained to p(0) = 1, with the 2^i factor fused by integer addition on the
 *          exponent bits (no final multiply). Max rel error 4.3e-4 (measured); EXACT powers
 *          of two at integer x (so identities like exp(0) = 1 survive). Never returns inf:
 *          x >= 1023 clamps to the largest representable power; x <= -1022.5 underflows to ~0.
 *
 *          Clamp proof: x in [-1022.5, 1023] => r = nearbyint(x) in [-1022, 1023]
 *          (ties-to-even). p(f) in (0.70, 1.42) has biased exponent 1022 or 1023, so the
 *          fused biased exponent exp(p) + r spans [0, 2046] — never negative, never 2047:
 *          at r = 1023 the clamp forces f <= 0 where p <= p(0) = 1, giving at most
 *          1023 + 1023 = 2046 (the largest finite exponent, value 2^1023).
 * <!-- ************************************************************************************** -->
 */
inline double fast_exp2(double x) noexcept {
#ifdef AFTERGLOW_FAST_MATH
    x = std::fmax(x, -1022.5);
    x = std::fmin(x, 1023.0);

    double r = std::nearbyint(x); // single instruction (frintn / roundsd)
    double f = x - r;             // f in [-0.5, 0.5]
    int64_t i = static_cast<int64_t>(r);

    // Degree-3 minimax (Remez) for 2^f on [-0.5, 0.5], constrained to p(0) = 1;
    // max rel err 4.315e-4, ~2x better than the Taylor coefficients at identical cost.
    constexpr double c1 = 0x1.643b8cd1fb0ddp-1;
    constexpr double c2 = 0x1.efadb16a62f7ap-3;
    constexpr double c3 = 0x1.55989d930a31ap-5;

    double p = 1.0 + f * (c1 + f * (c2 + f * c3));

    // 2^i fused by integer add on the exponent field of p (see clamp proof above).
    return std::bit_cast<double>(std::bit_cast<int64_t>(p) + (i << 52));
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
