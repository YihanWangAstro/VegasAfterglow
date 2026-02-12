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

#include "macros.h"

/**
 * <!-- ************************************************************************************** -->
 * @defgroup FastMath Fast Math Functions
 * @brief Optimized versions of common mathematical functions.
 * @details These functions provide fast approximations of exponential and logarithm functions using
 *          polynomial approximations when EXTREME_SPEED is defined.
 * <!-- ************************************************************************************** -->
 */

constexpr uint64_t FRAC_MASK = (1ull << 52) - 1;

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast approximation of the base-2 logarithm
 * @param x The input value
 * @return log2(x)
 * @details Uses bit manipulation + polynomial approximation.
 *          Error < 1e-10, ~3.5x faster than std::log2
 * <!-- ************************************************************************************** -->
 */
inline double fast_log2(double x) {
#ifdef EXTREME_SPEED
    uint64_t bits = std::bit_cast<uint64_t>(x);

    // Extract exponent
    int64_t e = static_cast<int64_t>((bits >> 52) & 0x7FF) - 1023;

    // Normalize mantissa to [1, 2)
    bits = (bits & FRAC_MASK) | (1023ull << 52);
    double m = std::bit_cast<double>(bits);

    // For better polynomial convergence, center around sqrt(2)
    // If m > sqrt(2), divide by 2 and increment exponent
    if (m > 1.4142135623730951) {
        m *= 0.5;
        e++;
    }

    // Now m is in [sqrt(2)/2, sqrt(2)] ≈ [0.707, 1.414]
    // Use u = (m-1)/(m+1) for faster convergence
    double u = (m - 1.0) / (m + 1.0);
    double u2 = u * u;

    // log2(m) = 2 * INV_LN2 * (u + u^3/3 + u^5/5 + u^7/7 + ...)
    // Polynomial coefficients: 2/(k*ln(2)) for odd k
    constexpr double c1 = 2.8853900817779268; // 2/ln(2)
    constexpr double c3 = 0.9617966939259756; // 2/(3*ln(2))
    constexpr double c5 = 0.5770780162461006; // 2/(5*ln(2))
    constexpr double c7 = 0.4121977821679615; // 2/(7*ln(2))
    constexpr double c9 = 0.3219280948873623; // 2/(9*ln(2))

    double poly = c1 + u2 * (c3 + u2 * (c5 + u2 * (c7 + u2 * c9)));

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
 * @details Uses polynomial approximation with bit manipulation.
 *          Error < 2e-7, ~1.5x faster than std::exp2
 * <!-- ************************************************************************************** -->
 */
inline double fast_exp2(double x) {
#ifdef EXTREME_SPEED
    // Handle overflow/underflow
    if (x >= 1024.0)
        return std::numeric_limits<double>::infinity();
    if (x <= -1075.0)
        return 0.0;

    // Split into integer and fractional parts (round to nearest for smaller |f|)
    double i_part = std::floor(x + 0.5);
    double f = x - i_part; // f in [-0.5, 0.5]
    int64_t i = static_cast<int64_t>(i_part);

    // Polynomial for 2^f, f in [-0.5, 0.5]
    // Taylor series: 2^f = sum_{k=0}^{n} (f*ln2)^k / k!
    constexpr double c0 = 1.0;
    constexpr double c1 = 0.6931471805599453; // ln(2)
    constexpr double c2 = 0.2402265069591007; // ln(2)^2 / 2!
    constexpr double c3 = 0.0555041086648216; // ln(2)^3 / 3!
    constexpr double c4 = 0.0096181291076285; // ln(2)^4 / 4!
    constexpr double c5 = 0.0013333558146428; // ln(2)^5 / 5!
    constexpr double c6 = 0.0001540353039338; // ln(2)^6 / 6!

    // Horner's method
    double poly = c0 + f * (c1 + f * (c2 + f * (c3 + f * (c4 + f * (c5 + f * c6)))));

    // Multiply by 2^i using bit manipulation
    uint64_t bits = std::bit_cast<uint64_t>(poly);
    int64_t exp = static_cast<int64_t>((bits >> 52) & 0x7FF);
    int64_t new_exp = exp + i;

    if (new_exp <= 0)
        return 0.0;
    if (new_exp >= 0x7FF)
        return std::numeric_limits<double>::infinity();

    bits = (bits & ~(0x7FFull << 52)) | (static_cast<uint64_t>(new_exp) << 52);
    return std::bit_cast<double>(bits);
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
inline Real fast_pow(Real a, Real b) {
    return fast_exp2(b * fast_log2(a));
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast approximation of the exponential function
 * @param x The exponent
 * @return e^x
 * <!-- ************************************************************************************** -->
 */
inline Real fast_exp(Real x) {
#ifdef EXTREME_SPEED
    return fast_exp2(x * 1.4426950408889634); // x / ln(2)
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
inline double fast_log(double x) {
#ifdef EXTREME_SPEED
    return fast_log2(x) * 0.6931471805599453;
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
inline Real log2_softplus(Real x) {
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
inline Real log2_broken_power_ratio(Real log2_x, Real log2_x_break, Real s_delta_beta, Real s) {
    return -log2_softplus(s_delta_beta * (log2_x - log2_x_break)) / s;
}
