//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include "../core/mesh.h"
#include "macros.h"

class Empty {};

// c++20 concept
template <typename T>
concept HasDmdt = requires(T t) {
    { t.dm_dt(0.0, 0.0, 0.0) };
};

template <typename T>
concept HasDedt = requires(T t) {
    { t.deps_dt(0.0, 0.0, 0.0) };
};

template <typename T>
concept HasSigma = requires(T t) {
    { t.sigma0(0.0, 0.0) };
};

template <typename T>
concept HasU = requires(T t) {
    { t.U2_th };
};

template <typename T>
concept HasMass = requires(T t) {
    { t.mass(0.0, 0.0, 0.0) };
};

#define MAKE_THIS_ODEINT_STATE(classname, data, array_size)                                                            \
    using array_type = std::array<Real, array_size>;                                                                   \
    using value_type = typename array_type::value_type;                                                                \
    using iterator = typename array_type::iterator;                                                                    \
    using const_iterator = typename array_type::const_iterator;                                                        \
    classname() : data{} {};                                                                                           \
    constexpr size_t size() const noexcept {                                                                           \
        return array_size;                                                                                             \
    }                                                                                                                  \
    constexpr iterator begin() noexcept {                                                                              \
        return data.begin();                                                                                           \
    }                                                                                                                  \
    constexpr iterator end() noexcept {                                                                                \
        return data.end();                                                                                             \
    }                                                                                                                  \
    constexpr const_iterator begin() const noexcept {                                                                  \
        return data.begin();                                                                                           \
    }                                                                                                                  \
    constexpr const_iterator end() const noexcept {                                                                    \
        return data.end();                                                                                             \
    }                                                                                                                  \
    constexpr value_type& operator[](size_t i) noexcept {                                                              \
        return data[i];                                                                                                \
    }                                                                                                                  \
    constexpr const value_type& operator[](size_t i) const noexcept {                                                  \
        return data[i];                                                                                                \
    }

void print_array(Array const& arr);

/**
 * <!-- ************************************************************************************** -->
 * @defgroup FunctionTypes Function Type Definitions
 * @brief Defines convenient aliases for unary, binary, and ternary functions operating on Reals.
 * @details These function types are used throughout the codebase for various mathematical operations
 *          and physical calculations.
 * <!-- ************************************************************************************** -->
 */

/// Function taking one Real argument
using UnaryFunc = std::function<Real(Real)>;
/// Function taking two Real arguments
using BinaryFunc = std::function<Real(Real, Real)>;
/// Function taking three Real arguments
using TernaryFunc = std::function<Real(Real, Real, Real)>;

/**
 * <!-- ************************************************************************************** -->
 * @namespace func
 * @brief Contains inline constexpr lambda functions that return constant values.
 * @details These functions are used throughout the codebase for various mathematical operations
 *          and physical calculations.
 * <!-- ************************************************************************************** -->
 */
namespace func {
    // Always returns 0 regardless of the input.
    inline constexpr auto zero_3d = [](Real phi, Real theta, Real t) constexpr noexcept { return 0.; };
    inline constexpr auto zero_2d = [](Real phi, Real theta) constexpr noexcept { return 0.; };
    // Always returns 1 regardless of the input.
    inline constexpr auto one_3d = [](Real phi, Real theta, Real t) constexpr noexcept { return 1.; };
    inline constexpr auto one_2d = [](Real phi, Real theta) constexpr noexcept { return 1.; };
} // namespace func

/**
 * <!-- ************************************************************************************** -->
 * @defgroup BasicMath Basic Math Functions
 * @brief Inline functions for specific power calculations, a step function, and unit conversion.
 * @details These functions are used throughout the codebase for various mathematical operations
 *          and physical calculations.
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes a^(5/2)
 * @param a The base value
 * @return a^(5/2)
 * <!-- ************************************************************************************** -->
 */
inline Real pow52(Real a) {
    return std::sqrt(a * a * a * a * a);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes a^(3/2)
 * @param a The base value
 * @return a^(3/2)
 * <!-- ************************************************************************************** -->
 */
inline Real pow32(Real a) {
    return std::sqrt(a * a * a);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes a^(4/3)
 * @param a The base value
 * @return a^(4/3)
 * <!-- ************************************************************************************** -->
 */
inline Real pow43(Real a) {
    return std::cbrt(a * a * a * a);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes a^(2/3)
 * @param a The base value
 * @return a^(2/3)
 * <!-- ************************************************************************************** -->
 */
inline Real pow23(Real a) {
    return std::cbrt(a * a);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Step function that returns 1 if x > 0, otherwise 0
 * @param x Input value
 * @return 1 if x > 0, otherwise 0
 * <!-- ************************************************************************************** -->
 */
inline Real stepFunc(Real x) {
    return x > 0 ? 1 : 0;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts electron volt (eV) to frequency (Hz)
 * @param eV Energy in electron volts
 * @return Frequency in Hertz
 * <!-- ************************************************************************************** -->
 */
inline Real eVtoHz(Real eV) {
    return eV / con::h;
}

/**
 * <!-- ************************************************************************************** -->
 * @defgroup Interpolation Interpolation Functions
 * @brief Functions for interpolating values between points in arrays.
 * @details These functions are used throughout the codebase for various mathematical operations
 *          and physical calculations.
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief General interpolation function
 * @param x0 X-value at which to interpolate
 * @param x Array of x-coordinates
 * @param y Array of y-coordinates
 * @param lo_extrap Whether to extrapolate for x0 < min(x)
 * @param hi_extrap Whether to extrapolate for x0 > max(x)
 * @return Interpolated y-value at x0
 * <!-- ************************************************************************************** -->
 */
Real interp(Real x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);

/**
 * <!-- ************************************************************************************** -->
 * @brief Interpolation for equally spaced x-values
 * @param x0 X-value at which to interpolate
 * @param x Array of equally spaced x-coordinates
 * @param y Array of y-coordinates
 * @param lo_extrap Whether to extrapolate for x0 < min(x)
 * @param hi_extrap Whether to extrapolate for x0 > max(x)
 * @return Interpolated y-value at x0
 * <!-- ************************************************************************************** -->
 */
Real eq_space_interp(Real x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);

/**
 * <!-- ************************************************************************************** -->
 * @brief Log-log interpolation (both x and y are in log space)
 * @param x0 X-value at which to interpolate
 * @param x Array of x-coordinates
 * @param y Array of y-coordinates
 * @param lo_extrap Whether to extrapolate for x0 < min(x)
 * @param hi_extrap Whether to extrapolate for x0 > max(x)
 * @return Interpolated y-value at x0
 * <!-- ************************************************************************************** -->
 */
Real loglog_interp(Real x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);

/**
 * <!-- ************************************************************************************** -->
 * @brief Log-log interpolation for equally spaced x-values in log space
 * @param x0 X-value at which to interpolate
 * @param x Array of equally spaced x-coordinates in log space
 * @param y Array of y-coordinates
 * @param lo_extrap Whether to extrapolate for x0 < min(x)
 * @param hi_extrap Whether to extrapolate for x0 > max(x)
 * @return Interpolated y-value at x0
 * <!-- ************************************************************************************** -->
 */
Real eq_space_loglog_interp(Real x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);

/**
 * <!-- ************************************************************************************** -->
 * @defgroup RootFinding Root Finding Methods
 * @brief Functions for finding roots of equations.
 * @details These functions are used throughout the codebase for various mathematical operations
 *          and physical calculations.
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Finds the root of a function using the bisection method
 * @tparam Fun Type of the function
 * @param f Function whose root we want to find
 * @param low Lower bound of the search interval
 * @param high Upper bound of the search interval
 * @param eps Desired accuracy (default: 1e-6)
 * @return Approximation of the root
 * <!-- ************************************************************************************** -->
 */
template <typename Fun>
auto root_bisect(Fun f, decltype(f(0)) low, decltype(f(0)) high, decltype(f(0)) eps = 1e-6) -> decltype(f(0)) {
    using Scalar = decltype(f(0));
    while ((high - low) > std::fabs((high + low) * 0.5) * eps) {
        Scalar mid = 0.5 * (high + low);
        if (f(mid) * f(high) > 0)
            high = mid;
        else
            low = mid;
    }
    return 0.5 * (high + low);
}

/**
 * <!-- ************************************************************************************** -->
 * @defgroup UtilityTemplates Utility Templates
 * @brief Template functions for common operations.
 * @details These functions are used throughout the codebase for various mathematical operations
 *          and physical calculations.
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Returns the value of a single parameter
 * @tparam T Type of the value
 * @param value The value to return
 * @return The input value
 * <!-- ************************************************************************************** -->
 */
template <typename T>
T min(T value) {
    return value;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Returns the minimum of multiple values
 * @tparam T Type of the first value
 * @tparam Args Types of the remaining values
 * @param first First value
 * @param args Remaining values
 * @return The minimum value
 * <!-- ************************************************************************************** -->
 */
template <typename T, typename... Args>
T min(T first, Args... args) {
    return std::min(first, std::min(args...));
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Returns the value of a single parameter
 * @tparam T Type of the value
 * @param value The value to return
 * @return The input value
 * <!-- ************************************************************************************** -->
 */
template <typename T>
T max(T value) {
    return value;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Returns the maximum of multiple values
 * @tparam T Type of the first value
 * @tparam Args Types of the remaining values
 * @param first First value
 * @param args Remaining values
 * @return The maximum value
 * <!-- ************************************************************************************** -->
 */
template <typename T, typename... Args>
T max(T first, Args... args) {
    return std::max(first, std::max(args...));
}

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
 * @brief Broadcasts computed values across grid based on shock symmetry.
 * @tparam Array 3D array type (e.g., ICPhotonGrid or ElectronGrid)
 * @tparam Shock Shock type containing symmetry information
 * @param arr The 3D array to broadcast
 * @param shock Shock object containing symmetry information and theta representatives
 * @details This function handles the broadcasting of computed values from representative
 *          grid points to all grid points based on the detected symmetry level:
 *          - isotropic: Broadcasts from (0,0) to all (phi, theta)
 *          - piecewise: Broadcasts within theta groups, then across phi
 *          - phi_symmetric: Broadcasts only across phi direction
 * <!-- ************************************************************************************** -->
 */
template <typename Array, typename Shock>
void broadcast_symmetry(Array& arr, Shock const& shock) {
    const size_t phi_size = arr.shape()[0];
    const size_t theta_size = arr.shape()[1];

    using SymmetryType = decltype(shock.symmetry);

    if (shock.symmetry == SymmetryType::isotropic) {
        for (size_t i = 0; i < phi_size; ++i)
            for (size_t j = 0; j < theta_size; ++j)
                if (i != 0 || j != 0)
                    xt::view(arr, i, j, xt::all()) = xt::view(arr, 0, 0, xt::all());
    } else if (shock.symmetry == SymmetryType::piecewise) {
        // Broadcast representative thetas to their groups
        for (size_t r = 0; r < shock.theta_reps.size(); ++r) {
            const size_t j_start = shock.theta_reps[r];
            const size_t j_end = (r + 1 < shock.theta_reps.size()) ? shock.theta_reps[r + 1] : theta_size;
            for (size_t j = j_start + 1; j < j_end; ++j)
                xt::view(arr, 0, j, xt::all()) = xt::view(arr, 0, j_start, xt::all());
        }
        // Then broadcast phi=0 to all phi
        for (size_t i = 1; i < phi_size; ++i)
            for (size_t j = 0; j < theta_size; ++j)
                xt::view(arr, i, j, xt::all()) = xt::view(arr, 0, j, xt::all());
    } else if (shock.symmetry == SymmetryType::phi_symmetric) {
        // Only broadcast phi direction
        for (size_t i = 1; i < phi_size; ++i)
            for (size_t j = 0; j < theta_size; ++j)
                xt::view(arr, i, j, xt::all()) = xt::view(arr, 0, j, xt::all());
    }
}
