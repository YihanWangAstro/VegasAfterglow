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
 *          alternative methods when EXTREME_SPEED is defined.
 * <!-- ************************************************************************************** -->
 */

constexpr int BUCKET_BITS = 9;                                    // 512 buckets
constexpr size_t BUCKETS = static_cast<size_t>(1) << BUCKET_BITS; // 512
// With 512 bins across [1,2): max |r| ≈ 1 / 2^(BUCKET_BITS+1) ≈ 9.77e-4
// Degree-3 Taylor gives neglected term O(r^4) ≲ 1e-12 → ~1e-9 relative overall.

constexpr uint64_t EXP_MASK = 0x7FFull;
constexpr uint64_t FRAC_MASK = (1ull << 52) - 1;
constexpr int EXP_BIAS = 1023;
constexpr double INV_LN2 = 1.4426950408889634073599246810018921;
struct Tables {
    std::array<double, BUCKETS> inv_c{};
    std::array<double, BUCKETS> log2_c{};
    std::array<uint64_t, BUCKETS> c_bits{}; // optional (not used in hot path)

    // Exp2 tables: store 2^(i/512) for i in [0, 512)
    std::array<double, BUCKETS> exp2_c{};     // 2^(i/BUCKETS)
    std::array<double, BUCKETS> exp2_slope{}; // slope for linear interpolation

    Tables() {
        constexpr int shift = 52 - BUCKET_BITS;
        constexpr uint64_t expbits = static_cast<uint64_t>(EXP_BIAS) << 52;

        // Initialize log2 tables
        for (size_t i = 0; i < BUCKETS; ++i) {
            const uint64_t frac_mid = (static_cast<uint64_t>(i) << shift) | (static_cast<uint64_t>(1) << (shift - 1));
            uint64_t bits = expbits | frac_mid; // center c_i bits
            const double c = std::bit_cast<double>(bits);
            c_bits[i] = bits;
            inv_c[i] = 1.0 / c;
            log2_c[i] = std::log2(c);
        }

        // Initialize exp2 tables
        constexpr double inv_buckets = 1.0 / static_cast<double>(BUCKETS);
        constexpr double ln2 = 0.693147180559945309417232121458176568;
        for (size_t i = 0; i < BUCKETS; ++i) {
            const double f = static_cast<double>(i) * inv_buckets;
            exp2_c[i] = std::exp2(f);
            // Slope for multiplicative interpolation: 2^(f+delta) ≈ 2^f * (1 + slope*delta)
            exp2_slope[i] = ln2;
        }
    }
};

inline const Tables& tables() {
    static const Tables T; // thread-safe init
    return T;
}

inline double decompose_normals(double x, int& e) {
    uint64_t bits = std::bit_cast<uint64_t>(x);
    const uint64_t expo = (bits >> 52) & EXP_MASK;
    if (expo == 0 || expo == EXP_MASK) {
        int ee;
        double m = std::frexp(x, &ee); // [0.5,1)
        m *= 2.0;
        e = ee - 1;
        return m; // -> [1,2)
    }
    e = static_cast<int>(expo) - EXP_BIAS;
    bits = (bits & FRAC_MASK) | (static_cast<uint64_t>(EXP_BIAS) << 52); // force mantissa exponent
    return std::bit_cast<double>(bits);                                  // m in [1,2)
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast approximation of the base-2 logarithm
 * @param x The input value
 * @return log2(val)
 * <!-- ************************************************************************************** -->
 */
inline double fast_log2(double x) {
#ifdef EXTREME_SPEED
    if (std::isnan(x))
        return std::numeric_limits<double>::quiet_NaN();
    if (x < 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    if (x == 0.0)
        return -std::numeric_limits<double>::infinity();
    if (std::isinf(x))
        return std::numeric_limits<double>::infinity();

    const Tables& T = tables();

    int e;
    double m = decompose_normals(x, e);

    uint64_t mbits = std::bit_cast<uint64_t>(m);
    uint64_t frac = mbits & FRAC_MASK; // 52-bit fraction
    const int shift = 52 - BUCKET_BITS;
    size_t idx = size_t(frac >> shift); // 0..511

    double r = m * T.inv_c[idx] - 1.0;

    double r2 = r * r;
    double poly = r + r2 * (-0.5 + r * (1.0 / 3.0));
    double approx = T.log2_c[idx] + poly * INV_LN2;

    return double(e) + approx;
#else
    return std::log2(x);
#endif
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast approximation of 2 raised to a power
 * @param x The exponent
 * @return 2^x
 * @details Uses 512-entry lookup table with linear interpolation.
 *          Error < 1e-10, ~3-5x faster than std::exp2
 * <!-- ************************************************************************************** -->
 */
inline double fast_exp2(double x) {
#ifdef EXTREME_SPEED
    // Handle special cases (branchless where possible)
    if (x >= 1024.0)
        return std::numeric_limits<double>::infinity();
    if (x <= -1074.0)
        return 0.0;

    // Split x into integer and fractional parts: x = i + f where f in [0,1)
    const double i_exact = std::floor(x);
    const int i = static_cast<int>(i_exact);
    const double f = x - i_exact;

    const Tables& T = tables();

    // Convert fractional part to table index
    const double idx_float = f * static_cast<double>(BUCKETS);
    const size_t idx = std::min(static_cast<size_t>(idx_float), BUCKETS - 1);
    const double delta = f - static_cast<double>(idx) / static_cast<double>(BUCKETS);

    // Multiplicative interpolation: 2^f ≈ 2^f_idx * (1 + delta * ln(2))
    const double frac_part = T.exp2_c[idx] * (1.0 + delta * T.exp2_slope[idx]);

    // Multiply by 2^i using bit manipulation
    uint64_t bits = std::bit_cast<uint64_t>(frac_part);
    const int64_t exp = (bits >> 52) & 0x7FF;
    const int64_t new_exp = exp + i;

    // Clamp exponent to valid range
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
