//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include "../core/mesh.h"
#include "fast-math.h"
#include "macros.h"
#include "traits.h"

class Empty {};

struct SpectralSegment {
    Real slope;      // Power law index
    Real log2_lower; // log2(lower bound)
    Real log2_const; // log2(norm) - slope * log2(lower)
};

template <int N>
struct BrokenPowerLaw {
    void clear() { size_ = 0; }

    void first_segment(Real norm, Real lower, Real slope) {
        size_ = 0;
        const Real log2_lower = fast_log2(lower);
        const Real log2_val = fast_log2(norm);
        auto& s = data_[size_++];
        s.slope = slope;
        s.log2_lower = log2_lower;
        s.log2_const = log2_val - slope * log2_lower;
    }

    void add_segment(Real lower, Real slope) {
        const Real log2_lower = fast_log2(lower);
        const int prev = size_ - 1;
        const Real log2_val = data_[prev].log2_const + data_[prev].slope * log2_lower;
        auto& s = data_[size_++];
        s.slope = slope;
        s.log2_lower = log2_lower;
        s.log2_const = log2_val - slope * log2_lower;
    }

    [[nodiscard]] Real eval(Real x) const {
        const Real log2_x = fast_log2(x);
        for (int i = size_ - 1; i > 0; --i) {
            if (log2_x >= data_[i].log2_lower) {
                return fast_exp2(data_[i].log2_const + data_[i].slope * log2_x);
            }
        }
        return (size_ > 0) ? fast_exp2(data_[0].log2_const + data_[0].slope * log2_x) : 0.0;
    }

    [[nodiscard]] Real log2_eval(Real log2_x) const {
        for (int i = size_ - 1; i > 0; --i) {
            if (log2_x >= data_[i].log2_lower) {
                return data_[i].log2_const + data_[i].slope * log2_x;
            }
        }
        return data_[0].log2_const + data_[0].slope * log2_x;
    }

  private:
    int size_{0};
    std::array<SpectralSegment, N> data_{};
};

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
    inline constexpr auto zero_3d = [](Real /*phi*/, Real /*theta*/, Real /*t*/) constexpr noexcept { return 0.; };
    inline constexpr auto zero_2d = [](Real /*phi*/, Real /*theta*/) constexpr noexcept { return 0.; };
    // Always returns 1 regardless of the input.
    inline constexpr auto one_3d = [](Real /*phi*/, Real /*theta*/, Real /*t*/) constexpr noexcept { return 1.; };
    inline constexpr auto one_2d = [](Real /*phi*/, Real /*theta*/) constexpr noexcept { return 1.; };
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
inline Real pow52(Real a) noexcept {
    return std::sqrt(a * a * a * a * a);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes a^(3/2)
 * @param a The base value
 * @return a^(3/2)
 * <!-- ************************************************************************************** -->
 */
inline Real pow32(Real a) noexcept {
    return std::sqrt(a * a * a);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes a^(4/3)
 * @param a The base value
 * @return a^(4/3)
 * <!-- ************************************************************************************** -->
 */
inline Real pow43(Real a) noexcept {
    return std::cbrt(a * a * a * a);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes a^(2/3)
 * @param a The base value
 * @return a^(2/3)
 * <!-- ************************************************************************************** -->
 */
inline Real pow23(Real a) noexcept {
    return std::cbrt(a * a);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Step function that returns 1 if x > 0, otherwise 0
 * @param x Input value
 * @return 1 if x > 0, otherwise 0
 * <!-- ************************************************************************************** -->
 */
inline constexpr Real stepFunc(Real x) noexcept {
    return x > 0 ? 1 : 0;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts electron volt (eV) to frequency (Hz)
 * @param eV Energy in electron volts
 * @return Frequency in Hertz
 * <!-- ************************************************************************************** -->
 */
inline constexpr Real eVtoHz(Real eV) noexcept {
    return eV / con::h;
}

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
constexpr T min(T value) {
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
constexpr T min(T first, Args... args) {
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
constexpr T max(T value) {
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
constexpr T max(T first, Args... args) {
    return std::max(first, std::max(args...));
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Broadcasts computed values across grid based on coord symmetry.
 * @tparam Array 3D array type (e.g., ICPhotonGrid or ElectronGrid)
 * @param arr The 3D array to broadcast
 * @param coord Coordinate object containing symmetry information and theta representatives
 * @details This function handles the broadcasting of computed values from representative
 *          grid points to all grid points based on the detected symmetry level:
 *          - isotropic: Broadcasts from (0,0) to all (phi, theta)
 *          - piecewise: Broadcasts within theta groups, then across phi
 *          - phi_symmetric: Broadcasts only across phi direction
 * <!-- ************************************************************************************** -->
 */
template <typename Array>
void broadcast_symmetry(Array& arr, Coord const& coord) {
    const size_t phi_size = arr.shape()[0];
    const size_t theta_size = arr.shape()[1];

    if (coord.symmetry == Symmetry::isotropic) {
        for (size_t i = 0; i < phi_size; ++i)
            for (size_t j = 0; j < theta_size; ++j)
                if (i != 0 || j != 0)
                    xt::view(arr, i, j, xt::all()) = xt::view(arr, 0, 0, xt::all());
    } else if (coord.symmetry == Symmetry::piecewise) {
        // Broadcast representative thetas to their groups
        for (size_t r = 0; r < coord.theta_reps.size(); ++r) {
            const size_t j_start = coord.theta_reps[r];
            const size_t j_end = (r + 1 < coord.theta_reps.size()) ? coord.theta_reps[r + 1] : theta_size;
            for (size_t j = j_start + 1; j < j_end; ++j)
                xt::view(arr, 0, j, xt::all()) = xt::view(arr, 0, j_start, xt::all());
        }
        // Then broadcast phi=0 to all phi
        for (size_t i = 1; i < phi_size; ++i)
            for (size_t j = 0; j < theta_size; ++j)
                xt::view(arr, i, j, xt::all()) = xt::view(arr, 0, j, xt::all());
    } else if (coord.symmetry == Symmetry::phi_symmetric) {
        // Only broadcast phi direction
        for (size_t i = 1; i < phi_size; ++i)
            for (size_t j = 0; j < theta_size; ++j)
                xt::view(arr, i, j, xt::all()) = xt::view(arr, 0, j, xt::all());
    }
}
