//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

namespace afterglow {

#define AFTERGLOW_REQUIRE(condition, message)                                                                          \
    do {                                                                                                               \
        if (!(condition)) [[unlikely]] {                                                                               \
            throw std::invalid_argument(message);                                                                      \
        }                                                                                                              \
    } while (0)

#define AFTERGLOW_ENSURE(condition, message)                                                                           \
    do {                                                                                                               \
        if (!(condition)) [[unlikely]] {                                                                               \
            throw std::logic_error(message);                                                                           \
        }                                                                                                              \
    } while (0)

    // ---- physics-parameter validation helpers ---------------------------------
    //
    // Used at pybind constructor entry points to reject obviously-bad inputs
    // (NaN, inf, negative-where-positive-required) before they propagate into
    // unit conversion and the C++ core. Stringify the parameter name via `#x`
    // so error messages name the offending field.

    inline bool _finite_pos(double x) noexcept {
        return std::isfinite(x) && x > 0;
    }
    inline bool _finite_nonneg(double x) noexcept {
        return std::isfinite(x) && x >= 0;
    }

#define AFTERGLOW_REQUIRE_FINITE_POS(x)                                                                                \
    AFTERGLOW_REQUIRE(::afterglow::_finite_pos(x),                                                                     \
                      std::string(#x) + " must be finite and > 0, got " + std::to_string(x))

#define AFTERGLOW_REQUIRE_FINITE_NONNEG(x)                                                                             \
    AFTERGLOW_REQUIRE(::afterglow::_finite_nonneg(x),                                                                  \
                      std::string(#x) + " must be finite and >= 0, got " + std::to_string(x))

// Inclusive on `hi`, exclusive on `lo`. Use for ratios in (0, 1] and angles in (0, pi/2].
#define AFTERGLOW_REQUIRE_RANGE_OI(x, lo, hi)                                                                          \
    AFTERGLOW_REQUIRE(std::isfinite(x) && (x) > (lo) && (x) <= (hi),                                                   \
                      std::string(#x) + " must be in (" + std::to_string(lo) + ", " + std::to_string(hi) + "], got " + \
                          std::to_string(x))

// Strict lower bound. Use for Gamma0 > 1, p > 1.
#define AFTERGLOW_REQUIRE_GREATER_THAN(x, lo)                                                                          \
    AFTERGLOW_REQUIRE(std::isfinite(x) && (x) > (lo),                                                                  \
                      std::string(#x) + " must be finite and > " + std::to_string(lo) + ", got " + std::to_string(x))

} // namespace afterglow
