//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "utilities.h"

#include <cmath>
#include <iostream>
#include <numeric>

#include "macros.h"

//========================================================================================================
//                                  Utility Functions
//========================================================================================================

/**
 * <!-- ************************************************************************************** -->
 * @brief Prints the elements of an array to standard output.
 * @param arr The array to print
 ******************************************************************************************************************* -->
 */
void print_array(Array const& arr) {
    for (auto const& a : arr) {
        std::cout << a << " ";
    }
    std::cout << std::endl;
}

//========================================================================================================
//                                  Point Interpolation Functions (Foundation)
//========================================================================================================

/**
 * <!-- ************************************************************************************** -->
 * @brief Point-wise linear interpolation between two points.
 * @param x0 First x-coordinate
 * @param x1 Second x-coordinate
 * @param y0 First y-coordinate
 * @param y1 Second y-coordinate
 * @param xi The x-value at which to interpolate
 * @return The interpolated y-value
 * <!-- ************************************************************************************** -->
 */
Real point_interp(Real x0, Real x1, Real y0, Real y1, Real xi) {
    if (x0 == x1)
        return y0;
    const Real slope = (y1 - y0) / (x1 - x0);
    return y0 + slope * (xi - x0);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Point-wise logarithmic interpolation between two points.
 * @details Performs linear interpolation in log-log space.
 * @param x0 First x-coordinate
 * @param x1 Second x-coordinate
 * @param y0 First y-coordinate
 * @param y1 Second y-coordinate
 * @param xi The x-value at which to interpolate
 * @return The interpolated y-value
 * <!-- ************************************************************************************** -->
 */
Real point_loglog_interp(Real x0, Real x1, Real y0, Real y1, Real xi) {
    if (y0 == 0 || y1 == 0)
        return 0;
    if (x0 == x1)
        return y0;
    const Real log_x0 = std::log(x0);
    const Real log_x1 = std::log(x1);
    const Real log_y0 = std::log(y0);
    const Real log_y1 = std::log(y1);
    const Real slope = (log_y1 - log_y0) / (log_x1 - log_x0);
    return std::exp(log_y0 + slope * (std::log(xi) - log_x0));
}

//========================================================================================================
//                                  General Interpolation Functions (Linear)
//========================================================================================================

/**
 * <!-- ************************************************************************************** -->
 * @brief Linear interpolation for arbitrary x-values.
 * @details Finds the appropriate interval for interpolation and uses point_interp.
 * @param x0 The x-value at which to interpolate
 * @param x Array of x-coordinates (must be monotonically increasing)
 * @param y Array of y-coordinates
 * @param lo_extrap Whether to extrapolate below the minimum x-value
 * @param hi_extrap Whether to extrapolate above the maximum x-value
 * @return The interpolated y-value
 * <!-- ************************************************************************************** -->
 */
Real interp(Real x0, Array const& x, Array const& y, bool lo_extrap, bool hi_extrap) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
        std::cout << "incorrect array size for interpolation!\n";
        return 0;
    }
    const auto x_back = x(x.size() - 1);
    const auto y_back = y(y.size() - 1);

    if (x0 < x(0)) {
        return (!lo_extrap || x(0) == x0) ? y(0) : point_interp(x(0), x(1), y(0), y(1), x0);
    } else if (x0 > x_back) {
        return (!hi_extrap || x_back == x0) ? y_back
                                            : point_interp(x(x.size() - 2), x_back, y(y.size() - 2), y_back, x0);
    } else {
        const auto it = std::ranges::lower_bound(x, x0);
        const size_t idx = it - x.begin();
        if (*it == x0)
            return y(idx); // Exact match
        return point_interp(x(idx - 1), x(idx), y(idx - 1), y(idx), x0);
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Linear interpolation for equally spaced x-values.
 * @details Optimized version that avoids search when x-values are equally spaced.
 * @param xi The x-value at which to interpolate
 * @param x Array of x-coordinates (must be equally spaced)
 * @param y Array of y-coordinates
 * @param lo_extrap Whether to extrapolate below the minimum x-value
 * @param hi_extrap Whether to extrapolate above the maximum x-value
 * @return The interpolated y-value
 * <!-- ************************************************************************************** -->
 */
Real eq_space_interp(Real xi, Array const& x, Array const& y, bool lo_extrap, bool hi_extrap) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
        std::cout << "incorrect array size for interpolation!\n";
        return 0;
    }

    const auto x_back = x[x.size() - 1];
    const auto y_back = y[y.size() - 1];

    if (xi <= x[0])
        return (!lo_extrap || x[0] == xi) ? y[0] : point_interp(x[0], x[1], y[0], y[1], xi);
    else if (xi >= x_back)
        return (!hi_extrap || x_back == xi) ? y_back
                                            : point_interp(x[x.size() - 2], x_back, y[y.size() - 2], y_back, xi);
    else {
        const Real dx = x[1] - x[0];
        const size_t idx = static_cast<size_t>((xi - x[0]) / dx + 1);
        if (xi == x[idx])
            return y[idx];
        return point_interp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}

//========================================================================================================
//                                  General Interpolation Functions (Logarithmic)
//========================================================================================================

/**
 * <!-- ************************************************************************************** -->
 * @brief Logarithmic interpolation for arbitrary x-values.
 * @details Performs interpolation in log-log space for arbitrary x-values.
 * @param xi The x-value at which to interpolate
 * @param x Array of x-coordinates (must be monotonically increasing)
 * @param y Array of y-coordinates
 * @param lo_extrap Whether to extrapolate below the minimum x-value
 * @param hi_extrap Whether to extrapolate above the maximum x-value
 * @return The interpolated y-value
 * <!-- ************************************************************************************** -->
 */
Real loglog_interp(Real xi, const Array& x, const Array& y, bool lo_extrap, bool hi_extrap) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
        std::cout << "incorrect array size for interpolation!\n";
        return 0;
    }
    const auto x_back = x[x.size() - 1];
    const auto y_back = y[y.size() - 1];

    if (xi <= x[0]) {
        return (!lo_extrap || x[0] == xi) ? y[0] : point_loglog_interp(x[0], x[1], y[0], y[1], xi);
    } else if (xi >= x_back) {
        return (!hi_extrap || x_back == xi) ? y_back
                                            : point_loglog_interp(x[x.size() - 2], x_back, y[y.size() - 2], y_back, xi);
    } else {
        const auto it = std::ranges::lower_bound(x, xi);
        const size_t idx = it - x.begin();
        if (*it == xi)
            return y[idx]; // Exact match
        return point_loglog_interp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Logarithmic interpolation for equally spaced x-values in log space.
 * @details Optimized version that avoids search when x-values are equally spaced in log space.
 * @param xi The x-value at which to interpolate
 * @param x Array of x-coordinates (must be equally spaced in log space)
 * @param y Array of y-coordinates
 * @param lo_extrap Whether to extrapolate below the minimum x-value
 * @param hi_extrap Whether to extrapolate above the maximum x-value
 * @return The interpolated y-value
 * <!-- ************************************************************************************** -->
 */
Real eq_space_loglog_interp(Real xi, const Array& x, const Array& y, bool lo_extrap, bool hi_extrap) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
        std::cout << "incorrect array size for interpolation!\n";
        return 0;
    }
    const auto x_back = x[x.size() - 1];
    const auto y_back = y[y.size() - 1];

    if (xi <= x[0]) {
        // std::cout << "here!" << (!lo_extrap || x[0] == xi) ? y[0] : point_loglog_interp(x[0], x[1], y[0], y[1], xi);
        return (!lo_extrap || x[0] == xi) ? y[0] : point_loglog_interp(x[0], x[1], y[0], y[1], xi);
    } else if (xi >= x_back) {
        return (!hi_extrap || x_back == xi) ? y_back
                                            : point_loglog_interp(x[x.size() - 2], x_back, y[y.size() - 2], y_back, xi);
    } else {
        const Real log_x0 = std::log(x[0]);
        const Real dx = std::log(x[1]) - log_x0;
        const size_t idx = static_cast<size_t>((std::log(xi) - log_x0) / dx + 1);

        if (xi == x[idx])
            return y[idx]; // Exact match
        return point_loglog_interp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}
