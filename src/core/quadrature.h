//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include "mesh.h"

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes interval widths dg(j) = grid(j+1) - grid(j) from a grid of boundaries.
 * <!-- ************************************************************************************** -->
 */
inline void compute_bin_widths(Array const& grid, Array& dg) {
    const size_t m = grid.size();
    if (m >= 2) {
        dg = Array::from_shape({m - 1});
        for (size_t j = 0; j + 1 < m; ++j)
            dg(j) = grid(j + 1) - grid(j);
    } else {
        dg = Array::from_shape({0});
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts a boundary grid to a midpoint grid with bin widths (midpoint rule, O(h^2)).
 *
 * Given N boundary points, produces N-1 center points (geometric mean of adjacent boundaries)
 * and N-1 bin widths. Works with any grid (uniform, adaptive, logspace).
 * <!-- ************************************************************************************** -->
 */
inline void compute_midpoint_grid(Array const& boundaries, Array& centers, Array& weights) {
    const size_t n = boundaries.size();
    if (n < 2) {
        centers = Array::from_shape({0});
        weights = Array::from_shape({0});
        return;
    }
    const size_t m = n - 1;
    centers = Array::from_shape({m});
    weights = Array::from_shape({m});
    for (size_t j = 0; j < m; ++j) {
        centers(j) = std::sqrt(boundaries(j) * boundaries(j + 1));
        weights(j) = boundaries(j + 1) - boundaries(j);
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes left-rectangle quadrature weights (O(h)).
 *
 * w_j = grid_{j+1} - grid_j for j=0..N-2, w_{N-1} = 0.
 * <!-- ************************************************************************************** -->
 */
inline void compute_rectangle_weights(Array const& grid, Array& weights) {
    const size_t n = grid.size();
    weights = Array::from_shape({n});
    if (n < 2) {
        if (n == 1)
            weights(0) = 0;
        return;
    }
    for (size_t j = 0; j + 1 < n; ++j) {
        weights(j) = grid(j + 1) - grid(j);
    }
    weights(n - 1) = 0;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes trapezoidal quadrature weights for a non-uniform grid.
 *
 * w_0 = 0.5*(grid_1 - grid_0), w_j = 0.5*(grid_{j+1} - grid_{j-1}), w_N = 0.5*(grid_N - grid_{N-1})
 * <!-- ************************************************************************************** -->
 */
inline void compute_trapezoidal_weights(Array const& grid, Array& weights) {
    const size_t n = grid.size();
    weights = Array::from_shape({n});
    if (n < 2) {
        if (n == 1) {
            weights(0) = 0;
        }
        return;
    }
    weights(0) = 0.5 * (grid(1) - grid(0));
    for (size_t j = 1; j + 1 < n; ++j) {
        weights(j) = 0.5 * (grid(j + 1) - grid(j - 1));
    }
    weights(n - 1) = 0.5 * (grid(n - 1) - grid(n - 2));
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes composite Simpson's quadrature weights for a non-uniform grid.
 *
 * Pairs of intervals use non-uniform Simpson's 1/3 rule. If the number of
 * intervals is odd, the last interval falls back to trapezoidal.
 * <!-- ************************************************************************************** -->
 */
inline void compute_simpson_weights(Array const& grid, Array& weights) {
    const size_t n = grid.size();
    weights = xt::zeros<Real>({n});
    if (n < 2)
        return;

    // Process pairs of intervals with non-uniform Simpson's rule
    size_t j = 0;
    for (; j + 2 < n; j += 2) {
        const Real h1 = grid(j + 1) - grid(j);
        const Real h2 = grid(j + 2) - grid(j + 1);
        const Real H = h1 + h2;
        weights(j) += H / 6.0 * (2.0 * h1 - h2) / h1;
        weights(j + 1) += H * H * H / (6.0 * h1 * h2);
        weights(j + 2) += H / 6.0 * (2.0 * h2 - h1) / h2;
    }

    // Leftover single interval: trapezoidal fallback
    if (j + 1 < n) {
        const Real h = grid(j + 1) - grid(j);
        weights(j) += 0.5 * h;
        weights(j + 1) += 0.5 * h;
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes Boole's rule (Newton-Cotes 4th order, O(h^6)) weights for a uniform log-spaced grid.
 *
 * Operates in log-space where the grid is uniform, then applies the Jacobian (gamma_i)
 * to convert back to linear-space weights. Groups of 4 intervals use Boole's coefficients
 * [7, 32, 12, 32, 7] * 2h/45. Leftovers use Simpson's 3/8, 1/3, or trapezoidal.
 * <!-- ************************************************************************************** -->
 */
inline void compute_boole_weights(Array const& grid, Array& weights) {
    const size_t n = grid.size();
    weights = xt::zeros<Real>({n});
    if (n < 2)
        return;

    // Uniform spacing in log space
    const Real h = std::log(grid(1) / grid(0));

    // Process groups of 4 intervals (5 points) with Boole's rule
    const Real cb = 2.0 * h / 45.0;
    size_t j = 0;
    for (; j + 4 < n; j += 4) {
        weights(j) += cb * 7;
        weights(j + 1) += cb * 32;
        weights(j + 2) += cb * 12;
        weights(j + 3) += cb * 32;
        weights(j + 4) += cb * 7;
    }

    // Leftover intervals
    const size_t remaining = n - 1 - j;
    if (remaining == 3) {
        // Simpson's 3/8 rule: (3h/8) * [1, 3, 3, 1]
        const Real c38 = 3.0 * h / 8.0;
        weights(j) += c38;
        weights(j + 1) += c38 * 3;
        weights(j + 2) += c38 * 3;
        weights(j + 3) += c38;
    } else if (remaining == 2) {
        // Simpson's 1/3 rule: (h/3) * [1, 4, 1]
        const Real c13 = h / 3.0;
        weights(j) += c13;
        weights(j + 1) += c13 * 4;
        weights(j + 2) += c13;
    } else if (remaining == 1) {
        weights(j) += 0.5 * h;
        weights(j + 1) += 0.5 * h;
    }

    // Jacobian: dγ = γ · d(ln γ), convert log-space weights to linear-space
    for (size_t i = 0; i < n; ++i) {
        weights(i) *= grid(i);
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes 7-point Newton-Cotes (O(h^8)) weights for a uniform log-spaced grid.
 *
 * Groups of 6 intervals use coefficients [41, 216, 27, 272, 27, 216, 41] * h/140.
 * Leftovers cascade: 5→6-pt, 4→Boole's, 3→Simpson 3/8, 2→Simpson 1/3, 1→trapezoidal.
 * All weights are positive (highest NC order with this property is 8-point).
 * <!-- ************************************************************************************** -->
 */
inline void compute_nc7_weights(Array const& grid, Array& weights) {
    const size_t n = grid.size();
    weights = xt::zeros<Real>({n});
    if (n < 2)
        return;

    const Real h = std::log(grid(1) / grid(0));

    // Groups of 6 intervals (7 points): h/140 * [41, 216, 27, 272, 27, 216, 41]
    const Real c7 = h / 140.0;
    size_t j = 0;
    for (; j + 6 < n; j += 6) {
        weights(j) += c7 * 41;
        weights(j + 1) += c7 * 216;
        weights(j + 2) += c7 * 27;
        weights(j + 3) += c7 * 272;
        weights(j + 4) += c7 * 27;
        weights(j + 5) += c7 * 216;
        weights(j + 6) += c7 * 41;
    }

    const size_t remaining = n - 1 - j;
    if (remaining == 5) {
        // 6-point: 5h/288 * [19, 75, 50, 50, 75, 19]
        const Real c6 = 5.0 * h / 288.0;
        weights(j) += c6 * 19;
        weights(j + 1) += c6 * 75;
        weights(j + 2) += c6 * 50;
        weights(j + 3) += c6 * 50;
        weights(j + 4) += c6 * 75;
        weights(j + 5) += c6 * 19;
    } else if (remaining == 4) {
        const Real cb = 2.0 * h / 45.0;
        weights(j) += cb * 7;
        weights(j + 1) += cb * 32;
        weights(j + 2) += cb * 12;
        weights(j + 3) += cb * 32;
        weights(j + 4) += cb * 7;
    } else if (remaining == 3) {
        const Real c38 = 3.0 * h / 8.0;
        weights(j) += c38;
        weights(j + 1) += c38 * 3;
        weights(j + 2) += c38 * 3;
        weights(j + 3) += c38;
    } else if (remaining == 2) {
        const Real c13 = h / 3.0;
        weights(j) += c13;
        weights(j + 1) += c13 * 4;
        weights(j + 2) += c13;
    } else if (remaining == 1) {
        weights(j) += 0.5 * h;
        weights(j + 1) += 0.5 * h;
    }

    for (size_t i = 0; i < n; ++i) {
        weights(i) *= grid(i);
    }
}
