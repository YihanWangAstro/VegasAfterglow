#include <boost/test/unit_test.hpp>
#include <cmath>

#include "core/quadrature.h"
#include "util/macros.h"

BOOST_AUTO_TEST_SUITE(Quadrature)

// ============================================================================
//  compute_bin_widths
// ============================================================================

BOOST_AUTO_TEST_CASE(bin_widths_uniform) {
    Array grid = {1.0, 2.0, 3.0, 4.0, 5.0};
    Array dg;
    compute_bin_widths(grid, dg);
    BOOST_CHECK_EQUAL(dg.size(), 4u);
    for (size_t i = 0; i < dg.size(); ++i) {
        BOOST_CHECK_CLOSE(dg(i), 1.0, 1e-10);
    }
}

BOOST_AUTO_TEST_CASE(bin_widths_nonuniform) {
    Array grid = {1.0, 3.0, 10.0};
    Array dg;
    compute_bin_widths(grid, dg);
    BOOST_CHECK_EQUAL(dg.size(), 2u);
    BOOST_CHECK_CLOSE(dg(0), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(dg(1), 7.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(bin_widths_single_point) {
    Array grid = {5.0};
    Array dg;
    compute_bin_widths(grid, dg);
    BOOST_CHECK_EQUAL(dg.size(), 0u);
}

BOOST_AUTO_TEST_CASE(bin_widths_empty) {
    Array grid = Array::from_shape({0});
    Array dg;
    compute_bin_widths(grid, dg);
    BOOST_CHECK_EQUAL(dg.size(), 0u);
}

// ============================================================================
//  compute_midpoint_grid
// ============================================================================

BOOST_AUTO_TEST_CASE(midpoint_grid_basic) {
    Array boundaries = {1.0, 4.0, 16.0};
    Array centers, weights;
    compute_midpoint_grid(boundaries, centers, weights);
    BOOST_CHECK_EQUAL(centers.size(), 2u);
    BOOST_CHECK_EQUAL(weights.size(), 2u);
    // Centers are geometric mean: sqrt(1*4)=2, sqrt(4*16)=8
    BOOST_CHECK_CLOSE(centers(0), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(centers(1), 8.0, 1e-10);
    // Weights are bin widths: 4-1=3, 16-4=12
    BOOST_CHECK_CLOSE(weights(0), 3.0, 1e-10);
    BOOST_CHECK_CLOSE(weights(1), 12.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(midpoint_grid_single_bin) {
    Array boundaries = {2.0, 8.0};
    Array centers, weights;
    compute_midpoint_grid(boundaries, centers, weights);
    BOOST_CHECK_EQUAL(centers.size(), 1u);
    BOOST_CHECK_CLOSE(centers(0), 4.0, 1e-10);
    BOOST_CHECK_CLOSE(weights(0), 6.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(midpoint_grid_too_few) {
    Array boundaries = {5.0};
    Array centers, weights;
    compute_midpoint_grid(boundaries, centers, weights);
    BOOST_CHECK_EQUAL(centers.size(), 0u);
    BOOST_CHECK_EQUAL(weights.size(), 0u);
}

// ============================================================================
//  compute_rectangle_weights
// ============================================================================

BOOST_AUTO_TEST_CASE(rectangle_weights_uniform) {
    Array grid = {0.0, 1.0, 2.0, 3.0};
    Array weights;
    compute_rectangle_weights(grid, weights);
    BOOST_CHECK_EQUAL(weights.size(), 4u);
    BOOST_CHECK_CLOSE(weights(0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(weights(1), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(weights(2), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(weights(3), 0.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(rectangle_weights_single_point) {
    Array grid = {5.0};
    Array weights;
    compute_rectangle_weights(grid, weights);
    BOOST_CHECK_EQUAL(weights.size(), 1u);
    BOOST_CHECK_CLOSE(weights(0), 0.0, 1e-10);
}

// ============================================================================
//  compute_trapezoidal_weights
// ============================================================================

BOOST_AUTO_TEST_CASE(trapezoidal_weights_uniform) {
    // Uniform grid: endpoints get h/2, interior get h
    Array grid = {0.0, 1.0, 2.0, 3.0, 4.0};
    Array weights;
    compute_trapezoidal_weights(grid, weights);
    BOOST_CHECK_EQUAL(weights.size(), 5u);
    BOOST_CHECK_CLOSE(weights(0), 0.5, 1e-10);
    BOOST_CHECK_CLOSE(weights(1), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(weights(2), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(weights(3), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(weights(4), 0.5, 1e-10);
}

BOOST_AUTO_TEST_CASE(trapezoidal_weights_two_points) {
    Array grid = {0.0, 2.0};
    Array weights;
    compute_trapezoidal_weights(grid, weights);
    BOOST_CHECK_EQUAL(weights.size(), 2u);
    BOOST_CHECK_CLOSE(weights(0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(weights(1), 1.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(trapezoidal_integrates_linear) {
    // Trapezoidal should integrate a linear function f(x)=x exactly
    Array grid = {0.0, 1.0, 2.0, 3.0, 4.0};
    Array weights;
    compute_trapezoidal_weights(grid, weights);
    Real integral = 0;
    for (size_t i = 0; i < grid.size(); ++i) {
        integral += grid(i) * weights(i);
    }
    // Exact integral of x from 0 to 4 = 8
    BOOST_CHECK_CLOSE(integral, 8.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(trapezoidal_single_point) {
    Array grid = {3.0};
    Array weights;
    compute_trapezoidal_weights(grid, weights);
    BOOST_CHECK_EQUAL(weights.size(), 1u);
    BOOST_CHECK_CLOSE(weights(0), 0.0, 1e-10);
}

// ============================================================================
//  compute_simpson_weights
// ============================================================================

BOOST_AUTO_TEST_CASE(simpson_integrates_cubic) {
    // Simpson's rule should exactly integrate up to cubic polynomials
    // f(x) = x^3 on [0, 2], exact integral = 4
    Array grid = {0.0, 1.0, 2.0};
    Array weights;
    compute_simpson_weights(grid, weights);
    Real integral = 0;
    for (size_t i = 0; i < grid.size(); ++i) {
        Real x = grid(i);
        integral += x * x * x * weights(i);
    }
    BOOST_CHECK_CLOSE(integral, 4.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(simpson_weights_positive) {
    Array grid = {0.0, 1.0, 2.0, 3.0, 4.0};
    Array weights;
    compute_simpson_weights(grid, weights);
    for (size_t i = 0; i < weights.size(); ++i) {
        BOOST_CHECK_GE(weights(i), 0.0);
    }
}

BOOST_AUTO_TEST_CASE(simpson_odd_intervals_fallback) {
    // 4 points = 3 intervals (odd): 2 intervals Simpson + 1 trapezoidal
    Array grid = {0.0, 1.0, 2.0, 3.0};
    Array weights;
    compute_simpson_weights(grid, weights);
    // Integrate f(x)=1 -> should give 3
    Real integral = 0;
    for (size_t i = 0; i < weights.size(); ++i) {
        integral += weights(i);
    }
    BOOST_CHECK_CLOSE(integral, 3.0, 1e-10);
}

// ============================================================================
//  compute_boole_weights
// ============================================================================

BOOST_AUTO_TEST_CASE(boole_weights_positive) {
    // Boole's rule weights should all be non-negative
    const size_t n = 17;
    Array grid = Array::from_shape({n});
    // Log-spaced grid (uniform in log)
    for (size_t i = 0; i < n; ++i) {
        grid(i) = std::exp(0.1 * static_cast<Real>(i) + 1.0);
    }
    Array weights;
    compute_boole_weights(grid, weights);
    BOOST_CHECK_EQUAL(weights.size(), n);
    for (size_t i = 0; i < n; ++i) {
        BOOST_CHECK_GE(weights(i), 0.0);
    }
}

BOOST_AUTO_TEST_CASE(boole_integrates_power_law) {
    // Integrate f(gamma) = gamma^2 over log-spaced grid
    // Analytic: integral gamma^2 dgamma from a to b = (b^3 - a^3)/3
    const size_t n = 21;
    const Real log_a = 0.0;
    const Real log_b = 2.0;
    const Real h = (log_b - log_a) / (n - 1);
    Array grid = Array::from_shape({n});
    for (size_t i = 0; i < n; ++i) {
        grid(i) = std::exp(log_a + h * i);
    }
    Array weights;
    compute_boole_weights(grid, weights);
    Real integral = 0;
    for (size_t i = 0; i < n; ++i) {
        integral += grid(i) * grid(i) * weights(i);
    }
    const Real a = std::exp(log_a);
    const Real b = std::exp(log_b);
    const Real exact = (b * b * b - a * a * a) / 3.0;
    BOOST_CHECK_CLOSE(integral, exact, 1.0); // Within 1%
}

// ============================================================================
//  compute_nc7_weights
// ============================================================================

BOOST_AUTO_TEST_CASE(nc7_weights_positive) {
    // NC7 weights should all be non-negative (guaranteed by property)
    const size_t n = 25;
    Array grid = Array::from_shape({n});
    for (size_t i = 0; i < n; ++i) {
        grid(i) = std::exp(0.1 * static_cast<Real>(i) + 1.0);
    }
    Array weights;
    compute_nc7_weights(grid, weights);
    BOOST_CHECK_EQUAL(weights.size(), n);
    for (size_t i = 0; i < n; ++i) {
        BOOST_CHECK_GE(weights(i), 0.0);
    }
}

BOOST_AUTO_TEST_CASE(nc7_integrates_power_law) {
    // Higher-order NC7 should be more accurate than Boole for smooth functions
    const size_t n = 25;
    const Real log_a = 0.0;
    const Real log_b = 2.0;
    const Real h = (log_b - log_a) / (n - 1);
    Array grid = Array::from_shape({n});
    for (size_t i = 0; i < n; ++i) {
        grid(i) = std::exp(log_a + h * i);
    }
    Array weights;
    compute_nc7_weights(grid, weights);
    Real integral = 0;
    for (size_t i = 0; i < n; ++i) {
        integral += grid(i) * grid(i) * weights(i);
    }
    const Real a = std::exp(log_a);
    const Real b = std::exp(log_b);
    const Real exact = (b * b * b - a * a * a) / 3.0;
    BOOST_CHECK_CLOSE(integral, exact, 0.1); // Within 0.1%
}

BOOST_AUTO_TEST_CASE(nc7_leftover_handling) {
    // Test with various sizes that exercise different leftover branches
    for (size_t n : {3u, 4u, 5u, 6u, 7u, 8u, 13u}) {
        Array grid = Array::from_shape({n});
        for (size_t i = 0; i < n; ++i) {
            grid(i) = std::exp(0.2 * static_cast<Real>(i) + 1.0);
        }
        Array weights;
        compute_nc7_weights(grid, weights);
        BOOST_CHECK_EQUAL(weights.size(), n);
        // All weights should be non-negative
        for (size_t i = 0; i < n; ++i) {
            BOOST_CHECK_GE(weights(i), 0.0);
        }
        // Sum of weights * 1 (constant func) should approximate the interval
        Real sum = 0;
        for (size_t i = 0; i < n; ++i) {
            sum += weights(i);
        }
        const Real exact = grid(n - 1) - grid(0);
        BOOST_CHECK_CLOSE(sum, exact, 5.0); // Within 5%
    }
}

// ============================================================================
//  Weight sum consistency: all methods should agree on integral of f=1
// ============================================================================

BOOST_AUTO_TEST_CASE(weight_methods_agree_on_constant) {
    // All quadrature methods integrating f(x)=1 should give the same result
    Array grid = {1.0, 2.0, 3.0, 4.0, 5.0};
    const Real exact = 4.0; // 5 - 1

    Array w_trap, w_simp;
    compute_trapezoidal_weights(grid, w_trap);
    compute_simpson_weights(grid, w_simp);

    Real sum_trap = 0, sum_simp = 0;
    for (size_t i = 0; i < grid.size(); ++i) {
        sum_trap += w_trap(i);
        sum_simp += w_simp(i);
    }
    BOOST_CHECK_CLOSE(sum_trap, exact, 1e-10);
    BOOST_CHECK_CLOSE(sum_simp, exact, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
