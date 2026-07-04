#include <boost/test/unit_test.hpp>
#include <cmath>
#include <limits>
#include <numbers>

#include "util/fast-math.h"

BOOST_AUTO_TEST_SUITE(FastMath)

// log2_softplus tests

// At x = 0, log2(1 + 2^0) = log2(2) = 1 exactly
BOOST_AUTO_TEST_CASE(log2_softplus_zero) {
    // log2(1 + 2^0) = log2(2) = 1
    BOOST_CHECK_CLOSE(log2_softplus(0.0), 1.0, 1e-10);
}

// Large-x asymptote: log2(1 + 2^x) -> x, checked at x = 50 and 100
BOOST_AUTO_TEST_CASE(log2_softplus_large_positive) {
    // For large x, log2(1 + 2^x) ~ x
    BOOST_CHECK_CLOSE(log2_softplus(50.0), 50.0, 1e-10);
    BOOST_CHECK_CLOSE(log2_softplus(100.0), 100.0, 1e-10);
}

// Large-negative asymptote: log2(1 + 2^x) -> 0, checked at x = -50 and -100
BOOST_AUTO_TEST_CASE(log2_softplus_large_negative) {
    // For large negative x, log2(1 + 2^x) ~ 0
    BOOST_CHECK_SMALL(log2_softplus(-50.0), 1e-10);
    BOOST_CHECK_SMALL(log2_softplus(-100.0), 1e-10);
}

// log2_softplus is non-decreasing over x in [-30, 30], sampled in unit steps
BOOST_AUTO_TEST_CASE(log2_softplus_monotonic) {
    Real prev = log2_softplus(-30.0);
    for (Real x = -29.0; x <= 30.0; x += 1.0) {
        Real val = log2_softplus(x);
        BOOST_CHECK_GE(val, prev);
        prev = val;
    }
}

// Strictly increasing across the x = 20 branch point: value(19) < value(20) < value(21)
BOOST_AUTO_TEST_CASE(log2_softplus_boundary_20) {
    // Branch points at ±20 — values just inside and outside should be close
    Real at_19 = log2_softplus(19.0);
    Real at_20 = log2_softplus(20.0);
    Real at_21 = log2_softplus(21.0);
    // Monotonic across the branch
    BOOST_CHECK_LT(at_19, at_20);
    BOOST_CHECK_LT(at_20, at_21);
}

// Crossing the x = 20 branch changes the value by < 1%; near x = -20 the result stays finite and non-negative
BOOST_AUTO_TEST_CASE(log2_softplus_continuity_at_branch) {
    // No discontinuity at x=20 and x=-20
    const Real eps = 0.01;
    Real below_20 = log2_softplus(20.0 - eps);
    Real above_20 = log2_softplus(20.0 + eps);
    BOOST_CHECK_CLOSE(below_20, above_20, 1.0); // within 1%

    // Near x=-20 the function transitions from exact to asymptotic;
    // both values are extremely small — just check they're non-negative and finite
    Real below_neg20 = log2_softplus(-20.0 - eps);
    Real above_neg20 = log2_softplus(-20.0 + eps);
    BOOST_CHECK_GE(below_neg20, 0.0);
    BOOST_CHECK_GE(above_neg20, 0.0);
    BOOST_CHECK(std::isfinite(below_neg20));
    BOOST_CHECK(std::isfinite(above_neg20));
}

// log2_broken_power_ratio tests

// Far below the break the slope correction vanishes (|result| < 0.01), preserving the original power law
BOOST_AUTO_TEST_CASE(log2_broken_power_ratio_far_below_break) {
    // Far below break, correction should be ~0
    Real result = log2_broken_power_ratio(-100.0, 0.0, 2.0, 1.0);
    BOOST_CHECK_SMALL(std::fabs(result), 0.01);
}

// Far above the break the correction is negative and grows more negative with distance (steeper slope)
BOOST_AUTO_TEST_CASE(log2_broken_power_ratio_far_above_break) {
    // Far above break, correction is proportional to distance
    Real r1 = log2_broken_power_ratio(50.0, 0.0, 2.0, 1.0);
    Real r2 = log2_broken_power_ratio(100.0, 0.0, 2.0, 1.0);
    // More negative further above
    BOOST_CHECK_LT(r2, r1);
    BOOST_CHECK_LT(r1, 0);
}

// Exactly at the break the correction is finite and negative (never positive)
BOOST_AUTO_TEST_CASE(log2_broken_power_ratio_at_break) {
    // At the break point, check it's a reasonable intermediate value
    Real result = log2_broken_power_ratio(0.0, 0.0, 2.0, 1.0);
    BOOST_CHECK(std::isfinite(result));
    BOOST_CHECK_LT(result, 0); // correction is always non-positive
}

// The correction is monotonically non-increasing in log2_x over [-20, 40], sampled in unit steps
BOOST_AUTO_TEST_CASE(log2_broken_power_ratio_monotonic) {
    // Correction is monotonically decreasing with increasing x
    const Real log2_break = 10.0;
    const Real s_delta = 3.0;
    const Real s = 2.0;
    Real prev = log2_broken_power_ratio(-20.0, log2_break, s_delta, s);
    for (Real x = -19.0; x <= 40.0; x += 1.0) {
        Real val = log2_broken_power_ratio(x, log2_break, s_delta, s);
        BOOST_CHECK_LE(val, prev + 1e-12);
        prev = val;
    }
}

// With sharpness s = 100 the transition is step-like: correction ~0 one unit below the break, clearly negative one unit above
BOOST_AUTO_TEST_CASE(log2_broken_power_ratio_extreme_sharpness) {
    // Very sharp transition (s=100) approaches step function
    Real below = log2_broken_power_ratio(-1.0, 0.0, 2.0, 100.0);
    Real above = log2_broken_power_ratio(1.0, 0.0, 2.0, 100.0);
    // Below break: correction near 0
    BOOST_CHECK_SMALL(std::fabs(below), 0.1);
    // Above break: correction significant
    BOOST_CHECK_LT(above, -0.01);
}

#ifndef AFTERGLOW_FAST_MATH

// Without AFTERGLOW_FAST_MATH the fast_* functions ARE their std counterparts —
// pin that contract bit-exactly.
BOOST_AUTO_TEST_CASE(fast_functions_match_std) {
    const double values[] = {0.001, 0.1, 1.0, 2.5, 10.0, 100.0, 1e6};
    for (double x : values) {
        BOOST_CHECK_EQUAL(fast_log2(x), std::log2(x));
        BOOST_CHECK_EQUAL(fast_log(x), std::log(x));
    }
    const double exp_values[] = {-10.0, -1.0, 0.0, 1.0, 5.0, 10.0};
    for (double x : exp_values) {
        BOOST_CHECK_EQUAL(fast_exp2(x), std::exp2(x));
        BOOST_CHECK_EQUAL(fast_exp(x), std::exp(x));
    }
    BOOST_CHECK_EQUAL(fast_pow(2.0, 3.0), std::exp2(3.0 * std::log2(2.0)));
    BOOST_CHECK_EQUAL(fast_log2(1e-300), std::log2(1e-300));
    BOOST_CHECK_EQUAL(fast_exp2(1024.0), std::exp2(1024.0));
    BOOST_CHECK_EQUAL(fast_exp2(-1100.0), std::exp2(-1100.0));
}

#else

// With AFTERGLOW_FAST_MATH the kernels are minimax polynomials; pin their
// accuracy bounds on dense sweeps of the full input range the code can produce.

// fast_log2 abs error <= 6e-8 (fitted bound 5.2e-8 + margin) over the full
// positive-normal range: every mantissa region x every magnitude decade.
BOOST_AUTO_TEST_CASE(fast_log2_accuracy_sweep) {
    double max_err = 0.0;
    for (int e = -1020; e <= 1020; e += 3) {
        const double base = std::exp2(static_cast<double>(e));
        for (int m = 0; m < 64; ++m) {
            const double x = base * (1.0 + m / 64.0);
            max_err = std::max(max_err, std::fabs(fast_log2(x) - std::log2(x)));
        }
    }
    BOOST_CHECK_LT(max_err, 6e-8);
}

// fast_exp2 rel error <= 1.3e-7 (fitted bound 1.13e-7 + margin) over the full
// non-overflowing exponent range, and EXACT powers of two at integer x.
BOOST_AUTO_TEST_CASE(fast_exp2_accuracy_sweep) {
    double max_rel = 0.0;
    for (int i = 0; i < 200000; ++i) {
        const double x = -1021.0 + i * (2043.0 / 199999.0);
        const double ref = std::exp2(x);
        max_rel = std::max(max_rel, std::fabs(fast_exp2(x) - ref) / ref);
    }
    BOOST_CHECK_LT(max_rel, 1.3e-7);
    for (int n = -1021; n <= 1023; n += 7) {
        BOOST_CHECK_EQUAL(fast_exp2(static_cast<double>(n)), std::exp2(static_cast<double>(n)));
    }
}

// fast_pow / fast_exp / fast_log inherit the kernel bounds: rel error <= 1e-6
// for exponents |b| <= 5 over 12 decades of base.
BOOST_AUTO_TEST_CASE(fast_pow_accuracy_sweep) {
    double max_rel = 0.0;
    for (int ia = -6; ia <= 6; ++ia) {
        const double a = std::pow(10.0, ia) * 3.7;
        for (int ib = -10; ib <= 10; ++ib) {
            const double b = ib / 2.0;
            const double ref = std::pow(a, b);
            max_rel = std::max(max_rel, std::fabs(fast_pow(a, b) - ref) / ref);
        }
    }
    BOOST_CHECK_LT(max_rel, 1e-6);
}

// Domain semantics match std where it matters: non-positive, subnormal, inf and
// NaN inputs delegate to libm; exp2 clamps instead of overflowing (by design).
BOOST_AUTO_TEST_CASE(fast_math_domain_semantics) {
    BOOST_CHECK_EQUAL(fast_log2(0.0), -std::numeric_limits<double>::infinity());
    BOOST_CHECK(std::isnan(fast_log2(-1.0)));
    BOOST_CHECK(std::isnan(fast_log2(std::numeric_limits<double>::quiet_NaN())));
    BOOST_CHECK_EQUAL(fast_log2(std::numeric_limits<double>::infinity()), std::numeric_limits<double>::infinity());
    BOOST_CHECK_EQUAL(fast_log2(1e-310), std::log2(1e-310)); // subnormal delegates

    BOOST_CHECK(std::isnan(fast_exp2(std::numeric_limits<double>::quiet_NaN())));
    // Never-inf clamp (documented): huge exponents saturate at 2^1023, deep
    // underflow lands at ~0 instead of exactly 0.
    BOOST_CHECK_EQUAL(fast_exp2(1024.0), std::exp2(1023.0));
    BOOST_CHECK_EQUAL(fast_exp2(std::numeric_limits<double>::infinity()), std::exp2(1023.0));
    BOOST_CHECK_LT(fast_exp2(-1100.0), 1e-307);
    BOOST_CHECK_GE(fast_exp2(-1100.0), 0.0);
}

#endif

BOOST_AUTO_TEST_SUITE_END()
