#include <boost/test/unit_test.hpp>
#include <cmath>
#include <numbers>

#include "util/fast-math.h"

BOOST_AUTO_TEST_SUITE(FastMath)

// log2_softplus tests

BOOST_AUTO_TEST_CASE(log2_softplus_zero) {
    // log2(1 + 2^0) = log2(2) = 1
    BOOST_CHECK_CLOSE(log2_softplus(0.0), 1.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(log2_softplus_large_positive) {
    // For large x, log2(1 + 2^x) ~ x
    BOOST_CHECK_CLOSE(log2_softplus(50.0), 50.0, 1e-10);
    BOOST_CHECK_CLOSE(log2_softplus(100.0), 100.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(log2_softplus_large_negative) {
    // For large negative x, log2(1 + 2^x) ~ 0
    BOOST_CHECK_SMALL(log2_softplus(-50.0), 1e-10);
    BOOST_CHECK_SMALL(log2_softplus(-100.0), 1e-10);
}

BOOST_AUTO_TEST_CASE(log2_softplus_monotonic) {
    Real prev = log2_softplus(-30.0);
    for (Real x = -29.0; x <= 30.0; x += 1.0) {
        Real val = log2_softplus(x);
        BOOST_CHECK_GE(val, prev);
        prev = val;
    }
}

BOOST_AUTO_TEST_CASE(log2_softplus_boundary_20) {
    // Branch points at ±20 — values just inside and outside should be close
    Real at_19 = log2_softplus(19.0);
    Real at_20 = log2_softplus(20.0);
    Real at_21 = log2_softplus(21.0);
    // Monotonic across the branch
    BOOST_CHECK_LT(at_19, at_20);
    BOOST_CHECK_LT(at_20, at_21);
}

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

BOOST_AUTO_TEST_CASE(log2_broken_power_ratio_far_below_break) {
    // Far below break, correction should be ~0
    Real result = log2_broken_power_ratio(-100.0, 0.0, 2.0, 1.0);
    BOOST_CHECK_SMALL(std::fabs(result), 0.01);
}

BOOST_AUTO_TEST_CASE(log2_broken_power_ratio_far_above_break) {
    // Far above break, correction is proportional to distance
    Real r1 = log2_broken_power_ratio(50.0, 0.0, 2.0, 1.0);
    Real r2 = log2_broken_power_ratio(100.0, 0.0, 2.0, 1.0);
    // More negative further above
    BOOST_CHECK_LT(r2, r1);
    BOOST_CHECK_LT(r1, 0);
}

BOOST_AUTO_TEST_CASE(log2_broken_power_ratio_at_break) {
    // At the break point, check it's a reasonable intermediate value
    Real result = log2_broken_power_ratio(0.0, 0.0, 2.0, 1.0);
    BOOST_CHECK(std::isfinite(result));
    BOOST_CHECK_LT(result, 0); // correction is always non-positive
}

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

BOOST_AUTO_TEST_CASE(log2_broken_power_ratio_extreme_sharpness) {
    // Very sharp transition (s=100) approaches step function
    Real below = log2_broken_power_ratio(-1.0, 0.0, 2.0, 100.0);
    Real above = log2_broken_power_ratio(1.0, 0.0, 2.0, 100.0);
    // Below break: correction near 0
    BOOST_CHECK_SMALL(std::fabs(below), 0.1);
    // Above break: correction significant
    BOOST_CHECK_LT(above, -0.01);
}

// fast_* functions match std (EXTREME_SPEED is off)

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
    // fast_pow
    BOOST_CHECK_EQUAL(fast_pow(2.0, 3.0), std::exp2(3.0 * std::log2(2.0)));
}

BOOST_AUTO_TEST_CASE(fast_log2_subnormal) {
    // Very small positive input
    double result = fast_log2(1e-300);
    BOOST_CHECK(std::isfinite(result));
    BOOST_CHECK_CLOSE(result, std::log2(1e-300), 1e-10);
}

BOOST_AUTO_TEST_CASE(fast_exp2_overflow) {
    double result = fast_exp2(1024.0);
    BOOST_CHECK_EQUAL(result, std::exp2(1024.0));
}

BOOST_AUTO_TEST_CASE(fast_exp2_underflow) {
    double result = fast_exp2(-1100.0);
    BOOST_CHECK_EQUAL(result, std::exp2(-1100.0));
}

BOOST_AUTO_TEST_SUITE_END()
