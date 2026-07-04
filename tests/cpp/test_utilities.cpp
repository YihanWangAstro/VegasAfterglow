#include <boost/test/unit_test.hpp>
#include <cmath>

#include "util/utilities.h"

BOOST_AUTO_TEST_SUITE(Utilities)

// Power functions

// pow52(a) equals a^(5/2) for integer and fractional bases, matching std::pow
BOOST_AUTO_TEST_CASE(pow52_values) {
    BOOST_CHECK_CLOSE(pow52(4.0), 32.0, 1e-10);
    BOOST_CHECK_CLOSE(pow52(1.0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(pow52(9.0), std::pow(9.0, 2.5), 1e-10);
    BOOST_CHECK_CLOSE(pow52(0.25), std::pow(0.25, 2.5), 1e-10);
}

// pow32(a) equals a^(3/2): pow32(4) = 8, pow32(9) = 27
BOOST_AUTO_TEST_CASE(pow32_values) {
    BOOST_CHECK_CLOSE(pow32(4.0), 8.0, 1e-10);
    BOOST_CHECK_CLOSE(pow32(1.0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(pow32(9.0), 27.0, 1e-10);
}

// pow43(a) equals a^(4/3): pow43(8) = 16, and matches std::pow for other bases
BOOST_AUTO_TEST_CASE(pow43_values) {
    BOOST_CHECK_CLOSE(pow43(8.0), 16.0, 1e-10);
    BOOST_CHECK_CLOSE(pow43(1.0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(pow43(27.0), std::pow(27.0, 4.0 / 3.0), 1e-10);
}

// pow23(a) equals a^(2/3): pow23(8) = 4, pow23(27) = 9
BOOST_AUTO_TEST_CASE(pow23_values) {
    BOOST_CHECK_CLOSE(pow23(8.0), 4.0, 1e-10);
    BOOST_CHECK_CLOSE(pow23(1.0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(pow23(27.0), 9.0, 1e-10);
}

// All fractional-power helpers return 0 at a = 0
BOOST_AUTO_TEST_CASE(pow_functions_zero) {
    BOOST_CHECK_SMALL(pow52(0.0), 1e-15);
    BOOST_CHECK_SMALL(pow32(0.0), 1e-15);
    BOOST_CHECK_SMALL(pow43(0.0), 1e-15);
    BOOST_CHECK_SMALL(pow23(0.0), 1e-15);
}

// All fractional-power helpers return 1 at a = 1
BOOST_AUTO_TEST_CASE(pow_functions_one) {
    BOOST_CHECK_CLOSE(pow52(1.0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(pow32(1.0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(pow43(1.0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(pow23(1.0), 1.0, 1e-10);
}

// Fractional-power helpers stay finite (no overflow) for a large input a = 1e10
// despite forming a^5/a^4 internally before the root
BOOST_AUTO_TEST_CASE(pow_functions_large) {
    Real x = 1e10;
    BOOST_CHECK(std::isfinite(pow52(x)));
    BOOST_CHECK(std::isfinite(pow32(x)));
    BOOST_CHECK(std::isfinite(pow43(x)));
    BOOST_CHECK(std::isfinite(pow23(x)));
}

// Step function

// stepFunc is a Heaviside step returning 1 only for strictly positive x, with stepFunc(0) = 0
BOOST_AUTO_TEST_CASE(step_func) {
    BOOST_CHECK_EQUAL(stepFunc(-1.0), 0.0);
    BOOST_CHECK_EQUAL(stepFunc(0.0), 0.0);
    BOOST_CHECK_EQUAL(stepFunc(1e-15), 1.0);
    BOOST_CHECK_EQUAL(stepFunc(100.0), 1.0);
}

// stepFunc keeps the strict x > 0 threshold even for tiny magnitudes (+/-1e-300)
BOOST_AUTO_TEST_CASE(step_func_near_zero) {
    BOOST_CHECK_EQUAL(stepFunc(1e-300), 1.0);
    BOOST_CHECK_EQUAL(stepFunc(-1e-300), 0.0);
}

// Unit conversion

// eVtoHz implements nu = E/h: pins 1 eV to the physical value 2.418e14 Hz within 1%
BOOST_AUTO_TEST_CASE(eV_to_Hz) {
    Real hz = eVtoHz(1.0 * unit::eV);
    // Result is in code units; convert to physical Hz for comparison
    BOOST_CHECK_CLOSE(hz / unit::Hz, 2.418e14, 1.0); // within 1%
}

// eVtoHz maps zero energy to exactly zero frequency
BOOST_AUTO_TEST_CASE(eV_to_Hz_zero) {
    BOOST_CHECK_EQUAL(eVtoHz(0.0), 0.0);
}

// Root finding

// Bisection finds the root x = 2 of x^2 - 4 on the bracket [0, 4]
BOOST_AUTO_TEST_CASE(root_bisect_quadratic) {
    auto f = [](Real x) { return x * x - 4.0; };
    Real root = root_bisect(f, 0.0, 4.0);
    BOOST_CHECK_CLOSE(root, 2.0, 1e-4);
}

// Bisection finds the root x = 3 of a linear function on [0, 10]
BOOST_AUTO_TEST_CASE(root_bisect_linear) {
    auto f = [](Real x) { return x - 3.0; };
    Real root = root_bisect(f, 0.0, 10.0);
    BOOST_CHECK_CLOSE(root, 3.0, 1e-4);
}

// Bisection converges to pi for sin(x) bracketed on [3, 4]
BOOST_AUTO_TEST_CASE(root_bisect_trig) {
    auto f = [](Real x) { return std::sin(x); };
    Real root = root_bisect(f, 3.0, 4.0);
    BOOST_CHECK_CLOSE(root, std::numbers::pi, 1e-4);
}

// Bisection still converges (to |root| < 1e-3) when the root sits exactly at the lower bracket endpoint
BOOST_AUTO_TEST_CASE(root_bisect_root_at_boundary) {
    auto f = [](Real x) { return x; };
    Real root = root_bisect(f, 0.0, 10.0);
    BOOST_CHECK_SMALL(root, 1e-3);
}

// Bisection returns the root x = 5 when the initial bracket is already narrower than its tolerance
BOOST_AUTO_TEST_CASE(root_bisect_narrow_interval) {
    auto f = [](Real x) { return x - 5.0; };
    Real root = root_bisect(f, 4.9999999, 5.0000001);
    BOOST_CHECK_CLOSE(root, 5.0, 1e-4);
}

// Variadic min/max

// Variadic min returns a single argument unchanged, the smallest of five, and handles equal arguments
BOOST_AUTO_TEST_CASE(variadic_min) {
    BOOST_CHECK_EQUAL(min(42.0), 42.0);
    BOOST_CHECK_EQUAL(min(5.0, 3.0, 7.0, 1.0, 9.0), 1.0);
    BOOST_CHECK_EQUAL(min(1.0, 1.0), 1.0);
}

// Variadic max returns a single argument unchanged, the largest of five, and handles equal arguments
BOOST_AUTO_TEST_CASE(variadic_max) {
    BOOST_CHECK_EQUAL(max(42.0), 42.0);
    BOOST_CHECK_EQUAL(max(5.0, 3.0, 7.0, 1.0, 9.0), 9.0);
    BOOST_CHECK_EQUAL(max(1.0, 1.0), 1.0);
}

// BrokenPowerLaw
//
// BOOST_CHECK_CLOSE tolerances are in PERCENT. The log2-space evaluation path
// goes through fast_log2/fast_exp2: exact under the default build, minimax
// polynomials (rel err ~1e-7) under AFTERGLOW_FAST_MATH.
#ifdef AFTERGLOW_FAST_MATH
constexpr Real BPL_TOL = 1e-4; // percent; kernel bound 1e-7 rel = 1e-5 % + margin
#else
constexpr Real BPL_TOL = 1e-8;
#endif

// A single segment evaluates as a pure power law 10 * x^{-2} at and above the normalization point x = 1
BOOST_AUTO_TEST_CASE(broken_power_law_single_segment) {
    BrokenPowerLaw<3> bpl;
    bpl.first_segment(10.0, 1.0, -2.0); // 10 * x^{-2}
    BOOST_CHECK_CLOSE(bpl.eval(1.0), 10.0, BPL_TOL);
    BOOST_CHECK_CLOSE(bpl.eval(2.0), 10.0 / 4.0, BPL_TOL);
    BOOST_CHECK_CLOSE(bpl.eval(10.0), 0.1, BPL_TOL);
}

// Adjacent segments join continuously: the value at the break x = 10 is 0.1
// and is approached from both sides
BOOST_AUTO_TEST_CASE(broken_power_law_continuity) {
    BrokenPowerLaw<3> bpl;
    bpl.first_segment(1.0, 1.0, -1.0); // x^{-1} below break
    bpl.add_segment(10.0, -2.0);       // steeper above break
    // At the break point, both segments should give the same value
    Real at_break = bpl.eval(10.0);
    // Continuity: value at break should be 1/10 = 0.1
    BOOST_CHECK_CLOSE(at_break, 0.1, 1e-4);
    // Values just below and above break should be close to the break value
    BOOST_CHECK_CLOSE(bpl.eval(9.999), at_break, 0.1);
    BOOST_CHECK_CLOSE(bpl.eval(10.001), at_break, 0.1);
}

// eval(x) agrees with exp2(log2_eval(log2(x))) across all segments, so the linear and log2 paths are consistent
BOOST_AUTO_TEST_CASE(broken_power_law_eval_matches_log2) {
    BrokenPowerLaw<3> bpl;
    bpl.first_segment(5.0, 2.0, -1.5);
    bpl.add_segment(100.0, -3.0);

    Real test_values[] = {3.0, 10.0, 50.0, 200.0, 1000.0};
    for (Real x : test_values) {
        Real from_eval = bpl.eval(x);
        Real from_log2 = std::exp2(bpl.log2_eval(std::log2(x)));
        BOOST_CHECK_CLOSE(from_eval, from_log2, BPL_TOL);
    }
}

// A rising/flat/falling three-segment law evaluates correctly inside each segment and at both break points
BOOST_AUTO_TEST_CASE(broken_power_law_three_segments) {
    BrokenPowerLaw<3> bpl;
    bpl.first_segment(1.0, 1.0, 1.0); // x^1 below 10
    bpl.add_segment(10.0, 0.0);       // flat between 10 and 100
    bpl.add_segment(100.0, -1.0);     // x^{-1} above 100

    // In first segment: f(5) = 5
    BOOST_CHECK_CLOSE(bpl.eval(5.0), 5.0, BPL_TOL);
    // At first break: f(10) = 10
    BOOST_CHECK_CLOSE(bpl.eval(10.0), 10.0, BPL_TOL);
    // In flat region: f(50) = 10
    BOOST_CHECK_CLOSE(bpl.eval(50.0), 10.0, BPL_TOL);
    // At second break: f(100) = 10
    BOOST_CHECK_CLOSE(bpl.eval(100.0), 10.0, BPL_TOL);
    // Above: f(200) = 10 * (200/100)^{-1} = 5
    BOOST_CHECK_CLOSE(bpl.eval(200.0), 5.0, BPL_TOL);
}

// clear() removes all segments, after which eval() returns 0
BOOST_AUTO_TEST_CASE(broken_power_law_clear) {
    BrokenPowerLaw<3> bpl;
    bpl.first_segment(10.0, 1.0, -2.0);
    bpl.clear();
    BOOST_CHECK_EQUAL(bpl.eval(5.0), 0.0);
}

// A very steep slope evaluates accurately without overflow: x^100 gives eval(2) = 2^100
BOOST_AUTO_TEST_CASE(broken_power_law_very_steep) {
    BrokenPowerLaw<3> bpl;
    bpl.first_segment(1.0, 1.0, 100.0); // x^100
    Real val = bpl.eval(2.0);
    BOOST_CHECK_CLOSE(val, std::pow(2.0, 100.0), 1e-4);
}

// Literal operators

// The _r user-defined literal yields Real values exactly equal to the plain floating-point literals
BOOST_AUTO_TEST_CASE(literal_operator_r) {
    Real a = 1.5_r;
    Real b = 42_r;
    BOOST_CHECK_EQUAL(a, 1.5);
    BOOST_CHECK_EQUAL(b, 42.0);
}

BOOST_AUTO_TEST_SUITE_END()
