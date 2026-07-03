#include <boost/test/unit_test.hpp>
#include <cmath>

#include "radiation/power-law-syn.h"
#include "radiation/synchrotron.h"
#include "util/macros.h"

// Forward declare the regime helper (defined in synchrotron.cpp)
size_t determine_regime(Real a, Real c, Real m);

namespace {

    constexpr Real kP = 2.3;
    constexpr Real kI_nu_max = 1e-20;
    constexpr Real kNu_M = 1e22;

    PowerLawSyn make_photon(Real nu_a, Real nu_m, Real nu_c, Real nu_M = kNu_M, Real I_nu_max = kI_nu_max,
                            Real p = kP) {
        PowerLawSyn ph;
        ph.p = p;
        ph.I_nu_max = I_nu_max;
        ph.nu_m = nu_m;
        ph.nu_c = nu_c;
        ph.nu_a = nu_a;
        ph.nu_M = nu_M;
        ph.regime = determine_regime(nu_a, nu_c, nu_m);
        ph.build();
        return ph;
    }

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(PowerLawSynTests)

// ============================================================================
//  build() and cached values
// ============================================================================

// build() caches log2_I_nu_max as a nonzero value equal to log2(I_nu_max)
BOOST_AUTO_TEST_CASE(build_sets_log2_values) {
    auto ph = make_photon(1e10, 1e14, 1e16);
    // build() sets log2_I_nu_max; log2_nu_c is set by the factory (not build)
    BOOST_CHECK_NE(ph.log2_I_nu_max, 0.0);
    BOOST_CHECK_CLOSE(ph.log2_I_nu_max, std::log2(ph.I_nu_max), 1e-6);
}

// ============================================================================
//  compute_I_nu positive and finite
// ============================================================================

// Intensity is positive and finite at every frequency from nu_a up to well above nu_c
BOOST_AUTO_TEST_CASE(compute_I_nu_positive) {
    auto ph = make_photon(1e10, 1e14, 1e16);
    for (Real nu : {1e10, 1e12, 1e14, 1e15, 1e16, 1e18, 1e20}) {
        Real I = ph.compute_I_nu(nu);
        BOOST_CHECK_GT(I, 0.0);
        BOOST_CHECK(std::isfinite(I));
    }
}

// ============================================================================
//  Exponential cutoff above nu_M
// ============================================================================

// Exponential cutoff: intensity three decades above nu_M is suppressed by more than 10 orders of magnitude
BOOST_AUTO_TEST_CASE(compute_I_nu_cutoff) {
    auto ph = make_photon(1e10, 1e14, 1e16);
    Real I_at_M = ph.compute_I_nu(ph.nu_M);
    Real I_far = ph.compute_I_nu(1e25);
    BOOST_CHECK_LT(I_far, I_at_M * 1e-10);
}

// ============================================================================
//  compute_I_nu matches compute_log2_I_nu
// ============================================================================

// Linear-space compute_I_nu agrees with exp2(compute_log2_I_nu(log2(nu))) to 0.1% across the spectrum
BOOST_AUTO_TEST_CASE(compute_I_nu_matches_log2) {
    auto ph = make_photon(1e10, 1e14, 1e16);
    for (Real nu : {1e10, 1e12, 1e14, 1e16, 1e18}) {
        Real I_linear = ph.compute_I_nu(nu);
        Real I_from_log2 = std::exp2(ph.compute_log2_I_nu(std::log2(nu)));
        BOOST_CHECK_CLOSE(I_linear, I_from_log2, 0.1);
    }
}

// ============================================================================
//  All 6 regimes produce valid spectra
// ============================================================================

// Regime 1: a <= m <= c (slow cooling, weak absorption) is classified as 1; intensity positive at nu_m and nu_c
BOOST_AUTO_TEST_CASE(regime_1) {
    auto ph = make_photon(1e8, 1e12, 1e16);
    BOOST_CHECK_EQUAL(ph.regime, 1u);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_m), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_c), 0.0);
}

// Regime 2: m <= a <= c is classified as 2; intensity positive at nu_a and nu_c
BOOST_AUTO_TEST_CASE(regime_2) {
    auto ph = make_photon(1e12, 1e8, 1e16);
    BOOST_CHECK_EQUAL(ph.regime, 2u);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_a), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_c), 0.0);
}

// Regime 3: a <= c <= m (fast cooling, weak absorption) is classified as 3; intensity positive at nu_c and nu_m
BOOST_AUTO_TEST_CASE(regime_3) {
    auto ph = make_photon(1e8, 1e16, 1e12);
    BOOST_CHECK_EQUAL(ph.regime, 3u);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_c), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_m), 0.0);
}

// Regime 4: c <= a <= m is classified as 4; intensity positive at nu_a and nu_m
BOOST_AUTO_TEST_CASE(regime_4) {
    auto ph = make_photon(1e12, 1e16, 1e8);
    BOOST_CHECK_EQUAL(ph.regime, 4u);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_a), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_m), 0.0);
}

// Regime 5: m <= c <= a (strong absorption) is classified as 5; intensity positive at nu_a
BOOST_AUTO_TEST_CASE(regime_5) {
    auto ph = make_photon(1e16, 1e8, 1e12);
    BOOST_CHECK_EQUAL(ph.regime, 5u);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_a), 0.0);
}

// Regime 6: c <= m <= a (strong absorption) is classified as 6; intensity positive at nu_a
BOOST_AUTO_TEST_CASE(regime_6) {
    auto ph = make_photon(1e16, 1e12, 1e8);
    BOOST_CHECK_EQUAL(ph.regime, 6u);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_a), 0.0);
}

// ============================================================================
//  Spectrum continuity
// ============================================================================

// Piecewise spectrum has no jumps at the breaks: adjacent points on a 100-point log grid
// from nu_a to nu_M differ by less than a factor of 4
BOOST_AUTO_TEST_CASE(spectrum_continuity) {
    auto ph = make_photon(1e10, 1e14, 1e16);
    constexpr int N = 100;
    const Real log_min = std::log10(ph.nu_a);
    const Real log_max = std::log10(ph.nu_M);
    const Real dlog = (log_max - log_min) / (N - 1);

    Real prev_I = ph.compute_I_nu(std::pow(10.0, log_min));
    for (int i = 1; i < N; ++i) {
        Real nu = std::pow(10.0, log_min + i * dlog);
        Real I = ph.compute_I_nu(nu);
        if (prev_I > 0 && I > 0) {
            Real ratio = I / prev_I;
            BOOST_CHECK_LT(ratio, 4.0);
            BOOST_CHECK_GT(ratio, 1.0 / 4.0);
        }
        prev_I = I;
    }
}

// ============================================================================
//  Edge cases
// ============================================================================

// I_nu_max = 0 yields an effectively zero spectrum (|I_nu| < 1e-300) at all frequencies
BOOST_AUTO_TEST_CASE(build_zero_I_nu_max) {
    auto ph = make_photon(1e10, 1e14, 1e16, kNu_M, 0.0);
    for (Real nu : {1e12, 1e14, 1e16}) {
        BOOST_CHECK_SMALL(ph.compute_I_nu(nu), 1e-300);
    }
}

// Degenerate breaks nu_m = nu_c still yield a finite, positive intensity at the coincident break
BOOST_AUTO_TEST_CASE(build_equal_frequencies) {
    // nu_m = nu_c should not crash
    auto ph = make_photon(1e10, 1e14, 1e14);
    BOOST_CHECK(std::isfinite(ph.compute_I_nu(1e14)));
    BOOST_CHECK_GT(ph.compute_I_nu(1e14), 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
