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

BOOST_AUTO_TEST_CASE(build_sets_log2_values) {
    auto ph = make_photon(1e10, 1e14, 1e16);
    // build() sets log2_I_nu_max; log2_nu_c is set by the factory (not build)
    BOOST_CHECK_NE(ph.log2_I_nu_max, 0.0);
    BOOST_CHECK_CLOSE(ph.log2_I_nu_max, std::log2(ph.I_nu_max), 1e-6);
}

// ============================================================================
//  compute_I_nu positive and finite
// ============================================================================

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

BOOST_AUTO_TEST_CASE(compute_I_nu_cutoff) {
    auto ph = make_photon(1e10, 1e14, 1e16);
    Real I_at_M = ph.compute_I_nu(ph.nu_M);
    Real I_far = ph.compute_I_nu(1e25);
    BOOST_CHECK_LT(I_far, I_at_M * 1e-10);
}

// ============================================================================
//  compute_I_nu matches compute_log2_I_nu
// ============================================================================

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

BOOST_AUTO_TEST_CASE(regime_1) {
    auto ph = make_photon(1e8, 1e12, 1e16);
    BOOST_CHECK_EQUAL(ph.regime, 1u);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_m), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_c), 0.0);
}

BOOST_AUTO_TEST_CASE(regime_2) {
    auto ph = make_photon(1e12, 1e8, 1e16);
    BOOST_CHECK_EQUAL(ph.regime, 2u);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_a), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_c), 0.0);
}

BOOST_AUTO_TEST_CASE(regime_3) {
    auto ph = make_photon(1e8, 1e16, 1e12);
    BOOST_CHECK_EQUAL(ph.regime, 3u);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_c), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_m), 0.0);
}

BOOST_AUTO_TEST_CASE(regime_4) {
    auto ph = make_photon(1e12, 1e16, 1e8);
    BOOST_CHECK_EQUAL(ph.regime, 4u);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_a), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_m), 0.0);
}

BOOST_AUTO_TEST_CASE(regime_5) {
    auto ph = make_photon(1e16, 1e8, 1e12);
    BOOST_CHECK_EQUAL(ph.regime, 5u);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_a), 0.0);
}

BOOST_AUTO_TEST_CASE(regime_6) {
    auto ph = make_photon(1e16, 1e12, 1e8);
    BOOST_CHECK_EQUAL(ph.regime, 6u);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_a), 0.0);
}

// ============================================================================
//  Spectrum continuity
// ============================================================================

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

BOOST_AUTO_TEST_CASE(build_zero_I_nu_max) {
    auto ph = make_photon(1e10, 1e14, 1e16, kNu_M, 0.0);
    for (Real nu : {1e12, 1e14, 1e16}) {
        BOOST_CHECK_SMALL(ph.compute_I_nu(nu), 1e-300);
    }
}

BOOST_AUTO_TEST_CASE(build_equal_frequencies) {
    // nu_m = nu_c should not crash
    auto ph = make_photon(1e10, 1e14, 1e14);
    BOOST_CHECK(std::isfinite(ph.compute_I_nu(1e14)));
    BOOST_CHECK_GT(ph.compute_I_nu(1e14), 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
