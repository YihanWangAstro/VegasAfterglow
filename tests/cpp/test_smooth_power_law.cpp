#include <boost/test/unit_test.hpp>
#include <cmath>

#include "radiation/smooth-power-law-syn.h"
#include "radiation/synchrotron.h"
#include "util/macros.h"

// Forward declare the regime helper (defined in synchrotron.cpp, not in any header)
size_t determine_regime(Real a, Real c, Real m);

namespace {

    // Common test parameters
    constexpr Real kP = 2.3;
    constexpr Real kI_nu_max = 1e-20;
    constexpr Real kNu_M = 1e22;

    // Helper: build a SmoothPowerLawSyn for a given set of characteristic frequencies.
    // Arguments are nu_a, nu_m, nu_c (the regime is computed automatically).
    SmoothPowerLawSyn make_photon(Real nu_a, Real nu_m, Real nu_c, Real nu_M = kNu_M, Real I_nu_max = kI_nu_max,
                                  Real p = kP) {
        SmoothPowerLawSyn ph;
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

BOOST_AUTO_TEST_SUITE(SmoothPowerLaw)

// ---------------------------------------------------------------------------
// 1. build_sets_cached_values
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(build_sets_cached_values) {
    // After build(), the log2 cached values should be nonzero (for nonzero inputs)
    auto ph = make_photon(1e10, 1e14, 1e16);

    BOOST_CHECK_NE(ph.log2_nu_m, 0.0);
    BOOST_CHECK_NE(ph.log2_nu_c, 0.0);
    BOOST_CHECK_NE(ph.log2_nu_a, 0.0);
    BOOST_CHECK_NE(ph.log2_nu_M, 0.0);
    BOOST_CHECK_NE(ph.log2_I_nu_max, 0.0);

    // Verify log2 values are consistent with the linear values
    BOOST_CHECK_CLOSE(ph.log2_nu_m, std::log2(ph.nu_m), 1e-6);
    BOOST_CHECK_CLOSE(ph.log2_nu_c, std::log2(ph.nu_c), 1e-6);
    BOOST_CHECK_CLOSE(ph.log2_nu_a, std::log2(ph.nu_a), 1e-6);
    BOOST_CHECK_CLOSE(ph.log2_nu_M, std::log2(ph.nu_M), 1e-6);
    BOOST_CHECK_CLOSE(ph.log2_I_nu_max, std::log2(ph.I_nu_max), 1e-6);
}

// ---------------------------------------------------------------------------
// 2. compute_I_nu_positive
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_I_nu_positive) {
    // I_nu should be positive for frequencies within the optically thin range [nu_a, nu_M]
    auto ph = make_photon(1e10, 1e14, 1e16);

    Real test_freqs[] = {1e10, 1e11, 1e13, 1e14, 1e15, 1e16, 1e18, 1e20};
    for (Real nu : test_freqs) {
        Real I = ph.compute_I_nu(nu);
        BOOST_CHECK_GT(I, 0.0);
        BOOST_CHECK(std::isfinite(I));
    }
}

// ---------------------------------------------------------------------------
// 3. compute_I_nu_zero_outside
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_I_nu_zero_outside) {
    // Far above nu_M the exponential cutoff should suppress the spectrum to ~0
    auto ph = make_photon(1e10, 1e14, 1e16);

    Real I_at_nu_M = ph.compute_I_nu(ph.nu_M);
    Real I_far_above = ph.compute_I_nu(1e25); // well above nu_M = 1e22

    // The far-above value should be negligibly small compared to the value at nu_M
    BOOST_CHECK_LT(I_far_above, I_at_nu_M * 1e-10);
}

// ---------------------------------------------------------------------------
// 4. compute_I_nu_matches_log2
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_I_nu_matches_log2) {
    // compute_I_nu(nu) should equal exp2(compute_log2_I_nu(log2(nu)))
    auto ph = make_photon(1e10, 1e14, 1e16);

    Real test_freqs[] = {1e10, 1e12, 1e14, 1e15, 1e16, 1e18, 1e20};
    for (Real nu : test_freqs) {
        Real I_linear = ph.compute_I_nu(nu);
        Real I_from_log2 = std::exp2(ph.compute_log2_I_nu(std::log2(nu)));
        // Allow 0.01% relative tolerance for fast_log2/fast_exp2 round-trip
        BOOST_CHECK_CLOSE(I_linear, I_from_log2, 0.01);
    }
}

// ---------------------------------------------------------------------------
// 5. spectrum_continuity
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(spectrum_continuity) {
    // No jumps > factor of 4 between adjacent log-spaced points across the full range
    auto ph = make_photon(1e10, 1e14, 1e16);

    constexpr int N = 100;
    const Real log_nu_min = std::log10(ph.nu_a);
    const Real log_nu_max = std::log10(ph.nu_M);
    const Real dlog = (log_nu_max - log_nu_min) / (N - 1);

    Real prev_I = ph.compute_I_nu(std::pow(10.0, log_nu_min));
    for (int i = 1; i < N; ++i) {
        Real nu = std::pow(10.0, log_nu_min + i * dlog);
        Real I = ph.compute_I_nu(nu);
        if (prev_I > 0 && I > 0) {
            Real ratio = I / prev_I;
            BOOST_CHECK_LT(ratio, 4.0);
            BOOST_CHECK_GT(ratio, 1.0 / 4.0);
        }
        prev_I = I;
    }
}

// ---------------------------------------------------------------------------
// 6. spectrum_peak_near_expected
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(spectrum_peak_near_expected) {
    // For optically thin slow cooling (nu_a << nu_m < nu_c), peak should be near nu_m
    auto ph = make_photon(1e8, 1e14, 1e18);

    constexpr int N = 200;
    const Real log_nu_min = std::log10(1e10);
    const Real log_nu_max = std::log10(1e20);
    const Real dlog = (log_nu_max - log_nu_min) / (N - 1);

    Real max_I = 0;
    Real nu_peak = 0;
    for (int i = 0; i < N; ++i) {
        Real nu = std::pow(10.0, log_nu_min + i * dlog);
        Real I = ph.compute_I_nu(nu);
        if (I > max_I) {
            max_I = I;
            nu_peak = nu;
        }
    }

    // Peak should be within a factor of ~10 of nu_m for optically thin slow cooling
    // (the smooth broken power law can shift the peak somewhat)
    BOOST_CHECK_GT(nu_peak, ph.nu_m / 10.0);
    BOOST_CHECK_LT(nu_peak, ph.nu_m * 10.0);
}

// ---------------------------------------------------------------------------
// 7. regime_1_slow_cooling — nu_a < nu_m < nu_c (a < m < c)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(regime_1_slow_cooling) {
    auto ph = make_photon(1e8, 1e12, 1e16);

    BOOST_CHECK_EQUAL(ph.regime, 1u);

    // Spectrum should be positive at characteristic frequencies
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_a), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_m), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_c), 0.0);

    // Below nu_a the spectrum should be rising (self-absorbed segment)
    Real I_low = ph.compute_I_nu(ph.nu_a * 0.1);
    Real I_at_a = ph.compute_I_nu(ph.nu_a);
    BOOST_CHECK_LT(I_low, I_at_a);
}

// ---------------------------------------------------------------------------
// 8. regime_2_slow_cooling — nu_m < nu_a < nu_c (m < a < c)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(regime_2_slow_cooling) {
    auto ph = make_photon(1e12, 1e8, 1e16);

    BOOST_CHECK_EQUAL(ph.regime, 2u);

    // Self-absorption suppresses spectrum below nu_a
    Real I_below_a = ph.compute_I_nu(ph.nu_a * 0.01);
    Real I_at_a = ph.compute_I_nu(ph.nu_a);
    BOOST_CHECK_LT(I_below_a, I_at_a);

    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_c), 0.0);
}

// ---------------------------------------------------------------------------
// 9. regime_3 — nu_a < nu_c < nu_m (a < c < m, fast cooling)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(regime_3) {
    auto ph = make_photon(1e8, 1e16, 1e12);

    BOOST_CHECK_EQUAL(ph.regime, 3u);

    // Fast cooling: spectrum positive across range
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_c), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_m), 0.0);

    // High-frequency tail should decay
    Real I_at_m = ph.compute_I_nu(ph.nu_m);
    Real I_above_m = ph.compute_I_nu(ph.nu_m * 100.0);
    BOOST_CHECK_LT(I_above_m, I_at_m);
}

// ---------------------------------------------------------------------------
// 10. regime_4 — nu_c < nu_a < nu_m (c < a < m)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(regime_4) {
    auto ph = make_photon(1e12, 1e16, 1e8);

    BOOST_CHECK_EQUAL(ph.regime, 4u);

    // Spectrum should be positive and well-behaved
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_c), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_a), 0.0);
    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_m), 0.0);
    BOOST_CHECK(std::isfinite(ph.compute_I_nu(1e14)));
}

// ---------------------------------------------------------------------------
// 11. regime_5 — nu_m < nu_c < nu_a (m < c < a)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(regime_5) {
    auto ph = make_photon(1e16, 1e8, 1e12);

    BOOST_CHECK_EQUAL(ph.regime, 5u);

    // Deep self-absorption: below nu_a should be suppressed
    Real I_below = ph.compute_I_nu(ph.nu_a * 0.01);
    Real I_at = ph.compute_I_nu(ph.nu_a);
    BOOST_CHECK_LT(I_below, I_at);

    BOOST_CHECK_GT(ph.compute_I_nu(ph.nu_m), 0.0);
}

// ---------------------------------------------------------------------------
// 12. regime_6 — nu_c < nu_m < nu_a (c < m < a)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(regime_6) {
    auto ph = make_photon(1e16, 1e12, 1e8);

    BOOST_CHECK_EQUAL(ph.regime, 6u);

    // Deep self-absorption regime: spectrum below nu_a is suppressed
    Real I_below = ph.compute_I_nu(ph.nu_a * 0.01);
    Real I_at = ph.compute_I_nu(ph.nu_a);
    BOOST_CHECK_LT(I_below, I_at);

    BOOST_CHECK(std::isfinite(ph.compute_I_nu(ph.nu_c)));
    BOOST_CHECK(std::isfinite(ph.compute_I_nu(ph.nu_m)));
}

// ---------------------------------------------------------------------------
// 13. build_zero_I_nu_max
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(build_zero_I_nu_max) {
    // When I_nu_max = 0 the spectrum should be identically zero everywhere
    auto ph = make_photon(1e10, 1e14, 1e16, kNu_M, 0.0);

    Real test_freqs[] = {1e8, 1e12, 1e14, 1e16, 1e20};
    for (Real nu : test_freqs) {
        Real I = ph.compute_I_nu(nu);
        BOOST_CHECK_SMALL(I, 1e-300);
    }
}

// ---------------------------------------------------------------------------
// 14. build_equal_frequencies — degenerate break nu_m = nu_c
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(build_equal_frequencies) {
    // nu_m = nu_c should not crash and should produce a valid spectrum
    auto ph = make_photon(1e10, 1e14, 1e14);

    BOOST_CHECK(std::isfinite(ph.compute_I_nu(1e12)));
    BOOST_CHECK(std::isfinite(ph.compute_I_nu(1e14)));
    BOOST_CHECK(std::isfinite(ph.compute_I_nu(1e16)));
    BOOST_CHECK_GT(ph.compute_I_nu(1e14), 0.0);
}

// ---------------------------------------------------------------------------
// 15. compute_I_nu_at_zero — nu = 0 should not crash
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_I_nu_at_zero) {
    auto ph = make_photon(1e10, 1e14, 1e16);

    Real I = ph.compute_I_nu(0.0);
    // nu=0 -> log2(0)=-inf -> may produce NaN; just verify no crash
    // The function is not designed for nu=0 input; check it doesn't throw
    (void)I; // suppress unused warning; we only care it doesn't crash
}

// ---------------------------------------------------------------------------
// 16. compute_I_nu_at_infinity — very high frequency returns ~0
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_I_nu_at_infinity) {
    auto ph = make_photon(1e10, 1e14, 1e16);

    Real I = ph.compute_I_nu(1e30);
    // Exponential cutoff at nu >> nu_M should kill the spectrum
    BOOST_CHECK_SMALL(I, 1e-100);
}

// ---------------------------------------------------------------------------
// 17. spectrum_nu_m_equals_nu_c — smooth, no discontinuity
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(spectrum_nu_m_equals_nu_c) {
    // With nu_m = nu_c, the spectrum should still be continuous
    auto ph = make_photon(1e10, 1e15, 1e15);

    constexpr int N = 50;
    const Real log_nu_min = std::log10(1e12);
    const Real log_nu_max = std::log10(1e18);
    const Real dlog = (log_nu_max - log_nu_min) / (N - 1);

    Real prev_I = ph.compute_I_nu(std::pow(10.0, log_nu_min));
    for (int i = 1; i < N; ++i) {
        Real nu = std::pow(10.0, log_nu_min + i * dlog);
        Real I = ph.compute_I_nu(nu);
        if (prev_I > 0 && I > 0) {
            Real ratio = I / prev_I;
            // No discontinuity: ratio should stay within a generous factor-of-4 window
            BOOST_CHECK_LT(ratio, 4.0);
            BOOST_CHECK_GT(ratio, 1.0 / 4.0);
        }
        prev_I = I;
    }
}

// ---------------------------------------------------------------------------
// 18. spectrum_nu_a_very_small — effectively no self-absorption
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(spectrum_nu_a_very_small) {
    // nu_a << nu_m: self-absorption should have negligible effect on optically thin part
    auto ph_thin = make_photon(1e2, 1e14, 1e16);
    auto ph_ref = make_photon(1e10, 1e14, 1e16);

    // Above nu_m, the two spectra should be nearly identical
    Real test_freqs[] = {1e14, 1e15, 1e16, 1e18};
    for (Real nu : test_freqs) {
        Real I_thin = ph_thin.compute_I_nu(nu);
        Real I_ref = ph_ref.compute_I_nu(nu);
        // Should agree within 10% above the optically thin break
        BOOST_CHECK_CLOSE(I_thin, I_ref, 10.0);
    }
}

// ---------------------------------------------------------------------------
// 19. spectrum_nu_a_very_large — deeply self-absorbed
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(spectrum_nu_a_very_large) {
    // nu_a >> nu_m, nu_c: the spectrum below nu_a is heavily self-absorbed
    auto ph = make_photon(1e18, 1e12, 1e14);

    // Below nu_a, the self-absorbed spectrum should be much weaker than at nu_a
    Real I_at_a = ph.compute_I_nu(ph.nu_a);
    Real I_below = ph.compute_I_nu(ph.nu_a * 0.001);
    BOOST_CHECK_LT(I_below, I_at_a);

    // The spectrum at nu_a should still be positive and finite
    BOOST_CHECK_GT(I_at_a, 0.0);
    BOOST_CHECK(std::isfinite(I_at_a));

    // Above nu_a the spectrum should decay (power-law tail + exponential cutoff)
    Real I_above = ph.compute_I_nu(ph.nu_a * 100.0);
    BOOST_CHECK_LT(I_above, I_at_a);
}

BOOST_AUTO_TEST_SUITE_END()
