#include <boost/test/unit_test.hpp>
#include <cmath>

// clang-format off
#include "radiation/synchrotron.h"        // Must come first (provides SynElectrons, Array)
#include "radiation/synchrotron-kernel.h"
// clang-format on
#include "util/macros.h"

// Forward declare the regime helper (defined in synchrotron.cpp)
size_t determine_regime(Real a, Real c, Real m);

BOOST_AUTO_TEST_SUITE(SynchrotronKernel)

// ============================================================================
//  synchrotron_F — analytic kernel approximation
// ============================================================================

BOOST_AUTO_TEST_CASE(synchrotron_F_zero) {
    BOOST_CHECK_EQUAL(synchrotron_F(0.0), 0.0);
}

BOOST_AUTO_TEST_CASE(synchrotron_F_negative) {
    BOOST_CHECK_EQUAL(synchrotron_F(-1.0), 0.0);
}

BOOST_AUTO_TEST_CASE(synchrotron_F_very_large) {
    // x > 500 returns 0
    BOOST_CHECK_EQUAL(synchrotron_F(501.0), 0.0);
}

BOOST_AUTO_TEST_CASE(synchrotron_F_positive_for_positive_x) {
    // F(x) should be positive for 0 < x < 500
    for (Real x : {0.001, 0.01, 0.1, 0.286, 1.0, 5.0, 10.0, 50.0, 100.0}) {
        BOOST_CHECK_GT(synchrotron_F(x), 0.0);
        BOOST_CHECK(std::isfinite(synchrotron_F(x)));
    }
}

BOOST_AUTO_TEST_CASE(synchrotron_F_peak_near_0_29) {
    // Peak of F(x) is near x ~ 0.286 with F_peak ~ 0.918
    // Sample to find approximate peak
    Real max_F = 0;
    Real x_peak = 0;
    for (Real x = 0.01; x < 2.0; x *= 1.05) {
        Real F = synchrotron_F(x);
        if (F > max_F) {
            max_F = F;
            x_peak = x;
        }
    }
    BOOST_CHECK_CLOSE(x_peak, 0.286, 15.0); // Peak within 15% of 0.286
    BOOST_CHECK_CLOSE(max_F, 0.918, 10.0);  // Peak value within 10% of 0.918
}

BOOST_AUTO_TEST_CASE(synchrotron_F_asymptote_small_x) {
    // F(x) ~ 2.15 * x^{1/3} for small x
    for (Real x : {1e-5, 1e-4, 1e-3}) {
        Real F = synchrotron_F(x);
        Real asympt = 2.15 * std::cbrt(x);
        BOOST_CHECK_CLOSE(F, asympt, 5.0); // Within 5%
    }
}

BOOST_AUTO_TEST_CASE(synchrotron_F_decays_exponentially) {
    // For large x, F should decay rapidly
    Real F_10 = synchrotron_F(10.0);
    Real F_50 = synchrotron_F(50.0);
    Real F_100 = synchrotron_F(100.0);
    BOOST_CHECK_GT(F_10, F_50);
    BOOST_CHECK_GT(F_50, F_100);
}

// ============================================================================
//  synchrotron_thin_I_nu
// ============================================================================

BOOST_AUTO_TEST_CASE(synchrotron_thin_I_nu_positive) {
    SynElectrons elec;
    elec.gamma_m = 100.0;
    elec.gamma_c = 1000.0;
    elec.gamma_a = 10.0;
    elec.gamma_M = 1e6;
    elec.p = 2.3;
    elec.N_e = 1e50;
    elec.column_den = 1e30;
    elec.Y_c = 0;
    elec.regime = 1;

    constexpr Real sin_angle_ave = con::pi / 4;
    const Real nu_L = 3 * con::e / (4 * con::pi * con::me * con::c);
    const Real P_coeff_4pi = std::sqrt(3.0) * con::e2 * sin_angle_ave / (3.0 * con::c);

    Real nu = nu_L * 100.0 * 100.0; // Near nu_m
    Real I = synchrotron_thin_I_nu(nu, elec, nu_L, P_coeff_4pi);
    BOOST_CHECK(std::isfinite(I));
    BOOST_CHECK_GT(I, 0.0);
}

BOOST_AUTO_TEST_CASE(synchrotron_thin_I_nu_zero_B) {
    SynElectrons elec;
    elec.gamma_m = 100.0;
    elec.gamma_c = 1000.0;
    elec.gamma_a = 10.0;
    elec.gamma_M = 1e6;
    elec.p = 2.3;
    elec.column_den = 1e30;
    elec.regime = 1;

    Real I = synchrotron_thin_I_nu(1e14, elec, 0.0, 1.0);
    BOOST_CHECK_EQUAL(I, 0.0);
}

// ============================================================================
//  SynKernelPhoton — build and query
// ============================================================================

namespace {

    SynKernelPhoton make_kernel_photon(Real nu_a, Real nu_m, Real nu_c, Real nu_M = 1e22, Real I_nu_max = 1e-20,
                                       Real p = 2.3) {
        SynKernelPhoton ph;
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

BOOST_AUTO_TEST_CASE(kernel_photon_build_valid) {
    auto ph = make_kernel_photon(1e10, 1e14, 1e16);
    // After build(), compute_I_nu should give positive results
    Real I = ph.compute_I_nu(1e14);
    BOOST_CHECK(std::isfinite(I));
    BOOST_CHECK_GT(I, 0.0);
}

BOOST_AUTO_TEST_CASE(kernel_photon_I_nu_positive) {
    auto ph = make_kernel_photon(1e10, 1e14, 1e16);
    for (Real nu : {1e10, 1e12, 1e14, 1e16, 1e18}) {
        Real I = ph.compute_I_nu(nu);
        BOOST_CHECK(std::isfinite(I));
        BOOST_CHECK_GT(I, 0.0);
    }
}

BOOST_AUTO_TEST_CASE(kernel_photon_log2_consistency) {
    auto ph = make_kernel_photon(1e10, 1e14, 1e16);
    for (Real nu : {1e12, 1e14, 1e16}) {
        Real I_lin = ph.compute_I_nu(nu);
        Real I_log = std::exp2(ph.compute_log2_I_nu(std::log2(nu)));
        BOOST_CHECK_CLOSE(I_lin, I_log, 0.1);
    }
}

BOOST_AUTO_TEST_CASE(kernel_photon_continuity) {
    auto ph = make_kernel_photon(1e10, 1e14, 1e16);
    constexpr int N = 60;
    const Real log_min = std::log10(1e10);
    const Real log_max = std::log10(1e20);
    const Real dlog = (log_max - log_min) / (N - 1);

    Real prev_I = ph.compute_I_nu(std::pow(10.0, log_min));
    for (int i = 1; i < N; ++i) {
        Real nu = std::pow(10.0, log_min + i * dlog);
        Real I = ph.compute_I_nu(nu);
        if (prev_I > 0 && I > 0) {
            Real ratio = I / prev_I;
            BOOST_CHECK_LT(ratio, 8.0);
            BOOST_CHECK_GT(ratio, 1.0 / 8.0);
        }
        prev_I = I;
    }
}

BOOST_AUTO_TEST_CASE(kernel_photon_zero_I_nu_max) {
    auto ph = make_kernel_photon(1e10, 1e14, 1e16, 1e22, 0.0);
    // I_nu_max=0 -> degenerate; just verify no crash
    Real I = ph.compute_I_nu(1e14);
    (void)I;
}

BOOST_AUTO_TEST_CASE(kernel_photon_all_regimes) {
    // Just verify all 6 regimes build without crash and produce finite output
    struct TestCase {
        Real nu_a, nu_m, nu_c;
    };
    TestCase cases[] = {
        {1e8, 1e12, 1e16}, // regime 1
        {1e12, 1e8, 1e16}, // regime 2
        {1e8, 1e16, 1e12}, // regime 3
        {1e12, 1e16, 1e8}, // regime 4
        {1e16, 1e8, 1e12}, // regime 5
        {1e16, 1e12, 1e8}, // regime 6
    };
    for (auto& tc : cases) {
        auto ph = make_kernel_photon(tc.nu_a, tc.nu_m, tc.nu_c);
        Real I = ph.compute_I_nu(std::sqrt(tc.nu_m * tc.nu_c));
        BOOST_CHECK(std::isfinite(I));
    }
}

BOOST_AUTO_TEST_SUITE_END()
