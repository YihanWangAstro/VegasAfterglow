#include <boost/test/unit_test.hpp>
#include <cmath>

#include "radiation/inverse-compton.h"
#include "radiation/synchrotron.h"
#include "util/macros.h"

// Forward declaration for internal function defined in synchrotron.cpp
size_t determine_regime(Real a, Real c, Real m);

BOOST_AUTO_TEST_SUITE(InverseCompton)

// ============================================================================
//  InverseComptonY — default constructor
// ============================================================================

// Default-constructed InverseComptonY has Y_T = 0 and Y(gamma) = 0 at every gamma.
BOOST_AUTO_TEST_CASE(default_constructor) {
    InverseComptonY Ys;
    BOOST_CHECK_CLOSE(Ys.Y_T, 0.0, 1e-10);
    BOOST_CHECK_CLOSE(Ys.gamma_spectrum(10.0), 0.0, 1e-10);
    BOOST_CHECK_CLOSE(Ys.gamma_spectrum(1000.0), 0.0, 1e-10);
}

// ============================================================================
//  InverseComptonY — Thomson regime (is_KN = false)
// ============================================================================

// Thomson regime (is_KN = false): constructor stores Y_T and Y(gamma) equals the constant Y_T at low gamma.
BOOST_AUTO_TEST_CASE(thomson_regime_constant_Y) {
    // Thomson regime: Y(gamma) = Y_T = constant for all gamma
    const Real gamma_m = 100.0;
    const Real gamma_c = 1000.0;
    const Real p = 2.3;
    const Real B = 1.0 * unit::Gauss;
    const Real Y_T = 0.5;

    InverseComptonY Ys(gamma_m, gamma_c, p, B, Y_T, false);
    BOOST_CHECK_CLOSE(Ys.Y_T, Y_T, 1e-10);

    // In Thomson regime (regime 0), Y should be constant = Y_T at low gamma
    Real Y_at_1 = Ys.gamma_spectrum(1.0);
    BOOST_CHECK_CLOSE(Y_at_1, Y_T, 1e-6);
}

// Thomson-regime Y(gamma) is finite and non-negative for gamma from 1 to 1e3.
BOOST_AUTO_TEST_CASE(thomson_regime_positive) {
    InverseComptonY Ys(100.0, 1000.0, 2.3, 1.0 * unit::Gauss, 1.0, false);
    for (Real gamma : {1.0, 10.0, 100.0, 1000.0}) {
        Real Y = Ys.gamma_spectrum(gamma);
        BOOST_CHECK(std::isfinite(Y));
        BOOST_CHECK_GE(Y, 0.0);
    }
}

// ============================================================================
//  InverseComptonY — KN regime (is_KN = true)
// ============================================================================

// KN-regime construction (is_KN = true, fast cooling) succeeds and yields a finite Y(gamma).
BOOST_AUTO_TEST_CASE(kn_regime_builds) {
    // KN regime should build without crash
    InverseComptonY Ys(100.0, 50.0, 2.3, 1.0 * unit::Gauss, 2.0, true);
    BOOST_CHECK(std::isfinite(Ys.gamma_spectrum(100.0)));
}

// KN-regime Y(gamma) stays finite and non-negative from gamma = 1 up to 1e6.
BOOST_AUTO_TEST_CASE(kn_regime_Y_finite) {
    // In KN regime, Y(gamma) should be finite and non-negative for all gamma
    InverseComptonY Ys(100.0, 50.0, 2.3, 1.0 * unit::Gauss, 5.0, true);
    for (Real gamma : {1.0, 10.0, 100.0, 1e4, 1e6}) {
        Real Y = Ys.gamma_spectrum(gamma);
        BOOST_CHECK(std::isfinite(Y));
        BOOST_CHECK_GE(Y, 0.0);
    }
}

// ============================================================================
//  InverseComptonY — nu_spectrum
// ============================================================================

// nu_spectrum (nu mapped to gamma via compute_syn_gamma, then gamma_spectrum) returns a finite, non-negative Y.
BOOST_AUTO_TEST_CASE(nu_spectrum_positive) {
    InverseComptonY Ys(100.0, 1000.0, 2.3, 1.0 * unit::Gauss, 1.0, false);
    Real B = 1.0 * unit::Gauss;
    // nu_spectrum converts nu -> gamma via compute_syn_gamma, then calls gamma_spectrum
    Real nu = compute_syn_freq(100.0, B);
    Real Y = Ys.nu_spectrum(nu);
    BOOST_CHECK(std::isfinite(Y));
    BOOST_CHECK_GE(Y, 0.0);
}

// ============================================================================
//  InverseComptonY — update_cooling_breaks
// ============================================================================

// update_cooling_breaks stores the new Y_T and leaves gamma_spectrum finite after the update.
BOOST_AUTO_TEST_CASE(update_cooling_breaks_no_crash) {
    InverseComptonY Ys(100.0, 50.0, 2.3, 1.0 * unit::Gauss, 2.0, true);
    // Update with new gamma_c and Y_T
    Ys.update_cooling_breaks(200.0, 3.0);
    BOOST_CHECK(std::isfinite(Ys.gamma_spectrum(100.0)));
    BOOST_CHECK_CLOSE(Ys.Y_T, 3.0, 1e-10);
}

// ============================================================================
//  compton_correction — KN cross section ratio
// ============================================================================

// Thomson limit: at photon energy x = h*nu/(me*c^2) = 1e-3, sigma/sigma_T equals 1 within 1%.
BOOST_AUTO_TEST_CASE(compton_correction_low_energy) {
    // For low-energy photons (x << 1): sigma/sigma_T ~ 1
    Real nu_low = 1e-3 * con::me * con::c2 / con::h;
    Real ratio = compton_correction(nu_low);
    BOOST_CHECK_CLOSE(ratio, 1.0, 1.0); // Within 1%
}

// Deep Klein-Nishina limit: at x = 1e3 the cross-section is strongly suppressed, 0 < sigma/sigma_T < 0.1.
BOOST_AUTO_TEST_CASE(compton_correction_high_energy) {
    // For high-energy photons (x >> 1): sigma/sigma_T << 1
    Real nu_high = 1e3 * con::me * con::c2 / con::h;
    Real ratio = compton_correction(nu_high);
    BOOST_CHECK_LT(ratio, 0.1);
    BOOST_CHECK_GT(ratio, 0.0);
}

// sigma/sigma_T is monotonically non-increasing with photon energy from x = 1e-2 to 100.
BOOST_AUTO_TEST_CASE(compton_correction_monotonic) {
    // sigma/sigma_T should decrease with increasing photon energy
    Real prev_ratio = compton_correction(1e-2 * con::me * con::c2 / con::h);
    for (Real x_fac : {0.1, 1.0, 10.0, 100.0}) {
        Real nu = x_fac * con::me * con::c2 / con::h;
        Real ratio = compton_correction(nu);
        BOOST_CHECK_LE(ratio, prev_ratio + 1e-10);
        prev_ratio = ratio;
    }
}

// nu = 0 returns exactly sigma/sigma_T = 0.
BOOST_AUTO_TEST_CASE(compton_correction_zero) {
    Real ratio = compton_correction(0.0);
    BOOST_CHECK_EQUAL(ratio, 0.0);
}

// Negative nu is guarded: returns exactly 0.
BOOST_AUTO_TEST_CASE(compton_correction_negative) {
    Real ratio = compton_correction(-1.0);
    BOOST_CHECK_EQUAL(ratio, 0.0);
}

// ============================================================================
//  eta_rad_Thomson
// ============================================================================

// Fast cooling (gamma_c < gamma_m): radiative efficiency eta is exactly 1.
BOOST_AUTO_TEST_CASE(eta_rad_thomson_fast_cooling) {
    // gamma_c < gamma_m: returns 1
    BOOST_CHECK_CLOSE(eta_rad_Thomson(100.0, 50.0, 2.3), 1.0, 1e-10);
}

// Slow cooling (gamma_c > gamma_m): eta equals (gamma_c/gamma_m)^(2-p).
BOOST_AUTO_TEST_CASE(eta_rad_thomson_slow_cooling) {
    // gamma_c > gamma_m: returns (gamma_c/gamma_m)^(2-p)
    const Real gamma_m = 100.0;
    const Real gamma_c = 1000.0;
    const Real p = 2.3;
    Real expected = fast_pow(gamma_c / gamma_m, 2 - p);
    BOOST_CHECK_CLOSE(eta_rad_Thomson(gamma_m, gamma_c, p), expected, 1e-6);
}

// Slow cooling with p = 2: eta = (gamma_c/gamma_m)^0 = 1 exactly.
BOOST_AUTO_TEST_CASE(eta_rad_thomson_p_equals_2) {
    // p=2: (gamma_c/gamma_m)^0 = 1 for slow cooling
    BOOST_CHECK_CLOSE(eta_rad_Thomson(100.0, 1000.0, 2.0), 1.0, 1e-10);
}

// ============================================================================
//  compute_Thomson_Y
// ============================================================================

// Thomson Y solving Y(1+Y) = eta*eps_e/eps_B is finite and strictly positive for slow-cooling parameters.
BOOST_AUTO_TEST_CASE(compute_Thomson_Y_positive) {
    RadParams rad;
    rad.eps_e = 0.1;
    rad.eps_B = 0.01;
    rad.p = 2.3;
    Real Y = compute_Thomson_Y(rad, 100.0, 1000.0);
    BOOST_CHECK(std::isfinite(Y));
    BOOST_CHECK_GT(Y, 0.0);
}

// ============================================================================
//  inverse_compton_correction
// ============================================================================

// With Y_c = 0 and Y(nu) = 0, the correction (1+Y_c)/(1+Y(nu)) is exactly 1 below nu_c.
BOOST_AUTO_TEST_CASE(inverse_compton_correction_below_nu_c) {
    // Below cooling frequency, correction should be ~1 (no IC effect)
    SmoothPowerLawSyn ph;
    ph.Y_c = 0;
    ph.Ys = InverseComptonY(); // default: Y=0 everywhere
    ph.nu_c = 1e16;

    Real corr = inverse_compton_correction(ph, 1e14); // below nu_c
    // (1 + 0) / (1 + 0) = 1
    BOOST_CHECK_CLOSE(corr, 1.0, 1e-10);
}

// With Y_c = 1 and Thomson-regime Ys, the correction (1+Y_c)/(1+Y(nu)) above nu_c is finite and positive.
BOOST_AUTO_TEST_CASE(inverse_compton_correction_with_Y) {
    // With Y_c > 0, correction = (1+Y_c) / (1+Y(nu))
    // For Thomson regime (Y = const = Y_T), correction = (1+Y_T)/(1+Y_T) = 1
    SmoothPowerLawSyn ph;
    ph.Y_c = 1.0;
    ph.Ys = InverseComptonY(100.0, 1000.0, 2.3, 1.0 * unit::Gauss, 1.0, false);
    ph.nu_c = 1e16;

    Real corr = inverse_compton_correction(ph, 1e18); // above nu_c
    BOOST_CHECK(std::isfinite(corr));
    BOOST_CHECK_GT(corr, 0.0);
}

// ============================================================================
//  ICPhoton template — basic construction and query
// ============================================================================

// First compute_I_nu query triggers lazy IC spectrum generation and returns a finite, non-negative intensity.
BOOST_AUTO_TEST_CASE(ICPhoton_construction) {
    // Build electron distribution
    SynElectrons elec;
    elec.gamma_m = 100.0;
    elec.gamma_c = 1000.0;
    elec.gamma_a = 10.0;
    elec.gamma_M = 1e6;
    elec.p = 2.3;
    elec.N_e = 1e50;
    elec.column_den = 1e30;
    elec.Y_c = 0;
    elec.Ys = InverseComptonY();
    elec.regime = determine_regime(10.0, 1000.0, 100.0);

    // Build photon distribution
    SmoothPowerLawSyn ph;
    ph.p = 2.3;
    ph.I_nu_max = 1e-20;
    ph.nu_m = 1e14;
    ph.nu_c = 1e16;
    ph.nu_a = 1e10;
    ph.nu_M = 1e22;
    ph.regime = determine_regime(1e10, 1e16, 1e14);
    ph.build();

    // Create ICPhoton
    ICPhoton<SynElectrons, SmoothPowerLawSyn> ic(elec, ph, false);

    // Query — first call triggers lazy generation
    Real I = ic.compute_I_nu(1e20);
    BOOST_CHECK(std::isfinite(I));
    BOOST_CHECK_GE(I, 0.0);
}

// compute_I_nu(nu) and exp2(compute_log2_I_nu(log2(nu))) agree within 1% when both are positive.
BOOST_AUTO_TEST_CASE(ICPhoton_log2_consistency) {
    SynElectrons elec;
    elec.gamma_m = 100.0;
    elec.gamma_c = 1000.0;
    elec.gamma_a = 10.0;
    elec.gamma_M = 1e6;
    elec.p = 2.3;
    elec.N_e = 1e50;
    elec.column_den = 1e30;
    elec.Y_c = 0;
    elec.Ys = InverseComptonY();
    elec.regime = determine_regime(10.0, 1000.0, 100.0);

    SmoothPowerLawSyn ph;
    ph.p = 2.3;
    ph.I_nu_max = 1e-20;
    ph.nu_m = 1e14;
    ph.nu_c = 1e16;
    ph.nu_a = 1e10;
    ph.nu_M = 1e22;
    ph.regime = determine_regime(1e10, 1e16, 1e14);
    ph.build();

    ICPhoton<SynElectrons, SmoothPowerLawSyn> ic(elec, ph, false);

    Real nu = 1e20;
    Real I_lin = ic.compute_I_nu(nu);
    Real I_log = std::exp2(ic.compute_log2_I_nu(std::log2(nu)));
    if (I_lin > 0 && I_log > 0) {
        BOOST_CHECK_CLOSE(I_lin, I_log, 1.0);
    }
}

// With Klein-Nishina corrections enabled, compute_I_nu still returns a finite, non-negative intensity.
BOOST_AUTO_TEST_CASE(ICPhoton_with_KN) {
    SynElectrons elec;
    elec.gamma_m = 100.0;
    elec.gamma_c = 1000.0;
    elec.gamma_a = 10.0;
    elec.gamma_M = 1e6;
    elec.p = 2.3;
    elec.N_e = 1e50;
    elec.column_den = 1e30;
    elec.Y_c = 0;
    elec.Ys = InverseComptonY();
    elec.regime = determine_regime(10.0, 1000.0, 100.0);

    SmoothPowerLawSyn ph;
    ph.p = 2.3;
    ph.I_nu_max = 1e-20;
    ph.nu_m = 1e14;
    ph.nu_c = 1e16;
    ph.nu_a = 1e10;
    ph.nu_M = 1e22;
    ph.regime = determine_regime(1e10, 1e16, 1e14);
    ph.build();

    ICPhoton<SynElectrons, SmoothPowerLawSyn> ic(elec, ph, true);

    Real I = ic.compute_I_nu(1e20);
    BOOST_CHECK(std::isfinite(I));
    BOOST_CHECK_GE(I, 0.0);
}

// ============================================================================
//  compton_correction — quantitative check against the exact KN cross-section
// ============================================================================

// Log-space LUT matches the exact total Klein-Nishina cross-section ratio at x = 0.1, 1, 10
// within 0.5% (measured LUT + fast_log2 error < 2e-4 relative).
BOOST_AUTO_TEST_CASE(compton_correction_matches_exact_KN) {
    // compton_correction(nu) evaluates sigma_KN/sigma_T at x = h*nu/(me*c^2)
    // via a log-space lookup table. Compare against the exact total
    // Klein-Nishina cross-section ratio computed directly:
    //   sigma/sigma_T = 3/4 * [ (1+x)/x^3 * (2x(1+x)/(1+2x) - ln(1+2x))
    //                           + ln(1+2x)/(2x) - (1+3x)/(1+2x)^2 ]
    auto sigma_kn_exact = [](Real x) {
        const Real l = std::log1p(2.0 * x);
        return 0.75 * ((1.0 + x) / (x * x * x) * (2.0 * x * (1.0 + x) / (1.0 + 2.0 * x) - l) + 0.5 * l / x -
                       (1.0 + 3.0 * x) / ((1.0 + 2.0 * x) * (1.0 + 2.0 * x)));
    };

    // Measured LUT + fast_log2 error is < 2e-4 relative at all three points;
    // assert within 0.5% for headroom.
    for (Real x : {0.1, 1.0, 10.0}) {
        Real nu = x * con::me * con::c2 / con::h;
        BOOST_CHECK_CLOSE(compton_correction(nu), sigma_kn_exact(x), 0.5);
    }
}

// ============================================================================
//  compute_Thomson_Y — quantitative limits of the Y(1+Y) = b root
// ============================================================================

// Fast-cooling Y is the exact root of Y(1+Y) = b = eps_e/eps_B (to 1e-6 relative), and matches
// the asymptotic limits Y ~ b for b << 1 and Y ~ sqrt(b) for b >> 1 within 10%.
BOOST_AUTO_TEST_CASE(compute_Thomson_Y_limits) {
    // The implementation solves Y(1+Y) = b with b = eta * eps_e/eps_B.
    // Fast cooling (gamma_c < gamma_m) gives eta = 1, so b = eps_e/eps_B.
    RadParams rad;
    rad.p = 2.3;
    const Real gamma_m = 1000.0;
    const Real gamma_c = 100.0;

    // b << 1: linear limit Y ~ b (measured: Y = 9.990e-4 for b = 1e-3)
    rad.eps_e = 1e-4;
    rad.eps_B = 0.1;
    const Real b_small = rad.eps_e / rad.eps_B;
    Real Y_small = compute_Thomson_Y(rad, gamma_m, gamma_c);
    BOOST_CHECK_CLOSE(Y_small, b_small, 10.0);
    // Exact root of Y(1+Y) = b to 1e-6 relative
    BOOST_CHECK_CLOSE(Y_small * (1.0 + Y_small), b_small, 1e-4);

    // b >> 1: square-root limit Y ~ sqrt(b) (measured: Y = 99.50 for b = 1e4)
    rad.eps_e = 0.1;
    rad.eps_B = 1e-5;
    const Real b_large = rad.eps_e / rad.eps_B;
    Real Y_large = compute_Thomson_Y(rad, gamma_m, gamma_c);
    BOOST_CHECK_CLOSE(Y_large, std::sqrt(b_large), 10.0);
    BOOST_CHECK_CLOSE(Y_large * (1.0 + Y_large), b_large, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
