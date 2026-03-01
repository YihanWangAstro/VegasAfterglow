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

BOOST_AUTO_TEST_CASE(default_constructor) {
    InverseComptonY Ys;
    BOOST_CHECK_CLOSE(Ys.Y_T, 0.0, 1e-10);
    BOOST_CHECK_CLOSE(Ys.gamma_spectrum(10.0), 0.0, 1e-10);
    BOOST_CHECK_CLOSE(Ys.gamma_spectrum(1000.0), 0.0, 1e-10);
}

// ============================================================================
//  InverseComptonY — Thomson regime (is_KN = false)
// ============================================================================

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

BOOST_AUTO_TEST_CASE(kn_regime_builds) {
    // KN regime should build without crash
    InverseComptonY Ys(100.0, 50.0, 2.3, 1.0 * unit::Gauss, 2.0, true);
    BOOST_CHECK(std::isfinite(Ys.gamma_spectrum(100.0)));
}

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

BOOST_AUTO_TEST_CASE(compton_correction_low_energy) {
    // For low-energy photons (x << 1): sigma/sigma_T ~ 1
    Real nu_low = 1e-3 * con::me * con::c2 / con::h;
    Real ratio = compton_correction(nu_low);
    BOOST_CHECK_CLOSE(ratio, 1.0, 1.0); // Within 1%
}

BOOST_AUTO_TEST_CASE(compton_correction_high_energy) {
    // For high-energy photons (x >> 1): sigma/sigma_T << 1
    Real nu_high = 1e3 * con::me * con::c2 / con::h;
    Real ratio = compton_correction(nu_high);
    BOOST_CHECK_LT(ratio, 0.1);
    BOOST_CHECK_GT(ratio, 0.0);
}

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

BOOST_AUTO_TEST_CASE(compton_correction_zero) {
    Real ratio = compton_correction(0.0);
    BOOST_CHECK_EQUAL(ratio, 0.0);
}

BOOST_AUTO_TEST_CASE(compton_correction_negative) {
    Real ratio = compton_correction(-1.0);
    BOOST_CHECK_EQUAL(ratio, 0.0);
}

// ============================================================================
//  eta_rad_Thomson
// ============================================================================

BOOST_AUTO_TEST_CASE(eta_rad_thomson_fast_cooling) {
    // gamma_c < gamma_m: returns 1
    BOOST_CHECK_CLOSE(eta_rad_Thomson(100.0, 50.0, 2.3), 1.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(eta_rad_thomson_slow_cooling) {
    // gamma_c > gamma_m: returns (gamma_c/gamma_m)^(2-p)
    const Real gamma_m = 100.0;
    const Real gamma_c = 1000.0;
    const Real p = 2.3;
    Real expected = fast_pow(gamma_c / gamma_m, 2 - p);
    BOOST_CHECK_CLOSE(eta_rad_Thomson(gamma_m, gamma_c, p), expected, 1e-6);
}

BOOST_AUTO_TEST_CASE(eta_rad_thomson_p_equals_2) {
    // p=2: (gamma_c/gamma_m)^0 = 1 for slow cooling
    BOOST_CHECK_CLOSE(eta_rad_Thomson(100.0, 1000.0, 2.0), 1.0, 1e-10);
}

// ============================================================================
//  compute_CMB_Y
// ============================================================================

BOOST_AUTO_TEST_CASE(compute_CMB_Y_positive) {
    Real Y = compute_CMB_Y(1.0 * unit::Gauss, 1.0);
    BOOST_CHECK(std::isfinite(Y));
    BOOST_CHECK_GT(Y, 0.0);
}

BOOST_AUTO_TEST_CASE(compute_CMB_Y_zero_B) {
    BOOST_CHECK_EQUAL(compute_CMB_Y(0.0, 1.0), 0.0);
}

BOOST_AUTO_TEST_CASE(compute_CMB_Y_zero_z) {
    BOOST_CHECK_EQUAL(compute_CMB_Y(1.0 * unit::Gauss, 0.0), 0.0);
}

BOOST_AUTO_TEST_CASE(compute_CMB_Y_increases_with_z) {
    // Higher redshift -> more CMB energy density -> larger Y
    Real Y_z1 = compute_CMB_Y(1.0 * unit::Gauss, 1.0);
    Real Y_z5 = compute_CMB_Y(1.0 * unit::Gauss, 5.0);
    BOOST_CHECK_GT(Y_z5, Y_z1);
}

BOOST_AUTO_TEST_CASE(compute_CMB_Y_decreases_with_B) {
    // Stronger B -> more magnetic energy density -> smaller Y_CMB
    Real Y_B1 = compute_CMB_Y(1.0 * unit::Gauss, 1.0);
    Real Y_B10 = compute_CMB_Y(10.0 * unit::Gauss, 1.0);
    BOOST_CHECK_GT(Y_B1, Y_B10);
}

// ============================================================================
//  compute_Thomson_Y
// ============================================================================

BOOST_AUTO_TEST_CASE(compute_Thomson_Y_positive) {
    RadParams rad;
    rad.eps_e = 0.1;
    rad.eps_B = 0.01;
    rad.p = 2.3;
    rad.cmb_cooling = false;
    Real Y = compute_Thomson_Y(rad, 100.0, 1000.0);
    BOOST_CHECK(std::isfinite(Y));
    BOOST_CHECK_GT(Y, 0.0);
}

BOOST_AUTO_TEST_CASE(compute_Thomson_Y_with_CMB) {
    RadParams rad;
    rad.eps_e = 0.1;
    rad.eps_B = 0.01;
    rad.p = 2.3;
    rad.cmb_cooling = true;
    Real Y_with = compute_Thomson_Y(rad, 100.0, 1000.0, 1.0 * unit::Gauss, 1.0);

    rad.cmb_cooling = false;
    Real Y_without = compute_Thomson_Y(rad, 100.0, 1000.0, 1.0 * unit::Gauss, 1.0);

    BOOST_CHECK_GT(Y_with, Y_without);
}

// ============================================================================
//  inverse_compton_correction
// ============================================================================

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

BOOST_AUTO_TEST_SUITE_END()
