#include <boost/test/unit_test.hpp>
#include <cmath>
#include <limits>

#include "radiation/synchrotron.h"
#include "util/macros.h"

// Forward declarations for internal functions defined in synchrotron.cpp.
// order() is inline in the .cpp so it cannot be linked from here, but
// determine_regime() has external linkage and can be called directly.
size_t determine_regime(Real a, Real c, Real m);

// These helpers are also defined (non-inline) in synchrotron.cpp:
Real compute_single_elec_P_nu_max(Real B, Real p);
Real compute_syn_I_peak(Real B, Real p, Real column_den);
Real compute_syn_gamma_M(Real B, Real Y, Real p);
Real compute_syn_gamma_m(Real Gamma_th, Real gamma_M, Real eps_e, Real p, Real xi);
Real cool_after_crossing(Real gamma_x, Real gamma_m_x, Real gamma_m, Real dt_comv, Real B, Real Y);
Real compute_syn_gamma_a(Real B, Real I_syn_peak, Real gamma_m, Real gamma_c, Real gamma_M, Real p,
                         InverseComptonY const& Ys, Real Y_c);
Real cyclotron_correction(Real gamma_m, Real p);

BOOST_AUTO_TEST_SUITE(Synchrotron)

// ============================================================================
//  determine_regime — six canonical orderings
// ============================================================================

// Regime 1: a <= m <= c  (slow cooling, weak absorption)
BOOST_AUTO_TEST_CASE(determine_regime_1) {
    // a=1, c=10, m=5  =>  order(a,m,c) = order(1,5,10) = true  =>  regime 1
    BOOST_CHECK_EQUAL(determine_regime(1, 10, 5), 1u);
}

// Regime 2: m <= a <= c
BOOST_AUTO_TEST_CASE(determine_regime_2) {
    // a=5, c=10, m=1  =>  order(m,a,c) = order(1,5,10) = true  =>  regime 2
    BOOST_CHECK_EQUAL(determine_regime(5, 10, 1), 2u);
}

// Regime 3: a <= c <= m  (fast cooling, weak absorption)
BOOST_AUTO_TEST_CASE(determine_regime_3) {
    // a=1, c=5, m=10  =>  order(a,c,m) = order(1,5,10) = true  =>  regime 3
    BOOST_CHECK_EQUAL(determine_regime(1, 5, 10), 3u);
}

// Regime 4: c <= a <= m
BOOST_AUTO_TEST_CASE(determine_regime_4) {
    // a=5, c=1, m=10  =>  order(c,a,m) = order(1,5,10) = true  =>  regime 4
    BOOST_CHECK_EQUAL(determine_regime(5, 1, 10), 4u);
}

// Regime 5: m <= c <= a
BOOST_AUTO_TEST_CASE(determine_regime_5) {
    // a=10, c=5, m=1  =>  order(m,c,a) = order(1,5,10) = true  =>  regime 5
    BOOST_CHECK_EQUAL(determine_regime(10, 5, 1), 5u);
}

// Regime 6: c <= m <= a
BOOST_AUTO_TEST_CASE(determine_regime_6) {
    // a=10, c=1, m=5  =>  order(c,m,a) = order(1,5,10) = true  =>  regime 6
    BOOST_CHECK_EQUAL(determine_regime(10, 1, 5), 6u);
}

// ============================================================================
//  determine_regime — degenerate / boundary cases
// ============================================================================

// All three values equal: the first branch (a<=m<=c) matches
BOOST_AUTO_TEST_CASE(determine_regime_equal_values) {
    size_t r = determine_regime(5, 5, 5);
    // With a=c=m=5, order(a,m,c) = order(5,5,5) = true  =>  regime 1
    BOOST_CHECK_EQUAL(r, 1u);
}

// Two values equal, third distinct — test three sub-cases
BOOST_AUTO_TEST_CASE(determine_regime_two_equal) {
    // a == m < c  =>  order(a,m,c) still true  =>  regime 1
    BOOST_CHECK_EQUAL(determine_regime(3, 10, 3), 1u);

    // a == c < m  =>  order(a,c,m) = order(3,3,10) = true => regime 3
    //   (order(a,m,c) = order(3,10,3) = false first)
    BOOST_CHECK_EQUAL(determine_regime(3, 3, 10), 3u);

    // c == m < a  =>  order(m,c,a) = order(5,5,10) = true  =>  regime 5
    BOOST_CHECK_EQUAL(determine_regime(10, 5, 5), 5u);
}

// ============================================================================
//  compute_syn_freq / compute_syn_gamma — formulae and roundtrip
// ============================================================================

// Characteristic synchrotron frequency matches nu = 3*e/(4*pi*me*c) * B * gamma^2 exactly and is positive
BOOST_AUTO_TEST_CASE(compute_syn_freq_values) {
    const Real gamma = 100.0;
    const Real B = 1.0 * unit::Gauss;
    const Real nu = compute_syn_freq(gamma, B);

    // nu = 3*e/(4*pi*me*c) * B * gamma^2
    const Real expected = 3 * con::e / (4 * con::pi * con::me * con::c) * B * gamma * gamma;
    BOOST_CHECK_CLOSE(nu, expected, 1e-10);
    BOOST_CHECK_GT(nu, 0);
}

// Emitting Lorentz factor matches gamma = sqrt(4*pi*me*c/(3*e) * nu/B) exactly and is positive
BOOST_AUTO_TEST_CASE(compute_syn_gamma_values) {
    const Real B = 0.1 * unit::Gauss;
    const Real nu = 1e14 * unit::Hz;
    const Real gamma = compute_syn_gamma(nu, B);

    // gamma = sqrt(4*pi*me*c/(3*e) * nu/B)
    const Real expected = std::sqrt((4 * con::pi * con::me * con::c / (3 * con::e)) * (nu / B));
    BOOST_CHECK_CLOSE(gamma, expected, 1e-10);
    BOOST_CHECK_GT(gamma, 0);
}

// compute_syn_gamma inverts compute_syn_freq: gamma -> nu -> gamma recovers the input
BOOST_AUTO_TEST_CASE(compute_syn_freq_gamma_roundtrip) {
    const Real gamma_in = 300.0;
    const Real B = 0.5 * unit::Gauss;

    const Real nu = compute_syn_freq(gamma_in, B);
    const Real gamma_out = compute_syn_gamma(nu, B);

    BOOST_CHECK_CLOSE(gamma_out, gamma_in, 1e-8);
}

// ============================================================================
//  compute_syn_freq / compute_syn_gamma — edge cases
// ============================================================================

// Zero magnetic field yields exactly zero synchrotron frequency
BOOST_AUTO_TEST_CASE(compute_syn_freq_zero_B) {
    // B=0 should give zero frequency
    BOOST_CHECK_EQUAL(compute_syn_freq(100.0, 0.0), 0.0);
}

// For positive B and nu the emitting Lorentz factor is finite and positive
BOOST_AUTO_TEST_CASE(compute_syn_gamma_positive_B) {
    // For positive B and nu, gamma should be finite and positive
    const Real gamma = compute_syn_gamma(1e10 * unit::Hz, 1.0 * unit::Gauss);
    BOOST_CHECK(std::isfinite(gamma));
    BOOST_CHECK_GT(gamma, 0);
}

// ============================================================================
//  compute_gamma_c
// ============================================================================

// Without IC cooling (Y=0), gamma_c equals the quadratic root (gamma_bar + sqrt(gamma_bar^2 + 4))/2
// with gamma_bar = 6*pi*me*c/(sigmaT*B^2*t), and exceeds 1
BOOST_AUTO_TEST_CASE(compute_gamma_c_no_IC) {
    // Without IC cooling (Y=0), check the analytic formula
    const Real B = 1.0 * unit::Gauss;
    const Real t = 1.0 * unit::day;
    const Real gc = compute_gamma_c(t, B, 0.0);

    // gamma_bar = 6*pi*me*c / (sigmaT * B^2 * t)
    const Real gamma_bar = (6 * con::pi * con::me * con::c / con::sigmaT) / (B * B * t);
    const Real expected = (gamma_bar + std::sqrt(gamma_bar * gamma_bar + 4)) / 2;

    BOOST_CHECK_CLOSE(gc, expected, 1e-10);
    BOOST_CHECK_GT(gc, 1.0);
}

// IC cooling (Y>0) lowers gamma_c relative to the Y=0 value while keeping it above 1
BOOST_AUTO_TEST_CASE(compute_gamma_c_with_IC) {
    // IC cooling (Y > 0) should lower gamma_c compared to Y=0
    const Real B = 1.0 * unit::Gauss;
    const Real t = 1.0 * unit::day;

    const Real gc_no_IC = compute_gamma_c(t, B, 0.0);
    const Real gc_with_IC = compute_gamma_c(t, B, 1.0);

    BOOST_CHECK_LT(gc_with_IC, gc_no_IC);
    BOOST_CHECK_GT(gc_with_IC, 1.0);
}

// gamma_c never drops below 1 across a grid of B (0.01-100 G) and comoving time (1 s - 1 yr)
BOOST_AUTO_TEST_CASE(compute_gamma_c_positive) {
    // gamma_c should always be >= 1 for physical parameters
    for (Real B : {0.01 * unit::Gauss, 1.0 * unit::Gauss, 100.0 * unit::Gauss}) {
        for (Real t : {1.0 * unit::sec, 1.0 * unit::day, 1.0 * unit::yr}) {
            Real gc = compute_gamma_c(t, B, 0.0);
            BOOST_CHECK_GE(gc, 1.0);
        }
    }
}

// With B=0 there is no synchrotron cooling, so gamma_c is effectively unbounded (> 1e10)
BOOST_AUTO_TEST_CASE(compute_gamma_c_zero_B) {
    // B=0 means no synchrotron cooling: gamma_c should be very large
    const Real gc = compute_gamma_c(1.0 * unit::day, 0.0, 0.0);
    BOOST_CHECK_GT(gc, 1e10);
}

// In the long-time limit (t = 1e20 s) gamma_c approaches 1 from above: 1 <= gamma_c < 2
BOOST_AUTO_TEST_CASE(compute_gamma_c_huge_t) {
    // Very large comoving time -> gamma_c approaches 1 from above
    const Real gc = compute_gamma_c(1e20 * unit::sec, 1.0 * unit::Gauss, 0.0);
    BOOST_CHECK_GE(gc, 1.0);
    BOOST_CHECK_LT(gc, 2.0);
}

// ============================================================================
//  compute_single_elec_P_nu_max / compute_syn_I_peak
// ============================================================================

// Peak single-electron synchrotron power is finite and positive for B = 1 G
BOOST_AUTO_TEST_CASE(compute_single_elec_P_nu_max_positive) {
    const Real P = compute_single_elec_P_nu_max(1.0 * unit::Gauss, 2.3);
    BOOST_CHECK(std::isfinite(P));
    BOOST_CHECK_GT(P, 0);
}

// Peak synchrotron specific intensity is finite and positive for physical B, p, and column density
BOOST_AUTO_TEST_CASE(compute_syn_I_peak_positive) {
    const Real I = compute_syn_I_peak(1.0 * unit::Gauss, 2.3, 1e20);
    BOOST_CHECK(std::isfinite(I));
    BOOST_CHECK_GT(I, 0);
}

// ============================================================================
//  SynElectrons — spectrum and column density
// ============================================================================

// Helper to build a minimal SynElectrons struct for testing
static SynElectrons make_test_electrons(Real gamma_a, Real gamma_m, Real gamma_c, Real gamma_M, Real p) {
    SynElectrons e;
    e.gamma_m = gamma_m;
    e.gamma_c = gamma_c;
    e.gamma_a = gamma_a;
    e.gamma_M = gamma_M;
    e.p = p;
    e.N_e = 1e50;
    e.column_den = 1e30;
    e.Y_c = 0;
    e.Ys = InverseComptonY(); // default: no IC correction
    e.regime = determine_regime(gamma_a, gamma_c, gamma_m);
    return e;
}

// Slow-cooling ordering is classified as regime 1 and N(gamma) is finite and positive
// from below gamma_m up past gamma_c
BOOST_AUTO_TEST_CASE(syn_electrons_spectrum_slow_cooling) {
    // Slow cooling (regime 1): gamma_a < gamma_m < gamma_c
    auto e = make_test_electrons(10, 100, 1000, 1e8, 2.3);
    BOOST_CHECK_EQUAL(e.regime, 1u);

    // Evaluate at several points — must be finite and positive
    for (Real gamma : {50.0, 100.0, 500.0, 1000.0, 5000.0}) {
        Real N = e.compute_N_gamma(gamma);
        BOOST_CHECK(std::isfinite(N));
        BOOST_CHECK_GT(N, 0);
    }
}

// Fast-cooling ordering is classified as regime 3 and N(gamma) is finite and positive
// from below gamma_c up past gamma_m
BOOST_AUTO_TEST_CASE(syn_electrons_spectrum_fast_cooling) {
    // Fast cooling (regime 3): gamma_a < gamma_c < gamma_m
    auto e = make_test_electrons(5, 1000, 50, 1e8, 2.5);
    BOOST_CHECK_EQUAL(e.regime, 3u);

    for (Real gamma : {30.0, 50.0, 200.0, 1000.0, 5000.0}) {
        Real N = e.compute_N_gamma(gamma);
        BOOST_CHECK(std::isfinite(N));
        BOOST_CHECK_GT(N, 0);
    }
}

// Per-gamma electron column density is finite and positive across the distribution
BOOST_AUTO_TEST_CASE(syn_electrons_column_den) {
    auto e = make_test_electrons(10, 100, 1000, 1e8, 2.3);

    for (Real gamma : {50.0, 100.0, 500.0, 2000.0}) {
        Real col = e.compute_column_den(gamma);
        BOOST_CHECK(std::isfinite(col));
        BOOST_CHECK_GT(col, 0);
    }
}

// The regime flips 1 -> 3 as gamma_m crosses gamma_c, and N(gamma) stays finite and
// positive on both sides of the slow/fast-cooling transition
BOOST_AUTO_TEST_CASE(syn_electrons_regime_transitions) {
    // Near the boundary between regimes: gamma_m ~ gamma_c
    // Slightly slow-cooling
    auto e_slow = make_test_electrons(1, 100, 101, 1e8, 2.3);
    BOOST_CHECK_EQUAL(e_slow.regime, 1u);

    // Slightly fast-cooling
    auto e_fast = make_test_electrons(1, 101, 100, 1e8, 2.3);
    BOOST_CHECK_EQUAL(e_fast.regime, 3u);

    // Both should produce finite, positive spectra near the transition
    Real N_slow = e_slow.compute_N_gamma(100.5);
    Real N_fast = e_fast.compute_N_gamma(100.5);
    BOOST_CHECK(std::isfinite(N_slow));
    BOOST_CHECK(std::isfinite(N_fast));
    BOOST_CHECK_GT(N_slow, 0);
    BOOST_CHECK_GT(N_fast, 0);
}

// ============================================================================
//  Functions defined only in synchrotron.cpp (extern linkage, no header decl)
// ============================================================================

// --- compute_syn_gamma_M ---

// Burn-off limit matches gamma_M = sqrt(6*pi*e/(sigmaT*B*(1+Y))) at Y=0 and exceeds 1
BOOST_AUTO_TEST_CASE(compute_syn_gamma_M_positive) {
    const Real B = 1.0 * unit::Gauss;
    const Real gM = compute_syn_gamma_M(B, 0.0, 2.3);
    // gamma_M = sqrt(6*pi*e / (sigmaT * B * (1+Y)))
    const Real expected = std::sqrt(6 * con::pi * con::e / con::sigmaT / (B * 1.0));
    BOOST_CHECK_CLOSE(gM, expected, 1e-8);
    BOOST_CHECK_GT(gM, 1.0);
}

// IC cooling (Y=1) reduces gamma_M below the Y=0 value while keeping it above 1
BOOST_AUTO_TEST_CASE(compute_syn_gamma_M_with_IC) {
    const Real B = 1.0 * unit::Gauss;
    const Real gM_no_IC = compute_syn_gamma_M(B, 0.0, 2.3);
    const Real gM_with_IC = compute_syn_gamma_M(B, 1.0, 2.3);
    // IC cooling should reduce gamma_M
    BOOST_CHECK_LT(gM_with_IC, gM_no_IC);
    BOOST_CHECK_GT(gM_with_IC, 1.0);
}

// B=0 removes the synchrotron burn-off limit entirely: gamma_M is infinite
BOOST_AUTO_TEST_CASE(compute_syn_gamma_M_zero_B) {
    // B=0 -> no synchrotron limit -> infinity
    const Real gM = compute_syn_gamma_M(0.0, 0.0, 2.3);
    BOOST_CHECK(std::isinf(gM));
}

// --- compute_syn_gamma_m ---

// For p > 2, gamma_m - 1 equals the shock-heating partition (p-2)/(p-1) * eps_e*(Gamma_th-1)*(mp/me)/xi
BOOST_AUTO_TEST_CASE(compute_syn_gamma_m_basic_p_gt_2) {
    // p > 2: gamma_m - 1 = (p-2)/(p-1) * eps_e*(Gamma_th-1)*(mp/me)/xi
    const Real Gamma_th = 10.0;
    const Real gamma_M = 1e8;
    const Real eps_e = 0.1;
    const Real p = 2.5;
    const Real xi = 1.0;

    const Real gm = compute_syn_gamma_m(Gamma_th, gamma_M, eps_e, p, xi);
    const Real expected_minus_1 = (p - 2) / (p - 1) * eps_e * (Gamma_th - 1) * (con::mp / con::me) / xi;
    BOOST_CHECK_CLOSE(gm, expected_minus_1 + 1, 1e-8);
    BOOST_CHECK_GT(gm, 1.0);
}

// The p=2 branch (bisection root solve) yields a finite gamma_m above 1
BOOST_AUTO_TEST_CASE(compute_syn_gamma_m_p_equals_2) {
    // p=2 uses root_bisect internally
    const Real gm = compute_syn_gamma_m(10.0, 1e8, 0.1, 2.0, 1.0);
    BOOST_CHECK(std::isfinite(gm));
    BOOST_CHECK_GT(gm, 1.0);
}

// The p < 2 formula branch yields a finite gamma_m above 1
BOOST_AUTO_TEST_CASE(compute_syn_gamma_m_p_lt_2) {
    // p < 2: different formula branch
    const Real gm = compute_syn_gamma_m(10.0, 1e8, 0.1, 1.8, 1.0);
    BOOST_CHECK(std::isfinite(gm));
    BOOST_CHECK_GT(gm, 1.0);
}

// Vanishing eps_e drives gamma_m toward unity: result stays bounded in (1, 10)
BOOST_AUTO_TEST_CASE(compute_syn_gamma_m_low_eps_e) {
    // Very small eps_e -> gamma_m close to 1
    const Real gm = compute_syn_gamma_m(2.0, 1e8, 1e-10, 2.3, 1.0);
    BOOST_CHECK_GT(gm, 1.0);
    BOOST_CHECK_LT(gm, 10.0);
}

// --- cool_after_crossing ---

// When gamma_m == gamma_m_x the adiabatic factor is 1 and gamma_x is returned unchanged
BOOST_AUTO_TEST_CASE(cool_after_crossing_no_change) {
    // If gamma_m == gamma_m_x, f_ad = 1, result = gamma_x
    const Real result = cool_after_crossing(100.0, 50.0, 50.0, 1.0, 1.0, 0.0);
    BOOST_CHECK_CLOSE(result, 100.0, 1e-10);
}

// With gamma_m < gamma_m_x, adiabatic expansion cools the electron to
// (gamma_x - 1)*f_ad + 1 with f_ad = (gamma_m - 1)/(gamma_m_x - 1) < 1
BOOST_AUTO_TEST_CASE(cool_after_crossing_adiabatic_expansion) {
    // gamma_m < gamma_m_x means f_ad < 1: adiabatic cooling reduces the result
    const Real gamma_x = 100.0;
    const Real gamma_m_x = 50.0;
    const Real gamma_m = 30.0; // smaller than gamma_m_x
    const Real result = cool_after_crossing(gamma_x, gamma_m_x, gamma_m, 1.0, 1.0, 0.0);
    // f_ad = (30-1)/(50-1) = 29/49 < 1
    const Real f_ad = (gamma_m - 1.0) / (gamma_m_x - 1.0);
    const Real expected = (gamma_x - 1.0) * f_ad + 1.0;
    BOOST_CHECK_CLOSE(result, expected, 1e-10);
    BOOST_CHECK_LT(result, gamma_x);
}

// The post-crossing Lorentz factor never drops below 1
BOOST_AUTO_TEST_CASE(cool_after_crossing_positive) {
    // Result should always be >= 1
    const Real result = cool_after_crossing(5.0, 100.0, 50.0, 10.0, 1.0, 0.0);
    BOOST_CHECK_GE(result, 1.0);
}

// --- compute_syn_gamma_a ---

// Self-absorption Lorentz factor gamma_a is finite and >= 1 for slow-cooling parameters
BOOST_AUTO_TEST_CASE(compute_syn_gamma_a_positive) {
    // gamma_a should be positive and >= 1 for physical parameters
    const Real B = 1.0 * unit::Gauss;
    const Real I_peak = 1e-20;
    const Real gamma_m = 100.0;
    const Real gamma_c = 1000.0;
    const Real gamma_M = 1e8;
    const Real p = 2.3;

    const Real ga = compute_syn_gamma_a(B, I_peak, gamma_m, gamma_c, gamma_M, p, InverseComptonY{}, 0.0);
    BOOST_CHECK(std::isfinite(ga));
    BOOST_CHECK_GE(ga, 1.0);
}

// gamma_a stays finite and >= 1 in fast cooling (gamma_c < gamma_m)
BOOST_AUTO_TEST_CASE(compute_syn_gamma_a_fast_cooling) {
    // Fast cooling: gamma_c < gamma_m
    const Real B = 0.1 * unit::Gauss;
    const Real I_peak = 1e-22;
    const Real ga = compute_syn_gamma_a(B, I_peak, 1000.0, 50.0, 1e8, 2.3, InverseComptonY{}, 0.0);
    BOOST_CHECK(std::isfinite(ga));
    BOOST_CHECK_GE(ga, 1.0);
}

// Larger peak intensity means stronger self-absorption, so gamma_a increases monotonically with I_peak
BOOST_AUTO_TEST_CASE(compute_syn_gamma_a_large_I_peak) {
    // Large I_peak -> stronger self-absorption -> larger gamma_a
    const Real B = 1.0 * unit::Gauss;
    const Real ga_small = compute_syn_gamma_a(B, 1e-25, 100.0, 1000.0, 1e8, 2.3, InverseComptonY{}, 0.0);
    const Real ga_large = compute_syn_gamma_a(B, 1e-18, 100.0, 1000.0, 1e8, 2.3, InverseComptonY{}, 0.0);
    BOOST_CHECK_GT(ga_large, ga_small);
}

// --- cyclotron_correction ---

// For gamma_m >> 1 the cyclotron correction factor approaches 1 (within 0.2%)
BOOST_AUTO_TEST_CASE(cyclotron_correction_large_gamma_m) {
    // For gamma_m >> 1: f ~ 1
    const Real f = cyclotron_correction(1000.0, 2.3);
    BOOST_CHECK_CLOSE(f, 1.0, 0.2);
}

// As gamma_m -> 1 the correction factor vanishes: bounded in [0, 0.1)
BOOST_AUTO_TEST_CASE(cyclotron_correction_gamma_m_near_1) {
    // gamma_m ~ 1: f ~ 0
    const Real f = cyclotron_correction(1.01, 2.3);
    BOOST_CHECK_LT(f, 0.1);
    BOOST_CHECK_GE(f, 0.0);
}

// For p > 3 the correction equals ((gamma_m-1)/gamma_m)^((p-1)/2)
BOOST_AUTO_TEST_CASE(cyclotron_correction_p_gt_3) {
    // For p > 3: f = ((gamma_m-1)/gamma_m)^((p-1)/2)
    const Real gamma_m = 10.0;
    const Real p = 3.5;
    const Real f = cyclotron_correction(gamma_m, p);
    const Real base = (gamma_m - 1.0) / gamma_m;
    const Real expected = fast_pow(base, (p - 1.0) / 2.0);
    BOOST_CHECK_CLOSE(f, expected, 1e-6);
}

// For p <= 3 the correction equals (gamma_m-1)/gamma_m exactly, with no extra exponent
BOOST_AUTO_TEST_CASE(cyclotron_correction_p_le_3) {
    // For p <= 3: f = (gamma_m-1)/gamma_m (no extra power)
    const Real gamma_m = 10.0;
    const Real p = 2.5;
    const Real f = cyclotron_correction(gamma_m, p);
    const Real expected = (gamma_m - 1.0) / gamma_m;
    BOOST_CHECK_CLOSE(f, expected, 1e-10);
}

// ============================================================================
//  SynElectrons — quantitative distribution slopes
// ============================================================================

// Local log-slope of the electron number distribution: d log N / d log gamma
static Real electron_local_slope(SynElectrons const& e, Real gamma) {
    const Real N_lo = e.compute_N_gamma(gamma / 1.01);
    const Real N_hi = e.compute_N_gamma(gamma * 1.01);
    return (std::log(N_hi) - std::log(N_lo)) / (2.0 * std::log(1.01));
}

// Slow-cooling log-slopes: -p between gamma_m and gamma_c, steepening to -(p+1) above
// gamma_c (pins measured midpoint slopes -2.299 and -3.301 to within 0.05)
BOOST_AUTO_TEST_CASE(syn_electrons_slope_slow_cooling) {
    // Slow cooling with cooling break: N(gamma) ~ gamma^-p between gamma_m and
    // gamma_c, steepening to gamma^-(p+1) between gamma_c and gamma_M.
    // (Ys default-constructed and Y_c=0, so no IC correction enters.)
    const Real p = 2.3;
    auto e = make_test_electrons(1.0, 100.0, 1e5, 1e8, p);
    e.N_e = 1.0;
    BOOST_CHECK_EQUAL(e.regime, 1u);

    // Geometric midpoint of [gamma_m, gamma_c] (measured: -2.299)
    BOOST_CHECK_LT(std::abs(electron_local_slope(e, 3162.0) - (-p)), 0.05);
    // Geometric midpoint of [gamma_c, gamma_M] (measured: -3.301)
    BOOST_CHECK_LT(std::abs(electron_local_slope(e, 3.16e6) - (-(p + 1.0))), 0.05);
}

// Fast-cooling log-slopes: -2 between gamma_c and gamma_m, steepening to -(p+1) above
// gamma_m (pins measured midpoint slopes -1.983 and -3.317 to within 0.05)
BOOST_AUTO_TEST_CASE(syn_electrons_slope_fast_cooling) {
    // Fast cooling: N(gamma) ~ gamma^-2 between gamma_c and gamma_m,
    // steepening to gamma^-(p+1) between gamma_m and gamma_M.
    const Real p = 2.3;
    auto e = make_test_electrons(1.0, 1e5, 100.0, 1e8, p);
    e.N_e = 1.0;
    BOOST_CHECK_EQUAL(e.regime, 3u);

    // Geometric midpoint of [gamma_c, gamma_m] (measured: -1.983)
    BOOST_CHECK_LT(std::abs(electron_local_slope(e, 3162.0) - (-2.0)), 0.05);
    // Geometric midpoint of [gamma_m, gamma_M] (measured: -3.317)
    BOOST_CHECK_LT(std::abs(electron_local_slope(e, 3.16e6) - (-(p + 1.0))), 0.05);
}

// ============================================================================
//  SynElectrons — smoke tests for the absorption-dominated regimes
// ============================================================================

// Ordering gamma_m < gamma_a < gamma_c is classified as regime 2 and N(gamma) stays finite and positive
BOOST_AUTO_TEST_CASE(syn_electrons_regime_2_smoke) {
    // Regime 2: gamma_m < gamma_a < gamma_c
    auto e = make_test_electrons(1e3, 100.0, 1e5, 1e8, 2.3);
    BOOST_CHECK_EQUAL(e.regime, 2u);
    for (Real gamma : {50.0, 500.0, 1e4, 1e6}) {
        Real N = e.compute_N_gamma(gamma);
        BOOST_CHECK(std::isfinite(N));
        BOOST_CHECK_GT(N, 0.0);
    }
}

// Ordering gamma_c < gamma_a < gamma_m is classified as regime 4 and N(gamma) stays finite and positive
BOOST_AUTO_TEST_CASE(syn_electrons_regime_4_smoke) {
    // Regime 4: gamma_c < gamma_a < gamma_m
    auto e = make_test_electrons(1e3, 1e5, 100.0, 1e8, 2.3);
    BOOST_CHECK_EQUAL(e.regime, 4u);
    for (Real gamma : {50.0, 500.0, 1e4, 1e6}) {
        Real N = e.compute_N_gamma(gamma);
        BOOST_CHECK(std::isfinite(N));
        BOOST_CHECK_GT(N, 0.0);
    }
}

// Ordering gamma_m < gamma_c < gamma_a is classified as regime 5 and N(gamma) stays finite and positive
BOOST_AUTO_TEST_CASE(syn_electrons_regime_5_smoke) {
    // Regime 5: gamma_m < gamma_c < gamma_a
    auto e = make_test_electrons(1e5, 100.0, 1e3, 1e8, 2.3);
    BOOST_CHECK_EQUAL(e.regime, 5u);
    for (Real gamma : {50.0, 500.0, 1e4, 1e6}) {
        Real N = e.compute_N_gamma(gamma);
        BOOST_CHECK(std::isfinite(N));
        BOOST_CHECK_GT(N, 0.0);
    }
}

// Ordering gamma_c < gamma_m < gamma_a is classified as regime 6 and N(gamma) stays finite and positive
BOOST_AUTO_TEST_CASE(syn_electrons_regime_6_smoke) {
    // Regime 6: gamma_c < gamma_m < gamma_a
    auto e = make_test_electrons(1e5, 1e3, 100.0, 1e8, 2.3);
    BOOST_CHECK_EQUAL(e.regime, 6u);
    for (Real gamma : {50.0, 500.0, 1e4, 1e6}) {
        Real N = e.compute_N_gamma(gamma);
        BOOST_CHECK(std::isfinite(N));
        BOOST_CHECK_GT(N, 0.0);
    }
}

BOOST_AUTO_TEST_SUITE_END()
