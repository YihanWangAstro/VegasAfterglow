#include <boost/test/unit_test.hpp>
#include <cmath>

#include "environment/jet.h"
#include "environment/medium.h"
#include "util/macros.h"

BOOST_AUTO_TEST_SUITE(Environment)

// ---------------------------------------------------------------------------
// ISM tests
// ---------------------------------------------------------------------------

// ISM density is homogeneous: rho returns the same value regardless of phi, theta, or r.
BOOST_AUTO_TEST_CASE(ism_constant_density) {
    ISM ism(1.0);
    Real rho0 = ism.rho(0, 0, 1e15);
    // Density is the same regardless of phi, theta, r
    BOOST_CHECK_EQUAL(ism.rho(0.5, 0.3, 1e10), rho0);
    BOOST_CHECK_EQUAL(ism.rho(1.0, 1.0, 1e20), rho0);
    BOOST_CHECK_EQUAL(ism.rho(3.14, 1.57, 1e5), rho0);
}

// ISM mass density equals number density times proton mass: rho = n * m_p.
BOOST_AUTO_TEST_CASE(ism_density_value) {
    ISM ism(1.0);
    Real expected = 1.0 * con::mp;
    BOOST_CHECK_CLOSE(ism.rho(0, 0, 1e15), expected, 1e-10);
}

// Swept-up mass follows mass(r) = rho * r^3 / 3 for a constant-density medium.
BOOST_AUTO_TEST_CASE(ism_mass_cubic) {
    ISM ism(1.0);
    Real r = 1e16;
    Real expected = (1.0 * con::mp) * r * r * r / 3.0;
    BOOST_CHECK_CLOSE(ism.mass(r), expected, 1e-10);
}

// Zero number density gives exactly zero rho and zero swept-up mass.
BOOST_AUTO_TEST_CASE(ism_zero_density) {
    ISM ism(0.0);
    BOOST_CHECK_EQUAL(ism.rho(0, 0, 1e15), 0.0);
    BOOST_CHECK_EQUAL(ism.mass(1e16), 0.0);
}

// ---------------------------------------------------------------------------
// Wind tests
// ---------------------------------------------------------------------------

// With n_ism = 0 (no ISM floor) and n0 = inf (r02 = 0), wind density follows the pure profile rho = A / r^2.
BOOST_AUTO_TEST_CASE(wind_density_r_squared) {
    // For large r with no ISM floor, density ~ A/r^2
    Real A_star = 1.0;
    Wind wind(A_star, 0, con::inf);
    Real A = A_star * 5e11 * unit::g / unit::cm;

    Real r = 1e18;
    Real expected = A / (r * r); // r02 = 0 when n0 = inf
    BOOST_CHECK_CLOSE(wind.rho(0, 0, r), expected, 1e-6);
}

// Wind density never falls below the ISM floor n_ism * m_p at any radius.
BOOST_AUTO_TEST_CASE(wind_density_floor) {
    Real A_star = 0.1;
    Real n_ism = 1.0;
    Wind wind(A_star, n_ism);

    Real floor = n_ism * con::mp;
    // At any radius, density should be >= the ISM floor
    BOOST_CHECK_GE(wind.rho(0, 0, 1e10), floor);
    BOOST_CHECK_GE(wind.rho(0, 0, 1e15), floor);
    BOOST_CHECK_GE(wind.rho(0, 0, 1e20), floor);
    BOOST_CHECK_GE(wind.rho(0, 0, 1e30), floor);
}

// Pure wind (no ISM floor, r02 = 0): swept-up mass grows linearly, mass(r) = A * r.
BOOST_AUTO_TEST_CASE(wind_mass_linear) {
    // Wind(A_star, 0) with n0=inf gives r02=0, so mass ~ A*r for large r
    Real A_star = 1.0;
    Wind wind(A_star, 0, con::inf);
    Real A = A_star * 5e11 * unit::g / unit::cm;

    Real r = 1e18;
    Real expected = A * r; // no ISM term, r02=0
    BOOST_CHECK_CLOSE(wind.mass(r), expected, 1e-10);
}

// With both wind and ISM floor, mass(r) sums the two contributions: A * r + rho_ism * r^3 / 3.
BOOST_AUTO_TEST_CASE(wind_mass_with_floor) {
    // With both A_star and n_ism, mass includes A*r and rho_ism*r^3/3 terms
    Real A_star = 1.0;
    Real n_ism = 1.0;
    Wind wind(A_star, n_ism, con::inf);
    Real A = A_star * 5e11 * unit::g / unit::cm;
    Real rho_ism = n_ism * con::mp;

    Real r = 1e18;
    Real expected = A * r + rho_ism * r * r * r / 3.0;
    BOOST_CHECK_CLOSE(wind.mass(r), expected, 1e-10);
}

// A_star = 0 reduces the wind to a uniform ISM: rho equals n_ism * m_p at all radii.
BOOST_AUTO_TEST_CASE(wind_zero_A_star) {
    // Wind(0, n_ism) density equals ISM floor only
    Real n_ism = 1.0;
    Wind wind(0.0, n_ism);
    Real expected = n_ism * con::mp;
    BOOST_CHECK_CLOSE(wind.rho(0, 0, 1e15), expected, 1e-10);
    BOOST_CHECK_CLOSE(wind.rho(0, 0, 1e20), expected, 1e-10);
}

// A finite small-radius wind density n0 gives r02 > 0, so rho stays finite and positive as r -> 0.
BOOST_AUTO_TEST_CASE(wind_very_small_r) {
    // When n0 is finite, r02 > 0 prevents divergence at r -> 0
    Real A_star = 1.0;
    Real n0 = 1e6; // finite floor density
    Wind wind(A_star, 0, n0);
    Real rho_at_zero = wind.rho(0, 0, 0.0);
    // Should be finite (A / r02) since r02 > 0
    BOOST_CHECK(std::isfinite(rho_at_zero));
    BOOST_CHECK_GT(rho_at_zero, 0.0);
}

// At very large r the A/r^2 term is negligible and rho converges to the ISM floor (within 0.001%).
BOOST_AUTO_TEST_CASE(wind_very_large_r) {
    // At very large r, wind density -> rho_ism floor
    Real A_star = 1.0;
    Real n_ism = 1.0;
    Wind wind(A_star, n_ism);
    Real floor = n_ism * con::mp;

    // At r = 1e30, A/r^2 is negligible
    Real rho_far = wind.rho(0, 0, 1e30);
    BOOST_CHECK_CLOSE(rho_far, floor, 1e-3); // within 0.001%
}

// ---------------------------------------------------------------------------
// TophatJet tests
// ---------------------------------------------------------------------------

// Inside the core (theta < theta_c) a tophat jet has eps_k > 0 and Gamma0 > 1.
BOOST_AUTO_TEST_CASE(tophat_jet_inside_core) {
    Real theta_c = 0.1;
    Real E_iso = 1e52;
    Real Gamma0 = 300.0;
    TophatJet jet(theta_c, E_iso, Gamma0);

    // Inside core: eps_k > 0 and Gamma0 > 1
    BOOST_CHECK_GT(jet.eps_k(0, 0.05), 0.0);
    BOOST_CHECK_GT(jet.Gamma0(0, 0.05), 1.0);
}

// Outside the core a tophat jet carries nothing: eps_k = 0 and Gamma0 = 1 exactly.
BOOST_AUTO_TEST_CASE(tophat_jet_outside_core) {
    Real theta_c = 0.1;
    TophatJet jet(theta_c, 1e52, 300.0);

    // Outside core: eps_k = 0, Gamma0 = 1
    BOOST_CHECK_EQUAL(jet.eps_k(0, 0.2), 0.0);
    BOOST_CHECK_EQUAL(jet.Gamma0(0, 0.2), 1.0);
}

// On-axis energy per solid angle of a tophat jet equals E_iso / (4*pi).
BOOST_AUTO_TEST_CASE(tophat_jet_energy_value) {
    Real E_iso = 1e52;
    TophatJet jet(0.1, E_iso, 300.0);

    Real expected = E_iso / (4 * con::pi);
    BOOST_CHECK_CLOSE(jet.eps_k(0, 0.0), expected, 1e-10);
}

// The core edge is exclusive: theta == theta_c yields eps_k = 0 and Gamma0 = 1, just inside yields jet values.
BOOST_AUTO_TEST_CASE(tophat_jet_at_boundary) {
    Real theta_c = 0.1;
    TophatJet jet(theta_c, 1e52, 300.0);

    // theta == theta_c is NOT less than theta_c, so returns 0/1
    BOOST_CHECK_EQUAL(jet.eps_k(0, 0.1), 0.0);
    BOOST_CHECK_EQUAL(jet.Gamma0(0, 0.1), 1.0);

    // Just inside
    BOOST_CHECK_GT(jet.eps_k(0, 0.0999), 0.0);
    BOOST_CHECK_GT(jet.Gamma0(0, 0.0999), 1.0);
}

// The jet axis (theta = 0) is inside the core: eps_k = E_iso/(4*pi) and Gamma0 equals the input value.
BOOST_AUTO_TEST_CASE(tophat_jet_theta_zero) {
    TophatJet jet(0.1, 1e52, 300.0);
    // theta = 0 is the center, always inside
    BOOST_CHECK_CLOSE(jet.eps_k(0, 0.0), 1e52 / (4 * con::pi), 1e-10);
    BOOST_CHECK_CLOSE(jet.Gamma0(0, 0.0), 300.0, 1e-10);
}

// The equator (theta = pi/2) lies outside the core: eps_k = 0 and Gamma0 = 1.
BOOST_AUTO_TEST_CASE(tophat_jet_theta_pi_half) {
    TophatJet jet(0.1, 1e52, 300.0);
    // theta = pi/2 is far outside core
    BOOST_CHECK_EQUAL(jet.eps_k(0, con::pi / 2), 0.0);
    BOOST_CHECK_EQUAL(jet.Gamma0(0, con::pi / 2), 1.0);
}

// ---------------------------------------------------------------------------
// GaussianJet tests
// ---------------------------------------------------------------------------

// A Gaussian jet peaks on-axis: eps_k(0) = E_iso/(4*pi) and Gamma0(0) recovers the input Gamma0.
BOOST_AUTO_TEST_CASE(gaussian_jet_peak) {
    Real E_iso = 1e52;
    Real Gamma0 = 300.0;
    GaussianJet jet(0.1, E_iso, Gamma0);

    // At theta=0, eps_k = E_iso/(4*pi) * exp(0) = E_iso/(4*pi)
    Real expected_eps = E_iso / (4 * con::pi);
    BOOST_CHECK_CLOSE(jet.eps_k(0, 0.0), expected_eps, 1e-10);
    // At theta=0, Gamma0 = (Gamma0-1)*exp(0) + 1 = Gamma0
    BOOST_CHECK_CLOSE(jet.Gamma0(0, 0.0), Gamma0, 1e-10);
}

// Gaussian jet eps_k decreases monotonically with increasing theta.
BOOST_AUTO_TEST_CASE(gaussian_jet_falloff) {
    GaussianJet jet(0.1, 1e52, 300.0);

    // eps_k decreases with increasing theta
    Real eps_0 = jet.eps_k(0, 0.0);
    Real eps_05 = jet.eps_k(0, 0.05);
    Real eps_1 = jet.eps_k(0, 0.1);
    Real eps_3 = jet.eps_k(0, 0.3);

    BOOST_CHECK_GT(eps_0, eps_05);
    BOOST_CHECK_GT(eps_05, eps_1);
    BOOST_CHECK_GT(eps_1, eps_3);
}

// On-axis peak values hold for a narrower, more energetic jet (theta_c = 0.05, E_iso = 1e53, Gamma0 = 500).
BOOST_AUTO_TEST_CASE(gaussian_jet_theta_zero) {
    Real E_iso = 1e53;
    Real Gamma0 = 500.0;
    GaussianJet jet(0.05, E_iso, Gamma0);

    // Peak values at theta=0
    BOOST_CHECK_CLOSE(jet.eps_k(0, 0.0), E_iso / (4 * con::pi), 1e-10);
    BOOST_CHECK_CLOSE(jet.Gamma0(0, 0.0), Gamma0, 1e-10);
}

// Far outside the core (theta = 20*theta_c), eps_k drops below 1e-50 of the peak and Gamma0 -> 1.
BOOST_AUTO_TEST_CASE(gaussian_jet_large_theta) {
    GaussianJet jet(0.05, 1e52, 300.0);
    Real eps_peak = jet.eps_k(0, 0.0);

    // theta >> theta_c: eps_k falls to negligible fraction of peak
    Real eps_far = jet.eps_k(0, 1.0); // 1.0 >> 0.05
    BOOST_CHECK_LT(eps_far / eps_peak, 1e-50);

    // Gamma0 approaches 1 at large theta
    Real gamma_far = jet.Gamma0(0, 1.0);
    BOOST_CHECK_CLOSE(gamma_far, 1.0, 1e-6);
}

// ---------------------------------------------------------------------------
// PowerLawJet tests
// ---------------------------------------------------------------------------

// On-axis the power-law profile 1/(1 + (theta/theta_c)^k) equals 1, so eps_k = E_iso/(4*pi) and Gamma0 = input.
BOOST_AUTO_TEST_CASE(power_law_jet_peak) {
    Real E_iso = 1e52;
    Real Gamma0 = 300.0;
    PowerLawJet jet(0.1, E_iso, Gamma0, 4.0, 4.0);

    // At theta=0, fast_pow(0/theta_c, k) = 0, so eps_k = eps_k_ / (1+0) = eps_k_
    Real expected_eps = E_iso / (4 * con::pi);
    BOOST_CHECK_CLOSE(jet.eps_k(0, 0.0), expected_eps, 1e-10);
    BOOST_CHECK_CLOSE(jet.Gamma0(0, 0.0), Gamma0, 1e-10);
}

// In the wing, eps_k/eps_peak follows 1/(1 + (theta/theta_c)^k_e), within 1% (fast_pow vs std::pow).
BOOST_AUTO_TEST_CASE(power_law_jet_wing) {
    Real theta_c = 0.1;
    Real k_e = 4.0;
    Real k_g = 4.0;
    PowerLawJet jet(theta_c, 1e52, 300.0, k_e, k_g);

    // For theta >> theta_c, eps_k ~ eps_k_ / (theta/theta_c)^k_e
    Real theta = 1.0; // 10 * theta_c
    Real eps = jet.eps_k(0, theta);
    Real eps_peak = jet.eps_k(0, 0.0);

    // Expected ratio: 1 / (1 + (10)^4) ~ 1/10001
    Real expected_ratio = 1.0 / (1.0 + std::pow(theta / theta_c, k_e));
    BOOST_CHECK_CLOSE(eps / eps_peak, expected_ratio, 1.0); // within 1% (fast_pow vs std::pow)
}

// ---------------------------------------------------------------------------
// two_component function
// ---------------------------------------------------------------------------

// two_component returns the core value for theta <= theta_c, the wing value for theta <= theta_w,
// and 0 beyond; both boundaries are inclusive.
BOOST_AUTO_TEST_CASE(two_component_jet) {
    auto f = math::two_component(0.1, 0.3, 100.0, 10.0);

    // Inside core
    BOOST_CHECK_EQUAL(f(0, 0.05), 100.0);
    // At core boundary (theta <= theta_c)
    BOOST_CHECK_EQUAL(f(0, 0.1), 100.0);
    // In wing
    BOOST_CHECK_EQUAL(f(0, 0.2), 10.0);
    // At wing boundary (theta <= theta_w)
    BOOST_CHECK_EQUAL(f(0, 0.3), 10.0);
    // Outside
    BOOST_CHECK_EQUAL(f(0, 0.5), 0.0);
}

// ---------------------------------------------------------------------------
// MediumVariant dispatch
// ---------------------------------------------------------------------------

// medium_rho dispatches through a MediumVariant holding an ISM and returns n * m_p.
BOOST_AUTO_TEST_CASE(medium_variant_dispatch) {
    MediumVariant mv = ISM(1.0);
    BOOST_CHECK_CLOSE(medium_rho(mv, 0, 0, 1e15), 1.0 * con::mp, 1e-10);
}

// medium_rho dispatch works for every variant alternative: ISM, Wind, and function-backed Medium.
BOOST_AUTO_TEST_CASE(medium_variant_all_types) {
    // ISM
    MediumVariant mv1 = ISM(1.0);
    BOOST_CHECK_CLOSE(medium_rho(mv1, 0, 0, 1e15), 1.0 * con::mp, 1e-10);

    // Wind
    MediumVariant mv2 = Wind(0.1);
    BOOST_CHECK_GT(medium_rho(mv2, 0, 0, 1e15), 0.0);

    // Medium with custom function via evn::ISM factory
    auto rho_func = evn::ISM(1.0);
    MediumVariant mv3 = Medium(rho_func, true);
    BOOST_CHECK_CLOSE(medium_rho(mv3, 0, 0, 1e15), 1.0 * con::mp, 1e-10);
}

// ---------------------------------------------------------------------------
// Ejecta tests
// ---------------------------------------------------------------------------

// Ejecta forwards user-supplied lambdas: eps_k and Gamma0 reproduce the custom angular profiles.
BOOST_AUTO_TEST_CASE(ejecta_custom_functions) {
    auto custom_eps = [](Real /*phi*/, Real theta) -> Real { return theta < 0.1 ? 1e50 : 0.0; };
    auto custom_gamma = [](Real /*phi*/, Real theta) -> Real { return theta < 0.1 ? 200.0 : 1.0; };

    Ejecta ej(custom_eps, custom_gamma);

    BOOST_CHECK_CLOSE(ej.eps_k(0, 0.05), 1e50, 1e-10);
    BOOST_CHECK_EQUAL(ej.eps_k(0, 0.2), 0.0);
    BOOST_CHECK_CLOSE(ej.Gamma0(0, 0.05), 200.0, 1e-10);
    BOOST_CHECK_EQUAL(ej.Gamma0(0, 0.2), 1.0);
}

// A default-constructed Ejecta is empty: eps_k = 0 and Gamma0 = 1 at every angle.
BOOST_AUTO_TEST_CASE(ejecta_default_constructor) {
    Ejecta ej;

    // Default: eps_k = 0, Gamma0 = 1 everywhere
    BOOST_CHECK_EQUAL(ej.eps_k(0, 0), 0.0);
    BOOST_CHECK_EQUAL(ej.eps_k(1.0, 0.5), 0.0);
    BOOST_CHECK_EQUAL(ej.Gamma0(0, 0), 1.0);
    BOOST_CHECK_EQUAL(ej.Gamma0(1.0, 0.5), 1.0);
}

// ---------------------------------------------------------------------------
// evn:: factory functions
// ---------------------------------------------------------------------------

// The evn::ISM factory returns a density function giving n * m_p independent of position.
BOOST_AUTO_TEST_CASE(evn_ism_factory) {
    auto rho = evn::ISM(1.0);
    Real expected = 1.0 * con::mp;

    // Returns correct density at any position
    BOOST_CHECK_CLOSE(rho(0, 0, 1e15), expected, 1e-10);
    BOOST_CHECK_CLOSE(rho(1.0, 0.5, 1e20), expected, 1e-10);
}

// The evn::wind factory density is positive everywhere and decreases with radius.
BOOST_AUTO_TEST_CASE(evn_wind_factory) {
    auto rho = evn::wind(1.0);

    // Density should be positive at any finite radius
    BOOST_CHECK_GT(rho(0, 0, 1e15), 0.0);
    BOOST_CHECK_GT(rho(0, 0, 1e18), 0.0);

    // Density decreases with radius (for wind profile)
    Real rho_near = rho(0, 0, 1e15);
    Real rho_far = rho(0, 0, 1e18);
    BOOST_CHECK_GT(rho_near, rho_far);
}

// ---------------------------------------------------------------------------
// Quantitative jet energetics
// ---------------------------------------------------------------------------

// Angle-integrated Gaussian jet energy (Simpson's rule) matches the analytic small-angle result
// E_iso/2 * theta_c^2 * [1 - exp(-(pi/2)^2/(2*theta_c^2))] within 1% (measured correction: 3.3e-3).
BOOST_AUTO_TEST_CASE(gaussian_jet_total_energy) {
    // GaussianJet profile is eps_k(theta) = E_iso/(4*pi) * exp(-theta^2/(2*theta_c^2)).
    // Total energy E = 2*pi * integral_0^{pi/2} eps_k(theta) sin(theta) dtheta.
    // With sin(theta) ~ theta (theta_c = 0.1 concentrates the integrand at
    // small angles) the analytic value is
    //   E = E_iso/2 * theta_c^2 * [1 - exp(-(pi/2)^2/(2*theta_c^2))].
    // The sin(theta) vs theta correction is -theta_c^2/3 ~ -0.33% (measured
    // relative difference: 3.3e-3), well inside the 1% tolerance.
    const Real theta_c = 0.1;
    const Real E_iso = 1e52;
    GaussianJet jet(theta_c, E_iso, 300.0);

    // Simpson's rule over [0, pi/2] with 200 intervals
    constexpr int N = 200;
    const Real h = (con::pi / 2) / N;
    Real sum = jet.eps_k(0, 0.0) * std::sin(0.0) + jet.eps_k(0, con::pi / 2) * std::sin(con::pi / 2);
    for (int i = 1; i < N; ++i) {
        const Real theta = i * h;
        sum += (i % 2 == 1 ? 4.0 : 2.0) * jet.eps_k(0, theta) * std::sin(theta);
    }
    const Real E_total = 2 * con::pi * sum * h / 3.0;

    const Real cut = 1.0 - std::exp(-(con::pi / 2) * (con::pi / 2) / (2 * theta_c * theta_c));
    const Real E_analytic = E_iso / 2.0 * theta_c * theta_c * cut;

    BOOST_CHECK_CLOSE(E_total, E_analytic, 1.0);
}

// Wing ratio eps_k(5*theta_c)/eps_k(10*theta_c) equals (1 + 10^k)/(1 + 5^k) within 0.5% and approaches
// the pure power-law value 2^k_e within 4% (measured deviations pinned in the body comment).
BOOST_AUTO_TEST_CASE(power_law_jet_wing_slope) {
    // PowerLawJet profile is eps_k(theta) = E_iso/(4*pi) / (1 + (theta/theta_c)^k_e),
    // so the wing ratio eps_k(5*theta_c)/eps_k(10*theta_c) equals
    // (1 + 10^k)/(1 + 5^k) exactly, approaching the pure power-law value 2^k
    // for large k_e (measured: exact form matches to machine precision;
    // 2^k differs by 2.9% at k_e=2 and 0.006% at k_e=6).
    const Real theta_c = 0.1;
    for (Real k_e : {2.0, 6.0}) {
        PowerLawJet jet(theta_c, 1e52, 300.0, k_e, 4.0);
        const Real ratio = jet.eps_k(0, 5 * theta_c) / jet.eps_k(0, 10 * theta_c);

        const Real exact_form = (1.0 + std::pow(10.0, k_e)) / (1.0 + std::pow(5.0, k_e));
        BOOST_CHECK_CLOSE(ratio, exact_form, 0.5);

        // Pure power-law prediction 2^k_e holds within a few percent
        BOOST_CHECK_CLOSE(ratio, std::pow(2.0, k_e), 4.0);
    }
}

BOOST_AUTO_TEST_SUITE_END()
