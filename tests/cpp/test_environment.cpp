#include <boost/test/unit_test.hpp>
#include <cmath>

#include "environment/jet.h"
#include "environment/medium.h"
#include "util/macros.h"

BOOST_AUTO_TEST_SUITE(Environment)

// ---------------------------------------------------------------------------
// ISM tests
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ism_constant_density) {
    ISM ism(1.0);
    Real rho0 = ism.rho(0, 0, 1e15);
    // Density is the same regardless of phi, theta, r
    BOOST_CHECK_EQUAL(ism.rho(0.5, 0.3, 1e10), rho0);
    BOOST_CHECK_EQUAL(ism.rho(1.0, 1.0, 1e20), rho0);
    BOOST_CHECK_EQUAL(ism.rho(3.14, 1.57, 1e5), rho0);
}

BOOST_AUTO_TEST_CASE(ism_density_value) {
    ISM ism(1.0);
    Real expected = 1.0 * con::mp;
    BOOST_CHECK_CLOSE(ism.rho(0, 0, 1e15), expected, 1e-10);
}

BOOST_AUTO_TEST_CASE(ism_mass_cubic) {
    ISM ism(1.0);
    Real r = 1e16;
    Real expected = (1.0 * con::mp) * r * r * r / 3.0;
    BOOST_CHECK_CLOSE(ism.mass(r), expected, 1e-10);
}

BOOST_AUTO_TEST_CASE(ism_zero_density) {
    ISM ism(0.0);
    BOOST_CHECK_EQUAL(ism.rho(0, 0, 1e15), 0.0);
    BOOST_CHECK_EQUAL(ism.mass(1e16), 0.0);
}

// ---------------------------------------------------------------------------
// Wind tests
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(wind_density_r_squared) {
    // For large r with no ISM floor, density ~ A/r^2
    Real A_star = 1.0;
    Wind wind(A_star, 0, con::inf);
    Real A = A_star * 5e11 * unit::g / unit::cm;

    Real r = 1e18;
    Real expected = A / (r * r); // r02 = 0 when n0 = inf
    BOOST_CHECK_CLOSE(wind.rho(0, 0, r), expected, 1e-6);
}

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

BOOST_AUTO_TEST_CASE(wind_mass_linear) {
    // Wind(A_star, 0) with n0=inf gives r02=0, so mass ~ A*r for large r
    Real A_star = 1.0;
    Wind wind(A_star, 0, con::inf);
    Real A = A_star * 5e11 * unit::g / unit::cm;

    Real r = 1e18;
    Real expected = A * r; // no ISM term, r02=0
    BOOST_CHECK_CLOSE(wind.mass(r), expected, 1e-10);
}

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

BOOST_AUTO_TEST_CASE(wind_zero_A_star) {
    // Wind(0, n_ism) density equals ISM floor only
    Real n_ism = 1.0;
    Wind wind(0.0, n_ism);
    Real expected = n_ism * con::mp;
    BOOST_CHECK_CLOSE(wind.rho(0, 0, 1e15), expected, 1e-10);
    BOOST_CHECK_CLOSE(wind.rho(0, 0, 1e20), expected, 1e-10);
}

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

BOOST_AUTO_TEST_CASE(tophat_jet_inside_core) {
    Real theta_c = 0.1;
    Real E_iso = 1e52;
    Real Gamma0 = 300.0;
    TophatJet jet(theta_c, E_iso, Gamma0);

    // Inside core: eps_k > 0 and Gamma0 > 1
    BOOST_CHECK_GT(jet.eps_k(0, 0.05), 0.0);
    BOOST_CHECK_GT(jet.Gamma0(0, 0.05), 1.0);
}

BOOST_AUTO_TEST_CASE(tophat_jet_outside_core) {
    Real theta_c = 0.1;
    TophatJet jet(theta_c, 1e52, 300.0);

    // Outside core: eps_k = 0, Gamma0 = 1
    BOOST_CHECK_EQUAL(jet.eps_k(0, 0.2), 0.0);
    BOOST_CHECK_EQUAL(jet.Gamma0(0, 0.2), 1.0);
}

BOOST_AUTO_TEST_CASE(tophat_jet_energy_value) {
    Real E_iso = 1e52;
    TophatJet jet(0.1, E_iso, 300.0);

    Real expected = E_iso / (4 * con::pi);
    BOOST_CHECK_CLOSE(jet.eps_k(0, 0.0), expected, 1e-10);
}

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

BOOST_AUTO_TEST_CASE(tophat_jet_theta_zero) {
    TophatJet jet(0.1, 1e52, 300.0);
    // theta = 0 is the center, always inside
    BOOST_CHECK_CLOSE(jet.eps_k(0, 0.0), 1e52 / (4 * con::pi), 1e-10);
    BOOST_CHECK_CLOSE(jet.Gamma0(0, 0.0), 300.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(tophat_jet_theta_pi_half) {
    TophatJet jet(0.1, 1e52, 300.0);
    // theta = pi/2 is far outside core
    BOOST_CHECK_EQUAL(jet.eps_k(0, con::pi / 2), 0.0);
    BOOST_CHECK_EQUAL(jet.Gamma0(0, con::pi / 2), 1.0);
}

// ---------------------------------------------------------------------------
// GaussianJet tests
// ---------------------------------------------------------------------------

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

BOOST_AUTO_TEST_CASE(gaussian_jet_theta_zero) {
    Real E_iso = 1e53;
    Real Gamma0 = 500.0;
    GaussianJet jet(0.05, E_iso, Gamma0);

    // Peak values at theta=0
    BOOST_CHECK_CLOSE(jet.eps_k(0, 0.0), E_iso / (4 * con::pi), 1e-10);
    BOOST_CHECK_CLOSE(jet.Gamma0(0, 0.0), Gamma0, 1e-10);
}

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

BOOST_AUTO_TEST_CASE(power_law_jet_peak) {
    Real E_iso = 1e52;
    Real Gamma0 = 300.0;
    PowerLawJet jet(0.1, E_iso, Gamma0, 4.0, 4.0);

    // At theta=0, fast_pow(0/theta_c, k) = 0, so eps_k = eps_k_ / (1+0) = eps_k_
    Real expected_eps = E_iso / (4 * con::pi);
    BOOST_CHECK_CLOSE(jet.eps_k(0, 0.0), expected_eps, 1e-10);
    BOOST_CHECK_CLOSE(jet.Gamma0(0, 0.0), Gamma0, 1e-10);
}

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

BOOST_AUTO_TEST_CASE(medium_variant_dispatch) {
    MediumVariant mv = ISM(1.0);
    BOOST_CHECK_CLOSE(medium_rho(mv, 0, 0, 1e15), 1.0 * con::mp, 1e-10);
}

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

BOOST_AUTO_TEST_CASE(ejecta_custom_functions) {
    auto custom_eps = [](Real /*phi*/, Real theta) -> Real { return theta < 0.1 ? 1e50 : 0.0; };
    auto custom_gamma = [](Real /*phi*/, Real theta) -> Real { return theta < 0.1 ? 200.0 : 1.0; };

    Ejecta ej(custom_eps, custom_gamma);

    BOOST_CHECK_CLOSE(ej.eps_k(0, 0.05), 1e50, 1e-10);
    BOOST_CHECK_EQUAL(ej.eps_k(0, 0.2), 0.0);
    BOOST_CHECK_CLOSE(ej.Gamma0(0, 0.05), 200.0, 1e-10);
    BOOST_CHECK_EQUAL(ej.Gamma0(0, 0.2), 1.0);
}

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

BOOST_AUTO_TEST_CASE(evn_ism_factory) {
    auto rho = evn::ISM(1.0);
    Real expected = 1.0 * con::mp;

    // Returns correct density at any position
    BOOST_CHECK_CLOSE(rho(0, 0, 1e15), expected, 1e-10);
    BOOST_CHECK_CLOSE(rho(1.0, 0.5, 1e20), expected, 1e-10);
}

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

BOOST_AUTO_TEST_SUITE_END()
