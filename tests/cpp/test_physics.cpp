#include <boost/test/unit_test.hpp>
#include <cmath>

#include "core/physics.h"
#include "util/macros.h"

BOOST_AUTO_TEST_SUITE(Physics)

// ========================== relativistic ==========================

BOOST_AUTO_TEST_CASE(gamma_to_beta_rest) {
    // gamma = 1 is at rest: beta should be exactly 0
    Real beta = physics::relativistic::gamma_to_beta(1.0);
    BOOST_CHECK_SMALL(beta, 1e-15);
}

BOOST_AUTO_TEST_CASE(gamma_to_beta_mildly_relativistic) {
    // gamma = 2 -> beta = sqrt(3)/2
    Real beta = physics::relativistic::gamma_to_beta(2.0);
    BOOST_CHECK_CLOSE(beta, std::sqrt(3.0) / 2.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(gamma_to_beta_ultra_relativistic) {
    // gamma = 1000 -> beta = sqrt(1 - 1/gamma^2) ~ 1 - 5e-7
    Real beta = physics::relativistic::gamma_to_beta(1000.0);
    Real expected = std::sqrt(1.0 - 1.0 / (1000.0 * 1000.0));
    BOOST_CHECK_CLOSE(beta, expected, 1e-10);
    // Also verify it is very close to 1
    BOOST_CHECK_GT(beta, 1.0 - 1e-5);
    BOOST_CHECK_LT(beta, 1.0);
}

BOOST_AUTO_TEST_CASE(gamma_to_beta_monotonic) {
    // Larger gamma should always give larger beta
    Real prev = physics::relativistic::gamma_to_beta(1.0);
    Real gammas[] = {1.01, 1.1, 2.0, 10.0, 100.0, 1000.0, 1e6};
    for (Real g : gammas) {
        Real beta = physics::relativistic::gamma_to_beta(g);
        BOOST_CHECK_GT(beta, prev);
        prev = beta;
    }
}

BOOST_AUTO_TEST_CASE(gamma_to_beta_near_one) {
    // gamma = 1 + 1e-10: should not produce NaN, and beta should be small but positive
    Real gamma = 1.0 + 1e-10;
    Real beta = physics::relativistic::gamma_to_beta(gamma);
    BOOST_CHECK(std::isfinite(beta));
    BOOST_CHECK_GE(beta, 0.0);
    BOOST_CHECK_LT(beta, 1.0);
}

BOOST_AUTO_TEST_CASE(gamma_to_beta_huge) {
    // gamma = 1e15: extreme ultrarelativistic limit
    Real gamma = 1e15;
    Real beta = physics::relativistic::gamma_to_beta(gamma);
    BOOST_CHECK(std::isfinite(beta));
    // beta should be extremely close to 1
    BOOST_CHECK_GT(beta, 1.0 - 1e-10);
    BOOST_CHECK_LE(beta, 1.0);
}

// ========================== thermo ==========================

BOOST_AUTO_TEST_CASE(adiabatic_idx_non_relativistic) {
    // gamma = 1 -> adiabatic index = 4/3 + 1/3 = 5/3
    Real idx = physics::thermo::adiabatic_idx(1.0);
    BOOST_CHECK_CLOSE(idx, 5.0 / 3.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(adiabatic_idx_ultra_relativistic) {
    // gamma = 1000 -> adiabatic index ~ 4/3 + 1/3000 ~ 4/3
    Real idx = physics::thermo::adiabatic_idx(1000.0);
    Real expected = 4.0 / 3.0 + 1.0 / 3000.0;
    BOOST_CHECK_CLOSE(idx, expected, 1e-10);
    // Should be very close to 4/3
    BOOST_CHECK_CLOSE(idx, 4.0 / 3.0, 0.1); // within 0.1%
}

BOOST_AUTO_TEST_CASE(adiabatic_idx_range) {
    // For any gamma >= 1, the adiabatic index should be in [4/3, 5/3]
    Real gammas[] = {1.0, 1.001, 1.5, 2.0, 10.0, 100.0, 1000.0, 1e6, 1e15};
    for (Real g : gammas) {
        Real idx = physics::thermo::adiabatic_idx(g);
        BOOST_CHECK_GE(idx, 4.0 / 3.0 - 1e-15);
        BOOST_CHECK_LE(idx, 5.0 / 3.0 + 1e-15);
    }
}

BOOST_AUTO_TEST_CASE(adiabatic_idx_near_one) {
    // gamma = 1 + 1e-12: should approach 5/3 without numerical issues
    Real gamma = 1.0 + 1e-12;
    Real idx = physics::thermo::adiabatic_idx(gamma);
    BOOST_CHECK(std::isfinite(idx));
    BOOST_CHECK_CLOSE(idx, 5.0 / 3.0, 1e-6); // very close to 5/3
}

// ========================== scales ==========================

// Typical GRB parameters used across scale tests
namespace {
    const Real E_iso = 1e52 * unit::erg;
    const Real n_ism = 1.0;
    const Real Gamma0 = 300.0;
    const Real engine_dura = 1.0 * unit::sec;
} // anonymous namespace

BOOST_AUTO_TEST_CASE(sedov_length_known) {
    // Sedov length = cbrt(E_iso / (4*pi/3 * n_ism * mp * c^2))
    Real expected = std::cbrt(E_iso / (4.0 * con::pi / 3.0 * n_ism * con::mp * con::c2));
    Real result = physics::scales::sedov_length(E_iso, n_ism);
    BOOST_CHECK_CLOSE(result, expected, 1e-10);
    // Should be positive and finite
    BOOST_CHECK_GT(result, 0.0);
    BOOST_CHECK(std::isfinite(result));
}

BOOST_AUTO_TEST_CASE(thin_shell_dec_radius) {
    // R_thin = cbrt(3 * E_iso / (4 * pi * mp * c^2 * n_ism * Gamma0^2))
    Real expected = std::cbrt(3.0 * E_iso / (4.0 * con::pi * con::mp * con::c2 * n_ism * Gamma0 * Gamma0));
    Real result = physics::scales::thin_shell_dec_radius(E_iso, n_ism, Gamma0);
    BOOST_CHECK_CLOSE(result, expected, 1e-10);
    BOOST_CHECK_GT(result, 0.0);
    BOOST_CHECK(std::isfinite(result));
}

BOOST_AUTO_TEST_CASE(thick_shell_dec_radius) {
    // R_thick = (3 * E_iso * T * c / (4 * pi * n_ism * mp * c^2))^(1/4)
    Real expected =
        std::sqrt(std::sqrt(3.0 * E_iso * engine_dura / n_ism * con::c / (4.0 * con::pi * con::mp * con::c2)));
    Real result = physics::scales::thick_shell_dec_radius(E_iso, n_ism, engine_dura);
    BOOST_CHECK_CLOSE(result, expected, 1e-10);
    BOOST_CHECK_GT(result, 0.0);
    BOOST_CHECK(std::isfinite(result));
}

BOOST_AUTO_TEST_CASE(dec_radius_is_max) {
    // dec_radius should be max(thin_shell, thick_shell)
    Real r_dec = physics::scales::dec_radius(E_iso, n_ism, Gamma0, engine_dura);
    Real r_thin = physics::scales::thin_shell_dec_radius(E_iso, n_ism, Gamma0);
    Real r_thick = physics::scales::thick_shell_dec_radius(E_iso, n_ism, engine_dura);
    BOOST_CHECK_GE(r_dec, r_thin);
    BOOST_CHECK_GE(r_dec, r_thick);
    // Should equal exactly one of them
    Real expected_max = std::max(r_thin, r_thick);
    BOOST_CHECK_CLOSE(r_dec, expected_max, 1e-10);
}

BOOST_AUTO_TEST_CASE(shell_spreading_radius) {
    // R_spread = Gamma0^2 * c * T
    Real expected = Gamma0 * Gamma0 * con::c * engine_dura;
    Real result = physics::scales::shell_spreading_radius(Gamma0, engine_dura);
    BOOST_CHECK_CLOSE(result, expected, 1e-10);
    BOOST_CHECK_GT(result, 0.0);
    BOOST_CHECK(std::isfinite(result));
}

BOOST_AUTO_TEST_CASE(RS_crossing_radius) {
    // R_RS_cross = (l^3 * c * T)^(1/4)
    Real l = physics::scales::sedov_length(E_iso, n_ism);
    Real expected = std::sqrt(std::sqrt(l * l * l * con::c * engine_dura));
    Real result = physics::scales::RS_crossing_radius(E_iso, n_ism, engine_dura);
    BOOST_CHECK_CLOSE(result, expected, 1e-10);
    BOOST_CHECK_GT(result, 0.0);
    BOOST_CHECK(std::isfinite(result));
}

BOOST_AUTO_TEST_CASE(RS_transition_radius) {
    // R_RS_trans = l^1.5 / sqrt(c * T) / Gamma0^2
    Real l = physics::scales::sedov_length(E_iso, n_ism);
    Real expected = std::pow(l, 1.5) / std::sqrt(con::c * engine_dura) / Gamma0 / Gamma0;
    Real result = physics::scales::RS_transition_radius(E_iso, n_ism, Gamma0, engine_dura);
    BOOST_CHECK_CLOSE(result, expected, 1e-10);
    BOOST_CHECK_GT(result, 0.0);
    BOOST_CHECK(std::isfinite(result));
}

BOOST_AUTO_TEST_CASE(shell_thickness_param) {
    // xi = sqrt(l / (c * T)) * Gamma0^(-4/3)
    Real l = physics::scales::sedov_length(E_iso, n_ism);
    Real shell_width = con::c * engine_dura;
    Real expected = std::sqrt(l / shell_width) * std::pow(Gamma0, -4.0 / 3.0);
    Real result = physics::scales::shell_thickness_param(E_iso, n_ism, Gamma0, engine_dura);
    BOOST_CHECK_CLOSE(result, expected, 1e-10);
    BOOST_CHECK_GT(result, 0.0);
    BOOST_CHECK(std::isfinite(result));
}

BOOST_AUTO_TEST_CASE(calc_engine_duration_roundtrip) {
    // Compute xi from known parameters, then recover engine_dura from calc_engine_duration
    Real xi = physics::scales::shell_thickness_param(E_iso, n_ism, Gamma0, engine_dura);
    Real recovered_T = physics::scales::calc_engine_duration(E_iso, n_ism, Gamma0, xi);
    BOOST_CHECK_CLOSE(recovered_T, engine_dura, 1e-8);
}

BOOST_AUTO_TEST_CASE(sedov_length_extreme_params) {
    // Very large E_iso and very small n_ism: should still produce finite positive result
    Real E_extreme = 1e60 * unit::erg;
    Real n_extreme = 1e-10;
    Real result = physics::scales::sedov_length(E_extreme, n_extreme);
    BOOST_CHECK(std::isfinite(result));
    BOOST_CHECK_GT(result, 0.0);
    // Larger energy and smaller density should give larger Sedov length than fiducial
    Real result_fid = physics::scales::sedov_length(E_iso, n_ism);
    BOOST_CHECK_GT(result, result_fid);
}

BOOST_AUTO_TEST_CASE(dec_radius_zero_Gamma) {
    // Gamma0 = 1 is a degenerate case (non-relativistic outflow)
    // thin_shell_dec_radius with Gamma0=1 should still give a finite positive result
    Real r_thin = physics::scales::thin_shell_dec_radius(E_iso, n_ism, 1.0);
    Real r_thick = physics::scales::thick_shell_dec_radius(E_iso, n_ism, engine_dura);
    Real r_dec = physics::scales::dec_radius(E_iso, n_ism, 1.0, engine_dura);
    BOOST_CHECK(std::isfinite(r_thin));
    BOOST_CHECK(std::isfinite(r_thick));
    BOOST_CHECK(std::isfinite(r_dec));
    BOOST_CHECK_GT(r_dec, 0.0);
    // With Gamma0=1, thin_shell radius should be larger than with Gamma0=300
    Real r_thin_300 = physics::scales::thin_shell_dec_radius(E_iso, n_ism, 300.0);
    BOOST_CHECK_GT(r_thin, r_thin_300);
}

BOOST_AUTO_TEST_SUITE_END()
