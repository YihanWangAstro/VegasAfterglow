#include <boost/test/unit_test.hpp>
#include <cmath>

#include "dynamics/shock-physics.h"
#include "util/macros.h"

BOOST_AUTO_TEST_SUITE(ShockPhysics)

// ---------------------------------------------------------------------------
// 1. compute_downstr_4vel: non-magnetized (sigma = 0)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_downstr_4vel_non_magnetized) {
    Real gamma_rel = 10.0;
    Real sigma = 0.0;
    Real u_down = compute_downstr_4vel(gamma_rel, sigma);
    // Must be positive and finite for a strong relativistic shock
    BOOST_CHECK(std::isfinite(u_down));
    BOOST_CHECK_GT(u_down, 0.0);
    // Non-magnetized ultrarelativistic: u_down should be order unity to ~few
    BOOST_CHECK_LT(u_down, gamma_rel);
}

// ---------------------------------------------------------------------------
// 2. compute_downstr_4vel: magnetized (sigma = 0.1)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_downstr_4vel_magnetized) {
    Real gamma_rel = 10.0;
    Real sigma = 0.1;
    Real u_down = compute_downstr_4vel(gamma_rel, sigma);
    BOOST_CHECK(std::isfinite(u_down));
    BOOST_CHECK_GT(u_down, 0.0);
    // Magnetization changes the downstream velocity compared to sigma=0
    Real u_down_0 = compute_downstr_4vel(gamma_rel, 0.0);
    BOOST_CHECK(std::fabs(u_down - u_down_0) > 1e-10);
}

// ---------------------------------------------------------------------------
// 3. compute_upstr_4vel: known values
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_upstr_4vel_values) {
    Real u_down = 1.5;
    Real gamma_rel = 5.0;
    Real u_up = compute_upstr_4vel(u_down, gamma_rel);
    // u_up should be larger than u_down for gamma_rel > 1
    BOOST_CHECK_GT(u_up, u_down);
    BOOST_CHECK(std::isfinite(u_up));
    // Verify the formula directly
    Real expected = std::sqrt((1 + u_down * u_down) * std::fabs(gamma_rel * gamma_rel - 1)) + u_down * gamma_rel;
    BOOST_CHECK_CLOSE(u_up, expected, 1e-12);
}

// ---------------------------------------------------------------------------
// 4. compute_4vel_jump: sigma=0, ultrarelativistic limit
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_4vel_jump_sigma0_limit) {
    // For ultrarelativistic unmagnetized shocks, the 4-velocity jump ratio -> ~4*gamma_rel
    // when u_down -> 0. For moderate gamma_rel, it should be close to 4*gamma_rel
    // but not exact because u_down is not exactly zero.
    Real gamma_rel = 100.0;
    Real ratio = compute_4vel_jump(gamma_rel, 0.0);
    BOOST_CHECK(std::isfinite(ratio));
    BOOST_CHECK_GT(ratio, 1.0);
}

// ---------------------------------------------------------------------------
// 5. compute_4vel_jump: non-relativistic (gamma_rel ~ 1.001)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_4vel_jump_non_relativistic) {
    Real gamma_rel = 1.001;
    Real ratio = compute_4vel_jump(gamma_rel, 0.0);
    BOOST_CHECK(std::isfinite(ratio));
    BOOST_CHECK_GT(ratio, 1.0);
    // For non-relativistic shocks, the compression ratio should be moderate
    BOOST_CHECK_LT(ratio, 100.0);
}

// ---------------------------------------------------------------------------
// 6. compute_sound_speed: always in (0, c) for gamma_rel > 1
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_sound_speed_range) {
    Real gammas[] = {1.001, 1.1, 2.0, 10.0, 100.0, 1000.0};
    for (Real g : gammas) {
        Real cs = compute_sound_speed(g);
        BOOST_CHECK_GT(cs, 0.0);
        BOOST_CHECK_LT(cs, con::c);
        BOOST_CHECK(std::isfinite(cs));
    }
}

// ---------------------------------------------------------------------------
// 7. compute_sound_speed: ultrarelativistic -> cs ~ c/sqrt(3)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_sound_speed_ultra_relativistic) {
    Real gamma_rel = 1000.0;
    Real cs = compute_sound_speed(gamma_rel);
    Real expected = con::c / std::sqrt(3.0);
    // Should be close to c/sqrt(3) in the ultrarelativistic limit
    BOOST_CHECK_CLOSE(cs, expected, 1.0); // within 1%
}

// ---------------------------------------------------------------------------
// 8. compute_effective_Gamma: Gamma=1 -> returns 1
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_effective_Gamma_at_rest) {
    Real adx = 5.0 / 3.0;
    Real result = compute_effective_Gamma(adx, 1.0);
    // (adx * 1 - adx + 1) / 1 = 1
    BOOST_CHECK_CLOSE(result, 1.0, 1e-12);
    // Also test with ad_idx = 4/3
    adx = 4.0 / 3.0;
    result = compute_effective_Gamma(adx, 1.0);
    BOOST_CHECK_CLOSE(result, 1.0, 1e-12);
}

// ---------------------------------------------------------------------------
// 9. compute_effective_Gamma_dGamma: matches finite difference
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_effective_Gamma_dGamma_finite_diff) {
    Real adx = 4.0 / 3.0;
    Real Gamma = 10.0;
    Real dG = 1e-6;
    Real f_plus = compute_effective_Gamma(adx, Gamma + dG);
    Real f_minus = compute_effective_Gamma(adx, Gamma);
    Real numerical_deriv = (f_plus - f_minus) / dG;
    Real analytic_deriv = compute_effective_Gamma_dGamma(adx, Gamma);
    BOOST_CHECK_CLOSE(numerical_deriv, analytic_deriv, 0.01); // 0.01% tolerance
}

// ---------------------------------------------------------------------------
// 10. compute_comv_weibel_B: eps_B=0 -> B=0
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_comv_weibel_B_zero) {
    Real B = compute_comv_weibel_B(0.0, 1e10);
    BOOST_CHECK_SMALL(B, 1e-30);
}

// ---------------------------------------------------------------------------
// 11. compute_comv_weibel_B: positive inputs -> positive B
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_comv_weibel_B_positive) {
    Real eps_B = 0.01;
    Real e_th = 1e5;
    Real B = compute_comv_weibel_B(eps_B, e_th);
    BOOST_CHECK_GT(B, 0.0);
    BOOST_CHECK(std::isfinite(B));
    // Verify formula: sqrt(8*pi*eps_B*e_th)
    Real expected = std::sqrt(8 * con::pi * eps_B * e_th);
    BOOST_CHECK_CLOSE(B, expected, 1e-12);
}

// ---------------------------------------------------------------------------
// 12. compute_dr_dt: beta=0 -> 0
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_dr_dt_beta_zero) {
    Real drdt = compute_dr_dt(0.0);
    BOOST_CHECK_SMALL(drdt, 1e-30);
}

// ---------------------------------------------------------------------------
// 13. compute_dr_dt: both overloads agree
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_dr_dt_overloads_agree) {
    Real Gamma = 10.0;
    Real beta = std::sqrt(1.0 - 1.0 / (Gamma * Gamma));
    Real u = Gamma * beta;
    Real drdt_beta = compute_dr_dt(beta);
    Real drdt_Gu = compute_dr_dt(Gamma, u);
    BOOST_CHECK_CLOSE(drdt_beta, drdt_Gu, 1e-8);
}

// ---------------------------------------------------------------------------
// 14. compute_dtheta_dt: positive spreading rate
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_dtheta_dt_positive) {
    Real theta_s = 0.1;
    Real theta = 0.1;
    Real Gamma = 10.0;
    Real beta = std::sqrt(1.0 - 1.0 / (Gamma * Gamma));
    Real r = 1e16;
    Real drdt = compute_dr_dt(beta);
    Real rate = compute_dtheta_dt(theta_s, theta, drdt, r, Gamma);
    BOOST_CHECK_GT(rate, 0.0);
    BOOST_CHECK(std::isfinite(rate));
}

// ---------------------------------------------------------------------------
// 15. compute_dtheta_dt: both overloads match
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_dtheta_dt_overloads_agree) {
    Real theta_s = 0.1;
    Real theta = 0.1;
    Real Gamma = 10.0;
    Real u2 = Gamma * Gamma - 1;
    Real u = std::sqrt(u2);
    Real beta = u / Gamma;
    Real r = 1e16;
    Real drdt = compute_dr_dt(beta);
    Real rate1 = compute_dtheta_dt(theta_s, theta, drdt, r, Gamma);
    Real rate2 = compute_dtheta_dt(theta_s, theta, drdt, r, Gamma, u, u2);
    BOOST_CHECK_CLOSE(rate1, rate2, 1e-12);
}

// ---------------------------------------------------------------------------
// 16. compute_dt_dt_comv: Gamma=1, beta=0 -> 1
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_dt_dt_comv_rest) {
    Real result = compute_dt_dt_comv(1.0, 0.0);
    BOOST_CHECK_CLOSE(result, 1.0, 1e-12);
}

// ---------------------------------------------------------------------------
// 17. compute_dt_dt_comv: always positive
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_dt_dt_comv_positive) {
    Real gammas[] = {1.0, 1.01, 2.0, 10.0, 100.0, 1000.0};
    for (Real g : gammas) {
        Real beta = (g > 1.0) ? std::sqrt(1.0 - 1.0 / (g * g)) : 0.0;
        Real result = compute_dt_dt_comv(g, beta);
        BOOST_CHECK_GT(result, 0.0);
        BOOST_CHECK(std::isfinite(result));
    }
}

// ---------------------------------------------------------------------------
// 18. compute_upstr_B: known sigma and rho
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_upstr_B_values) {
    Real rho = 1e-24;
    Real sigma = 0.1;
    Real B = compute_upstr_B(rho, sigma);
    BOOST_CHECK_GT(B, 0.0);
    BOOST_CHECK(std::isfinite(B));
    // Verify formula
    Real expected = std::sqrt(4 * con::pi * con::c2 * sigma * rho);
    BOOST_CHECK_CLOSE(B, expected, 1e-12);
}

// ---------------------------------------------------------------------------
// 19. compute_rel_Gamma: same frame -> rel_Gamma = 1
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_rel_Gamma_same_frame) {
    Real gammas[] = {1.0, 2.0, 10.0, 100.0, 1000.0};
    for (Real g : gammas) {
        Real rel = compute_rel_Gamma(g, g);
        BOOST_CHECK_CLOSE(rel, 1.0, 1e-6);
    }
}

// ---------------------------------------------------------------------------
// 20. compute_rel_Gamma: both overloads agree
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_rel_Gamma_overloads_agree) {
    Real g1 = 10.0;
    Real g2 = 5.0;
    Real b1 = std::sqrt(1.0 - 1.0 / (g1 * g1));
    Real b2 = std::sqrt(1.0 - 1.0 / (g2 * g2));
    Real rel1 = compute_rel_Gamma(g1, g2);
    Real rel2 = compute_rel_Gamma(g1, g2, b1, b2);
    BOOST_CHECK_CLOSE(rel1, rel2, 1e-8);
}

// ---------------------------------------------------------------------------
// 21. compute_rel_Gamma: symmetric
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_rel_Gamma_symmetric) {
    Real g1 = 10.0;
    Real g2 = 50.0;
    Real rel12 = compute_rel_Gamma(g1, g2);
    Real rel21 = compute_rel_Gamma(g2, g1);
    BOOST_CHECK_CLOSE(rel12, rel21, 1e-12);
}

// ---------------------------------------------------------------------------
// 22. compute_Gamma_from_relative: roundtrip
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_Gamma_from_relative_roundtrip) {
    Real g1 = 10.0;
    Real g2 = 50.0;
    Real gamma_rel = compute_rel_Gamma(g1, g2);
    // Given g2 (as gamma4) and gamma_rel, recover g1
    Real recovered = compute_Gamma_from_relative(g2, gamma_rel);
    BOOST_CHECK_CLOSE(recovered, g1, 0.1); // 0.1% tolerance
}

// ---------------------------------------------------------------------------
// 23. compute_shock_heating_rate: Gamma_rel=1 -> 0
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_shock_heating_rate_no_jump) {
    Real rate = compute_shock_heating_rate(1.0, 1e20);
    BOOST_CHECK_SMALL(rate, 1e-20);
}

// ---------------------------------------------------------------------------
// 24. compute_shock_heating_rate: positive for Gamma_rel > 1
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_shock_heating_rate_positive) {
    Real rate = compute_shock_heating_rate(10.0, 1e20);
    BOOST_CHECK_GT(rate, 0.0);
    BOOST_CHECK(std::isfinite(rate));
    // Check formula: mdot * (Gamma_rel - 1) * c2
    Real expected = 1e20 * (10.0 - 1.0) * con::c2;
    BOOST_CHECK_CLOSE(rate, expected, 1e-12);
}

// ---------------------------------------------------------------------------
// 25. compute_adiabatic_cooling_rate: expanding (drdt > 0) -> negative
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_adiabatic_cooling_rate_sign) {
    Real ad_idx = 4.0 / 3.0;
    Real r = 1e16;
    Real Gamma = 10.0;
    Real u = 1e10;       // internal energy density
    Real drdt = 1e8;     // positive: expanding
    Real dGammadt = 0.0; // no acceleration
    Real rate = compute_adiabatic_cooling_rate(ad_idx, r, Gamma, u, drdt, dGammadt);
    // -(ad_idx - 1) * (3*drdt/r - 0) * u should be negative
    BOOST_CHECK_LT(rate, 0.0);
    BOOST_CHECK(std::isfinite(rate));
}

// ---------------------------------------------------------------------------
// 26. compute_adiabatic_cooling_rate2: expanding -> negative
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_adiabatic_cooling_rate2_sign) {
    Real ad_idx = 4.0 / 3.0;
    Real r = 1e16;
    Real x = 0.1; // some positive opening angle
    Real u = 1e10;
    Real drdt = 1e8;
    Real dxdt = 0.0;
    Real rate = compute_adiabatic_cooling_rate2(ad_idx, r, x, u, drdt, dxdt);
    // -(ad_idx - 1) * (2*drdt/r + 0) * u should be negative
    BOOST_CHECK_LT(rate, 0.0);
    BOOST_CHECK(std::isfinite(rate));
}

// ---------------------------------------------------------------------------
// 27. compute_shell_spreading_rate: positive for physical inputs
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_shell_spreading_rate_positive) {
    Real gamma_rel = 10.0;
    Real dtdt_comv = 0.5;
    Real rate = compute_shell_spreading_rate(gamma_rel, dtdt_comv);
    BOOST_CHECK_GT(rate, 0.0);
    BOOST_CHECK(std::isfinite(rate));
}

// ---------------------------------------------------------------------------
// 28. compute_Gamma_therm: mass=0 -> 1
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_Gamma_therm_zero_mass) {
    Real result = compute_Gamma_therm(1e30, 0.0);
    BOOST_CHECK_CLOSE(result, 1.0, 1e-12);
}

// ---------------------------------------------------------------------------
// 29. compute_Gamma_therm: U_th=0 -> 1
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_Gamma_therm_zero_energy) {
    Real result = compute_Gamma_therm(0.0, 1e20);
    BOOST_CHECK_CLOSE(result, 1.0, 1e-12);
}

// ---------------------------------------------------------------------------
// 30. compute_Gamma_therm: U_th > 0 -> > 1
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_Gamma_therm_positive_energy) {
    Real mass = 1e20;
    Real U_th = 1e25;
    Real result = compute_Gamma_therm(U_th, mass);
    BOOST_CHECK_GT(result, 1.0);
    BOOST_CHECK(std::isfinite(result));
    // Verify formula: U_th / (mass * c2) + 1
    Real expected = U_th / (mass * con::c2) + 1;
    BOOST_CHECK_CLOSE(result, expected, 1e-12);
}

// ---------------------------------------------------------------------------
// 31. compute_Gamma_therm: limiter=true, small U_th -> 1
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_Gamma_therm_limiter) {
    Real mass = 1e20;
    // gamma_therm_cut = 1 + 1e-6, so U_th that gives Gamma_th < gamma_therm_cut
    // Gamma_th = U_th / (mass * c2) + 1 < 1 + 1e-6
    // U_th < 1e-6 * mass * c2
    Real U_th = 1e-8 * mass * con::c2; // gives Gamma_th = 1 + 1e-8 < gamma_therm_cut
    Real result_no_limiter = compute_Gamma_therm(U_th, mass, false);
    BOOST_CHECK_GT(result_no_limiter, 1.0);

    Real result_limiter = compute_Gamma_therm(U_th, mass, true);
    BOOST_CHECK_CLOSE(result_limiter, 1.0, 1e-12);
}

// ---------------------------------------------------------------------------
// 32. compute_compression: known Gamma values
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_compression_values) {
    Real Gamma_up = 100.0;
    Real Gamma_down = 10.0;
    Real sigma = 0.0;
    Real comp = compute_compression(Gamma_up, Gamma_down, sigma);
    BOOST_CHECK(std::isfinite(comp));
    BOOST_CHECK_GT(comp, 1.0);
}

// ---------------------------------------------------------------------------
// 33. compute_downstr_B: Weibel + compressed upstream B
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_downstr_B_values) {
    Real eps_B = 0.01;
    Real rho_upstr = 1e-24;
    Real B_upstr = 1e-6;
    Real Gamma_th = 10.0;
    Real comp_ratio = 4.0;
    Real B_down = compute_downstr_B(eps_B, rho_upstr, B_upstr, Gamma_th, comp_ratio);
    BOOST_CHECK_GT(B_down, 0.0);
    BOOST_CHECK(std::isfinite(B_down));
    // Should include both Weibel and compressed upstream contributions
    Real rho_downstr = rho_upstr * comp_ratio;
    Real e_th = (Gamma_th - 1) * rho_downstr * con::c2;
    Real B_weibel = compute_comv_weibel_B(eps_B, e_th);
    Real B_compressed = B_upstr * comp_ratio;
    Real expected = B_weibel + B_compressed;
    BOOST_CHECK_CLOSE(B_down, expected, 1e-12);
}

// ---------------------------------------------------------------------------
// 34. compute_radiative_efficiency: slow cooling (gamma_m < gamma_c)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_radiative_efficiency_slow_cooling) {
    // Use parameters that ensure gamma_m < gamma_c (slow cooling)
    // Large t_comv and small eps_B push gamma_c high, small Gamma_th keeps gamma_m low
    RadParams rad;
    rad.eps_e = 0.1;
    rad.eps_B = 1e-4;
    rad.p = 2.3;
    rad.xi_e = 1.0;
    Real t_comv = 1e10;
    Real Gamma_th = 2.0;
    Real u = 1e-10;
    Real eff = compute_radiative_efficiency(t_comv, Gamma_th, u, rad);
    // Slow cooling: efficiency < eps_e
    BOOST_CHECK_GE(eff, 0.0);
    BOOST_CHECK_LE(eff, rad.eps_e);
    BOOST_CHECK(std::isfinite(eff));
}

// ---------------------------------------------------------------------------
// 35. compute_radiative_efficiency: fast cooling (gamma_m > gamma_c) -> eps_e
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_radiative_efficiency_fast_cooling) {
    // Use parameters that ensure gamma_m > gamma_c (fast cooling)
    // Very high Gamma_th makes gamma_m large, large eps_B and u make gamma_c small
    RadParams rad;
    rad.eps_e = 0.1;
    rad.eps_B = 0.1;
    rad.p = 2.3;
    rad.xi_e = 1.0;
    Real t_comv = 1e5;
    Real Gamma_th = 1000.0;
    Real u = 1e10;
    Real eff = compute_radiative_efficiency(t_comv, Gamma_th, u, rad);
    BOOST_CHECK_CLOSE(eff, rad.eps_e, 1e-10);
}

// ---------------------------------------------------------------------------
// 36. simpson_logspace: power law integral
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(simpson_logspace_power_law) {
    // Integrate f(u) = exp(3u) from u_min = ln(r) - 18 to u_max = ln(r).
    // Analytic: integral of e^{3u} du = e^{3u}/3 evaluated at bounds
    // = r^3/3 - r^3*e^{-54}/3. Since e^{-54} ~ 0, result ~ r^3/3.
    Real r = 1e10;
    Real result = simpson_logspace([](Real u) { return std::exp(3 * u); }, r);
    Real expected = r * r * r / 3.0;
    BOOST_CHECK_CLOSE(result, expected, 5.0); // 5% tolerance (limited by log-space quadrature)
}

// ---------------------------------------------------------------------------
// 37. enclosed_mass: uniform density ISM
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(enclosed_mass_ism) {
    Real r = 1e15;
    Real rho_const = 1e-24;
    Real mass = enclosed_mass([=](Real /*r_*/) { return rho_const; }, r);
    Real expected = rho_const * r * r * r / 3.0;
    BOOST_CHECK_CLOSE(mass, expected, 5.0); // 5% tolerance (limited by log-space quadrature)
}

// ---------------------------------------------------------------------------
// 38. compute_dr_dt: beta very close to c
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_dr_dt_near_c) {
    Real beta = 1.0 - 1e-15;
    Real drdt = compute_dr_dt(beta);
    // beta*c / (1 - beta) = (1 - 1e-15) * c / 1e-15 ~ c * 1e15
    BOOST_CHECK(std::isfinite(drdt));
    BOOST_CHECK_GT(drdt, 0.0);
    // Should be very large
    BOOST_CHECK_GT(drdt, 1e14 * con::c);
}

// ---------------------------------------------------------------------------
// 39. compute_rel_Gamma: one frame at rest (gamma = 1)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_rel_Gamma_one_at_rest) {
    Real g_moving = 100.0;
    // If one frame is at rest (gamma=1), relative Gamma equals the other
    Real rel = compute_rel_Gamma(1.0, g_moving);
    BOOST_CHECK_CLOSE(rel, g_moving, 1e-6);
    // Symmetric check
    Real rel2 = compute_rel_Gamma(g_moving, 1.0);
    BOOST_CHECK_CLOSE(rel2, g_moving, 1e-6);
}

// ---------------------------------------------------------------------------
// 40. compute_Gamma_from_relative: gamma_rel=1 (no relative motion)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_Gamma_from_relative_edge) {
    Real gamma4 = 10.0;
    Real gamma_rel = 1.0;
    // No relative motion: the derived Gamma should equal gamma4
    Real result = compute_Gamma_from_relative(gamma4, gamma_rel);
    BOOST_CHECK_CLOSE(result, gamma4, 0.1);
}

// ---------------------------------------------------------------------------
// 41. compute_4vel_jump: u_down=0 fallback path
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_4vel_jump_u_down_zero) {
    // gamma_rel very close to 1 where u_down could be 0
    Real gamma_rel = 1.0;
    Real ratio = compute_4vel_jump(gamma_rel, 0.0);
    // When u_down == 0, fallback is 4 * gamma_rel
    Real expected = 4.0 * gamma_rel;
    BOOST_CHECK_CLOSE(ratio, expected, 1e-10);
}

// ---------------------------------------------------------------------------
// 42. compute_radiative_efficiency: eps_B approaching 0
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_radiative_efficiency_zero_B) {
    // Very small eps_B pushes gamma_c -> infinity (slow cooling regime)
    // gamma_m / gamma_c -> 0, so efficiency should be very small or zero
    RadParams rad;
    rad.eps_e = 0.1;
    rad.eps_B = 1e-30;
    rad.p = 2.3;
    rad.xi_e = 1.0;
    Real t_comv = 1e6;
    Real Gamma_th = 10.0;
    Real u = 1e-5;
    Real eff = compute_radiative_efficiency(t_comv, Gamma_th, u, rad);
    // With extremely small eps_B, gamma_c is huge, so gamma_m/gamma_c < 1e-2 -> returns 0
    BOOST_CHECK_SMALL(eff, 1e-10);
}

// ---------------------------------------------------------------------------
// 43. compute_radiative_efficiency: p=2 boundary case
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(compute_radiative_efficiency_p_equals_2) {
    // p = 2: the exponent (p-2) = 0 in the slow-cooling branch
    // Also (p-2)/(p-1) = 0, so gamma_m = 1
    // Since gamma_m/gamma_c < 1 and p <= 2, should go to fast cooling branch: return eps_e
    RadParams rad;
    rad.eps_e = 0.1;
    rad.eps_B = 0.01;
    rad.p = 2.0;
    rad.xi_e = 1.0;
    Real t_comv = 1e6;
    Real Gamma_th = 10.0;
    Real u = 1e-5;
    Real eff = compute_radiative_efficiency(t_comv, Gamma_th, u, rad);
    // p <= 2 always returns eps_e (falls through to fast cooling branch)
    BOOST_CHECK_CLOSE(eff, rad.eps_e, 1e-10);
}

// ---------------------------------------------------------------------------
// 44. simpson_logspace: constant function
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(simpson_logspace_constant) {
    // Integrate f(u) = K (constant) from u_min = ln(r) - 18 to u_max = ln(r)
    // Result = K * 18
    Real K = 42.0;
    Real r = 1e5;
    Real result = simpson_logspace([=](Real /*u*/) { return K; }, r);
    Real expected = K * 18.0;
    BOOST_CHECK_CLOSE(result, expected, 1e-8);
}

// ---------------------------------------------------------------------------
// 45. enclosed_thermal_energy: consistency check (thermal <= total swept)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(enclosed_thermal_energy_consistency) {
    Real r = 1e15;
    Real rho_const = 1e-24;
    Real Gamma = 10.0;
    Real ad_idx = 4.0 / 3.0;
    Real eps_e = 0.1;

    Real mass = enclosed_mass([=](Real /*r_*/) { return rho_const; }, r);
    Real U_th = enclosed_thermal_energy([=](Real /*r_*/) { return rho_const; }, r, Gamma, ad_idx, eps_e);

    // Thermal energy should be positive
    BOOST_CHECK_GT(U_th, 0.0);
    BOOST_CHECK(std::isfinite(U_th));

    // Thermal energy should be less than total available energy: (1 - eps_e) * (Gamma - 1) * c2 * mass
    Real max_energy = (1 - eps_e) * (Gamma - 1) * con::c2 * mass;
    BOOST_CHECK_LE(U_th, max_energy * 1.01); // allow 1% numerical tolerance
}

BOOST_AUTO_TEST_SUITE_END()
