//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once
#include "reverse-shock.hpp"
#include "shock.h"

template <typename Ejecta, typename Medium, typename Value>
FRShockEqn<Ejecta, Medium, Value>::FRShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi,
                                              std::array<Real, lanes> const& thetas, RadParams const& rad_fwd,
                                              RadParams const& /*rad_rvs*/)
    : medium(medium),
      ejecta(ejecta),
      rad_fwd(rad_fwd),
      phi(phi),
      theta0_lanes(thetas),
      Gamma4(at_theta0([&](Real th) { return ejecta.Gamma0(phi, th); })),
      deps0_dt_(at_theta0([&](Real th) { return ejecta.eps_k(phi, th) / ejecta.T0; })) {
    using std::sqrt;
    u4_ = sqrt(Gamma4 * Gamma4 - 1) * con::c;
    beta4_ = physics::relativistic::gamma_to_beta(Gamma4);
    cs4_ = compute_sound_speed(Gamma4);
    dm0_dt_ = deps0_dt_ / (Gamma4 * con::c2);
    if constexpr (HasSigma<Ejecta>) {
        dm0_dt_ = dm0_dt_ / at_theta0([&](Real th) { return 1 + ejecta.sigma0(phi, th); });
    }
    eps_rad_fwd_ = RadiativeEfficiency(rad_fwd);
}

template <typename Ejecta, typename Medium, typename Value>
Value FRShockEqn<Ejecta, Medium, Value>::injection_efficiency(State const& diff) const noexcept {
    using std::min;
    const Value f = min(diff.m4 / dm0_dt_, Value(1));
    return where(dm0_dt_ > 0 && diff.m4 > 0, f, 0.);
}

template <typename Ejecta, typename Medium, typename Value>
bool FRShockEqn<Ejecta, Medium, Value>::crossing_complete(Real m3, Real m4, Real t, size_t w) const noexcept {
    if (m3 < 0.999 * m4)
        return false;
    // Check if injection has stopped (same window as compute_dm4_dt)
    if (injection_window(t) > 1e-6)
        return false;
    if constexpr (ScalarState::mass_inject) {
        if (ejecta.dm_dt(phi, theta0_lanes[w], t) > 0)
            return false;
    }
    return true;
}

template <typename Ejecta, typename Medium, typename Value>
Value FRShockEqn<Ejecta, Medium, Value>::compute_dGamma_dt(State const& state, State const& diff, Real t,
                                                           Value const& Gamma34) const noexcept {
    const Value ad_idx2 = physics::thermo::adiabatic_idx(state.Gamma);
    const Value ad_idx3 = physics::thermo::adiabatic_idx(Gamma34);

    const Value Gamma_eff2 = compute_effective_Gamma(ad_idx2, state.Gamma);
    const Value Gamma_eff3 = compute_effective_Gamma(ad_idx3, state.Gamma);

    const Value dGamma_eff2_dGamma = compute_effective_Gamma_dGamma(ad_idx2, state.Gamma);
    const Value dGamma_eff3_dGamma = compute_effective_Gamma_dGamma(ad_idx3, state.Gamma);

    Value deps_dt = 0.;

    if constexpr (State::energy_inject) {
        deps_dt = at_theta0([&](Real th) { return ejecta.deps_dt(phi, th, t); });
    }

    const Value a = (state.Gamma - 1) * con::c2 * diff.m2 + (state.Gamma - Gamma4) * con::c2 * diff.m3 +
                    Gamma_eff2 * diff.U2_th + Gamma_eff3 * diff.U3_th - deps_dt;
    const Value b =
        (state.m2 + state.m3) * con::c2 + dGamma_eff2_dGamma * state.U2_th + dGamma_eff3_dGamma * state.U3_th;

    // b == 0 needs no separate test: it lands dGamma on +-inf or nan, caught here.
    using std::isfinite;
    const Value dGamma = -a / b;
    return where(isfinite(dGamma), dGamma, 0.);
}

template <typename Ejecta, typename Medium, typename Value>
Value FRShockEqn<Ejecta, Medium, Value>::compute_dU2_dt(State const& state, State const& diff,
                                                        Value const& rho) const noexcept {
    const Value e_th = (state.Gamma - 1) * 4 * state.Gamma * rho * con::c2;
    const Value eps_rad = eps_rad_fwd_(state.t_comv, state.Gamma, e_th);

    const Value ad_idx = physics::thermo::adiabatic_idx(state.Gamma);

    const Value shock_heating = compute_shock_heating_rate(state.Gamma, diff.m2);

    const Value adiabatic_cooling =
        compute_adiabatic_cooling_rate2(ad_idx, state.r, state.x4, state.U2_th, diff.r, diff.x4);

    return (1 - eps_rad) * shock_heating + adiabatic_cooling;
}

template <typename Ejecta, typename Medium, typename Value>
Value FRShockEqn<Ejecta, Medium, Value>::compute_dU3_dt(State const& state, State const& diff,
                                                        Value const& Gamma34) const noexcept {
    const Value ad_idx = physics::thermo::adiabatic_idx(Gamma34);
    const Value adiabatic_cooling =
        compute_adiabatic_cooling_rate2(ad_idx, state.r, state.x3, state.U3_th, diff.r, diff.x3);

    const Value shock_heating = compute_shock_heating_rate(Gamma34, diff.m3);
    return shock_heating + adiabatic_cooling;
}

template <typename Ejecta, typename Medium, typename Value>
Value FRShockEqn<Ejecta, Medium, Value>::compute_dx3_dt(State const& state, State const& diff, Value const& f,
                                                        Value const& Gamma34, Value const& sigma,
                                                        Value const& comp_ratio) const noexcept {
    using std::abs, std::max, std::min, std::sqrt;
    const Value sound_expansion = compute_shell_sound_expansion_rate(Gamma34, diff.t_comv);

    // Blend crossing/sound-expansion using same effective mass fraction as dm3
    const Value remaining = max(state.m4 - state.m3, Value(0));
    const Value crossing_w = f + (1.0 - f) * remaining / state.m4;

    // Penetration factor: the reverse shock eats the shell only while shocked-region
    // compression advances on the unshocked ejecta (-> 0 for magnetized Gamma34 -> 1,
    // >= 1 for sigma = 0).
    const Value penetration = state.Gamma * comp_ratio / Gamma4 - 1;

    // Sound expansion at every lane (post-crossing, no shell, or no penetration).
    const auto is_expansion = state.m4 <= 0 || crossing_w < 1e-6 || penetration <= 0;
    if (all_of(is_expansion)) {
        return sound_expansion;
    }

    const Value beta3 = physics::relativistic::gamma_to_beta(state.Gamma);
    // (beta4 - beta3) / (1 - beta3) via Lorentz-factor identities: the direct differences
    // catastrophically cancel at crossing onset where Gamma3 ~ Gamma4.
    const Value dx3dt = (Gamma4 - state.Gamma) * (Gamma4 + state.Gamma) * (1 + beta3) * con::c /
                        (Gamma4 * Gamma4 * (beta3 + beta4_) * penetration);
    Value crossing = abs(dx3dt * state.Gamma);

    // Weak-penetration regularization (magnetized shells only): the crossing rate is
    // capped at the fast magnetosonic speed of the upstream.
    const auto penetrates = penetration >= 1; // >= 1 whenever sigma = 0: skip the cap entirely
    if (!all_of(penetrates)) {
        const Value cs = compute_sound_speed(Gamma34);
        const Value va2 = sigma / (1 + sigma); // relativistic (v_A/c)^2
        const Value cs2 = cs * cs / (con::c * con::c);
        const Value v_ms = sqrt(va2 + cs2 * (1 - va2)) * con::c * diff.t_comv;
        crossing = where(penetrates, crossing, min(crossing, v_ms));
    }

    const Value blended = crossing_w * crossing + (1.0 - crossing_w) * sound_expansion;
    return where(is_expansion, sound_expansion, blended);
}

template <typename Ejecta, typename Medium, typename Value>
Value FRShockEqn<Ejecta, Medium, Value>::compute_dm3_dt(State const& state, State const& diff, Value const& f,
                                                        Value const& comp_ratio) const noexcept {
    using std::max, std::min;
    const Value remaining = max(state.m4 - state.m3, Value(0));

    // Zero at every lane: nothing left to shock and injection over.
    const auto is_done = state.m4 <= 0 || (remaining <= 0 && f < 1e-6);
    if (all_of(is_done)) {
        return 0.;
    }

    // Blend column density: m4 during injection, (m4-m3) after injection
    const Value eff_mass = f * state.m4 + (1.0 - f) * remaining;
    const Value dm3dt = eff_mass * comp_ratio / state.x4 * diff.x3;

    // During active injection, cap the rate so dm3 <= dm4
    const Value ratio = state.m3 / state.m4;
    const Value cap_w = smoothstep(0, 1.0, ratio);
    const Value with_cap = (1.0 - cap_w) * dm3dt + cap_w * min(dm3dt, diff.m4);
    const Value rate = where(f > 1e-6, with_cap, dm3dt);

    return where(is_done, 0., rate);
}

template <typename Ejecta, typename Medium, typename Value>
Value FRShockEqn<Ejecta, Medium, Value>::compute_dx4_dt(State const& /*state*/, State const& diff,
                                                        Value const& f) const noexcept {
    const Value sound_expansion = cs4_ * diff.t_comv;
    return where(f > 1e-6, f * u4_ + (1 - f) * sound_expansion, sound_expansion);
}

template <typename Ejecta, typename Medium, typename Value>
Value FRShockEqn<Ejecta, Medium, Value>::compute_deps4_dt(State const& /*state*/, State const& /*diff*/,
                                                          Real t) const noexcept {
    Value deps4_dt = 0.;

    const Real inject_w = injection_window(t);
    if (inject_w > 1e-6) {
        deps4_dt = inject_w * deps0_dt_;
    }

    if constexpr (State::energy_inject) {
        deps4_dt += at_theta0([&](Real th) { return ejecta.deps_dt(phi, th, t); });
    }

    return deps4_dt;
}

template <typename Ejecta, typename Medium, typename Value>
Value FRShockEqn<Ejecta, Medium, Value>::compute_dm4_dt(State const& /*state*/, State const& /*diff*/,
                                                        Real t) const noexcept {
    Value dm4_dt = 0.;

    const Real inject_w = injection_window(t);
    if (inject_w > 1e-6) {
        dm4_dt = inject_w * dm0_dt_;
    }

    if constexpr (State::mass_inject) {
        dm4_dt += at_theta0([&](Real th) { return ejecta.dm_dt(phi, th, t); });
    }

    return dm4_dt;
}

template <typename Ejecta, typename Medium, typename Value>
void FRShockEqn<Ejecta, Medium, Value>::operator()(State const& state_raw, State& diff, Real t) const noexcept {
    using std::clamp, std::max, std::sqrt;
    // Evaluate all rates on the projection of the state onto the physical domain:
    // Gamma3 in [1, Gamma4], m3 in [0, m4], x3 >= 0, U3_th >= 0. Adaptive RK stages can
    // overshoot these bounds at the seed scale of the region-3 variables, and rates from
    // unphysical values are not self-correcting. On the physical domain the projection
    // is the identity, so healthy solves are unaffected.
    State state = state_raw;
    state.Gamma = clamp(state.Gamma, Value(1), Gamma4);
    state.m3 = clamp(state.m3, Value(0), max(state.m4, Value(0)));
    state.x3 = max(state.x3, Value(0));
    state.U3_th = max(state.U3_th, Value(0));

    const Value Gamma = state.Gamma;
    const Value u3 = sqrt((Gamma - 1) * (Gamma + 1));

    diff.r = compute_dr_dt(Gamma, u3);
    diff.t_comv = Gamma + u3;

    const Value rho = lane_rho(medium, phi, state.theta, state.r);
    diff.m2 = state.r * state.r * rho * diff.r;

    diff.eps4 = compute_deps4_dt(state, diff, t);
    diff.m4 = compute_dm4_dt(state, diff, t);

    const Value Gamma34 = compute_rel_Gamma(Gamma4, state.Gamma);
    const Value sigma = compute_shell_sigma(state);
    // The jump conditions (sigma > 0 solves a cubic) stay scalar, evaluated per lane.
    const Value comp_ratio = lane::apply<Value>(compute_4vel_jump, Gamma34, sigma);

    const Value f = injection_efficiency(diff);
    diff.x4 = compute_dx4_dt(state, diff, f);
    diff.x3 = compute_dx3_dt(state, diff, f, Gamma34, sigma, comp_ratio);

    diff.m3 = compute_dm3_dt(state, diff, f, comp_ratio);

    diff.U2_th = compute_dU2_dt(state, diff, rho);
    diff.U3_th = compute_dU3_dt(state, diff, Gamma34);

    diff.Gamma = compute_dGamma_dt(state, diff, t, Gamma34);

    diff.theta = 0.;
}

inline Real compute_init_comv_shell_width(Real Gamma4, Real t0, Real T);

template <typename Ejecta, typename Medium, typename Value>
void FRShockEqn<Ejecta, Medium, Value>::save_cross_state(ScalarState const& state, size_t w) {
    V3_comv_x[w] = state.r * state.r * state.x3;

    const Real sigma4 = compute_shell_sigma(state, w);
    const Real comp_ratio34 = compute_compression(lane::get(Gamma4, w), state.Gamma, sigma4);
    const Real rho4 = state.m4 / (state.r * state.r * state.x4);
    rho3_x[w] = rho4 * comp_ratio34;

    const Real B4 = compute_upstr_B(rho4, sigma4);
    B3_ordered_x[w] = B4 * comp_ratio34;
}

template <typename Ejecta, typename Medium, typename Value>
void FRShockEqn<Ejecta, Medium, Value>::set_init_state_lane(ScalarState& state, Real t0, size_t w) const noexcept {
    const Real G4 = lane::get(Gamma4, w);
    const Real beta4 = physics::relativistic::gamma_to_beta(G4);

    // beta / (1 - beta) = beta * Gamma^2 * (1 + beta): avoids 1 - beta, which
    // loses ~2 Gamma^2 ulp of relative precision for ultrarelativistic shells.
    state.r = beta4 * con::c * t0 * G4 * G4 * (1 + beta4);
    state.t_comv = state.r / std::sqrt((G4 - 1) * (G4 + 1)) / con::c;
    state.theta = theta0_lanes[w];

    const Real dt = std::min(t0, ejecta.T0);
    state.eps4 = lane::get(deps0_dt_, w) * dt;
    state.m4 = lane::get(dm0_dt_, w) * dt;
    state.x4 = compute_init_comv_shell_width(G4, t0, ejecta.T0);

    auto rho_func = [&](Real r_) { return medium.rho(phi, theta0_lanes[w], r_); };

    state.m2 = enclosed_mass(rho_func, state.r);

    const Real m_jet_total = lane::get(dm0_dt_, w) * ejecta.T0;
    if (m_jet_total > 0 && state.m2 > 0) {
        state.Gamma = G4 / (1 + state.m2 / m_jet_total);
    } else {
        state.Gamma = G4;
    }

    Real ad_idx = physics::thermo::adiabatic_idx(state.Gamma);
    state.U2_th =
        enclosed_thermal_energy(rho_func, state.r, state.Gamma, ad_idx, rad_fwd.radiative ? rad_fwd.eps_e : 0.0);

    const Real Gamma34 = compute_rel_Gamma(G4, state.Gamma);
    if (Gamma34 > 1 && state.m4 > 0 && state.x4 > 0) {
        constexpr Real seed_frac = 1e-8;
        const Real sigma = compute_shell_sigma(state, w);
        const Real comp_ratio = compute_4vel_jump(Gamma34, sigma);
        state.x3 = state.x4 * seed_frac;
        state.m3 = state.m4 * comp_ratio * state.x3 / state.x4;
        state.U3_th = (Gamma34 - 1) * state.m3 * con::c2;
    } else {
        state.m3 = 0;
        state.U3_th = 0;
        state.x3 = 0;
    }
}

template <typename Ejecta, typename Medium, typename Value>
void FRShockEqn<Ejecta, Medium, Value>::set_init_state(State& state, Real t0) const noexcept {
    for (size_t w = 0; w < lanes; ++w) {
        ScalarState sw;
        set_init_state_lane(sw, t0, w);
        insert_lane(state, w, sw);
    }
}

template <typename Ejecta, typename Medium, typename Value>
Value FRShockEqn<Ejecta, Medium, Value>::compute_shell_sigma(State const& state) const noexcept {
    const Value sigma = state.eps4 / (Gamma4 * state.m4 * con::c2) - 1;
    return where(sigma > con::sigma_cut, sigma, 0.);
}

template <typename Ejecta, typename Medium, typename Value>
Real FRShockEqn<Ejecta, Medium, Value>::compute_shell_sigma(ScalarState const& state, size_t w) const noexcept {
    const Real sigma = state.eps4 / (lane::get(Gamma4, w) * state.m4 * con::c2) - 1;
    return (sigma > con::sigma_cut) ? sigma : 0;
}

//---------------------------------------------------------------------------------------------------------------------
// Helper functions
//---------------------------------------------------------------------------------------------------------------------

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the comoving shell width at initial radius.
 * @details Accounts for both pure injection phase and sound-speed shell expansion phase.
 * @param Gamma4 Lorentz factor of the unshocked ejecta
 * @param t0 Initial time
 * @param T Engine duration
 * @return The comoving shell width
 * <!-- ************************************************************************************** -->
 */
inline Real compute_init_comv_shell_width(Real Gamma4, Real t0, Real T) {
    const Real beta4 = physics::relativistic::gamma_to_beta(Gamma4);
    if (t0 < T) { // pure injection
        return Gamma4 * t0 * beta4 * con::c;
    } else { // injection + sound-speed shell expansion
        const Real cs = compute_sound_speed(Gamma4);
        return Gamma4 * T * beta4 * con::c + cs * (t0 - T) * Gamma4;
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Saves the state of lane w's reverse shock at a grid point.
 * <!-- ************************************************************************************** -->
 */
template <typename Eqn, typename State>
void save_rvs_shock_state(size_t i, size_t j, size_t k, Eqn const& eqn, State const& state, Shock& shock, size_t w) {
    if (k <= shock.injection_idx(i, j)) {
        const Real Gamma4 = lane::get(eqn.Gamma4, w);
        const Real sigma4 = eqn.compute_shell_sigma(state, w);

        const Real comp_ratio34 = compute_compression(Gamma4, state.Gamma, sigma4);
        const Real rho4 = state.m4 / (state.r * state.r * state.x4);
        const Real Gamma3_th = compute_Gamma_therm(state.U3_th, state.m3, true);

        const Real B4 = compute_upstr_B(rho4, sigma4);
        const Real B3 = compute_downstr_B(shock.rad.eps_B, rho4, B4, Gamma3_th, comp_ratio34);

        write_shock_state(shock, i, j, k, state.t_comv, state.r, state.theta, state.Gamma, Gamma3_th, B3, state.m3);
    } else {
        Real V3_comv = state.r * state.r * state.x3;
        const Real comp_ratio = eqn.V3_comv_x[w] / V3_comv;
        const Real Gamma3_th = compute_Gamma_therm(state.U3_th, state.m3);

        const Real B3 = compute_downstr_B(shock.rad.eps_B, eqn.rho3_x[w], eqn.B3_ordered_x[w], Gamma3_th, comp_ratio);

        write_shock_state(shock, i, j, k, state.t_comv, state.r, state.theta, state.Gamma, Gamma3_th, B3, state.m3);
    }
}

inline void reverse_shock_early_extrap(size_t i, size_t j, Shock& shock) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    size_t idx_cut = 0;
    for (; idx_cut < t_size; ++idx_cut) {
        if (shock.Gamma_th(i, j, idx_cut) > con::gamma_therm_cut) {
            break;
        }
    }

    Real gamma_slope = 0;
    Real B_slope = 0;
    Real N_p_slope = 0;
    const Real log2_r = fast_log2(shock.r(i, j, idx_cut));
    const Real log2_Gamma_th = fast_log2(shock.Gamma_th(i, j, idx_cut) - 1);
    const Real log2_B = fast_log2(shock.B(i, j, idx_cut));
    const Real log2_N_p = fast_log2(shock.N_p(i, j, idx_cut));

    constexpr size_t offset = 2;

    if (idx_cut == 0 || idx_cut >= t_size - offset || idx_cut >= shock.injection_idx(i, j)) {
        return;
    } else {
        gamma_slope = (fast_log2(shock.Gamma_th(i, j, idx_cut + offset) - 1) - log2_Gamma_th) /
                      (fast_log2(shock.r(i, j, idx_cut + offset)) - log2_r);

        B_slope = (fast_log2(shock.B(i, j, idx_cut + offset)) - log2_B) /
                  (fast_log2(shock.r(i, j, idx_cut + offset)) - log2_r);

        N_p_slope = (fast_log2(shock.N_p(i, j, idx_cut + offset)) - log2_N_p) /
                    (fast_log2(shock.r(i, j, idx_cut + offset)) - log2_r);
    }

    for (size_t k = 0; k < idx_cut; k++) {
        const Real dlog2_r = fast_log2(shock.r(i, j, k)) - log2_r;
        shock.Gamma_th(i, j, k) = 1 + fast_exp2(log2_Gamma_th + gamma_slope * dlog2_r);
        shock.B(i, j, k) = fast_exp2(log2_B + B_slope * dlog2_r);
        shock.N_p(i, j, k) = fast_exp2(log2_N_p + N_p_slope * dlog2_r);
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Bisects lane w's crossing-end time on the stepper's dense output.
 * @details On entry, crossing must be incomplete at t_lo and complete at t_hi (the ends of the
 *          last accepted step). On exit, state holds the solution at the returned time.
 * <!-- ************************************************************************************** -->
 */
template <typename Stepper, typename Eqn>
Real locate_crossing_time(Stepper& stepper, Eqn const& eqn, typename Eqn::State& state, Real t_lo, Real t_hi,
                          size_t w) {
    for (int iter = 0; iter < 100 && (t_hi - t_lo) > 1e-12 * t_hi; ++iter) {
        const Real t_mid = 0.5 * (t_lo + t_hi);
        stepper.calc_state(t_mid, state);
        if (eqn.crossing_complete(lane::get(state.m3, w), lane::get(state.m4, w), t_mid, w)) {
            t_hi = t_mid;
        } else {
            t_lo = t_mid;
        }
    }
    stepper.calc_state(t_hi, state);
    return t_hi;
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Integrates a batch of shock-pair rows in lockstep, storing both shocks on each row's
 *        t grid and handling each lane's reverse-shock crossing at its own time.
 * <!-- ************************************************************************************** -->
 */
template <typename Eqn>
void grid_solve_shock_pair(size_t i, std::array<size_t, Eqn::lanes> const& js, size_t n_live, Coord const& coord,
                           Shock& shock_fwd, Shock& shock_rvs, Eqn& eqn, Real rtol) {
    using namespace boost::numeric::odeint;
    using State = typename Eqn::State;
    constexpr size_t lanes = Eqn::lanes;

    const size_t t_size = coord.t.shape()[2];
    const auto [t_row, t0, t_end] = batch_time_span(coord, i, js, n_live, eqn, 0.01 * unit::sec);

    State state;
    eqn.set_init_state(state, t0);

    constexpr Real RS_Gamma_limit = 1.03;
    const auto live = mark_live_lanes<lanes>(state, n_live, RS_Gamma_limit, [&](size_t w) {
        set_stopping_shock(i, js[w], shock_fwd, extract_lane(state, w));
        set_stopping_shock(i, js[w], shock_rvs, extract_lane(state, w));
    });
    if (!lane::any_of(live)) {
        return;
    }

    // Magnetized shells solve the tolerance-sensitive sigma > 0 jump conditions:
    // tighten one decade so they deliver the same accuracy class as sigma = 0.
    const auto sigma0 = eqn.compute_shell_sigma(state);
    for (size_t w = 0; w < n_live; ++w) {
        if (live[w] && lane::get(sigma0, w) > 0) {
            rtol *= defaults::solver::magnetized_rtol_factor;
            break;
        }
    }

    auto stepper = make_dense_output(rtol, rtol, runge_kutta_dopri5<State, Real, State, Real, LaneAlgebra>());
    stepper.initialize(state, t0, 1e-9 * t0);

    std::array<size_t, lanes> k{};
    std::array<bool, lanes> crossing;
    crossing.fill(true);
    std::array<bool, lanes> inj_pending{};
    std::array<Real, lanes> t_cross{};
    Real t_step_start = t0;
    for (size_t steps = 0; stepper.current_time() <= t_end;) {
        stepper.do_step(eqn);
        if (++steps > defaults::solver::max_ode_steps) {
            std::fprintf(stderr, "Warning: reverse shock ODE exceeded %zu steps at (i=%zu, j=%zu), giving up\n",
                         defaults::solver::max_ode_steps, i, js[0]);
            break;
        }
        if (stepper.current_time() + stepper.current_time_step() == stepper.current_time()) {
            // dt has collapsed below one ulp of t: time can no longer advance. Fail fast.
            std::fprintf(stderr,
                         "Warning: reverse shock ODE stalled (dt below time ulp) at (i=%zu, j=%zu), giving up\n", i,
                         js[0]);
            break;
        }
        for (size_t w = 0; w < n_live; ++w) {
            if (live[w] && crossing[w] &&
                eqn.crossing_complete(lane::get(stepper.current_state().m3, w), lane::get(stepper.current_state().m4, w),
                                      stepper.current_time(), w)) {
                // Lane w's crossing ended within this step: freeze its snapshot at the
                // bisected event time rather than at the next stored grid time.
                t_cross[w] = locate_crossing_time(stepper, eqn, state, t_step_start, stepper.current_time(), w);
                eqn.save_cross_state(extract_lane(state, w), w);
                crossing[w] = false;
                inj_pending[w] = true;
            }
        }
        t_step_start = stepper.current_time();
        for (size_t w = 0; w < n_live; ++w) {
            if (!live[w]) {
                continue;
            }
            while (k[w] < t_size && stepper.current_time() > t_row[w][k[w]]) {
                stepper.calc_state(t_row[w][k[w]], state);
                if (inj_pending[w] && t_row[w][k[w]] >= t_cross[w]) {
                    shock_rvs.injection_idx(i, js[w]) = k[w] > 0 ? k[w] : 1;
                    inj_pending[w] = false;
                }
                const auto sw = extract_lane(state, w);
                save_fwd_shock_state(i, js[w], k[w], eqn, sw, shock_fwd);
                save_rvs_shock_state(i, js[w], k[w], eqn, sw, shock_rvs, w);
                ++k[w];
            }
        }
    }
    for (size_t w = 0; w < n_live; ++w) {
        if (live[w]) {
            reverse_shock_early_extrap(i, js[w], shock_rvs);
        }
    }
}

template <typename Ejecta, typename Medium>
ShockPair generate_shock_pair(Coord const& coord, Medium const& medium, Ejecta const& jet, RadParams const& rad_fwd,
                              RadParams const& rad_rvs, Real rtol) {
    auto [phi_size, theta_size, t_size] = coord.shape();
    const size_t phi_size_needed = coord.t.shape()[0];
    Shock f_shock(phi_size_needed, theta_size, t_size, rad_fwd);
    Shock r_shock(phi_size_needed, theta_size, t_size, rad_rvs);

    for (size_t i = 0; i < phi_size_needed; ++i) {
        for_each_theta_batch(coord, [&](auto const& js, auto const& thetas, size_t n_live) {
            constexpr size_t width = std::tuple_size_v<std::remove_cvref_t<decltype(js)>>;
            using Value = std::conditional_t<width == 1, Real, lane::Batch>;
            FRShockEqn<Ejecta, Medium, Value> eqn(medium, jet, coord.phi(i), thetas, rad_fwd, rad_rvs);
            grid_solve_shock_pair(i, js, n_live, coord, f_shock, r_shock, eqn, rtol);
        });

        if (coord.symmetry >= Symmetry::phi_symmetric) {
            f_shock.broadcast_groups(coord);
            r_shock.broadcast_groups(coord);
            return std::make_pair(std::move(f_shock), std::move(r_shock));
        }
    }
    return std::make_pair(std::move(f_shock), std::move(r_shock));
}
