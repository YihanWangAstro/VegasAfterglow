//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once
#include "reverse-shock.hpp"
#include "shock.h"

template <typename Eqn, typename State>
bool is_crossing(Eqn const& eqn, State const& state, Real t) {
    Real dmdt = 0;

    if constexpr (State::mass_inject) {
        dmdt = eqn.ejecta.dm_dt(eqn.phi, state.theta, t);
    }

    return (state.m3 < state.m4) || dmdt > 0 || t < eqn.ejecta.T0;
}

template <typename Ejecta, typename Medium>
FRShockEqn<Ejecta, Medium>::FRShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta,
                                       RadParams const& rad_fwd, RadParams const& rad_rvs)
    : medium(medium),
      ejecta(ejecta),
      rad_fwd(rad_fwd),
      rad_rvs(rad_rvs),
      phi(phi),
      theta0(theta),
      Gamma4(ejecta.Gamma0(phi, theta)),
      deps0_dt(ejecta.eps_k(phi, theta) / ejecta.T0),
      dm0_dt(deps0_dt / (Gamma4 * con::c2)),
      u4(std::sqrt(Gamma4 * Gamma4 - 1) * con::c) {
    if constexpr (HasSigma<Ejecta>) {
        dm0_dt /= 1 + ejecta.sigma0(phi, theta);
    }
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dGamma_dt(State const& state, State const& diff, Real t) const noexcept {
    const Real Gamma34 = compute_rel_Gamma(Gamma4, state.Gamma);
    const Real ad_idx2 = adiabatic_idx(state.Gamma);
    const Real ad_idx3 = adiabatic_idx(Gamma34);

    Real Gamma_eff2 = compute_effective_Gamma(ad_idx2, state.Gamma);
    Real Gamma_eff3 = compute_effective_Gamma(ad_idx3, state.Gamma);

    Real dGamma_eff2_dGamma = compute_effective_Gamma_dGamma(ad_idx2, state.Gamma);
    Real dGamma_eff3_dGamma = compute_effective_Gamma_dGamma(ad_idx3, state.Gamma);

    Real deps_dt = 0;

    if constexpr (State::energy_inject) {
        deps_dt = ejecta.deps_dt(phi, state.theta, t);
    }

    const Real a = (state.Gamma - 1) * con::c2 * diff.m2 + (state.Gamma - Gamma4) * con::c2 * diff.m3 +
             Gamma_eff2 * diff.U2_th + Gamma_eff3 * diff.U3_th - deps_dt;
    const Real b = (state.m2 + state.m3) * con::c2 + dGamma_eff2_dGamma * state.U2_th + dGamma_eff3_dGamma * state.U3_th;

    if (b == 0 || std::isnan(-a / b) || std::isinf(-a / b)) {
        return 0;
    } else {
        return -a / b;
    }
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dU2_dt(State const& state, State const& diff, Real t) const noexcept {
    const Real e_th = (state.Gamma - 1) * 4 * state.Gamma * medium.rho(phi, state.theta, state.r) * con::c2;
    // Real V_comv = state.r * state.r * state.r / (12 * state.Gamma * state.Gamma);
    // Real e_th = state.U2_th / V_comv;
    const Real eps_rad = compute_radiative_efficiency(state.t_comv, state.Gamma, e_th, rad_fwd);

    const Real ad_idx = adiabatic_idx(state.Gamma);

    const Real shock_heating = compute_shock_heating_rate(state.Gamma, diff.m2);

    const Real adiabatic_cooling = compute_adiabatic_cooling_rate2(ad_idx, state.r, state.x4, state.U2_th, diff.r, diff.x4);

    return (1 - eps_rad) * shock_heating + adiabatic_cooling;
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dU3_dt(State const& state, State const& diff, Real t) const noexcept {
    const Real Gamma34 = compute_rel_Gamma(this->Gamma4, state.Gamma);
    const Real ad_idx = adiabatic_idx(Gamma34);
    const Real adiabatic_cooling = compute_adiabatic_cooling_rate2(ad_idx, state.r, state.x3, state.U3_th, diff.r, diff.x3);

    if (state.m3 < state.m4 || diff.m4 > 0) { // reverse shock still crossing
        const Real shock_heating = compute_shock_heating_rate(Gamma34, diff.m3);
        constexpr Real eps_rad = 0; ////compute_radiative_efficiency(state.t_comv, state.Gamma, e_th, rad_rvs);
        return (1 - eps_rad) * shock_heating + adiabatic_cooling;
    } else {
        return adiabatic_cooling;
    }
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dx3_dt(State const& state, State const& diff, Real t) const noexcept {
    if ((state.m3 < state.m4 || diff.m4 > 0) && (state.Gamma != this->Gamma4)) {
        const Real sigma = compute_shell_sigma(state);
        const Real Gamma34 = compute_rel_Gamma(this->Gamma4, state.Gamma);
        const Real beta3 = gamma_to_beta(state.Gamma);
        const Real beta4 = gamma_to_beta(this->Gamma4);
        Real comp_ratio = compute_4vel_jump(Gamma34, sigma);
        Real dx3dt = (beta4 - beta3) * con::c / ((1 - beta3) * (state.Gamma * comp_ratio / this->Gamma4 - 1));

        return std::fabs(dx3dt * state.Gamma);
    } else {
        const Real Gamma34 = compute_rel_Gamma(this->Gamma4, state.Gamma);
        return compute_shell_spreading_rate(Gamma34, diff.t_comv);
    }
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dx4_dt(State const& state, State const& diff, Real t) const noexcept {
    if (diff.m4 > 0) {
        return u4;
    } else {
        return compute_shell_spreading_rate(this->Gamma4, diff.t_comv);
    }
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dm2_dt(State const& state, State const& diff, Real t) const noexcept {
    return state.r * state.r * medium.rho(phi, state.theta, state.r) * diff.r;
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dm3_dt(State const& state, State const& diff, Real t) const noexcept {
    if ((state.m3 < state.m4 || diff.m4 > 0) && (state.Gamma != this->Gamma4)) {
        const Real sigma = compute_shell_sigma(state);
        const Real Gamma34 = compute_rel_Gamma(this->Gamma4, state.Gamma);
        Real comp_ratio = compute_4vel_jump(Gamma34, sigma);
        Real column_den3 = state.m4 * comp_ratio / state.x4;
        Real dm3dt = column_den3 * diff.x3;

        if (state.m3 >= state.m4) {
            return std::min(dm3dt, diff.m4);
        } else {
            return dm3dt;
        }
    } else {
        return 0.;
    }
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_deps4_dt(State const& state, State const& diff, Real t) const noexcept {
    Real deps4_dt = 0;

    if (t < ejecta.T0) {
        deps4_dt = deps0_dt;
    }

    if constexpr (State::energy_inject) {
        deps4_dt += ejecta.deps_dt(phi, theta0, t);
    }

    return deps4_dt;
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dm4_dt(State const& state, State const& diff, Real t) const noexcept {
    Real dm4_dt = 0;

    if (t < ejecta.T0) {
        dm4_dt = dm0_dt;
    }

    if constexpr (State::mass_inject) {
        dm4_dt += ejecta.dm_dt(phi, theta0, t);
    }

    return dm4_dt;
}

template <typename Ejecta, typename Medium>
void FRShockEqn<Ejecta, Medium>::operator()(State const& state, State& diff, Real t) {
    const Real beta3 = gamma_to_beta(state.Gamma);

    diff.r = compute_dr_dt(beta3);
    diff.t_comv = compute_dt_dt_comv(state.Gamma, beta3);

    diff.m2 = compute_dm2_dt(state, diff, t);

    diff.eps4 = compute_deps4_dt(state, diff, t);
    diff.m4 = compute_dm4_dt(state, diff, t);

    diff.x4 = compute_dx4_dt(state, diff, t);
    diff.x3 = compute_dx3_dt(state, diff, t);

    diff.m3 = compute_dm3_dt(state, diff, t);

    diff.U2_th = compute_dU2_dt(state, diff, t);
    diff.U3_th = compute_dU3_dt(state, diff, t);

    diff.Gamma = compute_dGamma_dt(state, diff, t);

    diff.theta = 0;
}

inline Real compute_init_comv_shell_width(Real Gamma4, Real t0, Real T);

template <typename Ejecta, typename Medium>
void FRShockEqn<Ejecta, Medium>::save_cross_state(State const& state) {
    r_x = state.r;
    u_x = std::sqrt(state.Gamma * state.Gamma - 1);

    V3_comv_x = r_x * r_x * state.x3;

    const Real sigma4 = compute_shell_sigma(state);
    const Real comp_ratio34 = compute_compression(Gamma4, state.Gamma, sigma4);
    const Real rho4 = state.m4 / (state.r * state.r * state.x4);
    rho3_x = rho4 * comp_ratio34;

    const Real B4 = compute_upstr_B(rho4, sigma4);
    B3_ordered_x = B4 * comp_ratio34;
}

inline Real calculate_init_m3(Real Gamma4, Real Gamma3, Real m2, Real sigma) {
    const Real Gamma34 = compute_rel_Gamma(Gamma4, Gamma3);
    const Real ad_idx2 = adiabatic_idx(Gamma3);
    const Real ad_idx3 = adiabatic_idx(Gamma34);

    const Real Gamma_eff2 = compute_effective_Gamma(ad_idx2, Gamma3);
    const Real Gamma_eff3 = compute_effective_Gamma(ad_idx3, Gamma3);

    return -m2 * (Gamma3 - 1 + Gamma_eff2 * (Gamma3 - 1)) / (Gamma3 - Gamma4 + Gamma_eff3 * (Gamma34 - 1)) /
           (1 + sigma);
}

template <typename Ejecta, typename Medium>
void FRShockEqn<Ejecta, Medium>::set_init_state(State& state, Real t0) const noexcept {
    const Real beta4 = gamma_to_beta(Gamma4);

    state.r = beta4 * con::c * t0 / (1 - beta4);
    state.t_comv = state.r / std::sqrt(Gamma4 * Gamma4 - 1) / con::c;
    state.theta = theta0;

    const Real dt = std::min(t0, ejecta.T0);
    state.eps4 = deps0_dt * dt;
    state.m4 = dm0_dt * dt;
    state.x4 = compute_init_comv_shell_width(Gamma4, t0, ejecta.T0);

    // constexpr Real Gamma34 = 1;
    // compute_Gamma_from_relative(Gamma4, Gamma34);

    state.m2 = medium.rho(phi, theta0, state.r) * state.r * state.r * state.r / 3;

    //Real sigma4 = compute_shell_sigma(state);
    state.Gamma = Gamma4;

    Real ad_idx = adiabatic_idx(state.Gamma);
    state.U2_th = (state.Gamma - 1) * state.m2 * con::c2 / ad_idx;

    state.m3 = 0;
    //state.m4 * 1e-6;                         // calculate_init_m3(Gamma4, state.Gamma, state.m2, sigma4);
    state.U3_th = 0;
    //state.m3* con::c2 * 1e-6; //(Gamma34 - 1) * state.m3* con::c2;

    state.x3 = 0; //state.x4 * 1e-6;
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the power-law index for post-crossing four-velocity evolution.
 * @details The index transitions from g_low=1.5 for low relative Lorentz factors to g_high=3.5
 *          for high relative Lorentz factors (Blandford-McKee limit).
 * @param gamma_rel Relative Lorentz factor
 * @param k Medium power law index (default: 0)
 * @return The power-law index for velocity evolution
 * <!-- ************************************************************************************** -->
 */
inline Real get_post_cross_g(Real gamma_rel, Real k = 0) {
    constexpr Real g_low = 1.5;  // k is the medium power law index
    constexpr Real g_high = 3.5; // Blandford-McKee limit// TODO: need to be modified for non ISM medium
    const Real p = std::sqrt(std::sqrt(std::fabs(gamma_rel - 1)));
    return g_low + (g_high - g_low) * p / (1 + p);
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_shell_sigma(State const& state) const {
    const Real sigma = state.eps4 / (Gamma4 * state.m4 * con::c2) - 1;
    return (sigma > con::sigma_cut) ? sigma : 0;
}

//---------------------------------------------------------------------------------------------------------------------
// Helper functions
//---------------------------------------------------------------------------------------------------------------------

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the comoving shell width at initial radius.
 * @details Accounts for both pure injection phase and shell spreading phase.
 * @param Gamma4 Lorentz factor of the unshocked ejecta
 * @param t0 Initial time
 * @param T Engine duration
 * @return The comoving shell width
 * <!-- ************************************************************************************** -->
 */
inline Real compute_init_comv_shell_width(Real Gamma4, Real t0, Real T) {
    const Real beta4 = gamma_to_beta(Gamma4);
    if (t0 < T) { // pure injection
        return Gamma4 * t0 * beta4 * con::c;
    } else { // injection+shell spreading
        const Real cs = compute_sound_speed(Gamma4);
        return Gamma4 * T * beta4 * con::c + cs * (t0 - T) * Gamma4;
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Saves the state of reverse shocks at a grid point.
 * @details Updates shock properties for both shocks and checks if crossing is complete.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param k Grid index for time
 * @param eqn The reverse shock equation system
 * @param state Current state of the system
 * @param shock Reverse shock object to update
 * <!-- ************************************************************************************** -->
 */
template <typename Eqn, typename State>
void save_rvs_shock_state(size_t i, size_t j, int k, Eqn const& eqn, State const& state, Shock& shock) {
    if (k <= shock.injection_idx(i, j)) {
        const Real Gamma4 = eqn.Gamma4;
        const Real sigma4 = eqn.compute_shell_sigma(state);

        const Real comp_ratio34 = compute_compression(Gamma4, state.Gamma, sigma4);
        const Real rho4 = state.m4 / (state.r * state.r * state.x4);
        const Real Gamma3_th = compute_Gamma_therm(state.U3_th, state.m3, true);

        const Real B4 = compute_upstr_B(rho4, sigma4);
        const Real B3 = compute_downstr_B(shock.rad.eps_B, rho4, B4, Gamma3_th, comp_ratio34);

        write_shock_state(shock, i, j, k, state.t_comv, state.r, state.theta, state.Gamma, Gamma3_th, B3, state.m3);
    } else {
        Real V3_comv = state.r * state.r * state.x3;
        const Real comp_ratio = eqn.V3_comv_x / V3_comv;
        const Real Gamma3_th = compute_Gamma_therm(state.U3_th, state.m3);

        const Real B3 = compute_downstr_B(shock.rad.eps_B, eqn.rho3_x, eqn.B3_ordered_x, Gamma3_th, comp_ratio);

        write_shock_state(shock, i, j, k, state.t_comv, state.r, state.theta, state.Gamma, Gamma3_th, B3, state.m3);
    }
}

inline void reverse_shock_early_extrap(size_t i, size_t j, Shock& shock) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    size_t idx_cut = 0;
    for (; idx_cut < t_size; ++idx_cut) {
        if (shock.Gamma_th(i, j, idx_cut) > con::Gamma_cut) {
            break;
        }
    }

    Real gamma_slope = 0;
    Real B_slope = 0;
    Real N_p_slope = 0;
    const Real r_lg2 = fast_log2(shock.r(i, j, idx_cut));
    const Real Gamma_th_lg2 = fast_log2(shock.Gamma_th(i, j, idx_cut) - 1);
    const Real B_lg2 = fast_log2(shock.B(i, j, idx_cut));
    const Real N_p_lg2 = fast_log2(shock.N_p(i, j, idx_cut));

    constexpr size_t off_set = 2;

    if (idx_cut == 0 || idx_cut >= t_size - off_set || idx_cut >= shock.injection_idx(i, j)) {
        return;
    } else {
        gamma_slope = (fast_log2(shock.Gamma_th(i, j, idx_cut + off_set) - 1) - Gamma_th_lg2) /
                      (fast_log2(shock.r(i, j, idx_cut + off_set)) - r_lg2);

        B_slope = (fast_log2(shock.B(i, j, idx_cut + off_set)) - B_lg2) /
                  (fast_log2(shock.r(i, j, idx_cut + off_set)) - r_lg2);

        N_p_slope = (fast_log2(shock.N_p(i, j, idx_cut + off_set)) - N_p_lg2) /
                    (fast_log2(shock.r(i, j, idx_cut + off_set)) - r_lg2);
    }

    for (size_t k = 0; k < idx_cut; k++) {
        const Real dr_lg2 = fast_log2(shock.r(i, j, k)) - r_lg2;
        shock.Gamma_th(i, j, k) = 1 + fast_exp2(Gamma_th_lg2 + gamma_slope * dr_lg2);
        shock.B(i, j, k) = fast_exp2(B_lg2 + B_slope * dr_lg2);
        shock.N_p(i, j, k) = fast_exp2(N_p_lg2 + N_p_slope * dr_lg2);
    }
}
/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Solves the reverse/forward shock ODE at a grid point.
 * @details Manages the evolution of both shocks before and after crossing.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param t View of time points
 * @param shock_fwd Forward shock object
 * @param shock_rvs Reverse shock object
 * @param eqn Reverse shock equation system
 * @param rtol Relative tolerance for ODE solver
 * <!-- ************************************************************************************** -->
 */
template <typename Eqn, typename View>
void grid_solve_shock_pair(size_t i, size_t j, View const& t, Shock& shock_fwd, Shock& shock_rvs, Eqn& eqn,
                           Real rtol = 1e-6) {

    using namespace boost::numeric::odeint;

    typename Eqn::State state;
    Real t0 = 0.01 * unit::sec;
    eqn.set_init_state(state, t0);

    constexpr Real RS_Gamma_limit = 1.03;
    if (state.Gamma <= RS_Gamma_limit) {
        set_stopping_shock(i, j, shock_fwd, state);
        set_stopping_shock(i, j, shock_rvs, state);
        return;
    }

    auto stepper = make_dense_output(rtol, rtol, runge_kutta_dopri5<typename Eqn::State>());
    // auto stepper = bulirsch_stoer_dense_out<typename Eqn::State>{rtol, rtol};
    stepper.initialize(state, t0, 1e-9 * t0);

    size_t k = 0;
    for (; t(k) < t0; k++) {
        eqn.set_init_state(state, t(k));
        save_fwd_shock_state(i, j, k, eqn, state, shock_fwd);
        save_rvs_shock_state(i, j, k, eqn, state, shock_rvs);
    }

    bool reverse_shock_crossing = true;
    while (stepper.current_time() <= t.back()) {
        stepper.do_step(eqn);
        while (k < t.size() && stepper.current_time() > t(k)) {
            stepper.calc_state(t(k), state);
            if (reverse_shock_crossing && !is_crossing(eqn, state, t(k))) {
                shock_rvs.injection_idx(i, j) = k;
                reverse_shock_crossing = false;
                eqn.save_cross_state(state);
            }
            save_fwd_shock_state(i, j, k, eqn, state, shock_fwd);
            save_rvs_shock_state(i, j, k, eqn, state, shock_rvs);
            ++k;
        }
    }
    reverse_shock_early_extrap(i, j, shock_rvs);
}

template <typename Ejecta, typename Medium>
ShockPair generate_shock_pair(Coord const& coord, Medium const& medium, Ejecta const& jet, RadParams const& rad_fwd,
                              RadParams const& rad_rvs, Real rtol) {
    auto [phi_size, theta_size, t_size] = coord.shape();
    const size_t phi_size_needed = coord.t.shape()[0];
    Shock f_shock(phi_size_needed, theta_size, t_size, rad_fwd);
    Shock r_shock(phi_size_needed, theta_size, t_size, rad_rvs);
    for (size_t i = 0; i < phi_size_needed; ++i) {
        // Real theta_s =
        //     jet_spreading_edge(jet, medium, coord.phi(i), coord.theta.front(), coord.theta.back(), coord.t.front());
        for (size_t j = 0; j < theta_size; ++j) {
            auto eqn_r = FRShockEqn(medium, jet, coord.phi(i), coord.theta(j), rad_fwd, rad_rvs);
            // Solve the forward-reverse shock shell
            grid_solve_shock_pair(i, j, xt::view(coord.t, i, j, xt::all()), f_shock, r_shock, eqn_r, rtol);
        }
    }
    return std::make_pair(std::move(f_shock), std::move(r_shock));
}
