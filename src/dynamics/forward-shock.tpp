//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once
#include "forward-shock.hpp"
#include "simple-shock.hpp"

template <typename Ejecta, typename Medium>
ForwardShockEqn<Ejecta, Medium>::ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta,
                                                 RadParams const& rad_params, Real theta_s)
    : medium(medium),
      ejecta(ejecta),
      phi(phi),
      theta0(theta),
      rad(rad_params),
      dOmega0(1 - std::cos(theta)),
      theta_s(theta_s){
    m_jet0 = ejecta.eps_k(phi, theta0) / ejecta.Gamma0(phi, theta0) / con::c2;
    if constexpr (HasSigma<Ejecta>) {
        m_jet0 /= 1 + ejecta.sigma0(phi, theta0);
    }
}

template <typename Ejecta, typename Medium>
void ForwardShockEqn<Ejecta, Medium>::operator()(State const& state, State& diff, Real t) const noexcept {
    const Real beta = gamma_to_beta(state.Gamma);

    diff.r = compute_dr_dt(beta);
    diff.t_comv = compute_dt_dt_comv(state.Gamma, beta);

    if (ejecta.spreading && state.theta < 0.5 * con::pi) {
        diff.theta = compute_dtheta_dt(theta_s, state.theta, diff.r, state.r, state.Gamma);
    } else {
        diff.theta = 0;
    }

    if constexpr (State::mass_inject) {
        diff.m_jet = ejecta.dm_dt(phi, theta0, t);
    }

    if constexpr (State::energy_inject) {
        diff.eps_jet = ejecta.deps_dt(phi, theta0, t);
    }

    Real rho = medium.rho(phi, state.theta, state.r);
    diff.m2 = state.r * state.r * rho * diff.r;
    const Real e_th = (state.Gamma - 1) * 4 * state.Gamma * rho * con::c2;
    const Real eps_rad = compute_radiative_efficiency(state.t_comv, state.Gamma, e_th, rad);
    const Real ad_idx = adiabatic_idx(state.Gamma);
    diff.Gamma = compute_dGamma_dt(state, diff, ad_idx);
    diff.U2_th = compute_dU_dt(eps_rad, state, diff, ad_idx);
}

template <typename Ejecta, typename Medium>
Real ForwardShockEqn<Ejecta, Medium>::compute_dGamma_dt(State const& state, State const& diff,
                                                        Real ad_idx) const noexcept {
    Real dm_dt_swept = diff.m2;
    Real m_swept = state.m2;
    const Real Gamma2 = state.Gamma * state.Gamma;
    const Real Gamma_eff = (ad_idx * (Gamma2 - 1) + 1) / state.Gamma;
    Real dGamma_eff = (ad_idx * (Gamma2 + 1) - 1) / Gamma2;
    Real dlnVdt = 3 / state.r * diff.r; // only r term

    Real m_jet = this->m_jet0;
    Real U = state.U2_th; // Internal energy per unit solid angle

    if (ejecta.spreading) {
        const Real cos_theta = std::cos(state.theta);
        const Real sin_theta = std::sin(state.theta);
        const Real f_spread = (1 - cos_theta) / dOmega0;
        dm_dt_swept = dm_dt_swept * f_spread + m_swept / dOmega0 * sin_theta * diff.theta;
        m_swept *= f_spread;
        dlnVdt += sin_theta / (1 - cos_theta) * diff.theta;
        U *= f_spread;
    }

    Real a1 = -(state.Gamma - 1) * (Gamma_eff + 1) * con::c2 * dm_dt_swept;
    const Real a2 = (ad_idx - 1) * Gamma_eff * U * dlnVdt;

    if constexpr (State::energy_inject) {
        a1 += diff.eps_jet;
    }

    if constexpr (State::mass_inject) {
        a1 -= state.Gamma * diff.m_jet * con::c2;
        m_jet = state.m_jet;
    }

    const Real b1 = (m_jet + m_swept) * con::c2;
    const Real b2 = (dGamma_eff + Gamma_eff * (ad_idx - 1) / state.Gamma) * U;

    return (a1 + a2) / (b1 + b2);
}

template <typename Ejecta, typename Medium>
Real ForwardShockEqn<Ejecta, Medium>::compute_dU_dt(Real eps_rad, State const& state, State const& diff,
                                                    Real ad_idx) const noexcept {
    Real dm_dt_swept = diff.m2;
    const Real m_swept = state.m2;
    Real dlnVdt = 3 / state.r * diff.r - diff.Gamma / state.Gamma;
    if (ejecta.spreading) {
        const Real factor = std::sin(state.theta) / (1 - std::cos(state.theta)) * diff.theta;
        dm_dt_swept = dm_dt_swept + m_swept * factor;
        dlnVdt += factor;
        dlnVdt += factor / (ad_idx - 1);
    }

    return (1 - eps_rad) * (state.Gamma - 1) * con::c2 * dm_dt_swept - (ad_idx - 1) * dlnVdt * state.U2_th;
}

template <typename Ejecta, typename Medium>
void ForwardShockEqn<Ejecta, Medium>::set_init_state(State& state, Real t0) const noexcept {
    Real Gamma4 = ejecta.Gamma0(phi, theta0);

    const Real beta4 = gamma_to_beta(Gamma4);
    state.r = beta4 * con::c * t0 / (1 - beta4);

    state.t_comv = state.r / std::sqrt(Gamma4 * Gamma4 - 1) / con::c;

    state.theta = theta0;

    state.m2 = medium.rho(phi, theta0, state.r) * state.r * state.r * state.r / 3;

    state.Gamma = Gamma4;

    if constexpr (State::energy_inject) {
        state.eps_jet = ejecta.eps_k(phi, theta0);
    }

    if constexpr (State::mass_inject) {
        state.m_jet = m_jet0;
    }

    Real ad_idx = adiabatic_idx(state.Gamma);

    state.U2_th = (state.Gamma - 1) * state.m2 * con::c2 / ad_idx;
}

template <typename Eqn, typename State>
void save_fwd_shock_state(size_t i, size_t j, size_t k, Eqn const& eqn, State const& state, Shock& shock) {
    // Set constant parameters for the unshocked medium
    constexpr Real gamma1 = 1; // Lorentz factor of unshocked medium (at rest)
    constexpr Real sigma = 0;  // Magnetization of unshocked medium
    constexpr Real B_upstr = 0;

    const Real comp_ratio = compute_compression(gamma1, state.Gamma, sigma);
    const Real rho = eqn.medium.rho(eqn.phi, state.theta, state.r);

    Real U_th = 0;
    if constexpr (HasU<State>) {
        U_th = state.U2_th;
    } else {
        U_th = (state.Gamma - 1) * state.m2 * con::c2;
    }

    const Real Gamma_th = compute_Gamma_therm(U_th, state.m2);

    const Real B = compute_downstr_B(shock.rad.eps_B, rho, B_upstr, Gamma_th, comp_ratio);

    write_shock_state(shock, i, j, k, state.t_comv, state.r, state.theta, state.Gamma, Gamma_th, B, state.m2);
}

template <typename FwdEqn, typename View>
void grid_solve_fwd_shock(size_t i, size_t j, View const& t, Shock& shock, FwdEqn const& eqn, Real rtol) {
    using namespace boost::numeric::odeint;

    // Initialize state array
    typename FwdEqn::State state;

    // Get initial time and set up initial conditions
    Real t0 = std::min(t.front(), 1 * unit::sec);
    eqn.set_init_state(state, t0);

    // Early exit if the initial Lorentz factor is below cutoff
    if (state.Gamma <= con::Gamma_cut) {
        set_stopping_shock(i, j, shock, state);
        return;
    }

    // Set up the ODE solver with adaptive step size control
    auto stepper = make_dense_output(rtol, rtol, runge_kutta_dopri5<typename FwdEqn::State>());

    stepper.initialize(state, t0, 0.01 * t0);

    // Solve ODE and update the shock state at each requested time point
    for (size_t k = 0; stepper.current_time() <= t.back();) {
        // Advance solution by one adaptive step
        stepper.do_step(eqn);

        // Update shock state for all-time points that have been passed in this step
        while (k < t.size() && stepper.current_time() > t(k)) {
            stepper.calc_state(t(k), state);
            save_fwd_shock_state(i, j, k, eqn, state, shock);
            ++k;
        }
    }
}

template <typename Ejecta, typename Medium>
Shock generate_fwd_shock(Coord const& coord, Medium const& medium, Ejecta const& jet, RadParams const& rad_params,
                         Real rtol) {
    auto [phi_size, theta_size, t_size] = coord.shape(); // Unpack coordinate dimensions
    const size_t phi_size_needed = coord.t.shape()[0];
    Shock shock(phi_size_needed, theta_size, t_size, rad_params);

    for (size_t i = 0; i < phi_size_needed; ++i) {
        Real theta_s = 0;
        if (jet.spreading) {
            theta_s =
                jet_spreading_edge(jet, medium, coord.phi(i), coord.theta.front(), coord.theta.back(), coord.t.front());
        }
        for (size_t j = 0; j < theta_size; ++j) {
            auto eqn = ForwardShockEqn(medium, jet, coord.phi(i), coord.theta(j), rad_params, theta_s);
            // auto eqn = SimpleShockEqn(medium, jet, coord.phi(i), coord.theta(j), rad_params, theta_s);
            //          Solve the shock shell for this theta slice
            grid_solve_fwd_shock(i, j, xt::view(coord.t, i, j, xt::all()), shock, eqn, rtol);
        }
    }

    return shock;
}
