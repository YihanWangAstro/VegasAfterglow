//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once

#include "simple-shock.hpp"

template <typename Ejecta, typename Medium>
SimpleShockEqn<Ejecta, Medium>::SimpleShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta,
                                               RadParams const& rad_params, Real theta_s)
    : medium(medium),
      ejecta(ejecta),
      phi(phi),
      theta0(theta),
      rad(rad_params),
      dOmega0(1 - std::cos(theta0)),
      theta_s(theta_s) {
    m_jet0 = ejecta.eps_k(phi, theta0) / ejecta.Gamma0(phi, theta0) / con::c2;
    if constexpr (HasSigma<Ejecta>) {
        m_jet0 /= 1 + ejecta.sigma0(phi, theta0);
    }
}

template <typename Ejecta, typename Medium>
void SimpleShockEqn<Ejecta, Medium>::operator()(State const& state, State& diff, Real t) const noexcept {
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

    diff.Gamma = dGamma_dt(eps_rad, state, diff);
}

template <typename Ejecta, typename Medium>
Real SimpleShockEqn<Ejecta, Medium>::dGamma_dt(Real eps_rad, State const& state, State const& diff) const noexcept {
    Real m_swept = state.m2;
    Real m_jet = this->m_jet0;
    Real dm_dt_swept = diff.m2;

    if (ejecta.spreading) {
        const Real f_spread = (1 - std::cos(state.theta)) / dOmega0;
        dm_dt_swept = dm_dt_swept * f_spread + m_swept / dOmega0 * std::sin(state.theta) * diff.theta;
        m_swept *= f_spread;
    }

    double a1 = (1 - state.Gamma * state.Gamma) * dm_dt_swept;

    if constexpr (State::energy_inject) {
        a1 += diff.eps_jet / con::c2;
    }

    if constexpr (State::mass_inject) {
        a1 -= state.Gamma * diff.m_jet;
        m_jet = state.m_jet;
    }

    return a1 / (m_jet + eps_rad * m_swept + 2 * (1 - eps_rad) * state.Gamma * m_swept);
}

template <typename Ejecta, typename Medium>
void SimpleShockEqn<Ejecta, Medium>::set_init_state(State& state, Real t0) const noexcept {
    state.Gamma = ejecta.Gamma0(phi, theta0);

    const Real beta0 = gamma_to_beta(state.Gamma);
    state.r = beta0 * con::c * t0 / (1 - beta0);

    state.t_comv = state.r / std::sqrt(state.Gamma * state.Gamma - 1) / con::c;

    state.m2 = medium.rho(phi, theta0, state.r) * state.r * state.r * state.r / 3;

    state.theta = theta0;

    if constexpr (State::energy_inject) {
        state.eps_jet = ejecta.eps_k(phi, theta0);
    }

    if constexpr (State::mass_inject) {
        state.m_jet = m_jet0;
    }
}
