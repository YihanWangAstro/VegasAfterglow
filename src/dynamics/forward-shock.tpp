//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once
#include <array>

#include "forward-shock.hpp"

template <typename Ejecta, typename Medium, typename Value>
ForwardShockEqn<Ejecta, Medium, Value>::ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi,
                                                        std::array<Real, lanes> const& thetas,
                                                        RadParams const& rad_params, Real theta_s)
    : medium(medium),
      ejecta(ejecta),
      phi(phi),
      theta0_lanes(thetas),
      rad(rad_params),
      dOmega0_(at_theta0([](Real th) { return 1 - std::cos(th); })),
      theta_s_(theta_s) {
    m_jet0_ = at_theta0([&](Real th) {
        Real m_jet = ejecta.eps_k(phi, th) / ejecta.Gamma0(phi, th) / con::c2;
        if constexpr (HasSigma<Ejecta>) {
            m_jet /= 1 + ejecta.sigma0(phi, th);
        }
        return m_jet;
    });
    eps_rad_ = RadiativeEfficiency(rad);
}

template <typename Ejecta, typename Medium, typename Value>
void ForwardShockEqn<Ejecta, Medium, Value>::operator()(State const& state, State& diff, Real t) const noexcept {
    using std::cos, std::sin, std::sqrt;
    const Value Gamma = state.Gamma;
    const Value u2 = (Gamma - 1) * (Gamma + 1);
    const Value u = sqrt(u2);

    diff.r = compute_dr_dt(Gamma, u);
    diff.t_comv = Gamma + u;

    if (ejecta.spreading) {
        const Value dtheta = compute_dtheta_dt(theta_s_, state.theta, diff.r, state.r, Gamma, u, u2);
        diff.theta = where(state.theta < 0.5 * con::pi, dtheta, 0.);
    } else {
        diff.theta = 0.;
    }

    if constexpr (State::mass_inject) {
        diff.m_jet = at_theta0([&](Real th) { return ejecta.dm_dt(phi, th, t); });
    }

    if constexpr (State::energy_inject) {
        diff.eps_jet = at_theta0([&](Real th) { return ejecta.deps_dt(phi, th, t); });
    }

    const Value rho = lane_rho(medium, phi, state.theta, state.r);
    diff.m2 = state.r * state.r * rho * diff.r;
    const Value e_th = (Gamma - 1) * 4 * Gamma * rho * con::c2;
    const Value eps_rad = eps_rad_(state.t_comv, Gamma, e_th);
    const Value ad_idx = physics::thermo::adiabatic_idx(Gamma);
    Value sin_theta = 0.;
    Value cos_theta = 1.;
    if (ejecta.spreading) {
        sin_theta = sin(state.theta);
        cos_theta = cos(state.theta);
    }
    diff.Gamma = compute_dGamma_dt(state, diff, ad_idx, sin_theta, cos_theta);
    diff.U2_th = compute_dU_dt(eps_rad, state, diff, ad_idx, sin_theta, cos_theta);
}

template <typename Ejecta, typename Medium, typename Value>
Value ForwardShockEqn<Ejecta, Medium, Value>::compute_dGamma_dt(State const& state, State const& diff,
                                                                Value const& ad_idx, Value const& sin_theta,
                                                                Value const& cos_theta) const noexcept {
    Value dm_dt_swept = diff.m2;
    Value m_swept = state.m2;
    const Value Gamma2 = state.Gamma * state.Gamma;
    const Value Gamma_eff = (ad_idx * (Gamma2 - 1) + 1) / state.Gamma;
    Value dGamma_eff = (ad_idx * (Gamma2 + 1) - 1) / Gamma2;
    Value dlnVdt = 3 / state.r * diff.r; // only r term

    Value m_jet = this->m_jet0_;
    Value U = state.U2_th; // Internal energy per unit solid angle

    if (ejecta.spreading) {
        const Value f_spread = (1 - cos_theta) / dOmega0_;
        dm_dt_swept = dm_dt_swept * f_spread + m_swept / dOmega0_ * sin_theta * diff.theta;
        m_swept *= f_spread;
        dlnVdt += sin_theta / (1 - cos_theta) * diff.theta;
        U *= f_spread;
    }

    Value a1 = -(state.Gamma - 1) * (Gamma_eff + 1) * con::c2 * dm_dt_swept;
    const Value a2 = (ad_idx - 1) * Gamma_eff * U * dlnVdt;

    if constexpr (State::energy_inject) {
        a1 += diff.eps_jet;
    }

    if constexpr (State::mass_inject) {
        a1 -= state.Gamma * diff.m_jet * con::c2;
        m_jet = state.m_jet;
    }

    const Value b1 = (m_jet + m_swept) * con::c2;
    const Value b2 = (dGamma_eff + Gamma_eff * (ad_idx - 1) / state.Gamma) * U;

    return (a1 + a2) / (b1 + b2);
}

template <typename Ejecta, typename Medium, typename Value>
Value ForwardShockEqn<Ejecta, Medium, Value>::compute_dU_dt(Value const& eps_rad, State const& state, State const& diff,
                                                            Value const& ad_idx, Value const& sin_theta,
                                                            Value const& cos_theta) const noexcept {
    Value dm_dt_swept = diff.m2;
    const Value m_swept = state.m2;
    Value dlnVdt = 3 / state.r * diff.r - diff.Gamma / state.Gamma;
    if (ejecta.spreading) {
        const Value factor = sin_theta / (1 - cos_theta) * diff.theta;
        dm_dt_swept = dm_dt_swept + m_swept * factor;
        dlnVdt += factor;
        dlnVdt += factor / (ad_idx - 1);
    }

    return (1 - eps_rad) * (state.Gamma - 1) * con::c2 * dm_dt_swept - (ad_idx - 1) * dlnVdt * state.U2_th;
}

template <typename Ejecta, typename Medium, typename Value>
void ForwardShockEqn<Ejecta, Medium, Value>::set_init_state_lane(ScalarState& state, Real t0, size_t w) const noexcept {
    const Real theta = theta0_lanes[w];
    const Real Gamma4 = ejecta.Gamma0(phi, theta);

    const Real beta4 = physics::relativistic::gamma_to_beta(Gamma4);
    // beta / (1 - beta) = beta * Gamma^2 * (1 + beta): avoids 1 - beta, which
    // loses ~2 Gamma^2 ulp of relative precision for ultrarelativistic shells.
    state.r = beta4 * con::c * t0 * Gamma4 * Gamma4 * (1 + beta4);

    state.t_comv = state.r / std::sqrt((Gamma4 - 1) * (Gamma4 + 1)) / con::c;

    state.theta = theta;

    state.m2 = enclosed_mass_medium(medium, phi, theta, state.r);

    // The mass swept before t0 (a finite fraction of the deceleration mass in winds)
    // has already decelerated the shell below Gamma0.
    const auto init = initial_state_energy_conserving(
        Gamma4, lane::get(m_jet0_, w), state.m2, rad.radiative ? rad.eps_e : 0.0,
        enclosed_thermal_integral_medium(medium, phi, theta, state.r, physics::thermo::adiabatic_idx(Gamma4)));
    state.Gamma = init.Gamma;

    if constexpr (State::energy_inject) {
        state.eps_jet = ejecta.eps_k(phi, theta);
    }

    if constexpr (State::mass_inject) {
        state.m_jet = lane::get(m_jet0_, w);
    }

    state.U2_th = init.U_th;
}

template <typename Ejecta, typename Medium, typename Value>
void ForwardShockEqn<Ejecta, Medium, Value>::set_init_state(State& state, Real t0) const noexcept {
    for (size_t w = 0; w < lanes; ++w) {
        ScalarState sw;
        set_init_state_lane(sw, t0, w);
        insert_lane(state, w, sw);
    }
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

template <typename FwdEqn>
void grid_solve_fwd_shock(size_t i, std::array<size_t, FwdEqn::lanes> const& js, size_t n_live, Coord const& coord,
                          Shock& shock, FwdEqn const& eqn, Real rtol) {
    using namespace boost::numeric::odeint;
    using State = typename FwdEqn::State;
    constexpr size_t lanes = FwdEqn::lanes;

    const size_t t_size = coord.t.shape()[2];
    const auto [t_row, t0, t_end] = batch_time_span(coord, i, js, n_live, eqn, 0.1 * unit::sec);

    State state;
    eqn.set_init_state(state, t0);

    const auto live = mark_live_lanes<lanes>(
        state, n_live, con::Gamma_cut, [&](size_t w) { set_stopping_shock(i, js[w], shock, extract_lane(state, w)); });
    if (!lane::any_of(live)) {
        return;
    }

    auto stepper = make_dense_output(rtol, rtol, runge_kutta_dopri5<State, Real, State, Real, LaneAlgebra>());
    stepper.initialize(state, t0, 0.01 * t0);

    std::array<size_t, lanes> k{};
    for (size_t steps = 0; stepper.current_time() <= t_end;) {
        stepper.do_step(eqn);
        if (++steps > defaults::solver::max_ode_steps) {
            std::fprintf(stderr, "Warning: forward shock ODE exceeded %zu steps at (i=%zu, j=%zu), giving up\n",
                         defaults::solver::max_ode_steps, i, js[0]);
            return;
        }

        for (size_t w = 0; w < n_live; ++w) {
            if (!live[w]) {
                continue;
            }
            while (k[w] < t_size && stepper.current_time() > t_row[w][k[w]]) {
                stepper.calc_state(t_row[w][k[w]], state);
                save_fwd_shock_state(i, js[w], k[w], eqn, extract_lane(state, w), shock);
                ++k[w];
            }
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

        for_each_theta_batch(coord, [&](auto const& js, auto const& thetas, size_t n_live) {
            constexpr size_t width = std::tuple_size_v<std::remove_cvref_t<decltype(js)>>;
            using Value = std::conditional_t<width == 1, Real, lane::Batch>;
            ForwardShockEqn<Ejecta, Medium, Value> eqn(medium, jet, coord.phi(i), thetas, rad_params, theta_s);
            grid_solve_fwd_shock(i, js, n_live, coord, shock, eqn, rtol);
        });

        if (coord.symmetry >= Symmetry::phi_symmetric) {
            shock.broadcast_groups(coord);
            return shock;
        }
    }

    return shock;
}
