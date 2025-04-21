//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include "shock.h"
#include "utilities.h"
/********************************************************************************************************************
 * CLASS: RState
 * DESCRIPTION: This is a helper class to provide named access to the reverse shock ODE state array (required by ODE
 *              solver) components.
 ********************************************************************************************************************/
template <typename StateArray>
struct RState {
    using T = decltype(std::declval<StateArray>()[0]);
    RState() = delete;

    RState(StateArray& y)
        : width(y[0]), N3(y[1]), r(y[2]), t_com(y[3]), theta(y[4]), M_sw(y[5]), M_ej(y[6]), E_ej(y[7]) {}
    T& width;
    T& N3;
    T& r;
    T& t_com;
    T& theta;
    T& M_sw;
    T& M_ej;
    T& E_ej;
};

/********************************************************************************************************************
 * CLASS: FRShockEqn
 * DESCRIPTION: Represents the reverse shock (or forward-reverse shock) equation for a given Jet and medium.
 *              It defines a state vector (an array of 8 Reals) and overloads operator() to compute the
 *              derivatives of the state with respect to radius r. It also declares a helper function to compute
 *              the derivative of N3 (number per solid angle) with respect to t.
 ********************************************************************************************************************/
class FRShockEqn {
   public:
    // State vector for reverse shock variables [D_jet, N3, r, t_com, theta]
    // - D_jet: comoving shell width,
    // - N3: electron number per unit solid angle in region 3
    // - r: radius
    // - t_com: comoving time
    // - theta: jet opening angle
    // - M_sw: swept mass per solid angle
    // - M_ej: ejecta mass per solid angle
    // - E_ej: ejecta energy per solid angle
    using StateArray = std::array<Real, 8>;

    FRShockEqn(Medium const& medium, Ejecta const& jet, Real phi, Real theta, Real eps_e);

    Medium const& medium;  // Reference to the medium properties
    Ejecta const& ejecta;  // Reference to the jet properties
    Real const phi{0};     // Angular coordinate phi
    Real const theta0{0};  // Angular coordinate theta
    Real const eps_e{0};   // Electron energy fraction
    Real gamma4{1};        // initial Lorentz factor of the jet
    Real u_x{0};           // reverse shock crossed four velocity
    Real r_x{0};           // reverse shock crossed radius

    // Reverse shock ODE equation
    void operator()(StateArray const& y, StateArray& dydt, Real t);

    // Set the shock state when the reverse shock crosses the jet.
    template <typename State>
    void setCrossState(State const& state, Real B, Real t);

    // calculate the Gamma3 during the shock crossing phase.
    template <typename State>
    Real crossingGamma3(State const& state, Real t) const;

    // calculate the Gamma_43 post shock crossing.
    template <typename State>
    Real crossedGamma_rel(State const& state) const;

    // calculate the magnetic field post shock crossing.
    template <typename State>
    Real crossedB(State const& state) const;

    // calculate the Gamma3 post shock crossing.
    Real crossedGamma3(Real gamma_rel, Real r) const;

   private:
    Real N0{0};               // normalized total electron (for post crossing scaling calculation).
    Real adiabatic_const{1};  // normalized adiabatic constant where C = rho^idx/p.
    Real Emag_const{1};       // normalized magnetic energy constant where C = B^2/p.
    Real ad_idx0{4. / 3};     // adiabatic index at the shock crossing.
    bool crossed{false};
};

Real calc_gamma3(Real gamma4, Real N2, Real N3, Real sigma);

inline bool isInjecting(FRShockEqn const& eqn, Real t0) {
    Real dE_dt = eqn.ejecta.dEdtdOmega(eqn.phi, eqn.theta0, t0);
    Real dM_dt = eqn.ejecta.dMdtdOmega(eqn.phi, eqn.theta0, t0);
    return t0 < eqn.ejecta.T0 || dM_dt != 0;
}

/********************************************************************************************************************
 * FUNCTION: calcInitN3
 * DESCRIPTION: calculate the initial N3 at t0 (to ensure the energy conservation)
 ********************************************************************************************************************/
template <typename State>
inline auto calcInitN3(FRShockEqn const& eqn, State const& state, Real gamma4, Real t0) {
    Real E_iso = eqn.ejecta.dE0dOmega(eqn.phi, eqn.theta0) * 4 * con::pi;
    Real n1 = eqn.medium.rho(eqn.phi, eqn.theta0, state.r) / con::mp;
    Real l = SedovLength(E_iso, n1);
    Real sigma0 = eqn.ejecta.sigma0(eqn.phi, eqn.theta0);
    Real beta4 = gammaTobeta(gamma4);

    Real Delta0 = con::c * eqn.ejecta.T0;  // lab frame shell width
    Real r_x = std::sqrt(std::sqrt(Delta0 * l * l * l / 3) / (1 + sigma0) * 3. / 4);
    Real Rs = con::c * eqn.ejecta.T0 * beta4 / (1 - beta4);

    if (r_x > Rs) {
        constexpr Real rtol = 1e-3;
        Real r_x_new = (1 + 2 * rtol) * r_x;
        for (; std::abs((r_x - r_x_new) / r_x) > rtol;) {
            Delta0 = con::c * eqn.ejecta.T0 + (r_x_new - Rs) / (std::sqrt(3) * gamma4 * gamma4);
            r_x = r_x_new;
            r_x_new = std::sqrt(std::sqrt(Delta0 * l * l * l / 3) / (1 + sigma0) * 3. / 4);
        }
    }

    Real N4_tot = E_iso / (4 * con::pi * gamma4 * (1 + sigma0) * con::mp * con::c2);
    Real N40 = state.M_ej / con::mp;

    if (!isInjecting(eqn, t0)) {
        Real N3 = N4_tot * std::pow(state.r / r_x, 1.5);
        if (N3 >= N40) {
            return std::make_tuple(true, N40);
        } else {
            return std::make_tuple(false, N3);
        }
    } else {
        Real N3 = N4_tot * state.r / r_x;
        return std::make_tuple(false, std::min(N3, N40));
    }
}

// comoving shell width at r (including the spreading from r=0 to r=r0)
inline Real initShellWidth(Real gamma, Real t0, Real T) {
    Real beta = gammaTobeta(gamma);
    if (t0 < T) {  // pure injection
        return gamma * t0 * beta * con::c;
    } else {  // injection+shell spreading
        Real cs = con::c / std::sqrt(3);
        return gamma * T * beta * con::c + cs * (t0 - T) * gamma;
    }
}
// inline Real initShellWidth(Real gamma, Real T, Real r) { return r / gamma / std::sqrt(3); }
/********************************************************************************************************************
 * FUNCTION: setReverseInit
 * DESCRIPTION: Set the initial conditions for the reverse shock ODE.
 ********************************************************************************************************************/
template <typename State>
bool setReverseInit(FRShockEqn const& eqn, State& state, Real t0) {
    Real gamma4 = eqn.ejecta.Gamma0(eqn.phi, eqn.theta0);  // Obtain initial Gamma from the jet
    Real beta4 = gammaTobeta(gamma4);
    state.theta = eqn.theta0;
    state.width = initShellWidth(gamma4, t0, eqn.ejecta.T0);
    Real e_iso = eqn.ejecta.dE0dOmega(eqn.phi, eqn.theta0);
    Real sigma0 = eqn.ejecta.sigma0(eqn.phi, eqn.theta0);
    state.E_ej = e_iso * std::min(t0 / eqn.ejecta.T0, 1.);
    state.M_ej = state.E_ej / (gamma4 * (1 + sigma0) * con::c2);

    state.r = beta4 * con::c * t0 / (1 - beta4);
    state.t_com = state.r / std::sqrt(gamma4 * gamma4 - 1) / con::c;
    Real rho = eqn.medium.rho(eqn.phi, eqn.theta0, state.r);

    if (t0 < eqn.ejecta.T0) {
        // thick shell deceleration radius where gamma starts to drop propto t^{-1/4}
        Real r_dec =
            0.35 * std::sqrt(3 * e_iso * (1 - beta4) /
                             (beta4 * con::c * eqn.ejecta.T0 * (1 + sigma0) * rho * con::c2 * gamma4 * gamma4));
        // r0 is larger than the decelertation radius, so r=beta4 * con::c * t0 / (1 - beta4) is not appropriate
        if (state.r > r_dec) {
            Real t_dec = r_dec * (1 - beta4) / (beta4 * con::c);
            state.r = std::sqrt(4 * gamma4 * gamma4 * r_dec * (t0 - t_dec) + r_dec * r_dec);
            rho = eqn.medium.rho(eqn.phi, eqn.theta0, state.r);

            Real t_dec_com = r_dec / std::sqrt(gamma4 * gamma4 - 1) / con::c;
            state.t_com =
                (std::pow(state.r, 1.5) - std::pow(r_dec, 1.5)) / (1.5 * gamma4 * std::sqrt(r_dec)) + t_dec_com;
        }
    }

    state.M_sw = 1 / 3. * state.r * state.r * state.r * rho;
    auto [crossed, N3] = calcInitN3(eqn, state, gamma4, t0);
    state.N3 = N3;
    return crossed;

    // Real gamma3 = calc_gamma3(gamma4, state.M_sw / con::mp, N3, eqn.ejecta.sigma0(eqn.phi, eqn.theta0));
}
/********************************************************************************************************************
 * METHOD: FRShockEqn::setCrossState(State const& state, Real B, Real t)
 * DESCRIPTION: Set some constants at the shock crossed point for later scaling calculation and switch the ODE from
 *              shock crossing state to crossed state.
 ********************************************************************************************************************/
template <typename State>
void FRShockEqn::setCrossState(State const& state, Real B, Real t) {
    this->r_x = state.r;
    constexpr Real n3_norm = 1;  // normalized n3/n3_x = 1
    this->N0 = n3_norm * state.width * state.r * state.r;

    Real gamma3 = crossingGamma3(state, t);
    this->u_x = std::sqrt(gamma3 * gamma3 - 1);

    Real gamma_rel = relativeLorentz(gamma4, gamma3);
    Real gamma_electron = (gamma_rel - 1) * con::mp / con::me * eps_e + 1;  // electron Lorentz factor

    // use electron adiabatic index see section 3.2 https://arxiv.org/pdf/astro-ph/9910241
    this->ad_idx0 = adiabaticIndex(gamma_electron);
    Real p_norm = (ad_idx0 - 1) * (gamma_rel - 1) * n3_norm;      // p = (ad-1)(gamma-1)n m_p c^2
    this->adiabatic_const = std::pow(n3_norm, ad_idx0) / p_norm;  // p \propto n^ad
    this->Emag_const = p_norm / (B * B);
    this->crossed = true;
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::crossingGamma3(State const& state, Real t)
 * DESCRIPTION: Shock crossing gamma3 calculation.
 ********************************************************************************************************************/
template <typename State>
Real FRShockEqn::crossingGamma3(State const& state, Real t) const {
    Real N2 = state.M_sw / con::mp;
    Real sigma4 = state.E_ej / (gamma4 * state.M_ej * con::c2) - 1;
    return calc_gamma3(gamma4, N2, state.N3, sigma4);
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::crossedGamma_rel(State const& state)
 * DESCRIPTION: Shock crossing gamma_43 calculation.
 ********************************************************************************************************************/
template <typename State>
Real FRShockEqn::crossedGamma_rel(State const& state) const {
    Real n3_norm = N0 / (state.width * state.r * state.r);        // proton number conservation
    Real p3_norm = std::pow(n3_norm, ad_idx0) / adiabatic_const;  // adiabatic expansion
    return p3_norm / ((ad_idx0 - 1) * n3_norm) + 1;
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::crossedB(State const& state)
 * DESCRIPTION: Post shock crossing magnetic field calculation.
 ********************************************************************************************************************/
template <typename State>
Real FRShockEqn::crossedB(State const& state) const {
    Real n3_norm = N0 / (state.width * state.r * state.r);        // proton number conservation
    Real p3_norm = std::pow(n3_norm, ad_idx0) / adiabatic_const;  // adiabatic expansion
    return std::sqrt(p3_norm / Emag_const);
}

/********************************************************************************************************************
 * INLINE FUNCTION: updateCrossedReverseShock
 * DESCRIPTION: Updates the post shock crossing state based on current ODE solution at t.
 ********************************************************************************************************************/
template <typename Eqn, typename State>
inline void updateCrossedReverseShock(size_t i, size_t j, size_t k, size_t k0, Eqn const& eqn, State const& state,
                                      Shock& shock) {
    Real r0 = shock.r(i, j, k0);
    shock.t_com(i, j, k) = state.t_com;
    shock.r(i, j, k) = state.r;
    shock.theta(i, j, k) = state.theta;
    shock.Gamma_rel(i, j, k) = eqn.crossedGamma_rel(state);
    shock.Gamma(i, j, k) = eqn.crossedGamma3(shock.Gamma_rel(i, j, k), state.r);
    shock.column_num_den(i, j, k) = shock.column_num_den(i, j, k0) * (r0 * r0) / (state.r * state.r);
    shock.B(i, j, k) = eqn.crossedB(state);
}

/********************************************************************************************************************
 * FUNCTION: reverseShockExists
 * DESCRIPTION: Check if reverse shock can be generated.
 ********************************************************************************************************************/
/*bool reverseShockExists(FRShockEqn const& eqn, Real t0) {           // TODO: for general t case
    Real gamma = eqn.ejecta.Gamma(eqn.phi, eqn.theta0, t0);  // Initial Lorentz factor
    Real beta0 = gammaTobeta(gamma);
    Real r0 = beta0 * con::c * t0 / (1 - beta0);
    Real D_jet0 = initShellWidth(gamma, eqn.ejecta.duration, r0);
    Real n4 = calc_n4(eqn.ejecta.dEdOmega(eqn.phi, eqn.theta0, t0), gamma, r0, D_jet0, eqn.jet_sigma);
    Real n1 = eqn.medium.rho(eqn.phi, eqn.theta0, r0) / con::mp;
    return eqn.jet_sigma < 8. / 3 * (gamma * gamma - 1) * n1 / n4;
}*/

/********************************************************************************************************************
 * FUNCTION: set_f_state
 * DESCRIPTION: Set the initial condition for forward shock ODE solver based on the reverse shock state at crossed
 *              point.
 ********************************************************************************************************************/
template <typename FState, typename RState>
void set_f_state(FState& state_fwd, RState const& state_rvs, Real gamma2) {
    state_fwd.r = state_rvs.r;
    state_fwd.t_com = state_rvs.t_com;
    state_fwd.theta = state_rvs.theta;
    state_fwd.Gamma = gamma2;
    state_fwd.M_sw = state_rvs.M_sw;
    state_fwd.M_ej = state_rvs.M_ej;
    state_fwd.E_ej = state_rvs.E_ej;
    state_fwd.u = (gamma2 - 1) * state_fwd.M_sw * con::c2;
}

template <typename Eqn, typename State>
bool updateForwardReverseShock(size_t i, size_t j, int k, Eqn const& eqn_rvs, State const& state, Real t,
                               Shock& shock_fwd, Shock& shock_rvs) {
    Real N4 = state.M_ej / con::mp;
    Real N2 = state.M_sw / con::mp;

    Real n4 = N4 / (state.r * state.r * state.width);
    Real sigma4 = state.E_ej / (eqn_rvs.gamma4 * state.M_ej * con::c2) - 1;
    if (sigma4 < 0) sigma4 = 0;

    Real n1 = eqn_rvs.medium.rho(eqn_rvs.phi, state.theta, state.r) / con::mp;
    constexpr Real gamma1 = 1;
    constexpr Real sigma1 = 0;

    Real gamma3 = calc_gamma3(eqn_rvs.gamma4, N2, state.N3, sigma4);

    updateShockState(shock_fwd, i, j, k, state, gamma3, gamma1, N2, n1, sigma1);
    updateShockState(shock_rvs, i, j, k, state, gamma3, eqn_rvs.gamma4, state.N3, n4, sigma4);
    return state.N3 >= N4 && t > eqn_rvs.ejecta.T0;
}

bool initializeShocks(size_t i, size_t j, Real t0, auto& state_fwd, auto& state_rvs, auto const& eqn_fwd, auto& eqn_rvs,
                      Real& dt) {}

// Handle low Gamma scenario
bool handleLowGamma(size_t i, size_t j, Shock& shock_fwd, Shock& shock_rvs, Array const& t, auto const& state_fwd,
                    auto const& state_rvs, auto const& eqn_rvs) {
    if (state_fwd.Gamma < con::Gamma_cut) {
        setStoppingShock(i, j, shock_fwd, state_fwd);
        setStoppingShock(i, j, shock_rvs, state_rvs);
        return true;
    }
    return false;
}

// Solve reverse shock until it crosses
size_t solveReverseShockUntilCross(size_t i, size_t j, const Array& t, auto& stepper_rvs, auto& eqn_rvs, auto& y_rvs,
                                   auto& state_rvs, Shock& shock_fwd, Shock& shock_rvs) {
    size_t k0 = 0;
    Real t_back = t.back();
    bool crossed = false;
    size_t k = 0;

    while (!crossed && stepper_rvs.current_time() <= t_back) {
        stepper_rvs.do_step(eqn_rvs);
        while (k < t.size() && stepper_rvs.current_time() > t(k)) {
            stepper_rvs.calc_state(t(k), y_rvs);
            crossed = updateForwardReverseShock(i, j, k, eqn_rvs, state_rvs, t(k), shock_fwd, shock_rvs);
            if (crossed) {
                k0 = k;
                shock_rvs.injection_idx(i, j) = k0;
                eqn_rvs.setCrossState(state_rvs, shock_rvs.B(i, j, k), t(k));
                stepper_rvs.initialize(y_rvs, t(k0), stepper_rvs.current_time_step());
                break;
            }
            k++;
        }
    }

    return k0;
}

// Solve forward shock after crossing
void solveForwardShockAfterCross(size_t i, size_t j, const Array& t, size_t k0, auto& stepper_fwd, const auto& eqn_fwd,
                                 auto& y_fwd, auto& state_fwd, Shock& shock_fwd) {
    size_t k = k0 + 1;
    Real t_back = t.back();

    while (stepper_fwd.current_time() <= t_back) {
        stepper_fwd.do_step(eqn_fwd);
        while (k < t.size() && stepper_fwd.current_time() > t(k)) {
            stepper_fwd.calc_state(t(k), y_fwd);
            updateForwardShock(i, j, k, eqn_fwd, state_fwd, shock_fwd);
            k++;
        }
    }
}

// Solve reverse shock after crossing
void solveReverseShockAfterCross(size_t i, size_t j, const Array& t, size_t k0, auto& stepper_rvs, auto& eqn_rvs,
                                 auto& y_rvs, auto& state_rvs, Shock& shock_rvs) {
    size_t k = k0 + 1;
    Real t_back = t.back();

    while (stepper_rvs.current_time() <= t_back) {
        stepper_rvs.do_step(eqn_rvs);
        while (k < t.size() && stepper_rvs.current_time() > t(k)) {
            stepper_rvs.calc_state(t(k), y_rvs);
            updateCrossedReverseShock(i, j, k, k0, eqn_rvs, state_rvs, shock_rvs);
            k++;
        }
    }
}

template <typename StateArray>
struct FState;
/********************************************************************************************************************
 * FUNCTION: solveFRShell
 * DESCRIPTION: Solve the reverse/forward shock ODE at grid (phi[i], theta[j]) as a function of t (on-axis
 *observation time).
 ********************************************************************************************************************/
template <typename FwdEqn, typename RvsEqn>
void solveFRShell(size_t i, size_t j, Array const& t, Shock& shock_fwd, Shock& shock_rvs, FwdEqn const& eqn_fwd,
                  RvsEqn& eqn_rvs, Real rtol = 1e-6) {
    using namespace boost::numeric::odeint;
    auto stepper_fwd = make_dense_output(0, rtol, runge_kutta_dopri5<typename FwdEqn::StateArray>());

    typename RvsEqn::StateArray y_rvs;
    typename FwdEqn::StateArray y_fwd;
    FState state_fwd(y_fwd);
    RState state_rvs(y_rvs);
    Real t0 = t.front();

    shock_fwd.injection_idx(i, j) = t.size();
    shock_rvs.injection_idx(i, j) = t.size();

    Real t_dec = setForwardInit(eqn_fwd, state_fwd, t0);
    Real dt = t_dec / 100;
    bool crossed = setReverseInit(eqn_rvs, state_rvs, t0);

    if (handleLowGamma(i, j, shock_fwd, shock_rvs, t, state_fwd, state_rvs, eqn_rvs)) return;

    if (crossed) {
        solveForwardShell(i, j, t, shock_fwd, eqn_fwd, rtol);
        return;
    }

    auto stepper_rvs = make_dense_output(0, rtol, runge_kutta_dopri5<typename RvsEqn::StateArray>());
    stepper_rvs.initialize(y_rvs, t0, dt);

    size_t k0 = solveReverseShockUntilCross(i, j, t, stepper_rvs, eqn_rvs, y_rvs, state_rvs, shock_fwd, shock_rvs);

    set_f_state(state_fwd, state_rvs, shock_rvs.Gamma(i, j, k0));
    stepper_fwd.initialize(y_fwd, t(k0), stepper_rvs.current_time_step());

    solveForwardShockAfterCross(i, j, t, k0, stepper_fwd, eqn_fwd, y_fwd, state_fwd, shock_fwd);
    solveReverseShockAfterCross(i, j, t, k0, stepper_rvs, eqn_rvs, y_rvs, state_rvs, shock_rvs);
}
