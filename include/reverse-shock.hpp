//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _REVERSESHOCK_
#define _REVERSESHOCK_
#include "shock.h"
/********************************************************************************************************************
 * CLASS: FRShockEqn
 * DESCRIPTION: Represents the reverse shock (or forward-reverse shock) equation for a given Jet and Injector.
 *              It defines a state vector (an array of 5 Reals) and overloads operator() to compute the
 *              derivatives of the state with respect to radius r. It also declares a helper function to compute
 *              the derivative of N3 (number per solid angle) with respect to r.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
class FRShockEqn {
   public:
    using State = std::array<Real, 5>;  // State vector for reverse shock variables

    FRShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, Real phi, Real theta);

    Medium const& medium;     // Reference to the medium properties
    Jet const& jet;           // Reference to the jet properties
    Injector const& inject;   // Reference to the injector properties
    Real const phi{0};        // Angular coordinate phi
    Real const theta{0};      // Angular coordinate theta
    Real const jet_sigma{0};  // Jet magnetization parameter
    Real gamma4{1};           // Initial Gamma parameter from the jet

    // Overloaded operator() to compute the derivatives of the state vector with respect to time t.
    void operator()(State const& y, State& dydt, Real t);

   private:
    // Helper function: computes the derivative of N3 (number per solid angle) with respect to t.
    inline Real dN3dtPerOmega(Real r, Real n1, Real n4, Real gamma3);
};

/********************************************************************************************************************
 * INLINE FUNCTION: calc_gamma3
 * DESCRIPTION: Computes gamma3 for the reverse shock based on radius, upstream and downstream densities,
 *              the jet Lorentz factor (gamma4), and magnetization (sigma).
 ********************************************************************************************************************/
inline Real calc_gamma3(Real n1, Real n4, Real gamma4, Real sigma) {
    Real C = n4 / n1 * (1 + sigma);
    Real gamma3 =
        gamma4 * std::sqrt(((C - 2) - 2 * std::sqrt(C * (gamma4 * gamma4 - 1) + 1)) / (C - 4 * gamma4 * gamma4));
    return gamma3;
}

/********************************************************************************************************************
 * FUNCTION: setReverseInit
 * DESCRIPTION: Initializes the state vector for reverse shock evolution at the given radius t0.
 ********************************************************************************************************************/
template <typename ShockEqn>
void setReverseInit(ShockEqn& eqn, typename ShockEqn::State& state, Real t0) {
    Real gamma4 = eqn.jet.Gamma0(eqn.phi, eqn.theta, t0);  // Obtain initial Gamma from the jet
    Real beta0 = gammaTobeta(gamma4);
    Real r0 = beta0 * con::c * t0 / (1 - beta0);
    Real t_com0 = r0 / std::sqrt(gamma4 * gamma4 - 1) / con::c;
    Real D_jet0 = con::c * eqn.jet.duration;
    Real dN3dOmega = 0;  // Initialize number per unit solid angle to zero
    state = {0., dN3dOmega, r0, t_com0, D_jet0};
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::dN3dtPerOmega
 * DESCRIPTION: Computes the derivative of N3 (number per unit solid angle) with respect to time.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Real FRShockEqn<Jet, Injector>::dN3dtPerOmega(Real r, Real n1, Real n4, Real gamma3) {
    Real gamma34 = (gamma4 / gamma3 + gamma3 / gamma4) / 2;
    Real ratio_u = u_UpStr2u_DownStr(gamma34, this->jet_sigma);

    /* The following commented-out section shows an alternative computation.
    Real ad_idx2 = adiabatic_index(Gamma);
    Real ad_idx3 = adiabatic_index(Gamma34);
    Real u3s_ = u_down_str(Gamma34, this->sigma);
    Real u4s_ = u_up_str(u3s_, Gamma34);
    Real n2 = n_down_str(n1, Gamma, this->sigma);
    Real e2 = e_thermal_down_str(Gamma, n2);
    Real p2 = (ad_idx2 - 1) * e2;
    Real pB4 = calc_pB4(n4, this->sigma);
    Real pB3 = pB4 * ratio_u * ratio_u;
    Real f_a = fa(Gamma34, u3s_, this->sigma);
    Real f_b = ratio_u / ((ad_idx3 * Gamma34 + 1) / (ad_idx3 - 1));
    Real f_c = fc(p2, pB3);
    Real F = f_a * f_b * f_c;
    */
    Real n3 = n4 * ratio_u;
    Real dxdr = 1. / (gamma4 * std::sqrt((1 + this->jet_sigma) * n4 / n1) * std::fabs(1 - gamma4 * n4 / gamma3 / n3));
    return n3 * r * r * gamma3 * dxdr;
}

/********************************************************************************************************************
 * CONSTRUCTOR: FRShockEqn::FRShockEqn
 * DESCRIPTION: Initializes an FRShockEqn object with references to the medium, jet, and injector, and sets the
 *              angular coordinates, jet magnetization, and initial Gamma.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
FRShockEqn<Jet, Injector>::FRShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, Real phi,
                                      Real theta)
    : medium(medium),
      jet(jet),
      inject(inject),
      phi(phi),
      theta(theta),
      jet_sigma(jet.sigma0(phi, theta, 0)),
      gamma4(jet.Gamma0(phi, theta, 0)) {}

/********************************************************************************************************************
 * METHOD: FRShockEqn::operator()(State const& y, State& dydt, Real t)
 * DESCRIPTION: Computes the derivatives for the reverse shock evolution.
 *              The state vector for FRShockEqn is similar to that of ForwardShockEqn.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
void FRShockEqn<Jet, Injector>::operator()(State const& y, State& dydt, Real t) {
    // y[0] is left blank.
    // Real N_per_Omega = y[1];
    Real r = y[2];
    // Real t_com = y[3];
    Real D_jet_lab = y[4];

    Real n4 = calc_n4(jet.dEdOmega(phi, theta, t), gamma4, r, D_jet_lab, jet_sigma);
    Real n1 = medium.rho(r) / con::mp;

    Real gamma3 = calc_gamma3(n1, n4, gamma4, jet_sigma);
    Real beta3 = gammaTobeta(gamma3);
    Real beta4 = gammaTobeta(gamma4);
    dydt[0] = 0;
    dydt[1] = dN3dtPerOmega(r, n1, n4, gamma3);
    dydt[2] = drdt(beta3);
    dydt[3] = dtdt_CoMoving(gamma3, beta3);
    dydt[4] = dDdt_Jet(gamma4, beta4);
}

/********************************************************************************************************************
 * STRUCT: CrossState
 * DESCRIPTION: Represents the state variables at the crossing between forward and reverse shock phases.
 ********************************************************************************************************************/
struct CrossState {
    Real gamma_rel;
    Real r;
    Real column_num_den;
    Real B;
};

/********************************************************************************************************************
 * INLINE FUNCTION: Blandford_McKee
 * DESCRIPTION: Updates the shock state at a grid cell using the Blandford–McKee self-similar solution.
 ********************************************************************************************************************/
inline void Blandford_McKee(size_t i, size_t j, size_t k, Shock& shock, CrossState const& state_c, Real r, Real t_com) {
    Real const g = 2.;
    shock.t_com[i][j][k] = t_com;
    shock.r[i][j][k] = r;  // TODO need to mofiy !!
    // The following lines are alternative formulations (commented out):
    // shock.Gamma[j][k] = (state_c.gamma - 1) * std::pow(r / state_c.r, -g) + 1;
    // shock.n_p[j][k] = state_c.n3 * std::pow(r / state_c.r, -6 * (3 + g) / 7);
    // shock.e_th[j][k] = state_c.e3 * std::pow(r / state_c.r, -8 * (3 + g) / 7);
    // shock.width_eff[j][k] = state_c.D_eff * std::pow(r / state_c.r, (6. * (3. + g) - 14.) / 7.);
    shock.Gamma_rel[i][j][k] = (state_c.gamma_rel - 1) * std::pow(r / state_c.r, -g) + 1;
    shock.column_num_den[i][j][k] = state_c.column_num_den * std::pow(r / state_c.r, -2);
    shock.B[i][j][k] = state_c.B * std::pow(r / state_c.r, (6. * (3. + g) - 14.) / 7.);
}

/********************************************************************************************************************
 * FUNCTION: solveFRShell
 * DESCRIPTION: Solves the evolution of the forward–reverse shock shell over the radius array r, updating both
 *              the forward shock (f_shock) and reverse shock (r_shock) objects accordingly.
 ********************************************************************************************************************/
template <typename FShockEqn, typename RShockEqn>
void solveFRShell(size_t i, size_t j, Array const& t, Shock& f_shock, Shock& r_shock, FShockEqn& eqn_f,
                  RShockEqn& eqn_r, Real rtol = 1e-9) {
    using namespace boost::numeric::odeint;
    Real t0 = t[0];
    Real dt = (t[1] - t[0]) / 100;
    typename FShockEqn::State state;

    auto stepper = bulirsch_stoer_dense_out<typename FShockEqn::State>{0, rtol};
    //  Alternatively:
    // auto stepper = make_dense_output(0, rtol, runge_kutta_dopri5<typename FShockEqn::State>());

    setForwardInit(eqn_f, state, t0);
    if (state[0] <= con::Gamma_cut) {  // If the initial Lorentz factor is too low, exit early
        return;
    }
    Real gamma3 = eqn_f.gamma4;
    Real dN3dOmega = state[1];
    Real r = state[2];
    Real t_com = state[3];
    Real D_jet = state[4];

    bool RS_crossing = reverseShockExists(eqn_r, r, gamma3, t0, D_jet);
    bool crossed = false;

    if (RS_crossing) {
        setReverseInit(eqn_r, state, t0);
    }

    Real n1 = 0, n4 = 0;

    Real dN4dOmega =
        eqn_r.jet.dEdOmega(eqn_r.phi, eqn_r.theta, t0) / (eqn_r.gamma4 * con::mp * con::c2 * (1 + eqn_r.jet_sigma));

    CrossState state_c;
    stepper.initialize(state, t0, dt);
    // Integrate the shell over the radius array r
    Real t_back = t[t.size() - 1];
    for (int k = 0; stepper.current_time() <= t_back;) {
        RS_crossing ? stepper.do_step(eqn_r) : stepper.do_step(eqn_f);

        for (; stepper.current_time() > t[k] && k < t.size(); k++) {
            stepper.calc_state(t[k], state);
            r = state[2];
            t_com = state[3];
            D_jet = state[4];
            n1 = eqn_f.medium.rho(r) / con::mp;
            Real mass = eqn_f.medium.mass(r);
            if (RS_crossing) {
                n4 = calc_n4(eqn_r.jet.dEdOmega(eqn_r.phi, eqn_r.theta, t[k]), eqn_r.gamma4, r, D_jet, eqn_r.jet_sigma);
                gamma3 = calc_gamma3(n1, n4, eqn_r.gamma4, eqn_r.jet_sigma);
            } else {
                gamma3 = state[0];
            }

            Real dM2dOmega = mass / (4 * con::pi);
            updateShockState(f_shock, i, j, k, r, gamma3, t_com, dM2dOmega, n1, eqn_f.jet_sigma);

            if (!crossed && RS_crossing) {
                dN3dOmega = state[1];
                Real gamma34 = (eqn_r.gamma4 / gamma3 + gamma3 / eqn_r.gamma4) / 2;
                Real dM3dOmega = dN3dOmega * con::mp * con::c2;
                updateShockState(r_shock, i, j, k, r, gamma34, t_com, dM3dOmega, n4, eqn_r.jet_sigma);

                crossed = dN3dOmega >= dN4dOmega;
                if (crossed) {
                    RS_crossing = false;
                    state_c = {r_shock.Gamma_rel[i][j][k], r_shock.r[i][j][k], r_shock.column_num_den[i][j][k],
                               r_shock.B[i][j][k]};
                    Real u0 = (gamma3 - 1) * mass / (4 * con::pi) * con::c2;
                    state = {gamma3, u0, r, t_com, D_jet};
                    eqn_f.gamma4 = gamma3;
                    stepper.initialize(state, t[k], dt);
                }
            } else if (!crossed && !RS_crossing) {
                RS_crossing = reverseShockExists(eqn_r, r, gamma3, t[k], D_jet);
                if (RS_crossing) {
                    state = {1, 0., r, t_com, D_jet};
                    stepper.initialize(state, t[k], dt);
                }
            } else {
                Blandford_McKee(i, j, k, r_shock, state_c, r, t_com);
            }
        }
    }
}

// Determines whether a reverse shock exists based on current shock parameters.
template <typename ShockEqn>
bool reverseShockExists(ShockEqn const& eqn, Real r, Real gamma, Real t, Real D_jet) {
    Real n4 = calc_n4(eqn.jet.dEdOmega(eqn.phi, eqn.theta, t), gamma, r, D_jet, eqn.jet_sigma);
    Real n1 = eqn.medium.rho(r) / con::mp;
    return eqn.jet_sigma < 8. / 3 * gamma * gamma * n1 / n4;
}
#endif