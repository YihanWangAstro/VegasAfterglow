//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _SIMPLESHOCK_
#define _SIMPLESHOCK_
#include "shock.h"
/********************************************************************************************************************
 * CLASS: SimpleShockEqn
 * DESCRIPTION: Represents the forward shock equation for a given Jet and Injector. It defines a state vector
 *              (an array of 5 Reals) and overloads operator() to compute the derivatives of the state with
 *              respect to radius t. It also declares helper functions for the derivatives. Simple version from
 *              Huang et al. 2000
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
class SimpleShockEqn {
   public:
    using StateArray = std::array<Real, 5>;  // State vector: typically [Gamma, u, r, t_com, theta_jet]

    SimpleShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, Real phi, Real theta, Real eps_e);

    Medium const& medium;     // Reference to the medium properties
    Jet const& jet;           // Reference to the jet properties
    Injector const& inject;   // Reference to the injector properties
    Real const phi{0};        // Angular coordinate phi
    Real const theta0{0};     // Angular coordinate theta
    Real const eps_e{0};      // Electron energy fraction
    Real const jet_sigma{0};  // Jet magnetization parameter
    Real gamma0{1};           // Initial Lorentz factor (or a related parameter)

    // Overloaded operator() to compute the derivatives of the state vector with respect to radius r.
    void operator()(StateArray const& y, StateArray& dydt, Real t);

   private:
    // Helper function: computes the derivative of Gamma with respect to t.
    inline Real dGammadt(Real t, Real Gamma, Real r, Real theta, Real drdt, Real dthetadt, Real rho);
    Real const dM0{0};         // Initial mass per unit solid angle
    Real const inj_Gamma0{0};  // Initial Gamma from the injector
    Real const inj_sigma{0};   // Injector magnetization parameter
    Real const dOmega0{0};     // Initial solid angle
};

/********************************************************************************************************************
 * METHOD: SimpleShockEqn::operator()(State const& y, State& dydr, Real t)
 * DESCRIPTION: Computes the derivatives of the state variables with respect to radius r.
 *              The state vector components are:
 *                  y[0] - Gamma (Lorentz factor)
 *                  y[1] - u (internal energy per solid angle)
 *                  y[2] - r (radius)
 *                  y[3] - t_com (co-moving time) [unused here]
 *                  y[4] - theta_jet (jet opening angle)
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
void SimpleShockEqn<Jet, Injector>::operator()(StateArray const& y, StateArray& dydt, Real t) {
    constFState state(y);

    Real rho = medium.rho(state.r);        // Get medium density at radius r
    Real beta = gammaTobeta(state.Gamma);  // Convert Gamma to beta (velocity/c)
    Real uv = state.Gamma * beta;
    Real beta4 = gammaTobeta(gamma0);  // Convert gamma4 to beta

    dydt[2] = drdt(beta);  // Compute derivative of r with respect to t

    if (jet.spreading && state.theta < 0.5 * con::pi && uv * state.theta < 0.5) {
        dydt[4] = dtheta_dt(uv, dydt[2], state.r, state.Gamma);
    } else {
        dydt[4] = 0;
    }

    dydt[0] = dGammadt(t, state.Gamma, state.r, state.theta, dydt[2], dydt[4], rho);  // d(Gamma)/dt
    dydt[1] = 0;
    dydt[3] = dtdt_CoMoving(state.Gamma);  // d(t_com)/dt
}

/********************************************************************************************************************
 * CONSTRUCTOR: SimpleShockEqn::SimpleShockEqn
 * DESCRIPTION: Initializes a SimpleShockEqn object with references to the medium, jet, and injector,
 *              along with the angular coordinates and energy fraction.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
SimpleShockEqn<Jet, Injector>::SimpleShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, Real phi,
                                              Real theta, Real eps_e)
    : medium(medium),
      jet(jet),
      inject(inject),
      phi(phi),
      theta0(theta),
      eps_e(eps_e),
      jet_sigma(jet.sigma0(phi, theta, 0)),
      gamma0(jet.Gamma0(phi, theta, 0)),
      dM0(jet.dEdOmega(phi, theta, 0) / (gamma0 * (1 + jet_sigma) * con::c2)),
      inj_Gamma0(inject.Gamma0(phi, theta, 0)),
      inj_sigma(inject.sigma0(phi, theta, 0)),
      dOmega0(1 - std::cos(theta0)) {}

/********************************************************************************************************************
 * METHOD: SimpleShockEqn::dGammadt
 * DESCRIPTION: Computes the derivative of Gamma with respect to radius t.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Real SimpleShockEqn<Jet, Injector>::dGammadt(Real t, Real Gamma, Real r, Real theta, Real drdt, Real dthetadt,
                                             Real rho) {
    Real dm = medium.mass(r) / (4 * con::pi);  // Mass per unit solid angle from medium
    Real dmdt = r * r * rho * drdt;
    Real dm_inj = inject.dEdOmega(phi, theta0, t) / (inj_Gamma0 * (1 + inj_sigma) * con::c2);  // Injected mass
    Real L_inj = inject.dLdOmega(phi, theta0, t);  // Injected luminosity per unit solid angle

    if (jet.spreading) {
        Real f_spread = (1 - std::cos(theta)) / dOmega0;
        dmdt = dmdt * f_spread + dm / dOmega0 * std::sin(theta) * dthetadt;
        dm *= f_spread;
    }

    double a1 = (1 - Gamma * Gamma) * dmdt;
    double a2 = L_inj / con::c2 * (1 - Gamma / (inj_Gamma0 * (1 + inj_sigma)));
    return (a1 + a2) / (dM0 + dm_inj + eps_e * dm + 2 * (1 - eps_e) * Gamma * dm);
}

#endif