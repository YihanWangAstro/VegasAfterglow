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
 *              respect to radius r. It also declares helper functions for the derivatives. Simple version from
 *              Huang et al. 2000
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
class SimpleShockEqn {
   public:
    using State = std::array<Real, 5>;  // State vector: typically [Gamma, u, t_eng, t_com, D_jet]

    SimpleShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, Real phi, Real theta, Real eps_e);

    Medium const& medium;     // Reference to the medium properties
    Jet const& jet;           // Reference to the jet properties
    Injector const& inject;   // Reference to the injector properties
    Real const phi{0};        // Angular coordinate phi
    Real const theta{0};      // Angular coordinate theta
    Real const eps_e{0};      // Electron energy fraction
    Real const jet_sigma{0};  // Jet magnetization parameter
    Real gamma4{1};           // Initial Lorentz factor (or a related parameter)

    // Overloaded operator() to compute the derivatives of the state vector with respect to radius r.
    void operator()(State const& y, State& dydr, Real r);

   private:
    // Helper function: computes the derivative of Gamma with respect to r.
    inline Real dGammadr(Real r, Real Gamma, Real t_eng, Real rho);
    Real const dM0{0};  // Initial mass per unit solid angle
};

/********************************************************************************************************************
 * METHOD: SimpleShockEqn::operator()(State const& y, State& dydr, Real r)
 * DESCRIPTION: Computes the derivatives of the state variables with respect to radius r.
 *              The state vector components are:
 *                  y[0] - Gamma (Lorentz factor)
 *                  y[1] - u (internal energy-related variable)
 *                  y[2] - t_eng (engine time)
 *                  y[3] - t_com (co-moving time) [unused here]
 *                  y[4] - D_jet (jet shell width) [unused here]
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
void SimpleShockEqn<Jet, Injector>::operator()(State const& y, State& dydr, Real r) {
    Real Gamma = y[0];
    Real t_eng = y[2];  // engine time
    // Real t_com = y[3];  // co-moving time (unused)
    // Real D_jet = y[4];  // co-moving jet shell width (unused)

    Real rho = medium.rho(r);          // Get medium density at radius r
    Real beta = gammaTobeta(Gamma);    // Convert Gamma to beta (velocity/c)
    Real beta4 = gammaTobeta(gamma4);  // Convert gamma4 to beta

    dydr[2] = dtdr_Engine(beta);               // Compute derivative of engine time with respect to r
    dydr[0] = dGammadr(r, Gamma, t_eng, rho);  // d(Gamma)/dr
    dydr[1] = 0;
    dydr[3] = dtdr_CoMoving(Gamma, beta);  // d(t_com)/dr
    dydr[4] = dDdr_Jet(gamma4, beta4);     // d(D_jet)/dr
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
      theta(theta),
      eps_e(eps_e),
      jet_sigma(jet.sigma0(phi, theta, 0)),
      gamma4(jet.Gamma0(phi, theta, 0)),
      dM0(jet.dEdOmega(phi, theta, 0) / (gamma4 * (1 + jet_sigma) * con::c2)) {
    // dM0dOmega(jet.dE0dOmega(theta) / (jet.Gamma0(theta) * con::c2)) is commented out.
}

/********************************************************************************************************************
 * METHOD: SimpleShockEqn::dGammadr
 * DESCRIPTION: Computes the derivative of Gamma with respect to radius r.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Real SimpleShockEqn<Jet, Injector>::dGammadr(Real r, Real Gamma, Real t_eng, Real rho) {
    Real dm = medium.mass(r) / (4 * con::pi);  // Mass per unit solid angle from medium
    return -(Gamma * Gamma - 1) / (dM0 + eps_e * dm + 2 * (1 - eps_e) * Gamma * dm) * r * r * rho;
}

#endif