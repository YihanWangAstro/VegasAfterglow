//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _SHOCKDYNAMICS_
#define _SHOCKDYNAMICS_

#include <tuple>

#include "boost/numeric/odeint.hpp"
#include "jet.h"
#include "medium.h"
#include "mesh.h"
#include "physics.h"

/********************************************************************************************************************
 * CLASS: Shock
 * DESCRIPTION: Represents a shock structure that stores grid-based data for the shock evolution, including
 *              comoving time, radius, theta(for jet spreading), bulk Lorentz factor, relative Lorentz factor,
 *              magnetic field, and proton column density.
 *              It also stores constant energy fractions (eps_e and eps_B) and provides a method to return
 *              the grid dimensions.
 ********************************************************************************************************************/
class Shock {
   public:
    Shock(size_t phi_size, size_t theta_size, size_t t_size, Real eps_e, Real eps_B);
    Shock() = delete;

    MeshGrid3d t_com;           // comoving time
    MeshGrid3d r;               // radius
    MeshGrid3d theta;           // theta for jet spreading
    MeshGrid3d Gamma;           // bulk lorentz factor
    MeshGrid3d Gamma_rel;       // relative lorentz factor between down stream and up stream
    MeshGrid3d B;               // comoving magnetic field
    MeshGrid3d column_num_den;  // down stream proton column number density
    Real eps_e{0};              // electron energy fraction
    Real eps_B{0};              // magnetic energy fraction

    auto shape() const { return std::make_tuple(phi_size, theta_size, t_size); }  // Returns grid dimensions

   private:
    size_t const phi_size{0};    // Number of grid points in phi direction
    size_t const theta_size{0};  // Number of grid points in theta direction
    size_t const t_size{0};      // Number of grid points in time direction
};

/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Shock Generation Interfaces
 * DESCRIPTION: These function templates declare interfaces to generate forward shocks (2D and 3D) and
 *              forward/reverse shock pairs.
 ********************************************************************************************************************/
using ShockPair = std::pair<Shock, Shock>;

template <typename Jet, typename Injector>
Shock genForwardShock(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, Real eps_e,
                      Real eps_B, Real rtol = 1e-6, bool is_axisymmetric = true);

template <typename Jet, typename Injector>
ShockPair genFRShocks(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, Real eps_e,
                      Real eps_B, Real rtol = 1e-6, bool is_axisymmetric = true);

/********************************************************************************************************************
 * INLINE FUNCTIONS: Shock Utilities
 * DESCRIPTION: This section defines a set of inline functions used in shock calculations. These functions compute
 *              various physical quantities such as the comoving magnetic field (via the Weibel instability),
 *              thermal energy density, time derivatives, jet width derivative, downstream number density, fluid
 *              velocities, and update the shock state.
 ********************************************************************************************************************/

Real u_DownStr(Real gamma_rel, Real sigma);
Real u_UpStr2u_DownStr(Real gamma_rel, Real sigma);
void updateShockState(Shock& shock, size_t i, size_t j, size_t k, Real r, Real theta, Real Gamma, Real Gamma_rel,
                      Real t_com, Real dNdOmega_up, Real n_up_str, Real sigma);

inline Real soundSpeed(Real Gamma_rel) {
    Real ad_idx = adiabaticIndex(Gamma_rel);
    return std::sqrt(ad_idx * (ad_idx - 1) * (Gamma_rel - 1) / (1 + (Gamma_rel - 1) * ad_idx));
}
inline Real coMovingWeibelB(Real eps_B, Real e_thermal) { return std::sqrt(8 * con::pi * eps_B * e_thermal); }
inline Real drdt(Real beta) { return (beta * con::c) / (1 - beta); }
inline Real dtheta_dt(Real uv, Real drdt, Real r, Real Gamma) {
    return 0.5 / Gamma * drdt / r * sqrt((2 * uv * uv + 3) / (4 * uv * uv + 3));
}
inline Real dtdt_CoMoving(Real Gamma) { return 1 / (Gamma - std::sqrt(Gamma * Gamma - 1)); };
inline Real calc_pB4(Real n4, Real sigma) { return sigma * n4 * con::mp * con::c2 / 2; }
inline Real u_UpStr(Real u_down, Real gamma_rel) {
    return std::sqrt((1 + u_down * u_down) * (gamma_rel * gamma_rel - 1)) + u_down * gamma_rel;
}
inline Real n_DownStr(Real n_up_str, Real gamma_rel, Real sigma) {
    return n_up_str * u_UpStr2u_DownStr(gamma_rel, sigma);
}
inline Real relativeLorentz(Real gamma1, Real gamma2) {
    return gamma1 * gamma2 - std::sqrt((gamma1 * gamma1 - 1) * (gamma2 * gamma2 - 1));
}
inline Real e_ThermalDownStr(Real gamma_rel, Real n_down_str) {
    return n_down_str * (gamma_rel - 1) * con::mp * con::c2;
}
// D_jet co-moving shell width
inline Real dDdt_Jet(Real Gamma_rel, Real dtdt_com) {
    Real cs = soundSpeed(Gamma_rel);
    return cs * dtdt_com;
}

inline Real dN3dt(Real r, Real n1, Real n4, Real gamma3, Real drdt, Real gamma0, Real sigma) {
    Real gamma34 = relativeLorentz(gamma0, gamma3);
    Real ratio_u = u_UpStr2u_DownStr(gamma34, sigma);
    Real n3 = n4 * ratio_u;
    Real dxdr = 1. / (gamma0 * std::sqrt((1 + sigma) * n4 / n1) * (1 - gamma0 * n4 / (gamma3 * n3)));
    return n3 * r * r * gamma3 * dxdr * drdt;
}

inline Real calc_n4(Real dEdOmega, Real Gamma0, Real r, Real D_jet, Real sigma) {
    return dEdOmega / (Gamma0 * con::mp * con::c2 * r * r * D_jet) / (1 + sigma);
}

inline void setStoppingShock(size_t i, size_t j, Shock& shock, Array const& t, Real r0, Real theta0) {
    shock.t_com[i][j] = t;
    std::fill(shock.r[i][j].begin(), shock.r[i][j].end(), r0);
    std::fill(shock.theta[i][j].begin(), shock.theta[i][j].end(), theta0);
}

template <typename Jet, typename Injector>
class SimpleShockEqn;
template <typename Jet, typename Injector>
class ForwardShockEqn;
template <typename Jet, typename Injector>
class FRShockEqn;

/********************************************************************************************************************
 * FUNCTION: genForwardShock
 * DESCRIPTION: Generates a forward shock using the provided coordinates, medium, jet, injector, and energy
 *              fractions. It creates a Shock object and iterates over phi, theta values, solving the shock
 *              evolution for each theta slice.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Shock genForwardShock(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, Real eps_e,
                      Real eps_B, Real rtol, bool is_axisymmetric) {
    auto [phi_size, theta_size, t_size] = coord.shape();  // Unpack coordinate dimensions
    size_t phi_size_needed = is_axisymmetric ? 1 : phi_size;
    Shock f_shock(phi_size_needed, theta_size, t_size, eps_e, eps_B);  // Create Shock with 1 phi slice
    for (size_t i = 0; i < phi_size_needed; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            // Create a ForwardShockEqn for each theta slice (phi is set to 0)
            // auto eqn = ForwardShockEqn(medium, jet, inject, 0, coord.theta[j], eps_e);
            auto eqn = SimpleShockEqn(medium, jet, inject, 0, coord.theta[j], eps_e);
            //  Solve the shock shell for this theta slice
            solveForwardShell(i, j, coord.t, f_shock, eqn, rtol);
        }
    }

    return f_shock;
}

/********************************************************************************************************************
 * FUNCTION: genFRShocks
 * DESCRIPTION: Generates a pair of forward and reverse shocks using the provided coordinates, medium, jet,
 *              injector, and energy fractions. It creates two Shock objects (one forward, one reverse) and solves
 *              the shock shells for each phi, theta slice.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
ShockPair genFRShocks(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, Real eps_e,
                      Real eps_B, Real rtol, bool is_axisymmetric) {
    auto [phi_size, theta_size, t_size] = coord.shape();
    size_t phi_size_needed = is_axisymmetric ? 1 : phi_size;
    Shock f_shock(phi_size_needed, theta_size, t_size, eps_e, eps_B);  // Forward shock for 1 phi slice
    Shock r_shock(phi_size_needed, theta_size, t_size, eps_e, eps_B);  // Reverse shock for 1 phi slice
    for (size_t i = 0; i < phi_size_needed; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            // Create equations for forward and reverse shocks for each theta slice (phi is 0)
            // auto eqn_f = ForwardShockEqn(medium, jet, inject, 0, coord.theta[j], eps_e);
            auto eqn_f = SimpleShockEqn(medium, jet, inject, 0, coord.theta[j], eps_e);
            auto eqn_r = FRShockEqn(medium, jet, inject, 0, coord.theta[j], eps_e);
            // Solve the forward-reverse shock shell
            solveFRShell(i, j, coord.t, f_shock, r_shock, eqn_f, eqn_r, rtol);
        }
    }

    return std::make_pair(std::move(f_shock), std::move(r_shock));
}

#include "forward-shock.hpp"
#include "reverse-shock.hpp"
#include "simple-shock.hpp"
#endif