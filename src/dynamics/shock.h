//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <tuple>

#include "../core/mesh.h"
#include "../core/physics.h"
#include "../environment/medium.h"
#include "boost/numeric/odeint.hpp"

/**
 * <!-- ************************************************************************************** -->
 * @class Shock
 * @brief Represents a shock wave in an astrophysical environment.
 * @details The class stores physical properties of the shock across a 3D grid defined by azimuthal angle (phi),
 *          polar angle (theta), and time bins. Provides methods for shock calculations, including relativistic
 *          jump conditions, magnetic field calculations, and energy density computations.
 * <!-- ************************************************************************************** -->
 */
class Shock {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Constructs a Shock object with the given grid dimensions and energy fractions.
     * @details Initializes various 3D grids for storing physical properties of the shock, including comoving time,
     *          radius, Lorentz factors, magnetic fields, and downstream densities.
     * @param phi_size Number of grid points in the phi direction
     * @param theta_size Number of grid points in theta direction
     * @param t_size Number of grid points in time direction
     * @param rad_params Radiation parameters
     * <!-- ************************************************************************************** -->
     */
    Shock(size_t phi_size, size_t theta_size, size_t t_size, RadParams const& rad_params);

    Shock() noexcept = default;

    MeshGrid3d t_comv;       ///< Comoving time
    MeshGrid3d r;            ///< Radius
    MeshGrid3d theta;        ///< Theta for jet spreading
    MeshGrid3d Gamma;        ///< Bulk Lorentz factor
    MeshGrid3d Gamma_th;     ///< Downstream internal Lorentz factor
    MeshGrid3d B;            ///< Comoving magnetic field
    MeshGrid3d N_p;          ///< Downstream proton number per solid angle
    IndexGrid injection_idx; ///< Beyond which grid index there is no electron injection
    MaskGrid required;       ///< Grid points actually required for final flux calculation
    RadParams rad;           ///< Radiation parameters

    /// Returns grid dimensions as a tuple
    [[nodiscard]] auto shape() const { return std::make_tuple(phi_size, theta_size, t_size); }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Resizes all grid components of the Shock object to new dimensions.
     * @param phi_size New number of grid points in the phi direction
     * @param theta_size New number of grid points in theta direction
     * @param t_size New number of grid points in time direction
     * <!-- ************************************************************************************** -->
     */
    void resize(size_t phi_size, size_t theta_size, size_t t_size);

  private:
    size_t phi_size{0};   ///< Number of grid points in phi direction
    size_t theta_size{0}; ///< Number of grid points in theta direction
    size_t t_size{0};     ///< Number of grid points in time direction
};

/**
 * <!-- ************************************************************************************** -->
 * @defgroup ShockUtilities Shock Utilities
 * @brief Inline functions used in shock calculations.
 * @details This section defines a set of inline functions used in shock calculations. These functions compute
 *          various physical quantities such as the comoving magnetic field (via the Weibel instability),
 *          thermal energy density, time derivatives, jet width derivative, downstream number density, fluid
 *          velocities, and update the shock state.
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the downstream four-velocity for a given relative Lorentz factor and magnetization parameter.
 * @param gamma_rel Relative Lorentz factor between upstream and downstream regions
 * @param sigma Magnetization parameter (ratio of magnetic to rest-mass energy density)
 * @return The downstream four-velocity in the shock frame
 * @details The calculation handles both magnetized (sigma > 0) and non-magnetized (sigma = 0) cases using
 *          different equations based on jump conditions across the shock front.
 * <!-- ************************************************************************************** -->
 */
Real compute_downstr_4vel(Real gamma_rel, Real sigma);

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the upstream four-velocity from downstream four-velocity and relative Lorentz factor.
 * @param u_down Downstream four-velocity
 * @param gamma_rel Relative Lorentz factor
 * @return The upstream four-velocity in the shock frame
 * <!-- ************************************************************************************** -->
 */
inline Real compute_upstr_4vel(Real u_down, Real gamma_rel) {
    return std::sqrt((1 + u_down * u_down) * std::fabs(gamma_rel * gamma_rel - 1)) + u_down * gamma_rel;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the ratio of upstream to downstream four-velocity across the shock front.
 * @param gamma_rel Relative Lorentz factor between upstream and downstream regions
 * @param sigma_upstr Magnetization parameter (ratio of magnetic to rest-mass energy density)
 * @return The ratio of upstream to downstream four-velocity
 * @details This ratio is a key parameter in determining various shock properties, such as compression ratio
 *          and jump conditions for density, pressure, and magnetic field.
 * <!-- ************************************************************************************** -->
 */
inline Real compute_4vel_jump(Real gamma_rel, Real sigma_upstr) {
    const Real u_down_s_ = compute_downstr_4vel(gamma_rel, sigma_upstr);
    const Real u_up_s_ = compute_upstr_4vel(u_down_s_, gamma_rel);
    Real ratio_u = u_up_s_ / u_down_s_;
    if (u_down_s_ == 0.) {
        ratio_u = 4 * gamma_rel; // (g_hat*gamma_rel+1)/(g_hat-1)
    }
    return ratio_u;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the sound speed in the shocked medium based on the relative Lorentz factor.
 * @param Gamma_rel Relative Lorentz factor between upstream and downstream regions
 * @return The sound speed as a fraction of light speed
 * <!-- ************************************************************************************** -->
 */
inline Real compute_sound_speed(Real Gamma_rel) {
    const Real ad_idx = adiabatic_idx(Gamma_rel);
    return std::sqrt(std::fabs(ad_idx * (ad_idx - 1) * (Gamma_rel - 1) / (1 + (Gamma_rel - 1) * ad_idx))) * con::c;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the effective Lorentz factor accounting for the adiabatic index.
 * @param adx Adiabatic index of the medium
 * @param Gamma Bulk Lorentz factor
 * @return The effective Lorentz factor
 * <!-- ************************************************************************************** -->
 */
inline Real compute_effective_Gamma(Real adx, Real Gamma) {
    return (adx * Gamma * Gamma - adx + 1) / Gamma;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the derivative of the effective Lorentz factor with respect to the bulk Lorentz factor.
 * @param adx Adiabatic index of the medium
 * @param Gamma Bulk Lorentz factor
 * @return The derivative of the effective Lorentz factor
 * <!-- ************************************************************************************** -->
 */
inline Real compute_effective_Gamma_dGamma(Real adx, Real Gamma) {
    const Real Gamma2 = Gamma * Gamma;
    return (adx * Gamma2 + adx - 1) / Gamma2;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the comoving magnetic field using the Weibel instability mechanism.
 * @param eps_B Fraction of thermal energy in magnetic fields
 * @param e_thermal Thermal energy density
 * @return The comoving magnetic field strength
 * <!-- ************************************************************************************** -->
 */
inline Real compute_comv_weibel_B(Real eps_B, Real e_thermal) {
    return std::sqrt(8 * con::pi * eps_B * e_thermal);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the time derivative of radius (dr/dt) based on the shock velocity.
 * @param beta Shock velocity as a fraction of light speed
 * @return The rate of change of radius with respect to observer time
 * <!-- ************************************************************************************** -->
 */
inline Real compute_dr_dt(Real beta) {
    return (beta * con::c) / (1 - beta);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the time derivative of theta (dθ/dt) for jet spreading.
 * @param theta_s Typical spreading angle of the jet
 * @param theta Theta of the current grid point
 * @param drdt Time derivative of radius
 * @param r Current radius
 * @param Gamma Current bulk Lorentz factor
 * @return The rate of change of the half-opening angle
 * <!-- ************************************************************************************** -->
 */
inline Real compute_dtheta_dt(Real theta_s, Real theta, Real drdt, Real r, Real Gamma) {
    constexpr Real Q = 7;
    const Real u2 = Gamma * Gamma - 1;
    const Real u = std::sqrt(u2);
    const Real f = 1 / (1 + u * theta_s * Q);
    return drdt / (2 * Gamma * r) * std::sqrt((2 * u2 + 3) / (4 * u2 + 3)) * f;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the time derivative of comoving time (dt_comv/dt) based on the Lorentz factor.
 * @param Gamma Bulk Lorentz factor
 * @param beta Bulk velocity as a fraction of light speed
 * @return The rate of change of comoving time with respect to observer time
 * <!-- ************************************************************************************** -->
 */
inline Real compute_dt_dt_comv(Real Gamma, Real beta) {
    return 1 / (Gamma * (1 - beta));
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the upstream magnetic field.
 * @param rho_up Upstream density
 * @param sigma Magnetization parameter
 * @return The upstream magnetic field
 * <!-- ************************************************************************************** -->
 */
inline Real compute_upstr_B(Real rho_up, Real sigma) {
    return std::sqrt((4 * con::pi * con::c2) * sigma * rho_up);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the relative Lorentz factor between two frames with given Lorentz factors.
 * @param gamma1 First frame's Lorentz factor
 * @param gamma2 Second frame's Lorentz factor
 * @return The relative Lorentz factor between the two frames
 * <!-- ************************************************************************************** -->
 */
inline Real compute_rel_Gamma(Real gamma1, Real gamma2) {
    return gamma1 * gamma2 - std::sqrt(std::fabs((gamma1 * gamma1 - 1) * (gamma2 * gamma2 - 1)));
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the relative Lorentz factor between two frames with given Lorentz factors and velocities.
 * @param gamma1 First frame's Lorentz factor
 * @param gamma2 Second frame's Lorentz factor
 * @param beta1 First frame's velocity as a fraction of light speed
 * @param beta2 Second frame's velocity as a fraction of light speed
 * @return The relative Lorentz factor between the two frames
 * <!-- ************************************************************************************** -->
 */
inline Real compute_rel_Gamma(Real gamma1, Real gamma2, Real beta1, Real beta2) {
    return gamma1 * gamma2 * (1 - beta1 * beta2);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes a Lorentz factor from a reference Lorentz factor and relative Lorentz factor.
 * @param gamma4 Reference Lorentz factor (typically for region 4, unshocked ejecta)
 * @param gamma_rel Relative Lorentz factor
 * @return The derived Lorentz factor
 * <!-- ************************************************************************************** -->
 */
inline Real compute_Gamma_from_relative(Real gamma4, Real gamma_rel) {
    const Real b = -2 * gamma4 * gamma_rel;
    const Real c = gamma4 * gamma4 + gamma_rel * gamma_rel - 1;
    return (-b - std::sqrt(std::fabs(b * b - 4 * c))) / 2;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the shock heating rate.
 * @param Gamma_rel Relative Lorentz factor
 * @param mdot Mass accretion rate
 * @return The shock heating rate
 * <!-- ************************************************************************************** -->
 */
inline Real compute_shock_heating_rate(Real Gamma_rel, Real mdot) {
    return mdot * (Gamma_rel - 1) * con::c2;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the adiabatic cooling rate.
 * @param ad_idx Adiabatic index
 * @param r Radius
 * @param Gamma Lorentz factor
 * @param u Internal energy density
 * @param drdt Rate of change of radius
 * @param dGammadt Rate of change of the Lorentz factor
 * @return The adiabatic cooling rate
 * <!-- ************************************************************************************** -->
 */
inline Real compute_adiabatic_cooling_rate(Real ad_idx, Real r, Real Gamma, Real u, Real drdt, Real dGammadt) {
    return -(ad_idx - 1) * (3 * drdt / r - dGammadt / Gamma) * u;
}

inline Real compute_adiabatic_cooling_rate2(Real ad_idx, Real r, Real x, Real u, Real drdt, Real dxdt) {
    Real dlnvdt = 2 * drdt / r;
    if (x > 0) {
        dlnvdt += dxdt / x;
    }
    return -(ad_idx - 1) * dlnvdt * u;
}
/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the rate at which the shock shell spreads in the comoving frame.
 * @param Gamma_rel Relative Lorentz factor
 * @param dtdt_comv Rate of change of comoving time with respect to burst time
 * @return The shell spreading rate in the comoving frame
 * <!-- ************************************************************************************** -->
 */
inline Real compute_shell_spreading_rate(Real Gamma_rel, Real dtdt_comv) {
    const Real cs = compute_sound_speed(Gamma_rel);
    return cs * dtdt_comv;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the radiative efficiency based on the radiative constant, comoving time, Lorentz factor, and density.
 * @param t_comv Comoving time
 * @param Gamma_th Thermal Lorentz factor
 * @param u    internal energy density
 * @param rad  Radiation parameters
 * @return The radiative efficiency
 * <!-- ************************************************************************************** -->
 */
inline Real compute_radiative_efficiency(Real t_comv, Real Gamma_th, Real u, RadParams const& rad) { //
    const Real gamma_m = (rad.p - 2) / (rad.p - 1) * rad.eps_e * (Gamma_th - 1) * con::mp / con::me / rad.xi_e + 1;
    const Real gamma_c_bar = (6 * con::pi * con::me * con::c / con::sigmaT) / (rad.eps_B * u * t_comv);
    const Real gamma_c = 0.5 * (gamma_c_bar + std::sqrt(gamma_c_bar * gamma_c_bar + 4));

    const Real g_m_g_c = std::fabs(gamma_m / gamma_c); // gamma_m/gamma_c
    if (g_m_g_c < 1 && rad.p > 2) {                    // slow cooling
        return rad.eps_e * fast_pow(g_m_g_c, rad.p - 2);
    } else { // fast cooling or p<=2
        return rad.eps_e;
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the thermal Lorentz factor from the thermal energy and mass.
 * @param U_th Thermal energy
 * @param mass Mass
 * @param limiter
 * @return The thermal Lorentz factor
 * <!-- ************************************************************************************** -->
 */
inline Real compute_Gamma_therm(Real U_th, Real mass, bool limiter = false) {
    if (mass == 0) [[unlikely]] {
        return 1;
    } else [[likely]] {
        const Real Gamma_th = U_th / (mass * con::c2) + 1;
        if (limiter && Gamma_th < con::Gamma_cut) {
            return 1;
        } else {
            return Gamma_th;
        }
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Saves the current state of the shock.
 * @param shock Reference to the Shock object to be updated
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param k Grid index for time
 * @param t_comv Comoving time
 * @param r Radius
 * @param theta Angle
 * @param Gamma Lorentz factor
 * @param Gamma_th Thermal Lorentz factor
 * @param B Magnetic field strength
 * @param mass Mass
 * <!-- ************************************************************************************** -->
 */

inline void write_shock_state(Shock& shock, size_t i, size_t j, size_t k, Real t_comv, Real r, Real theta, Real Gamma,
                              Real Gamma_th, Real B, Real mass) {
    shock.t_comv(i, j, k) = t_comv;
    shock.r(i, j, k) = r;
    shock.theta(i, j, k) = theta;
    shock.Gamma(i, j, k) = Gamma;
    shock.Gamma_th(i, j, k) = Gamma_th;
    shock.B(i, j, k) = B;
    shock.N_p(i, j, k) = mass / con::mp;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the compression ratio across the shock.
 * @param Gamma_upstr Lorentz factor upstream
 * @param Gamma_downstr Lorentz factor downstream
 * @param sigma_upstr Upstream magnetization
 * @return The compression ratio
 * <!-- ************************************************************************************** -->
 */
inline Real compute_compression(Real Gamma_upstr, Real Gamma_downstr, Real sigma_upstr) {
    const Real Gamma_rel = compute_rel_Gamma(Gamma_upstr, Gamma_downstr);
    return compute_4vel_jump(Gamma_rel, sigma_upstr);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the downstream magnetic field strength.
 * @param eps_B Magnetic energy fraction
 * @param rho_upstr Upstream rest mass density
 * @param B_upstr Upstream magnetic field strength
 * @param Gamma_th Thermal Lorentz factor
 * @param comp_ratio Compression ratio
 * @return The downstream magnetic field strength
 * <!-- ************************************************************************************** -->
 */
inline Real compute_downstr_B(Real eps_B, Real rho_upstr, Real B_upstr, Real Gamma_th, Real comp_ratio) {
    const Real rho_downstr = rho_upstr * comp_ratio;

    const Real e_th = (Gamma_th - 1) * rho_downstr * con::c2;

    return compute_comv_weibel_B(eps_B, e_th) + B_upstr * comp_ratio;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Sets a stopping shock state when the Lorentz factor drops below a threshold.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param shock Reference to the Shock object to be updated
 * @param state0 Initial state to be used for some parameters
 * <!-- ************************************************************************************** -->
 */
template <typename State>
inline void set_stopping_shock(size_t i, size_t j, Shock& shock, State const& state0);

//========================================================================================================
//                                  template function implementation
//========================================================================================================
template <typename State>
inline void set_stopping_shock(size_t i, size_t j, Shock& shock, State const& state0) {
    xt::view(shock.t_comv, i, j, xt::all()) = state0.t_comv;
    xt::view(shock.r, i, j, xt::all()) = state0.r;
    xt::view(shock.theta, i, j, xt::all()) = state0.theta;
    xt::view(shock.Gamma, i, j, xt::all()) = 1;
    xt::view(shock.Gamma_th, i, j, xt::all()) = 1;
    xt::view(shock.B, i, j, xt::all()) = 0;
    xt::view(shock.N_p, i, j, xt::all()) = 0;
}

/*
template <typename Eqn>
Real compute_dec_time(Eqn const& eqn, Real t_max) {
    Real e_k = eqn.ejecta.eps_k(eqn.phi, eqn.theta0);
    Real gamma = eqn.ejecta.Gamma0(eqn.phi, eqn.theta0);
    Real beta = gamma_to_beta(gamma);

    Real m_shell = e_k / (gamma * con::c2);

    if constexpr (HasSigma<decltype(eqn.ejecta)>) {
        m_shell /= 1 + eqn.ejecta.sigma0(eqn.phi, eqn.theta0);
    }

    Real r_dec = beta * con::c * t_max / (1 - beta);
    for (size_t i = 0; i < 30; i++) {
        Real m_swept = eqn.medium.mass(eqn.phi, eqn.theta0, r_dec);
        if (m_swept < m_shell / gamma) {
            break;
        }
        r_dec /= 3;
    }
    Real t_dec = r_dec * (1 - beta) / (beta * con::c);
    return t_dec;
}*/
