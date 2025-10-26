//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <cmath>

#include "../util/macros.h"
/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the deceleration radius of the shock.
 * @details For given isotropic energy E_iso, ISM density n_ism, initial Lorentz factor Gamma0,
 *          and engine duration, the deceleration radius is the maximum of the thin shell and thick shell
 *          deceleration radii.
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @param engine_dura Engine duration
 * @return The deceleration radius
 * <!-- ************************************************************************************** -->
 */
Real dec_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the deceleration radius for the thin shell case.
 * @details Uses the formula: R_dec = [3E_iso / (4π n_ism mp c^2 Gamma0^2)]^(1/3)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @return The thin shell deceleration radius
 * <!-- ************************************************************************************** -->
 */
Real thin_shell_dec_radius(Real E_iso, Real n_ism, Real Gamma0);

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the deceleration radius for the thick shell case.
 * @details Uses the formula: R_dec = [3 E_iso engine_dura c / (4π n_ism mp c^2)]^(1/4)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param engine_dura Engine duration
 * @return The thick shell deceleration radius
 * <!-- ************************************************************************************** -->
 */
Real thick_shell_dec_radius(Real E_iso, Real n_ism, Real engine_dura);

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the radius at which shell spreading becomes significant.
 * @details Uses the formula: R_spread = Gamma0^2 * c * engine_dura
 * @param Gamma0 Initial Lorentz factor
 * @param engine_dura Engine duration
 * @return The shell spreading radius
 * <!-- ************************************************************************************** -->
 */
Real shell_spreading_radius(Real Gamma0, Real engine_dura);

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the radius at which the reverse shock transitions.
 * @details Based on the Sedov length, engine duration, and initial Lorentz factor.
 *          Uses the formula: R_RS = (SedovLength^(1.5)) / (sqrt(c * engine_dura) * Gamma0^2)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @param engine_dura Engine duration
 * @return The reverse shock transition radius
 * <!-- ************************************************************************************** -->
 */
Real RS_transition_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the dimensionless parameter (ξ) that characterizes the shell geometry.
 * @details This parameter helps determine whether the shell behaves as thick or thin.
 *          Uses the formula: ξ = sqrt(Sedov_length / shell_width) * Gamma0^(-4/3)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @param engine_dura Engine duration
 * @return The shell thickness parameter ξ
 * <!-- ************************************************************************************** -->
 */
Real shell_thickness_param(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);

/**
 * <!-- ************************************************************************************** -->
 * @brief Calculates the engine duration needed to achieve a specific shell thickness parameter.
 * @details Uses the formula: T_engine = Sedov_l / (ξ^2 * Gamma0^(8/3) * c)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @param xi Target shell thickness parameter
 * @return The required engine duration
 * <!-- ************************************************************************************** -->
 */
Real calc_engine_duration(Real E_iso, Real n_ism, Real Gamma0, Real xi);

/**
 * <!-- ************************************************************************************** -->
 * @defgroup GammaConversions Gamma Conversion and Adiabatic Index Functions
 * @brief Helper functions for Lorentz factor conversions and adiabatic index calculations
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts Lorentz factor (gamma) to a velocity fraction (beta)
 * @param gamma Lorentz factor
 * @return Velocity fraction (beta = v/c)
 * <!-- ************************************************************************************** -->
 */
inline Real gamma_to_beta(Real gamma) {
    return std::sqrt(gamma * gamma - 1) / gamma;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes adiabatic index as a function of the Lorentz factor
 * @param gamma Lorentz factor
 * @return Adiabatic index
 * <!-- ************************************************************************************** -->
 */
inline Real adiabatic_idx(Real gamma) {
    return 4.0 / 3.0 + 1 / (3 * gamma);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the Sedov length—a characteristic scale for blast wave deceleration
 * @param E_iso Isotropic equivalent energy
 * @param n_ism ISM number density
 * @return Sedov length
 * @details The Sedov length is a characteristic scale defined as the cube root of (E_iso / (ρc²)),
 *          where ρ is the ambient medium mass density
 * <!-- ************************************************************************************** -->
 */
inline Real sedov_length(Real E_iso, Real n_ism) {
    return std::cbrt(E_iso / (4 * con::pi / 3 * n_ism * con::mp * con::c2));
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Returns the radius at which the reverse shock crosses, defined as the thick shell deceleration radius
 * @param E_iso Isotropic equivalent energy
 * @param n_ism ISM number density
 * @param engine_dura Engine duration
 * @return Reverse shock crossing radius
 * <!-- ************************************************************************************** -->
 */
inline Real RS_crossing_radius(Real E_iso, Real n_ism, Real engine_dura) {
    const Real l = sedov_length(E_iso, n_ism);
    return std::sqrt(std::sqrt(l * l * l * con::c * engine_dura));
}

//========================================================================================================
//                                  template function implementation
//========================================================================================================

/**
 * <!-- ************************************************************************************** -->
 * @brief Parameters for radiation transport
 * @details Parameters for radiation transport
 * <!-- ************************************************************************************** -->
 */
struct RadParams {
    Real eps_e{0.1};  ///< Electron energy fraction
    Real eps_B{0.01}; ///< Magnetic field energy fraction
    Real p{2.3};      ///< Electron energy distribution index
    Real xi_e{1};     ///< Electron self-absorption parameter
};
