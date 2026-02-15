//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once

#include "../core/mesh.h"
#include "../core/physics.h"
#include "../dynamics/shock.h"
#include "inverse-compton.h"
#include "power-law-syn.h"
#include "smooth-power-law-syn.h"
/**
 * <!-- ************************************************************************************** -->
 * @struct SynElectrons
 * @brief Represents synchrotron-emitting electrons in the comoving frame along with their energy distribution
 *        and properties.
 * <!-- ************************************************************************************** -->
 */
struct SynElectrons {
    // All values in comoving frame
    Real gamma_m{0};    ///< Minimum electron Lorentz factor
    Real gamma_c{0};    ///< Cooling electron Lorentz factor
    Real gamma_a{0};    ///< Self-absorption Lorentz factor
    Real gamma_M{0};    ///< Maximum electron Lorentz factor
    Real p{2.3};        ///< Power-law index for the electron energy distribution
    Real N_e{0};        ///< shock electron number PER SOLID ANGLE
    Real column_den{0}; ///< Column number density
    Real Y_c{0};        ///< Inverse Compton Y parameter at cooling frequency
    size_t regime{0};   ///< Regime indicator (1-6, determines spectral shape)
    InverseComptonY Ys; ///< InverseComptonY parameters for this electron population

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the comoving electron number (PER SOLID ANGLE) spectrum at a specific Lorentz factor.
     * @details Includes corrections for inverse Compton cooling effects above the cooling Lorentz factor.
     * @param gamma Electron Lorentz factor
     * @return Electron number per solid angle at the specified Lorentz factor
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real compute_N_gamma(Real gamma) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the column number density of the electron distribution at a specific Lorentz factor.
     * @details Includes corrections for inverse Compton cooling effects above the cooling Lorentz factor.
     * @param gamma Electron Lorentz factor
     * @return Column number density at the specified Lorentz factor
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real compute_column_den(Real gamma) const;

  private:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the comoving electron energy spectrum at a given Lorentz factor.
     * @details Different spectral forms apply based on the current regime and relative to
     *          characteristic Lorentz factors (gamma_a, gamma_c, gamma_m, gamma_M).
     * @param gamma Electron Lorentz factor
     * @return The normalized electron energy spectrum value
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] inline Real compute_spectrum(Real gamma) const;
};

/**
 * <!-- ************************************************************************************** -->
 * @defgroup SynchrotronGrids Synchrotron Grid Type Aliases
 * @brief Defines multi-dimensional grid types for Synchrotron Photons and Electrons.
 * <!-- ************************************************************************************** -->
 */

using SynPhotons = SmoothPowerLawSyn;

/// Type alias for 3D grid of synchrotron photons
using SynPhotonGrid = xt::xtensor<SynPhotons, 3>;
/// Type alias for 3D grid of synchrotron electrons
using SynElectronGrid = xt::xtensor<SynElectrons, 3>;

/**
 * <!-- ************************************************************************************** -->
 * @defgroup SynchrotronFunctions Synchrotron Grid Creation and Generation
 * @brief Functions to create and generate grids for Synchrotron electrons and photons.
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates and returns a new electron grid based on shock parameters
 * @details Initializes all electron properties including Lorentz factors, column densities,
 *          and peak intensities for each grid cell.
 * @param shock The shock object containing physical properties
 * @return A new grid of synchrotron electrons
 * <!-- ************************************************************************************** -->
 */
SynElectronGrid generate_syn_electrons(Shock const& shock);

/**
 * <!-- ************************************************************************************** -->
 * @brief Populates an existing electron grid with values based on shock parameters
 * @details Modifies a grid supplied by the caller rather than creating a new one.
 * @param electrons The electron grid to populate
 * @param shock The shock object containing physical properties
 * <!-- ************************************************************************************** -->
 */
void generate_syn_electrons(SynElectronGrid& electrons, Shock const& shock);

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates and returns a new photon grid based on shock and electron grid
 * @details Computes characteristic frequencies and updates calculation constants for each grid cell.
 *          Returns the populated photon grid.
 * @param shock The shock object containing physical properties
 * @param electrons The electron grid providing energy distribution information
 * @return A new grid of synchrotron photons
 * <!-- ************************************************************************************** -->
 */
SynPhotonGrid generate_syn_photons(Shock const& shock, SynElectronGrid const& electrons);

/**
 * <!-- ************************************************************************************** -->
 * @brief Populates an existing photon grid with values based on shock and electron grid
 * @param photons The photon grid to populate
 * @param shock The shock object containing physical properties
 * @param electrons The electron grid providing energy distribution information
 * <!-- ************************************************************************************** -->
 */
void generate_syn_photons(SynPhotonGrid& photons, Shock const& shock, SynElectronGrid const& electrons);

/**
 * <!-- ************************************************************************************** -->
 * @defgroup SynchrotronUpdates Synchrotron Update and Parameter Calculation
 * @brief Functions for updating electron grids and calculating synchrotron parameters.
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Calculates a cooling Lorentz factor based on comoving time, magnetic field, and IC parameters
 * @details Accounts for synchrotron and inverse Compton cooling using an iterative approach
 *          to handle the Lorentz factor-dependent IC cooling.
 * @param t_comv Comoving time
 * @param B Magnetic field
 * @param Y Inverse Compton Y parameter
 * @return The cooling Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_gamma_c(Real t_comv, Real B, Real Y);

/**
 * <!-- ************************************************************************************** -->
 * @brief Determines the electron Lorentz factor at which the number density peaks.
 * @details Based on the relative ordering of absorption, minimum, and cooling Lorentz factors.
 * @param gamma_a Absorption Lorentz factor
 * @param gamma_m Minimum electron Lorentz factor
 * @param gamma_c Cooling electron Lorentz factor
 * @return Peak Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_gamma_peak(Real gamma_a, Real gamma_m, Real gamma_c);

/**
 * <!-- ************************************************************************************** -->
 * @brief Calculates synchrotron frequency for a given Lorentz factor and magnetic field
 * @param gamma Electron Lorentz factor
 * @param B Magnetic field
 * @return The synchrotron frequency
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_freq(Real gamma, Real B);

/**
 * <!-- ************************************************************************************** -->
 * @brief Calculates electron Lorentz factor corresponding to a given synchrotron frequency
 * @param nu Synchrotron frequency
 * @param B Magnetic field
 * @return The corresponding electron Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_gamma(Real nu, Real B);
