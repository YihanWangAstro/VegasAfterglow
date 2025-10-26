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
     * @return Column number density at the specified Lorentz factor
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real compute_N_gamma(Real gamma) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the column number density of the electron distribution.
     * @details Uses the electron energy spectrum to compute the column number density.
     * @return Column number density of the electron distribution
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
 * @struct SynPhotons
 * @brief Represents synchrotron photons in the comoving frame and provides spectral functions.
 * <!-- ************************************************************************************** -->
 */
struct SynPhotons {
    // All values in comoving frame
    Real I_nu_max{0}; ///< Maximum specific synchrotron power PER SOLID ANGLE
    Real nu_m{0};     ///< Characteristic frequency corresponding to gamma_m
    Real nu_c{0};     ///< Cooling frequency corresponding to gamma_c
    Real nu_a{0};     ///< Self-absorption frequency
    Real nu_M{0};     ///< Maximum photon frequency
    Real p{2.3};      ///< Power-law index for the electron energy distribution

    Real log2_I_nu_max{0}; ///< Log2 of I_nu_max (for computational efficiency)
    Real log2_nu_m{0};     ///< Log2 of nu_m
    Real log2_nu_c{0};     ///< Log2 of nu_c
    Real log2_nu_a{0};     ///< Log2 of nu_a
    Real log2_nu_M{0};     ///< Log2 of nu_M
    Real Y_c{0};           ///< Inverse Compton Y parameter at cooling frequency
    size_t regime{0};      ///< Regime indicator (1-6, determines spectral shape)
    InverseComptonY Ys;    ///< InverseComptonY parameters for this electron population

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the comoving synchrotron specific intensity.
     * @details Includes inverse Compton corrections for frequencies above the cooling frequency.
     * @param nu Frequency at which to compute the specific intensity
     * @return The synchrotron specific intensity at the specified frequency
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real compute_I_nu(Real nu) const; ///< Linear power PER SOLID ANGLE

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the base-2 logarithm of comoving synchrotron specific intensity at a given frequency.
     * @details Optimized for numerical computation by using logarithmic arithmetic.
     * @param log2_nu Base-2 logarithm of the frequency
     * @return Base-2 logarithm of synchrotron specific intensity
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real
    compute_log2_I_nu(Real log2_nu) const; ///<  Log2 specific intensity (for computational efficiency)

    /**
     * <!-- ************************************************************************************** -->
     * @brief Updates cached calculation constants used for efficiently computing synchrotron spectra.
     * @details Constants vary based on the electron regime (1-6) and involve different power laws.
     * <!-- ************************************************************************************** -->
     */
    void update_constant();

  private:
    // Cached calculation constants for spectral computations
    // Optimized calculation constants
    Real C1_{0}; ///< Cached spectral coefficient 1
    Real C2_{0}; ///< Cached spectral coefficient 2
    Real C3_{0}; ///< Cached spectral coefficient 3

    // Log2 of calculation constants for faster computation
    Real log2_C1_{0};
    Real log2_C2_{0};
    Real log2_C3_{0};
    Real log2_C4_{0};

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the synchrotron spectrum at a given frequency based on the electron regime.
     * @details Implements the broken power-law with exponential cutoff formulae for different regimes.
     * @param nu The frequency at which to compute the spectrum
     * @return The normalized synchrotron spectrum value
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] inline Real compute_spectrum(Real nu) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the base-2 logarithm of synchrotron spectrum at a given frequency.
     * @details Uses logarithmic arithmetic for numerical stability in different spectral regimes.
     * @param log2_nu Base-2 logarithm of the frequency
     * @return Base-2 logarithm of the synchrotron spectrum
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] inline Real compute_log2_spectrum(Real log2_nu) const;
};

/**
 * <!-- ************************************************************************************** -->
 * @defgroup SynchrotronGrids Synchrotron Grid Type Aliases
 * @brief Defines multi-dimensional grid types for Synchrotron Photons and Electrons.
 * <!-- ************************************************************************************** -->
 */

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
