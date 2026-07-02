//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once

#include "../dynamics/shock.h"
#include "inverse-compton.h"
#include "smooth-power-law-syn.h"
#include "syn-concepts.h"
/**
 * <!-- ************************************************************************************** -->
 * @struct SynElectrons
 * @brief Represents synchrotron-emitting electrons in the comoving frame along with their energy distribution
 *        and properties.
 * <!-- ************************************************************************************** -->
 */
struct SynElectrons {
    constexpr static size_t n_breaks{4}; ///< Number of break gamma
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
    [[nodiscard]] Real compute_N_gamma(Real gamma) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the column number density of the electron distribution at a specific Lorentz factor.
     * @details Includes corrections for inverse Compton cooling effects above the cooling Lorentz factor.
     * @param gamma Electron Lorentz factor
     * @return Column number density at the specified Lorentz factor
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real compute_column_den(Real gamma) const noexcept;

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
    [[nodiscard]] inline Real compute_spectrum(Real gamma) const noexcept;
};

static_assert(SynElectronModel<SynElectrons>,
              "SynElectrons must satisfy the SynElectronModel concept (syn-concepts.h)");

/**
 * <!-- ************************************************************************************** -->
 * @defgroup SynchrotronGrids Synchrotron Grid Type Aliases
 * @brief Defines multi-dimensional grid types for Synchrotron Photons and Electrons.
 * <!-- ************************************************************************************** -->
 */

/// The active synchrotron photon model — one-line switch. Any replacement must satisfy the
/// SynPhotonModel concept (syn-concepts.h); conformance is static_assert'ed in each model's
/// header, so the compiler enforces the contract instead of comment discipline.
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
SynElectronGrid generate_syn_electrons(Shock const& shock, Coord const& coord);

/**
 * <!-- ************************************************************************************** -->
 * @brief Populates an existing electron grid with values based on shock parameters
 * @details Modifies a grid supplied by the caller rather than creating a new one.
 * @param electrons The electron grid to populate
 * @param shock The shock object containing physical properties
 * @param coord The coordinate object containing symmetry information
 * <!-- ************************************************************************************** -->
 */
void generate_syn_electrons(SynElectronGrid& electrons, Shock const& shock, Coord const& coord);

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates and returns a new photon grid based on shock and electron grid
 * @details Computes characteristic frequencies and updates calculation constants for each grid cell.
 *          Returns the populated photon grid.
 * @param shock The shock object containing physical properties
 * @param electrons The electron grid providing energy distribution information
 * @param coord The coordinate object containing symmetry information
 * @return A new grid of synchrotron photons
 * <!-- ************************************************************************************** -->
 */
SynPhotonGrid generate_syn_photons(Shock const& shock, SynElectronGrid const& electrons, Coord const& coord);

/**
 * <!-- ************************************************************************************** -->
 * @brief Populates an existing photon grid with values based on shock and electron grid
 * @param photons The photon grid to populate
 * @param shock The shock object containing physical properties
 * @param electrons The electron grid providing energy distribution information
 * @param coord The coordinate object containing symmetry information
 * <!-- ************************************************************************************** -->
 */
void generate_syn_photons(SynPhotonGrid& photons, Shock const& shock, SynElectronGrid const& electrons,
                          Coord const& coord);

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
 * @brief Passive cooling of a Lorentz factor frozen at reverse-shock crossing.
 * @param gamma_x Lorentz factor at crossing
 * @param gamma_m_x Minimum Lorentz factor at crossing
 * @param gamma_m Current minimum Lorentz factor
 * @param dt_comv Comoving time elapsed since crossing
 * @param B Magnetic field
 * @param Y Inverse Compton Y parameter
 * @return The cooled Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real cool_after_crossing(Real gamma_x, Real gamma_m_x, Real gamma_m, Real dt_comv, Real B, Real Y);

/**
 * <!-- ************************************************************************************** -->
 * @brief Applies relic cooling to one grid cell if it lies beyond the electron-injection cutoff.
 * @details Single owner of the relic-cooling rule: beyond Shock::injection_idx no new electrons
 *          are shocked (the reverse shock has finished crossing), so gamma_c and gamma_M evolve
 *          by passive cooling of the population frozen at the crossing cell (injection_idx - 1).
 *          Forward shocks always inject (injection_idx == t_size), so this is a no-op for them.
 * @param electrons Electron grid; cell (i, j, k) is updated in place
 * @param shock Shock providing injection_idx, t_comv, and B
 * @return true if the cell is a relic cell (gamma_c/gamma_M were overwritten)
 * <!-- ************************************************************************************** -->
 */
template <SynElectronModel Electrons>
bool cool_relic_electrons(ElectronGrid<Electrons>& electrons, Shock const& shock, size_t i, size_t j, size_t k) {
    if (!shock.is_relic(i, j, k)) {
        return false;
    }
    const size_t k_inj = shock.injection_idx(i, j);
    auto& elec = electrons(i, j, k);
    auto const& inj = electrons(i, j, k_inj - 1);
    const Real t_com = shock.t_comv(i, j, k);
    const Real B = shock.B(i, j, k);
    const Real dt_comv = t_com - shock.t_comv(i, j, k_inj - 1);
    elec.gamma_c = cool_after_crossing(inj.gamma_c, inj.gamma_m, elec.gamma_m, t_com, B, 0);
    elec.gamma_M = cool_after_crossing(inj.gamma_M, inj.gamma_m, elec.gamma_m, dt_comv, B, 0);
    return true;
}

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
