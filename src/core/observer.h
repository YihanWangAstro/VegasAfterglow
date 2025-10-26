//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once

#include <iostream>
#include <stdexcept>
#include <thread>

#include "../dynamics/shock.h"
#include "../util/macros.h"
#include "../util/utilities.h"
#include "mesh.h"
#include "xtensor/core/xnoalias.hpp"

/**
 * <!-- ************************************************************************************** -->
 * @class Observer
 * @brief Represents an observer in the GRB afterglow simulation.
 * @details This class handles the calculation of observed quantities such as specific flux, integrated flux,
 *          and spectra. It accounts for relativistic effects (Doppler boosting), cosmological effects (redshift),
 *          and geometric effects (solid angle). The observer can be placed at any viewing angle relative to the
 *          jet axis.
 * <!-- ************************************************************************************** -->
 */
class Observer {
  public:
    /// Default constructor
    Observer() = default;

    /// Grid of observation times
    MeshGrid3d time;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Sets up the Observer for flux calculation.
     * @details Initializes the observation time and Doppler factor grids, as well as the emission surface.
     * @param coord Coordinate grid containing angular information
     * @param shock Shock object containing the evolution data
     * @param luminosity_dist Luminosity distance to the source
     * @param redshift Redshift of the source
     * <!-- ************************************************************************************** -->
     */
    void observe(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Sets up the Observer for flux calculation at specific observation times.
     * @details Similar to observe(), but also marks required grid points for the given observation times.
     * @param t_obs Array of observation times
     * @param coord Coordinate grid containing angular information
     * @param shock Shock object containing the evolution data (modified to mark required points)
     * @param luminosity_dist Luminosity distance to the source
     * @param redshift Redshift of the source
     * <!-- ************************************************************************************** -->
     */
    void observe_at(Array const& t_obs, Coord const& coord, Shock& shock, Real luminosity_dist, Real redshift);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the specific flux at a single observed frequency
     * @tparam PhotonGrid Types of photon grid objects
     * @param t_obs Array of observation times
     * @param nu_obs Observed frequency
     * @param photons Parameter pack of photon grid objects
     * @details Returns the specific flux (as an Array) for a single observed frequency (nu_obs) by computing the
     *          specific flux over the observation times.
     * @return Array of specific flux values at each observation time
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    Array specific_flux(Array const& t_obs, Real nu_obs, PhotonGrid& photons);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the specific flux at multiple observed frequencies
     * @tparam PhotonGrid Types of photon grid objects
     * @param t_obs Array of observation times
     * @param nu_obs Array of observed frequencies
     * @param photons Parameter pack of photon grid objects
     * @return 2D grid of specific flux values (frequency × time)
     * @details Returns the specific flux (as a MeshGrid) for multiple observed frequencies (nu_obs) by computing
     *          the specific flux for each frequency and assembling the results into a grid. This method accounts for
     *          relativistic beaming and cosmological effects.
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    MeshGrid specific_flux(Array const& t_obs, Array const& nu_obs, PhotonGrid& photons);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the specific flux at a single observed frequency for multiple observation times
     * @tparam PhotonGrid Types of photon grid objects
     * @param t_obs Array of observation times
     * @param nu_obs Observed frequency
     * @param photons Parameter pack of photon grid objects
     * @return Array of specific flux values at each observation time for a single observed frequency
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    Array specific_flux_series(Array const& t_obs, Array const& nu_obs, PhotonGrid& photons);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the integrated flux over a frequency band
     * @tparam PhotonGrid Types of photon grid objects
     * @param t_obs Array of observation times
     * @param band_freq Array of frequency band boundaries
     * @param photons Parameter pack of photon grid objects
     * @details Computes the integrated flux over a frequency band specified by band_freq.
     *          It converts band boundaries to center frequencies, computes the specific flux at each frequency,
     *          and integrates (sums) the flux contributions weighted by the frequency bin widths.
     * @return Array of integrated flux values at each observation time
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    Array flux(Array const& t_obs, Array const& band_freq, PhotonGrid& photons);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the spectrum at multiple observation times
     * @tparam PhotonGrid Types of photon grid objects
     * @param freqs Array of frequencies
     * @param t_obs Array of observation times
     * @param photons Parameter pack of photon grid objects
     * @details Returns the spectra (as a MeshGrid) for multiple observation times by computing the specific flux
     *          for each frequency and transposing the result to get freq × time format.
     * @return 2D grid of spectra (frequency × time)
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    MeshGrid spectra(Array const& freqs, Array const& t_obs, PhotonGrid& photons);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the spectrum at a single observation time
     * @tparam PhotonGrid Types of photon grid objects
     * @param freqs Array of frequencies
     * @param t_obs Single observation time
     * @param photons Parameter pack of photon grid objects
     * @details Returns the spectrum (as an Array) at a single observation time by computing the specific flux
     *          for each frequency in the given array.
     * @return Array containing the spectrum at the given time
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    Array spectrum(Array const& freqs, Real t_obs, PhotonGrid& photons);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Updates the required grid points for observation.
     * @details Identifies grid points needed for interpolation at requested observation times.
     * @param required Mask grid to mark required points (modified in-place)
     * @param t_obs Array of observation times
     * <!-- ************************************************************************************** -->
     */
    void update_required(MaskGrid& required, Array const& t_obs);

    MeshGrid3d lg2_t;           ///< Log2 of observation time grid
    MeshGrid3d lg2_doppler;     ///< Log2 of Doppler factor grid
    MeshGrid3d lg2_geom_factor; ///< Log2 of observe frame geometric factor (solid angle * r^2 * D^3)
  private:
    Real one_plus_z{1}; ///< 1 + redshift
    Real lumi_dist{1};  ///< Luminosity distance

    // Grid dimensions
    size_t jet_3d{0};       ///< Flag indicating if the jet is non-axis-symmetric (non-zero if true)
    size_t eff_phi_grid{1}; ///< Effective number of phi grid points
    size_t theta_grid{0};   ///< Number of theta grid points
    size_t t_grid{0};       ///< Number of time grid points

    /**
     * <!-- ************************************************************************************** -->
     * @brief Builds the time grid and related structures for observation.
     * @details Initializes grid dimensions based on the shock evolution data and sets up arrays for
     *          time, Doppler factor, and emission surface.
     * @param coord Coordinate grid containing angular information
     * @param shock Shock object containing the evolution data
     * @param luminosity_dist Luminosity distance to the source
     * @param redshift Redshift of the source
     * <!-- ************************************************************************************** -->
     */
    void build_time_grid(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the observation time grid and Doppler factor grid.
     * @details For each grid point, computes the Doppler factor based on the Lorentz factor and
     *          calculates the observed time taking redshift into account.
     * @param coord Coordinate grid containing angular information
     * @param shock Shock object containing the evolution data
     * <!-- ************************************************************************************** -->
     */
    void calc_t_obs(Coord const& coord, Shock const& shock);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the observe frame solid angle for each grid point.
     * @details Computes the observe frame solid angle as the product of the differential
     *          cosine of theta and either 2π (if the effective phi size is 1) or the differential phi value.
     * @param coord Coordinate grid containing angular information
     * @param shock Shock object containing the evolution data
     * <!-- ************************************************************************************** -->
     */
    void calc_solid_angle(Coord const& coord, Shock const& shock);

    /**
     * <!-- ************************************************************************************** -->
     * @struct InterpState
     * @brief Helper structure for logarithmic interpolation state
     * <!-- ************************************************************************************** -->
     */
    struct InterpState {
        Real slope{0};       ///< Slope for logarithmic interpolation
        Real lg2_L_nu_lo{0}; ///< Lower boundary of specific luminosity (log2 scale)
        Real lg2_L_nu_hi{0}; ///< Upper boundary of specific luminosity (log2 scale)
        Real last_lg2_nu{0}; ///< Last log2 frequency (for interpolation)
        size_t last_hi{0};   ///< Index for the upper boundary in the grid
    };

    /**
     * <!-- ************************************************************************************** -->
     * @brief Interpolates the luminosity using the observation time (t_obs) in logarithmic space
     * @param state The interpolation state
     * @param lg2_t_obs Observation time (in log2 scale)
     * @param lg2_t_lo Lower boundary of observation time (in log2 scale)
     * @return The interpolated luminosity value
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] static Real loglog_interpolate(InterpState const& state, Real lg2_t_obs, Real lg2_t_lo) noexcept;

    [[nodiscard]] static Real interpolate(InterpState const& state, Real lg2_t_obs, Real lg2_t_lo) noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Validates and sets the interpolation boundaries
     * @tparam PhotonGrid Types of photon grid objects
     * @param state The interpolation state to update
     * @param eff_i Effective phi grid index (accounts for jet symmetry)
     * @param i Phi grid index
     * @param j Theta grid index
     * @param k Time grid index
     * @param lg2_nu_src Log2 of the observed frequency
     * @param photons Parameter pack of photon grid objects
     * @details Attempts to set the lower and upper boundary values for logarithmic interpolation.
     *          It updates the internal boundary members for:
     *            - Logarithmic observation time (t_obs_lo, log_t_ratio)
     *            - Logarithmic luminosity (L_lo, L_hi)
     *          The boundaries are set using data from the provided grids and the photon grids.
     *          Returns true if both lower and upper boundaries are finite such that interpolation can proceed
     * @return True if both lower and upper boundaries are valid for interpolation, false otherwise
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    bool set_boundaries(InterpState& state, size_t eff_i, size_t i, size_t j, size_t k, Real lg2_nu_src,
                        PhotonGrid& photons) noexcept;
};

//========================================================================================================
//                                  template function implementation
//========================================================================================================

/**
 * <!-- ************************************************************************************** -->
 * @brief Helper function that advances an iterator through an array until the value exceeds the target
 * @param value Target value to exceed
 * @param arr Array to iterate through
 * @param it Iterator position (updated by this function)
 * @details Used for efficiently finding the appropriate position in a sorted array without binary search.
 * <!-- ************************************************************************************** -->
 */
inline void iterate_to(Real value, Array const& arr, size_t& it) noexcept {
    while (it < arr.size() && arr(it) < value) {
        it++;
    }
}

template <typename PhotonGrid>
bool Observer::set_boundaries(InterpState& state, size_t eff_i, size_t i, size_t j, size_t k, Real lg2_nu_src,
                              PhotonGrid& photons) noexcept {
    if (state.last_hi == k + 1 && state.last_lg2_nu == lg2_nu_src) {
        if (!std::isfinite(state.slope)) {
            return false;
        } else {
            return true;
        }
    }

    const Real lg2_t_ratio = lg2_t(i, j, k + 1) - lg2_t(i, j, k);

    // continuing from the previous boundary, shift the high boundary to lower.
    // Calling .I_nu()/.log_I_nu() could be expensive.
    if (state.last_hi != 0 && k == state.last_hi && lg2_nu_src == state.last_lg2_nu) {
        state.lg2_L_nu_lo = state.lg2_L_nu_hi;
    } else {
        Real lg2_nu_lo = lg2_nu_src - lg2_doppler(i, j, k);
        state.lg2_L_nu_lo = photons(eff_i, j, k).compute_log2_I_nu(lg2_nu_lo) + lg2_geom_factor(i, j, k);
    }

    Real lg2_nu_hi = lg2_nu_src - lg2_doppler(i, j, k + 1);
    state.lg2_L_nu_hi = photons(eff_i, j, k + 1).compute_log2_I_nu(lg2_nu_hi) + lg2_geom_factor(i, j, k + 1);

    state.slope = (state.lg2_L_nu_hi - state.lg2_L_nu_lo) / lg2_t_ratio;

    if (!std::isfinite(state.slope)) {
        return false;
    }

    state.last_hi = k + 1;
    state.last_lg2_nu = lg2_nu_src;
    return true;
}

template <typename PhotonGrid>
MeshGrid Observer::specific_flux(Array const& t_obs, Array const& nu_obs, PhotonGrid& photons) {
    const size_t t_obs_len = t_obs.size();
    const size_t nu_len = nu_obs.size();

    const Array lg2_t_obs = xt::log2(t_obs);
    const Array lg2_nu_src = xt::log2(nu_obs) + std::log2(one_plus_z);

    MeshGrid F_nu({nu_len, t_obs_len}, 0);
    xt::xtensor<Real, 4> lg2_F_nu_ij({eff_phi_grid, theta_grid, nu_len, t_obs_len}, -con::inf);

    for (size_t i = 0; i < eff_phi_grid; i++) {
        size_t eff_i = i * jet_3d;
        for (size_t j = 0; j < theta_grid; j++) {
            // Skip observation times that are below the grid's start time
            size_t t_idx = 0;
            iterate_to(time(i, j, 0), t_obs, t_idx);

            InterpState state;
            for (size_t k = 0; k < t_grid - 1 && t_idx < t_obs_len; k++) {
                if (const Real t_hi = time(i, j, k + 1); t_hi >= t_obs(t_idx)) {
                    const size_t idx_start = t_idx;
                    iterate_to(t_hi, t_obs, t_idx);
                    const size_t idx_end = t_idx;

                    for (size_t l = 0; l < nu_len; l++) {
                        if (set_boundaries(state, eff_i, i, j, k, lg2_nu_src[l], photons)) [[likely]] {
                            for (size_t idx = idx_start; idx < idx_end; idx++) {
                                //F_nu(l, idx) += loglog_interpolate(state, lg2_t_obs(idx), lg2_t(i, j, k));
                                lg2_F_nu_ij(i, j, l, idx) = interpolate(state, lg2_t_obs(idx), lg2_t(i, j, k));
                            }
                        }
                    }
                }
            }
        }
    }

    lg2_F_nu_ij = xt::exp2(lg2_F_nu_ij);
    for (size_t i = 0; i < eff_phi_grid; i++) {
        for (size_t j = 0; j < theta_grid; j++) {
            for (size_t l = 0; l < nu_len; l++) {
                for (size_t idx = 0; idx < t_obs_len; idx++) {
                    F_nu(l, idx) += lg2_F_nu_ij(i, j, l, idx);
                }
            }
        }
    }

    F_nu *= one_plus_z / (lumi_dist * lumi_dist);

    return F_nu;
}

template <typename PhotonGrid>
Array Observer::specific_flux_series(Array const& t_obs, Array const& nu_obs, PhotonGrid& photons) {
    const size_t t_obs_len = t_obs.size();

    if (nu_obs.size() != t_obs_len) {
        std::cout << "nu_obs and t_obs must have the same length" << std::endl;
    }

    const Array lg2_t_obs = xt::log2(t_obs);
    const Array lg2_nu_src = xt::log2(nu_obs) + std::log2(one_plus_z);

    Array F_nu = xt::zeros<Real>({t_obs_len});
    xt::xtensor<Real, 3> lg2_F_nu_ij({eff_phi_grid, theta_grid, t_obs_len}, -con::inf);

    for (size_t i = 0; i < eff_phi_grid; i++) {
        size_t eff_i = i * jet_3d;
        for (size_t j = 0; j < theta_grid; j++) {
            size_t t_idx = 0;
            iterate_to(time(i, j, 0), t_obs, t_idx);

            InterpState state;
            for (size_t k = 0; t_idx < t_obs_len && k < t_grid - 1;) {
                if (time(i, j, k + 1) < t_obs(t_idx)) {
                    k++;
                } else {
                    if (set_boundaries(state, eff_i, i, j, k, lg2_nu_src(t_idx), photons)) [[likely]] {
                        //F_nu(t_idx) += loglog_interpolate(state, lg2_t_obs(t_idx), lg2_t(i, j, k));
                        lg2_F_nu_ij(i, j, t_idx) = interpolate(state, lg2_t_obs(t_idx), lg2_t(i, j, k));
                    }
                    t_idx++;
                }
            }
        }
    }

    lg2_F_nu_ij = xt::exp2(lg2_F_nu_ij);

    for (size_t i = 0; i < eff_phi_grid; i++) {
        for (size_t j = 0; j < theta_grid; j++) {
            for (size_t t_idx = 0; t_idx < t_obs_len; t_idx++) {
                F_nu(t_idx) += lg2_F_nu_ij(i, j, t_idx);
            }
        }
    }

    // Normalize the flux by the factor (1+z)/(lumi_dist^2).
    F_nu *= one_plus_z / (lumi_dist * lumi_dist);

    return F_nu;
}

template <typename PhotonGrid>
Array Observer::specific_flux(Array const& t_obs, Real nu_obs, PhotonGrid& photons) {
    return xt::view(specific_flux(t_obs, Array({nu_obs}), photons), 0);
}

template <typename PhotonGrid>
Array Observer::spectrum(Array const& freqs, Real t_obs, PhotonGrid& photons) {
    return xt::view(spectra(freqs, Array({t_obs}), photons), 0);
}

template <typename PhotonGrid>
MeshGrid Observer::spectra(Array const& freqs, Array const& t_obs, PhotonGrid& photons) {
    return xt::transpose(specific_flux(t_obs, freqs, photons));
}

template <typename PhotonGrid>
Array Observer::flux(Array const& t_obs, Array const& band_freq, PhotonGrid& photons) {
    Array nu_obs = boundary_to_center_log(band_freq);
    MeshGrid F_nu = specific_flux(t_obs, nu_obs, photons);
    Array flux({t_obs.size()}, 0);
    for (size_t i = 0; i < nu_obs.size(); ++i) {
        const Real dnu = band_freq(i + 1) - band_freq(i);
        for (size_t j = 0; j < flux.size(); ++j) {
            flux(j) += dnu * F_nu(i, j);
        }
    }
    return flux;
}
