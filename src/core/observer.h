//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once

#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

#include "../dynamics/shock.h"
#include "../util/macros.h"
#include "mesh.h"
#include "quadrature.h"

/**
 * @struct SkyImageResult
 * @brief Result of a batched sky image computation: 3D surface brightness map with shared extent and pixel solid angle.
 */
struct SkyImageResult {
    MeshGrid3d image;           ///< [n_frames, npixel, npixel] surface brightness (code units)
    std::array<Real, 4> extent; ///< {x_min, x_max, y_min, y_max} angular extent (rad)
    Real pixel_solid_angle{0};  ///< Pixel solid angle (sr)
};

/// Internal cell grid for Gaussian splatting (sky position, luminosity, kernel width per cell).
struct SplatGrid {
    MeshGrid x, y, L, sigma;
    size_t n_phi{0}, n_theta{0};
    Real xmin{}, xmax{}, ymin{}, ymax{};

    SplatGrid() = default;
    SplatGrid(size_t n_phi, size_t n_theta);

    /// Reset all arrays and bounds for reuse across frames.
    void reset();

    /// Compute per-cell Gaussian σ from neighbor distances.
    void compute_sigma();

    /// Render the splat grid into a frame of the 3D result image.
    void render(MeshGrid3d& image, size_t frame, Real pix_size, bool azimuthal_avg, bool mirror_y) const;
};

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
     * @brief Computes resolved sky images at multiple observer times and a single frequency.
     * @tparam PhotonGrid Type of photon grid object
     * @param coord Coordinate grid (needed for phi, theta_view projection)
     * @param shock Shock object (needed for r, theta at interpolated time)
     * @param t_obs Array of observer times (internal units)
     * @param nu_obs Observer frequency (internal units)
     * @param photons Photon grid object
     * @param npixel Number of pixels per side
     * @param pix_size Pixel size in radians (FOV = pix_size × npixel)
     * @return SkyImageResult with [n_frames, npixel, npixel] surface brightness, shared extent and pixel_solid_angle
     */
    template <typename PhotonGrid>
    SkyImageResult sky_image(Coord const& coord, Shock const& shock, Array const& t_obs, Real nu_obs,
                             PhotonGrid& photons, size_t npixel, Real pix_size);

    MeshGrid3d lg2_t;           ///< Log2 of observation time grid
    MeshGrid3d lg2_doppler;     ///< Log2 of Doppler factor grid
    MeshGrid3d lg2_geom_factor; ///< Log2 of observe frame geometric factor (solid angle * r^2 * D^3)
    Real one_plus_z{1};         ///< 1 + redshift
  private:
    Real lumi_dist{1}; ///< Luminosity distance
    // Grid dimensions
    size_t jet_3d{0};           ///< Flag indicating if the jet is non-axis-symmetric (non-zero if true)
    size_t eff_phi_grid{1};     ///< Effective number of phi grid points
    size_t theta_grid{0};       ///< Number of theta grid points
    size_t t_grid{0};           ///< Number of time grid points
    bool jet_spreading_{false}; ///< Whether the jet has lateral spreading (theta varies along k)

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
     * @brief Fused EAT grid + geometric factor computation for non-spreading jets.
     * @details Computes time, lg2_t, lg2_doppler, and lg2_geom_factor in a single pass,
     *          avoiding redundant memory reads of shock.r and lg2_doppler.
     * @param coord Coordinate grid containing angular information
     * @param shock Shock object containing the evolution data
     * <!-- ************************************************************************************** -->
     */
    void calc_eat_non_spreading(Coord const& coord, Shock const& shock);

    /// Fills an existing splat grid (sky positions, luminosities, Gaussian widths) for a single frame.
    template <typename PhotonGrid>
    void build_splat_grid(SplatGrid& grid, Coord const& coord, Shock const& shock, Real t_obs, Real nu_obs,
                          PhotonGrid& photons, std::vector<Real> const& cos_phi_tab,
                          std::vector<Real> const& sin_phi_tab);

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
     * @brief Interpolates the specific luminosity in log-log space, returning the result in linear scale.
     * @param state The interpolation state
     * @param lg2_t_obs Observation time (in log2 scale)
     * @param lg2_t_lo Lower boundary of observation time (in log2 scale)
     * @return The interpolated specific luminosity (linear scale)
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] static Real loglog_interpolate(InterpState const& state, Real lg2_t_obs, Real lg2_t_lo) noexcept;

    /**
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
        return std::isfinite(state.slope);
    }

    const Real lg2_t_ratio = lg2_t(i, j, k + 1) - lg2_t(i, j, k);

    // Reuse previous hi boundary as new lo when advancing sequentially
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

    const Real* lg2t = lg2_t_obs.data();
    const Real* lg2nu = lg2_nu_src.data();
    Real* fnu_ptr = F_nu.data();

    // Precomputed log2-luminosity at each time boundary × frequency. [k * nu_len + l]
    Array boundary_vals({t_grid * nu_len}, 0);
    Array slope_buf({nu_len}, 0);
    Array lo_buf({nu_len}, 0);
    Real* slope = slope_buf.data();
    Real* lo = lo_buf.data();

    for (size_t i = 0; i < eff_phi_grid; i++) {
        const size_t eff_i = i * jet_3d;
        for (size_t j = 0; j < theta_grid; j++) {
            // Precompute boundary values for all time grid points.
            for (size_t k = 0; k < t_grid; k++) {
                const Real lg2_dop = lg2_doppler(i, j, k);
                const Real lg2_geom = lg2_geom_factor(i, j, k);
                Real* vals = boundary_vals.data() + k * nu_len;
                for (size_t l = 0; l < nu_len; l++) {
                    vals[l] = photons(eff_i, j, k).compute_log2_I_nu(lg2nu[l] - lg2_dop) + lg2_geom;
                }
            }

            size_t t_idx = 0;
            iterate_to(lg2_t(i, j, 0), lg2_t_obs, t_idx);

            for (size_t k = 0; k < t_grid - 1 && t_idx < t_obs_len; k++) {
                if (lg2_t(i, j, k + 1) < lg2_t_obs(t_idx))
                    continue;

                const size_t idx_start = t_idx;
                iterate_to(lg2_t(i, j, k + 1), lg2_t_obs, t_idx);

                const Real* lo_row = boundary_vals.data() + k * nu_len;
                const Real* hi_row = boundary_vals.data() + (k + 1) * nu_len;
                const Real t_lo_val = lg2_t(i, j, k);
                const Real inv_t_ratio = 1.0 / (lg2_t(i, j, k + 1) - t_lo_val);

                for (size_t l = 0; l < nu_len; l++) {
                    const Real s = (hi_row[l] - lo_row[l]) * inv_t_ratio;
                    const bool ok = std::isfinite(s);
                    slope[l] = ok ? s : 0.0;
                    lo[l] = ok ? lo_row[l] : -con::inf;
                }
                for (size_t idx = idx_start; idx < t_idx; idx++) {
                    const Real dlg2_t = lg2t[idx] - t_lo_val;
                    for (size_t l = 0; l < nu_len; l++) {
                        fnu_ptr[l * t_obs_len + idx] += fast_exp2(lo[l] + dlg2_t * slope[l]);
                    }
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
        throw std::invalid_argument("specific_flux_series requires nu_obs and t_obs to have the same length");
    }

    const Array lg2_t_obs = xt::log2(t_obs);
    const Array lg2_nu_src = xt::log2(nu_obs) + std::log2(one_plus_z);

    Array F_nu = xt::zeros<Real>({t_obs_len});

    const Real* lg2t = lg2_t_obs.data();
    const Real* lg2nu = lg2_nu_src.data();
    Real* fnu = F_nu.data();

    for (size_t i = 0; i < eff_phi_grid; i++) {
        const size_t eff_i = i * jet_3d;
        for (size_t j = 0; j < theta_grid; j++) {
            size_t idx = 0;
            iterate_to(lg2_t(i, j, 0), lg2_t_obs, idx);
            InterpState state;

            for (size_t k = 0; idx < t_obs_len && k < t_grid - 1;) {
                if (lg2_t(i, j, k + 1) < lg2t[idx]) {
                    k++;
                } else {
                    if (set_boundaries(state, eff_i, i, j, k, lg2nu[idx], photons)) [[likely]] {
                        fnu[idx] += loglog_interpolate(state, lg2t[idx], lg2_t(i, j, k));
                    }
                    idx++;
                }
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
    MeshGrid F_nu = specific_flux(t_obs, band_freq, photons);
    Array weights;
    compute_boole_weights(band_freq, weights);
    Array flux({t_obs.size()}, 0);
    for (size_t i = 0; i < band_freq.size(); ++i) {
        for (size_t j = 0; j < flux.size(); ++j) {
            flux(j) += F_nu(i, j) * weights(i);
        }
    }
    return flux;
}

//========================================================================================================
//                                  sky_image template implementation
//========================================================================================================

template <typename PhotonGrid>
void Observer::build_splat_grid(SplatGrid& grid, Coord const& coord, Shock const& shock, Real t_obs, Real nu_obs,
                                PhotonGrid& photons, std::vector<Real> const& cos_phi_tab,
                                std::vector<Real> const& sin_phi_tab) {
    const Real lg2_t_target = fast_log2(t_obs);
    const Real lg2_nu_src = fast_log2(nu_obs * one_plus_z);
    const Real inv_d_A = one_plus_z * one_plus_z / lumi_dist;
    const Real cos_obs = std::cos(coord.theta_view);
    const Real sin_obs = std::sin(coord.theta_view);
    const bool is_axisym = (eff_phi_grid == 1);
    const size_t copies_per_phys = grid.n_phi / eff_phi_grid;

    for (size_t i_phys = 0; i_phys < eff_phi_grid; ++i_phys) {
        const size_t eff_i = i_phys * jet_3d;
        const Real cp_fixed = is_axisym ? 0 : std::cos(coord.phi[i_phys] - coord.phi_view);
        const Real sp_fixed = is_axisym ? 0 : std::sin(coord.phi[i_phys] - coord.phi_view);

        for (size_t j = 0; j < theta_grid; ++j) {
            if (lg2_t_target < lg2_t(i_phys, j, 0) || lg2_t_target >= lg2_t(i_phys, j, t_grid - 1)) {
                continue;
            }

            // Find time bracket
            size_t k = 0;
            while (k + 2 < t_grid && lg2_t(i_phys, j, k + 1) < lg2_t_target) {
                ++k;
            }

            // Interpolate luminosity
            InterpState state_L;
            if (!set_boundaries(state_L, eff_i, i_phys, j, k, lg2_nu_src, photons)) {
                continue;
            }
            const Real L_val = loglog_interpolate(state_L, lg2_t_target, lg2_t(i_phys, j, k));
            if (!std::isfinite(L_val) || L_val <= 0) {
                continue;
            }

            // Interpolate position
            const Real w = (lg2_t_target - lg2_t(i_phys, j, k)) / (lg2_t(i_phys, j, k + 1) - lg2_t(i_phys, j, k));
            const Real r_ang = (shock.r(eff_i, j, k) + w * (shock.r(eff_i, j, k + 1) - shock.r(eff_i, j, k))) * inv_d_A;
            const Real theta_val = jet_spreading_ ? shock.theta(eff_i, j, k) +
                                                        w * (shock.theta(eff_i, j, k + 1) - shock.theta(eff_i, j, k))
                                                  : shock.theta(eff_i, j, 0);
            const Real ct = std::cos(theta_val);
            const Real st = std::sin(theta_val);
            const Real L_per_copy = L_val / copies_per_phys;

            // Project onto sky plane for each azimuthal copy
            for (size_t c = 0; c < copies_per_phys; ++c) {
                const Real cp = is_axisym ? cos_phi_tab[c] : cp_fixed;
                const Real sp = is_axisym ? sin_phi_tab[c] : sp_fixed;
                const Real x_ang = r_ang * (st * cp * cos_obs - ct * sin_obs);
                const Real y_ang = r_ang * st * sp;

                const size_t i_splat = i_phys * copies_per_phys + c;
                grid.x(i_splat, j) = x_ang;
                grid.y(i_splat, j) = y_ang;
                grid.L(i_splat, j) = L_per_copy;
                grid.xmin = std::min(grid.xmin, x_ang);
                grid.xmax = std::max(grid.xmax, x_ang);
                grid.ymin = std::min(grid.ymin, y_ang);
                grid.ymax = std::max(grid.ymax, y_ang);
            }
        }
    }

    // Convert luminosity to surface brightness
    grid.L *= one_plus_z / (lumi_dist * lumi_dist);
    grid.compute_sigma();
}

template <typename PhotonGrid>
SkyImageResult Observer::sky_image(Coord const& coord, Shock const& shock, Array const& t_obs, Real nu_obs,
                                   PhotonGrid& photons, size_t npixel, Real pix_size) {
    const size_t n_frames = t_obs.size();
    SkyImageResult result;
    result.image = MeshGrid3d({n_frames, npixel, npixel}, 0);

    if (npixel == 0 || t_grid < 2) {
        result.extent = {0, 0, 0, 0};
        return result;
    }

    // Compute extent and pixel_solid_angle once (depends only on pix_size and npixel)
    const Real half_fov = 0.5 * pix_size * static_cast<Real>(npixel);
    result.extent = {-half_fov, half_fov, -half_fov, half_fov};
    result.pixel_solid_angle = pix_size * pix_size;

    // 2D Gaussian splatting
    const size_t n_splat_phi = (eff_phi_grid == 1) ? 4 * theta_grid : eff_phi_grid;
    SplatGrid grid(n_splat_phi, theta_grid);

    const size_t copies_per_phys = n_splat_phi / eff_phi_grid;
    std::vector<Real> cos_phi_tab, sin_phi_tab;
    if (eff_phi_grid == 1) {
        cos_phi_tab.resize(copies_per_phys);
        sin_phi_tab.resize(copies_per_phys);
        for (size_t c = 0; c < copies_per_phys; ++c) {
            const Real phi = 2.0 * con::pi * c / n_splat_phi;
            cos_phi_tab[c] = std::cos(phi);
            sin_phi_tab[c] = std::sin(phi);
        }
    }

    const bool is_axisym = (jet_3d == 0);
    const bool azimuthal_avg = is_axisym && (coord.theta_view == 0);
    const bool mirror_y = is_axisym && !azimuthal_avg;

    for (size_t f = 0; f < n_frames; ++f) {
        grid.reset();
        build_splat_grid(grid, coord, shock, t_obs(f), nu_obs, photons, cos_phi_tab, sin_phi_tab);
        grid.render(result.image, f, pix_size, azimuthal_avg, mirror_y);
    }
    return result;
}
