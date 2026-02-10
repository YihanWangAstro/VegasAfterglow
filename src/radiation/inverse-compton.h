//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <algorithm>
#include <array>

#include "../core/mesh.h"
#include "../core/physics.h"
#include "../dynamics/shock.h"
#include "../util/macros.h"
#include "../util/profiler.h"
#include "../util/utilities.h"

struct SpectralSegment {
    Real gamma_max; // Upper energy bound of this segment
    Real norm;      // Normalization (A) at the reference energy
    Real slope;     // Power law index (s)
    Real ref;       // Reference energy (gamma_ref)

    inline Real eval(Real gamma) const { return norm * fast_pow(gamma / ref, slope); }
};

/**
 * <!-- ************************************************************************************** -->
 * @struct InverseComptonY
 * @brief Handles Inverse Compton Y parameter calculations and related threshold values.
 * <!-- ************************************************************************************** -->
 */
struct InverseComptonY {
    /**
     * <!-- ************************************************************************************** -->
     * @brief Initializes an InverseComptonY object with frequency thresholds, magnetic field and Y parameter.
     * @details Computes characteristic gamma values and corresponding frequencies, then determines cooling regime.
     * @param gamma_m Characteristic Lorentz factor for the minimum Lorentz factor
     * @param gamma_c Characteristic Lorentz factor for the cooling Lorentz factor
     * @param p Spectral index of electron distribution
     * @param B Magnetic field strength
     * @param Y_T Thomson scattering Y parameter
     * <!-- ************************************************************************************** -->
     */
    InverseComptonY(Real gamma_m, Real gamma_c, Real p, Real B, Real Y_T, bool is_KN) noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Default constructor that initializes all member variables to zero.
     * <!-- ************************************************************************************** -->
     */
    InverseComptonY() noexcept;

    // Member variables
    Real gamma_m_hat{1}; ///< Lorentz factor threshold for minimum energy electrons
    Real gamma_c_hat{1}; ///< Lorentz factor threshold for cooling electrons
    Real gamma_self{1};  ///< Lorentz factor threshold for gamma_self = gamma_self_hat
    Real gamma0{1};      ///< Lorentz factor threshold for Y(gamma0) = 1
    Real Y_T{0};         ///< Thomson scattering Y parameter
    size_t regime{0};    ///< Indicator for the operating regime (1=fast IC cooling, 2=slow IC cooling, 3=special case)

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter for a given Lorentz factor and spectral index.
     * @details Different scaling relations apply depending on the cooling regime and gamma value.
     * @param gamma Electron Lorentz factor
     * @return The effective Y parameter at the given gamma
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real gamma_spectrum(Real gamma) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter for a given synchrotron frequency.
     * @details Different scaling relations apply depending on the cooling regime and frequency value.
     * @param nu Synchrotron frequency
     * @return The effective Y parameter at the given frequency
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real nu_spectrum(Real nu) const;

    void update_cooling_breaks(Real gamma_c, Real Y_T) noexcept;

  private:
    Real gamma_m_{1};    ///< Characteristic Lorentz factor for the minimum Lorentz factor
    Real B_{0};          ///< Magnetic field strength
    Real p_{2.3};        ///< Spectral index of electron distribution
    Real gamma_self3{1}; ///< Precomputed gamma_self^3 for efficiency
    void update_gamma0(Real gamma_c) noexcept;

    Real compute_gamma_hat(Real gamma) const noexcept;

    int active_segment_count = 0;
    std::array<SpectralSegment, 5> segments; // Max 5 segments needed for complex regimes

    void build_segments() noexcept;
};

/**
 * @brief Fast KN correction lookup: returns sigma_KN/sigmaT at frequency nu.
 * @param nu Photon frequency.
 * @return Dimensionless correction factor to Thomson cross-section.
 */
Real compton_correction(Real nu);

/**
 * <!-- ************************************************************************************** -->
 * @struct ICPhoton
 * @tparam Electrons Type of the electron distribution
 * @tparam Photons Type of the photon distribution
 * @brief Represents a single inverse Compton (IC) photon.
 * @details Contains methods to compute the photon intensity I_nu and to generate an IC photon spectrum based
 *          on electron and synchrotron photon properties.
 * <!-- ************************************************************************************** -->
 */
template <typename Electrons, typename Photons>
struct ICPhoton {
  public:
    /// Default constructor
    ICPhoton() = default;

    ICPhoton(Electrons const& electrons, Photons const& photons, bool KN) noexcept;
    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns the photon-specific intensity.
     * @param nu The frequency at which to compute the specific intensity
     * @return The specific intensity at the given frequency
     * <!-- ************************************************************************************** -->
     */
    Real compute_I_nu(Real nu);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the base-2 logarithm of the photon-specific intensity at a given frequency.
     * @param log2_nu The base-2 logarithm of the frequency
     * @return The base-2 logarithm of the photon specific intensity at the given frequency
     * <!-- ************************************************************************************** -->
     */
    Real compute_log2_I_nu(Real log2_nu);

    Photons photons;

    Electrons electrons;

  private:
    void generate_spectrum();

    /**
     * <!-- ************************************************************************************** -->
     * @brief Grid parameter structure holding ranges and sizes for gamma and nu
     * <!-- ************************************************************************************** -->
     */
    struct GridParams {
        Real gamma_min;
        Real gamma_max;
        size_t gamma_size;
        Real nu_min;
        Real nu_max;
        Real nu_IC_min;
        Real nu_IC_max;
        size_t spectrum_resol;
        std::array<Real, 4> nu_breaks;
        std::array<Real, 12> nu_IC_breaks{};
        std::array<Real, 12> nu_IC_break_scores{};
    };

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes grid parameters based on electron and photon distributions
     * @return GridParams structure containing all computed parameters
     * <!-- ************************************************************************************** -->
     */
    GridParams compute_grid_params() const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Initializes frequency and intensity grids
     * @param params Grid parameters
     * <!-- ************************************************************************************** -->
     */
    void initialize_grids(GridParams const& params);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Preprocesses electron and photon data into arrays for integration
     * @param params Grid parameters
     * @param gamma Output array for gamma values
     * @param nu_sync Output array for synchrotron frequencies
     * @param dnu_sync Output array for synchrotron frequency bin widths
     * <!-- ************************************************************************************** -->
     */
    void preprocess_distributions(GridParams const& params, Array& gamma, Array& nu_sync, Array& dnu_sync) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the IC spectrum through integration over electron and photon distributions
     * @param params Grid parameters
     * @param gamma Array of gamma values
     * @param nu_sync Array of synchrotron frequencies
     * @param dnu_sync Array of synchrotron frequency bin widths
     * <!-- ************************************************************************************** -->
     */
    void compute_IC_spectrum(GridParams const& params, Array const& gamma, Array const& nu_sync, Array const& dnu_sync);

    Array log2_nu_IC;

    Array log2_I_nu_IC;

    Array interp_slope;

    static constexpr size_t gamma_grid_per_order{4};

    static constexpr Real nu_grid_per_order{3};

    static constexpr Real ic_grid_per_order{3};

    bool KN{false}; // Klein-Nishina flag

    bool generated{false};
};

/// Defines a 3D grid (using xt::xtensor) for storing ICPhoton objects.
template <typename Electrons, typename Photons>
using ICPhotonGrid = xt::xtensor<ICPhoton<Electrons, Photons>, 3>;

template <typename Electrons>
using ElectronGrid = xt::xtensor<Electrons, 3>;

template <typename Photons>
using PhotonGrid = xt::xtensor<Photons, 3>;

inline constexpr Real IC_x0 = 0.47140452079103166;
/**
 * <!-- ************************************************************************************** -->
 * @brief Creates and generates an IC photon grid from electron and photon distributions
 * @tparam Electrons Type of the electron distribution
 * @tparam Photons Type of the photon distribution
 * @param electrons The electron grid
 * @param photons The photon grid
 * @param KN flag for Klein-Nishina
 * @param conical flag for conical solution
 * @return A 3D grid of IC photons
 * <!-- ************************************************************************************** -->
 */
template <typename Electrons, typename Photons>
ICPhotonGrid<Electrons, Photons> generate_IC_photons(ElectronGrid<Electrons> const& electrons,
                                                     PhotonGrid<Photons> const& photons, bool KN,
                                                     Shock const& shock) noexcept;

/**
 * <!-- ************************************************************************************** -->
 * @brief Applies Thomson cooling to electrons based on photon distribution
 * @tparam Electrons Type of the electron distribution
 * @tparam Photons Type of the photon distribution
 * @param electrons The electron grid to be modified
 * @param photons The photon grid
 * @param shock The shock properties
 * <!-- ************************************************************************************** -->
 */
template <typename Electrons, typename Photons>
void Thomson_cooling(ElectronGrid<Electrons>& electrons, PhotonGrid<Photons>& photons, Shock const& shock);

/**
 * <!-- ************************************************************************************** -->
 * @brief Applies Klein-Nishina cooling to electrons based on photon distribution
 * @tparam Electrons Type of the electron distribution
 * @tparam Photons Type of the photon distribution
 * @param electrons The electron grid to be modified
 * @param photons The photon grid
 * @param shock The shock properties
 * <!-- ************************************************************************************** -->
 */
template <typename Electrons, typename Photons>
void KN_cooling(ElectronGrid<Electrons>& electrons, PhotonGrid<Photons>& photons, Shock const& shock);

//========================================================================================================
//                                  template function implementation
//========================================================================================================

template <typename Electrons, typename Photons>
ICPhoton<Electrons, Photons>::ICPhoton(Electrons const& electrons, Photons const& photons, bool KN) noexcept
    : photons(photons), electrons(electrons), KN(KN) {}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::generate_spectrum() {
    GridParams params = compute_grid_params();
    initialize_grids(params);

    Array gamma, nu_sync, dnu_sync;
    preprocess_distributions(params, gamma, nu_sync, dnu_sync);
    compute_IC_spectrum(params, gamma, nu_sync, dnu_sync);

    generated = true;
}

template <typename Electrons, typename Photons>
typename ICPhoton<Electrons, Photons>::GridParams ICPhoton<Electrons, Photons>::compute_grid_params() const {
    GridParams params;

    constexpr Real eps_tail = 1e-2;     // target suppression level for single-cutoff tail
    constexpr Real safety_margin = 2.0; // retain headroom beyond estimated tail
    constexpr Real min_tail_factor = 5.0;

    const Real tail_factor = std::max(-std::log(eps_tail), min_tail_factor);

    params.gamma_min = std::min(electrons.gamma_m, electrons.gamma_c);
    params.gamma_max = std::max(electrons.gamma_M * tail_factor, params.gamma_min * 1.01);
    params.gamma_size =
        static_cast<size_t>(std::max(std::log10(params.gamma_max / params.gamma_min), 3.) * gamma_grid_per_order);

    params.nu_min = std::min(photons.nu_a, photons.nu_m) / 10;
    params.nu_max = std::max(photons.nu_M * tail_factor, params.nu_min * 1.01);

    params.nu_breaks = {photons.nu_a, photons.nu_m, photons.nu_c, photons.nu_M};

    params.nu_IC_min = 4 * IC_x0 * params.nu_min * params.gamma_min * params.gamma_min;

    const Real nu_ic_base = 4 * IC_x0 * photons.nu_M * electrons.gamma_M * electrons.gamma_M;
    const Real nu_ic_single_cut = std::max(nu_ic_base * tail_factor * tail_factor, nu_ic_base * tail_factor);
    params.nu_IC_max = nu_ic_single_cut * safety_margin;

    // Candidate IC breaks from synchrotron/electron break combinations.
    const std::array<Real, 4> nu_seed_breaks = {photons.nu_a, photons.nu_m, photons.nu_c, photons.nu_M};
    const std::array<Real, 3> gamma_breaks = {electrons.gamma_m, electrons.gamma_c, electrons.gamma_M};

    static_assert(std::tuple_size_v<decltype(params.nu_IC_breaks)> == nu_seed_breaks.size() * gamma_breaks.size());

    size_t k = 0;
    for (Real nu_seed : nu_seed_breaks) {
        for (Real gb : gamma_breaks) {
            const Real nu_ic = 4 * IC_x0 * nu_seed * gb * gb;
            params.nu_IC_breaks[k] = nu_ic;

            Real score = 0;
            if ((nu_seed > 0) && std::isfinite(nu_seed) && (gb > 0) && std::isfinite(gb)) {
                const Real ne = std::max(std::abs(electrons.compute_column_den(gb)), Real(0));
                const Real i_seed = std::max(photons.compute_I_nu(nu_seed), Real(0));
                const Real inv_nu2 = 1.0 / (nu_seed * nu_seed);
                const Real kn_weight = KN ? std::max(compton_correction(gb * nu_seed), Real(0)) : 1;
                score = ne * i_seed * inv_nu2 * kn_weight;
            }
            params.nu_IC_break_scores[k] = (std::isfinite(score) && score >= 0) ? score : 0;
            ++k;
        }
    }

    params.spectrum_resol = static_cast<size_t>(std::max(5 * std::log10(params.nu_IC_max / (params.nu_IC_min)), 5.));

    return params;
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::initialize_grids(GridParams const& params) {
    const size_t safe_resol = std::max(params.spectrum_resol, static_cast<size_t>(2));

    Array nu_ic;
    build_adaptive_grid(std::log2(params.nu_IC_min), std::log2(params.nu_IC_max),
                        std::span<const Real>(params.nu_IC_breaks), std::span<const Real>(params.nu_IC_break_scores),
                        ic_grid_per_order, nu_ic);

    // Ensure minimum spectrum resolution to avoid out-of-bounds access.
    if (nu_ic.size() < 2) {
        log2_nu_IC = xt::linspace(std::log2(params.nu_IC_min), std::log2(params.nu_IC_max), safe_resol);
    } else {
        log2_nu_IC = xt::log2(nu_ic);
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::preprocess_distributions(GridParams const& params, Array& gamma, Array& nu_sync,
                                                            Array& dnu_sync) const {
    // Adaptive grid for nu_sync (CDF + interpolation absorbs grid rearrangement)
    build_adaptive_grid(std::log2(params.nu_min), std::log2(params.nu_max), std::span<const Real>(params.nu_breaks),
                        nu_grid_per_order, nu_sync, dnu_sync);

    // Uniform log grid for gamma (centers only; bin widths = gamma * width_factor, computed in compute_IC_spectrum)
    logspace_center(std::log2(params.gamma_min), std::log2(params.gamma_max), params.gamma_size, gamma);
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::compute_IC_spectrum(GridParams const& params, Array const& gamma,
                                                       Array const& nu_sync, Array const& dnu_sync) {
    const size_t spec_size = log2_nu_IC.size();
    log2_I_nu_IC = Array::from_shape({spec_size});

    Array I_nu_IC = xt::zeros<Real>({spec_size});
    Array cdf = Array::from_shape({nu_sync.size()});
    Array f_vals = Array::from_shape({nu_sync.size()});
    Array inv_dnu = Array::from_shape({nu_sync.size() - 1});
    Array inv_nu2 = Array::from_shape({nu_sync.size()});
    Array nu_IC = xt::exp2(log2_nu_IC);

    const int nu_last = static_cast<int>(nu_sync.size()) - 1;
    const int spec_last = static_cast<int>(spec_size) - 1;

    // Uniform log gamma grid: dgamma(i) = gamma(i) * width_factor (constant ratio)
    const Real r_gamma =
        std::exp2((std::log2(params.gamma_max) - std::log2(params.gamma_min)) / static_cast<Real>(params.gamma_size));
    const Real width_factor = (r_gamma - 1.0) / std::sqrt(r_gamma);

    const size_t nu_size = nu_sync.size();

    // Raw pointers for hot-path array access (xtensor operator() has stride overhead)
    const Real* nu_p = nu_sync.data();
    const Real* dnu_p = dnu_sync.data();
    Real* cdf_p = cdf.data();
    Real* fv_p = f_vals.data();
    Real* I_p = I_nu_IC.data();
    Real* inv_dnu_p = inv_dnu.data();
    Real* inv_nu2_p = inv_nu2.data();
    const Real* nu_ic_p = nu_IC.data();

    // Precompute reciprocal interval widths (avoids division in hot loop)
    for (size_t j = 0; j + 1 < nu_size; ++j)
        inv_dnu_p[j] = 1.0 / dnu_p[j];
    for (size_t j = 0; j < nu_size; ++j)
        inv_nu2_p[j] = 1.0 / (nu_p[j] * nu_p[j]);

    // CDF lookup at arbitrary nu_seed: partial-bin trapezoidal interpolation
    auto cdf_at = [&](int idx, Real nu_seed) {
        const Real frac = std::max((nu_seed - nu_p[idx]) * inv_dnu_p[idx], Real(0));
        const Real rem = 1.0 - frac;
        const Real f_seed = fv_p[idx] * rem + fv_p[idx + 1] * frac;
        return std::max(cdf_p[idx + 1] + 0.5 * (f_seed + fv_p[idx + 1]) * rem * dnu_p[idx], Real(0));
    };

    const Real nu_sync_max = nu_p[nu_last];

    // Per-gamma upper nu_IC bound from nu_seed <= nu_sync_max.
    // gamma is monotonic, so k_max is monotonic too; a rolling cursor beats per-gamma binary search.
    auto advance_k_max = [&](Real gi2, int& k_max_cursor) {
        const Real nu_ic_limit = 4 * IC_x0 * gi2 * nu_sync_max;
        while (k_max_cursor + 1 <= spec_last && nu_ic_p[k_max_cursor + 1] <= nu_ic_limit) {
            ++k_max_cursor;
        }
        return k_max_cursor;
    };

    // Scan CDF and accumulate I_nu_IC for one electron energy bin.
    auto accumulate = [&](Real weight, Real inv_4x0_gi2, int k_max) {
        int scan_idx = 0;
        for (int k = 0; k <= k_max; ++k) {
            const Real nu_seed = nu_ic_p[k] * inv_4x0_gi2;
            while (scan_idx < nu_last && nu_p[scan_idx + 1] <= nu_seed) {
                ++scan_idx;
            }
            if (scan_idx >= nu_last) {
                break; // CDF is zero beyond the grid
            }
            I_p[k] += weight * cdf_at(scan_idx, nu_seed);
        }
    };

    if (!KN) {
        int k_max_cursor = -1;
        fv_p[nu_last] = photons.compute_I_nu(nu_p[nu_last]) * inv_nu2_p[nu_last];
        cdf_p[nu_last] = 0;
        for (int j = nu_last - 1; j >= 0; --j) {
            fv_p[j] = photons.compute_I_nu(nu_p[j]) * inv_nu2_p[j];
            cdf_p[j] = cdf_p[j + 1] + 0.5 * (fv_p[j] + fv_p[j + 1]) * dnu_p[j];
        }

        for (size_t i = 0; i < gamma.size(); ++i) {
            const Real gi = gamma(i);
            const Real dN_e = electrons.compute_column_den(gi) * width_factor;
            if (dN_e <= 0) {
                continue;
            }
            const Real gi2 = gi * gi;
            const int k_max = advance_k_max(gi2, k_max_cursor);
            if (k_max < 0) {
                continue;
            }
            accumulate(dN_e / gi, 1 / (4 * IC_x0 * gi2), k_max);
        }
    } else {
        int k_max_cursor = -1;
        Array I_nu_sync = Array::from_shape({nu_sync.size()});
        Real* Is_p = I_nu_sync.data();
        for (int j = 0; j <= nu_last; ++j)
            Is_p[j] = photons.compute_I_nu(nu_p[j]);

        for (size_t i = 0; i < gamma.size(); ++i) {
            const Real gi = gamma(i);
            const Real dN_e = electrons.compute_column_den(gi) * width_factor;
            if (dN_e <= 0) {
                continue;
            }

            {
                const Real nu = nu_p[nu_last];
                fv_p[nu_last] = Is_p[nu_last] * compton_correction(gi * nu) * inv_nu2_p[nu_last];
            }
            cdf_p[nu_last] = 0;
            for (int j = nu_last - 1; j >= 0; --j) {
                const Real nu = nu_p[j];
                fv_p[j] = Is_p[j] * compton_correction(gi * nu) * inv_nu2_p[j];
                cdf_p[j] = cdf_p[j + 1] + 0.5 * (fv_p[j] + fv_p[j + 1]) * dnu_p[j];
            }

            const Real gi2 = gi * gi;
            const int k_max = advance_k_max(gi2, k_max_cursor);
            if (k_max < 0) {
                continue;
            }
            accumulate(dN_e / gi, 1 / (4 * IC_x0 * gi2), k_max);
        }
    }

    const Real log2_scale = fast_log2(0.25 * con::sigmaT);
    for (size_t i = 0; i < spec_size; ++i) {
        log2_I_nu_IC(i) = fast_log2(I_nu_IC(i)) + log2_nu_IC(i) + log2_scale;
    }

    interp_slope = Array::from_shape({spec_size - 1});
    for (size_t i = 0; i + 1 < spec_size; ++i) {
        const Real dlog2 = log2_nu_IC(i + 1) - log2_nu_IC(i);
        interp_slope(i) = (dlog2 != 0) ? (log2_I_nu_IC(i + 1) - log2_I_nu_IC(i)) / dlog2 : 0;
    }

    /*AFTERGLOW_PROFILE_COUNT(ic_spectra_generated);
    AFTERGLOW_PROFILE_COUNT_N(ic_nu_ic_bins_total, spec_size);
    AFTERGLOW_PROFILE_MAX(ic_nu_ic_bins_max, spec_size);
    AFTERGLOW_PROFILE_COUNT_N(ic_nu_ic_bins_uniform_total, params.spectrum_resol);
    AFTERGLOW_PROFILE_MAX(ic_nu_ic_bins_uniform_max, params.spectrum_resol);*/
}

template <typename Electrons, typename Photons>
Real ICPhoton<Electrons, Photons>::compute_I_nu(Real nu) {
    return std::exp2(compute_log2_I_nu(std::log2(nu)));
}

template <typename Electrons, typename Photons>
Real ICPhoton<Electrons, Photons>::compute_log2_I_nu(Real log2_nu) {
    if (generated == false) {
        generate_spectrum();
    }

    const size_t n = log2_nu_IC.size();
    if (n < 2) {
        return -con::inf;
    }

    // Short tables are common here; linear scan is usually faster than binary search.
    size_t idx = 0;
    while (idx + 1 < n && log2_nu_IC(idx + 1) <= log2_nu) {
        ++idx;
    }

    if (idx >= interp_slope.size()) {
        return -con::inf;
    } else [[likely]] {
        return log2_I_nu_IC(idx) + (log2_nu - log2_nu_IC(idx)) * interp_slope(idx);
    }
}

template <typename Electrons, typename Photons>
ICPhotonGrid<Electrons, Photons> generate_IC_photons(ElectronGrid<Electrons> const& electrons,
                                                     PhotonGrid<Photons> const& photons, bool KN,
                                                     Shock const& shock) noexcept {
    size_t phi_size = electrons.shape()[0];
    size_t theta_size = electrons.shape()[1];
    size_t t_size = electrons.shape()[2];

    ICPhotonGrid<Electrons, Photons> IC_ph({phi_size, theta_size, t_size});

    const size_t phi_compute = (shock.symmetry != Symmetry::structured) ? 1 : phi_size;

    // Compute only representative cells based on symmetry
    for (size_t i = 0; i < phi_compute; ++i) {
        for (size_t j : shock.theta_reps) {
            for (size_t k = 0; k < t_size; ++k) {
                IC_ph(i, j, k) = ICPhoton(electrons(i, j, k), photons(i, j, k), KN);
            }
        }
    }

    // Eagerly generate spectra for representative cells before broadcasting.
    // Avoids redundant spectrum generation in broadcast copies.
    if (shock.symmetry != Symmetry::structured) {
        for (size_t i = 0; i < phi_compute; ++i) {
            for (size_t j : shock.theta_reps) {
                for (size_t k = 0; k < t_size; ++k)
                    IC_ph(i, j, k).compute_log2_I_nu(0);
            }
        }
    }

    // Broadcast based on symmetry
    broadcast_symmetry(IC_ph, shock);

    return IC_ph;
}

inline Real eta_rad_Thomson(Real gamma_m, Real gamma_c, Real p) {
    if (gamma_c < gamma_m) {
        return 1;
    } else {
        return fast_pow(gamma_c / gamma_m, 2 - p);
    }
}

inline Real compute_Thomson_Y(RadParams const& rad, Real gamma_m, Real gamma_c) {
    Real eta_e = eta_rad_Thomson(gamma_m, gamma_c, rad.p);
    Real b = eta_e * rad.eps_e / rad.eps_B;
    return 0.5 * (std::sqrt(1. + 4. * b) - 1.);
}

Real compute_gamma_c(Real t_comv, Real B, Real Y);

Real cool_after_crossing(Real gamma_x, Real gamma_m_x, Real gamma_m, Real dt_comv, Real B, Real Y);

Real compute_syn_gamma_M(Real B, Real Y, Real p);

Real compute_syn_I_peak(Real B, Real p, Real column_den);

Real compute_syn_gamma_a(Real B, Real I_nu_peak, Real gamma_m, Real gamma_c, Real gamma_M, Real p);

size_t determine_regime(Real gamma_a, Real gamma_c, Real gamma_m);

Real compute_gamma_0(Real Y0, Real gamma_m, Real gamma_m_hat);

void update_gamma_c_Thomson(Real& gamma_c, InverseComptonY& Ys, RadParams const& rad, Real B, Real t_com, Real gamma_m,
                            Real gamma_c_last);

void update_gamma_c_KN(Real& gamma_c, InverseComptonY& Ys, RadParams const& rad, Real B, Real t_com, Real gamma_m,
                       Real gamma_c_last);

void update_gamma_M(Real& gamma_M, InverseComptonY const& Ys, Real p, Real B);

Real compute_syn_gamma(Real nu, Real B);

template <typename Electrons, typename Photons, typename Updater>
void IC_cooling(ElectronGrid<Electrons>& electrons, PhotonGrid<Photons>& photons, Shock const& shock,
                Updater&& update_gamma_c) {
    const size_t phi_size = electrons.shape()[0];
    const size_t t_size = electrons.shape()[2];

    const size_t phi_compute = (shock.symmetry != Symmetry::structured) ? 1 : phi_size;

    for (size_t i = 0; i < phi_compute; ++i) {
        for (size_t j : shock.theta_reps) {
            const size_t k_inj = shock.injection_idx(i, j);

            for (size_t k = 0; k < t_size; ++k) {
                const Real t_com = shock.t_comv(i, j, k);
                const Real B = shock.B(i, j, k);

                auto& elec = electrons(i, j, k);
                auto& Ys = elec.Ys;
                const Real p = elec.p;
                const Real gamma_c_last = electrons(i, j, k > 0 ? k - 1 : 0).gamma_c;

                update_gamma_c(elec.gamma_c, Ys, shock.rad, B, t_com, elec.gamma_m, gamma_c_last);

                update_gamma_M(elec.gamma_M, Ys, p, B);

                if (k >= k_inj) {
                    const auto& inj = electrons(i, j, k_inj - 1);
                    const Real dt_comv = t_com - shock.t_comv(i, j, k_inj - 1);
                    elec.gamma_c = cool_after_crossing(inj.gamma_c, inj.gamma_m, elec.gamma_m, t_com, B, 0);

                    elec.gamma_M = cool_after_crossing(inj.gamma_M, inj.gamma_m, elec.gamma_m, dt_comv, B, 0);
                }

                const Real I_nu_peak = compute_syn_I_peak(B, p, elec.column_den);
                elec.gamma_a = compute_syn_gamma_a(B, I_nu_peak, elec.gamma_m, elec.gamma_c, elec.gamma_M, p);
                elec.regime = determine_regime(elec.gamma_a, elec.gamma_c, elec.gamma_m);
                elec.Y_c = Ys.gamma_spectrum(elec.gamma_c);
            }
        }
    }
    broadcast_symmetry(electrons, shock);
    generate_syn_photons(photons, shock, electrons);
}

template <typename Electrons, typename Photons>
void Thomson_cooling(ElectronGrid<Electrons>& electrons, PhotonGrid<Photons>& photons, Shock const& shock) {
    IC_cooling(electrons, photons, shock, update_gamma_c_Thomson);
}

template <typename Electrons, typename Photons>
void KN_cooling(ElectronGrid<Electrons>& electrons, PhotonGrid<Photons>& photons, Shock const& shock) {
    IC_cooling(electrons, photons, shock, update_gamma_c_KN);
}

template <typename Photon>
Real inverse_compton_correction(Photon const& ph, Real nu) {
    return (1. + ph.Y_c) / (1 + ph.Ys.nu_spectrum(nu));
}
