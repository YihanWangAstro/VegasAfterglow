//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <algorithm>
#include <array>
#include <cstdio>
#include <limits>

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

        Real nu_min;
        Real nu_max;
        Real nu_IC_min;
        Real nu_IC_max;
        std::array<Real, 4> gamma_breaks{};
        std::array<Real, 4> nu_breaks;
        std::array<Real, 4> nu_break_weight{};
        std::array<Real, 4> gamma_break_weight{};
        std::array<Real, 16> nu_IC_breaks{};
        std::array<Real, 16> nu_IC_break_weight{};
    };

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes grid parameters based on electron and photon distributions
     * @return GridParams structure containing all computed parameters
     * <!-- ************************************************************************************** -->
     */
    GridParams compute_grid_params() const;

    /**
     * @brief Computes weighted break scores for nu_sync, gamma, and nu_ic adaptive grids.
     */
    void compute_break_weight(GridParams& params) const;

    void initialize_grids(GridParams const& params, Array& gamma, Array& dgamma, Array& nu_sync, Array& dnu_sync,
                          Array& nu_IC);

    void sample_distributions(Array const& gamma, Array const& dgamma, Array const& nu_sync, Array& dN_e_boost,
                              Array& I_nu_sync);

    void compute_IC_spectrum(Array const& gamma, Array const& dN_e_boost, Array const& nu_sync, Array const& dnu_sync,
                             Array const& nu_IC, Array const& I_nu_sync);

    Array log2_nu_IC;

    Array log2_I_nu_IC;

    Array interp_slope;

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
    Array gamma, dgamma, nu_sync, dnu_sync, nu_IC;
    initialize_grids(params, gamma, dgamma, nu_sync, dnu_sync, nu_IC);

    Array dN_e_boost, I_nu_sync;
    sample_distributions(gamma, dgamma, nu_sync, dN_e_boost, I_nu_sync);
    compute_IC_spectrum(gamma, dN_e_boost, nu_sync, dnu_sync, nu_IC, I_nu_sync);

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
    params.gamma_breaks = {electrons.gamma_a, electrons.gamma_m, electrons.gamma_c, electrons.gamma_M};

    params.nu_min = std::min(photons.nu_a, photons.nu_m) / 10;
    params.nu_max = std::max(photons.nu_M * tail_factor, params.nu_min * 1.01);

    params.nu_breaks = {photons.nu_a, photons.nu_m, photons.nu_c, photons.nu_M};

    params.nu_IC_min = 4 * IC_x0 * params.nu_min * params.gamma_min * params.gamma_min;

    const Real nu_ic_base = 4 * IC_x0 * photons.nu_M * electrons.gamma_M * electrons.gamma_M;
    const Real nu_ic_single_cut = std::max(nu_ic_base * tail_factor * tail_factor, nu_ic_base * tail_factor);
    params.nu_IC_max = nu_ic_single_cut * safety_margin;

    compute_break_weight(params);

    return params;
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::compute_break_weight(GridParams& params) const {
    params.nu_break_weight.fill(0);
    params.gamma_break_weight.fill(0);
    params.nu_IC_break_weight.fill(0);
    params.nu_IC_breaks.fill(std::numeric_limits<Real>::quiet_NaN());

    const std::array<Real, 4> nu_seed_breaks = params.nu_breaks;
    const std::array<Real, 4> gamma_breaks = params.gamma_breaks;

    static_assert(std::tuple_size_v<decltype(GridParams::nu_IC_breaks)> == nu_seed_breaks.size() * gamma_breaks.size());

    size_t k = 0;
    for (size_t i_seed_idx = 0; i_seed_idx < nu_seed_breaks.size(); ++i_seed_idx) {
        const Real nu_seed = nu_seed_breaks[i_seed_idx];
        const Real I_seed = photons.compute_I_nu(nu_seed);
        params.nu_break_weight[i_seed_idx] = I_seed / (nu_seed * nu_seed);

        for (size_t j_gamma_idx = 0; j_gamma_idx < gamma_breaks.size(); ++j_gamma_idx) {
            const Real gamma_b = gamma_breaks[j_gamma_idx];
            const Real nu_IC = 2 * IC_x0 * nu_seed * gamma_b * gamma_b;

            const Real dN_e = electrons.compute_column_den(gamma_b);
            const Real KN_correction = KN ? compton_correction(gamma_b * nu_seed) : Real(1);

            params.nu_IC_break_weight[k] = dN_e * I_seed * KN_correction;
            params.nu_IC_breaks[k] = nu_IC;
            ++k;

            if (i_seed_idx == 0) {
                params.gamma_break_weight[j_gamma_idx] = dN_e / (gamma_b * gamma_b);
            }
        }
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::initialize_grids(GridParams const& params, Array& gamma, Array& dgamma,
                                                    Array& nu_sync, Array& dnu_sync, Array& nu_IC) {
    constexpr Real gamma_grid_per_order{15};
    constexpr Real nu_grid_per_order{3};
    constexpr Real ic_grid_per_order{2};
    // IC output frequency grid
    adaptive_grid_with_breaks(std::log2(params.nu_IC_min), std::log2(params.nu_IC_max),
                              std::span<const Real>(params.nu_IC_breaks),
                              std::span<const Real>(params.nu_IC_break_weight), ic_grid_per_order, nu_IC, 1, 2, 3);
    log2_nu_IC = xt::log2(nu_IC);

    // Synchrotron seed photon frequency grid
    adaptive_grid_with_breaks(std::log2(params.nu_min), std::log2(params.nu_max),
                              std::span<const Real>(params.nu_breaks), std::span<const Real>(params.nu_break_weight),
                              nu_grid_per_order, nu_sync, 1, 0.5, 3);
    compute_bin_widths(nu_sync, dnu_sync);

    // Electron energy grid
    Array gamma_boundary;
    adaptive_grid_with_breaks(
        std::log2(params.gamma_min), std::log2(params.gamma_max), std::span<const Real>(params.gamma_breaks),
        std::span<const Real>(params.gamma_break_weight), gamma_grid_per_order, gamma_boundary, 1, 1, 3);
    compute_bin_widths(gamma_boundary, dgamma);

    gamma = Array::from_shape({gamma_boundary.size() - 1});
    boundary_to_center_log(gamma_boundary, gamma);
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::sample_distributions(Array const& gamma, Array const& dgamma, Array const& nu_sync,
                                                        Array& dN_e_boost, Array& I_nu_sync) {
    dN_e_boost = Array::from_shape({gamma.size()});
    for (size_t i = 0; i < gamma.size(); ++i) {
        dN_e_boost(i) = electrons.compute_column_den(gamma(i)) * dgamma(i) / (gamma(i) * gamma(i));
    }

    I_nu_sync = Array::from_shape({nu_sync.size()});
    for (size_t j = 0; j < nu_sync.size(); ++j) {
        I_nu_sync(j) = photons.compute_I_nu(nu_sync(j));
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::compute_IC_spectrum(Array const& gamma, Array const& dN_e, Array const& nu_sync,
                                                       Array const& dnu_sync, Array const& nu_IC,
                                                       Array const& I_nu_sync) {
    const size_t spec_size = log2_nu_IC.size();
    log2_I_nu_IC = Array::from_shape({spec_size});

    const size_t nu_size = nu_sync.size();
    const int g_size = static_cast<int>(std::min(gamma.size(), dN_e.size()));
    Array I_nu_IC = xt::zeros<Real>({spec_size});
    Array cdf = Array::from_shape({nu_size});
    Array f_vals = Array::from_shape({nu_size});
    Array inv_dnu = Array::from_shape({nu_size - 1});
    Array inv_nu2 = Array::from_shape({nu_size});

    const int nu_last = static_cast<int>(nu_size) - 1;
    const int spec_last = static_cast<int>(spec_size) - 1;
    const Real inv_4x0 = 1 / (4 * IC_x0);

    // Raw pointers for hot-path array access
    const Real* nu_p = nu_sync.data();
    const Real* dnu_p = dnu_sync.data();
    const Real* nu_ic_p = nu_IC.data();
    const Real* Is_p = I_nu_sync.data();
    const Real* g_p = gamma.data();
    const Real* dNe_p = dN_e.data();
    Real* cdf_p = cdf.data();
    Real* fv_p = f_vals.data();
    Real* I_p = I_nu_IC.data();
    Real* inv_dnu_p = inv_dnu.data();
    Real* inv_nu2_p = inv_nu2.data();

    for (size_t j = 0; j + 1 < nu_size; ++j) {
        inv_dnu_p[j] = 1.0 / dnu_p[j];
    }
    for (size_t j = 0; j < nu_size; ++j) {
        inv_nu2_p[j] = 1.0 / (nu_p[j] * nu_p[j]);
    }

    const Real nu_sync_max = nu_p[nu_last];

    auto accumulate = [&](Real weight, Real inv_4x0_gi2, Real gi2) {
        const Real limit = 4 * IC_x0 * gi2 * nu_sync_max;
        int scan = 0;
        for (int k = 0; k <= spec_last; ++k) {
            if (nu_ic_p[k] > limit) {
                break;
            }

            const Real nu_seed = nu_ic_p[k] * inv_4x0_gi2;
            while (scan < nu_last && nu_p[scan + 1] <= nu_seed) {
                ++scan;
            }
            if (scan >= nu_last) {
                break;
            }

            const Real frac = std::max((nu_seed - nu_p[scan]) * inv_dnu_p[scan], Real(0));
            const Real rem = 1.0 - frac;
            const Real f_seed = fv_p[scan] * rem + fv_p[scan + 1] * frac;
            const Real cdf_val =
                std::max(cdf_p[scan + 1] + 0.5 * (f_seed + fv_p[scan + 1]) * rem * dnu_p[scan], Real(0));
            I_p[k] += weight * cdf_val;
        }
    };

    auto build_cdf = [&](Real gamma_i) {
        fv_p[nu_last] =
            Is_p[nu_last] * (KN ? compton_correction(gamma_i * nu_p[nu_last]) : Real(1)) * inv_nu2_p[nu_last];
        cdf_p[nu_last] = 0;
        for (int j = nu_last - 1; j >= 0; --j) {
            fv_p[j] = Is_p[j] * (KN ? compton_correction(gamma_i * nu_p[j]) : Real(1)) * inv_nu2_p[j];
            cdf_p[j] = cdf_p[j + 1] + 0.5 * (fv_p[j] + fv_p[j + 1]) * dnu_p[j];
        }
    };

    // Thomson: CDF is gamma-independent, build once
    if (!KN)
        build_cdf(0);

    for (int i = 0; i < g_size; ++i) {
        const Real weight = dNe_p[i];
        if (weight <= 0) {
            continue;
        }
        const Real gamma_i = g_p[i];
        const Real gamma_i2 = gamma_i * gamma_i;
        if (KN) {
            build_cdf(gamma_i);
        }
        accumulate(weight, inv_4x0 / gamma_i2, gamma_i2);
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
}

template <typename Electrons, typename Photons>
Real ICPhoton<Electrons, Photons>::compute_I_nu(Real nu) {
    return std::exp2(compute_log2_I_nu(std::log2(nu)));
}

template <typename Electrons, typename Photons>
Real ICPhoton<Electrons, Photons>::compute_log2_I_nu(Real log2_nu) {
    if (!generated) {
        generate_spectrum();
    }

    const size_t n = log2_nu_IC.size();
    if (n < 2 || log2_nu > log2_nu_IC(n - 1)) {
        return -con::inf;
    }

    size_t idx = 0;
    while (idx + 2 < n && log2_nu_IC(idx + 1) <= log2_nu) {
        ++idx;
    }

    return log2_I_nu_IC(idx) + (log2_nu - log2_nu_IC(idx)) * interp_slope(idx);
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
