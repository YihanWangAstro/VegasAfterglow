//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <algorithm>
#include <array>
#include <limits>

#include "../core/mesh.h"
#include "../core/physics.h"
#include "../core/quadrature.h"
#include "../dynamics/shock.h"
#include "../util/macros.h"
#include "../util/utilities.h"

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
    Real gamma0{1};      ///< Lorentz factor threshold for Y(gamma0) = 1x
    Real Y_T{0};         ///< Thomson scattering Y parameter
    size_t regime{0};    ///< Indicator for the operating regime (1=fast IC cooling, 2=slow IC cooling, 3=special case)

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter for a given electron Lorentz factor.
     * @details Different scaling relations apply depending on the cooling regime and gamma value.
     * @param gamma Electron Lorentz factor
     * @return The effective Y parameter at the given gamma
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real gamma_spectrum(Real gamma) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter at the electron Lorentz factor corresponding to a given synchrotron frequency.
     * @details Converts frequency to electron Lorentz factor via the synchrotron relation, then evaluates Y(gamma).
     * @param nu Synchrotron frequency
     * @return The effective Y parameter at the corresponding electron Lorentz factor
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

    BrokenPowerLaw<5> segments_;

    void build_segments() noexcept;
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast KN correction lookup: returns sigma_KN/sigmaT at frequency nu.
 * @param nu Photon frequency.
 * @return Dimensionless correction factor to Thomson cross-section.
 * <!-- ************************************************************************************** -->
 */
Real compton_correction(Real nu);

/**
 * <!-- ************************************************************************************** -->
 * @struct ICPhoton
 * @tparam Electrons Type of the electron distribution
 * @tparam Photons Type of the photon distribution
 * @brief Represents a single inverse Compton (IC) photon.
 * @details Contains methods to compute the photon intensity I_nu and to generate an IC photon spectrum based
 *          on electron and seed photon properties.
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
     * <!-- ************************************************************************************** -->
     * @brief Computes weighted break scores for nu_seed, gamma, and nu_ic adaptive grids.
     * <!-- ************************************************************************************** -->
     */
    void compute_break_weight(GridParams& params) const;

    void initialize_grids(GridParams const& params, Array& gamma, Array& dgamma, Array& nu_seed, Array& dnu_seed,
                          Array& nu_IC);

    void sample_distributions(Array const& gamma, Array const& dgamma, Array const& nu_seed, Array& dN_e_boost,
                              Array& I_nu_seed);

    void compute_IC_spectrum(Array const& gamma, Array const& dN_e_boost, Array const& nu_seed, Array const& dnu_seed,
                             Array const& nu_IC, Array const& I_nu_seed);

    static void build_cdf_thomson(Array const& I_nu_seed, Array const& nu_seed, Array const& dnu_seed, Array& fv_buf,
                                  Array& cdf_buf);

    static void build_cdf_KN(Real gamma_i, Array const& I_nu_seed, Array const& nu_seed, Array const& dnu_seed,
                             Array& fv_buf, Array& cdf_buf);

    static void accumulate_IC(Real dN_e_boost, Real gamma_i2, Array const& nu_IC, Array const& nu_seed,
                              Array const& dnu_seed, Array const& fv_buf, Array const& cdf_buf, Array& I_buf,
                              size_t spec_size);

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
inline constexpr Real inv_4x0 = 1 / (4 * IC_x0);
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
    Array gamma, dgamma, nu_seed, dnu_seed, nu_IC;
    initialize_grids(params, gamma, dgamma, nu_seed, dnu_seed, nu_IC);

    Array dN_e_boost, I_nu_seed;
    sample_distributions(gamma, dgamma, nu_seed, dN_e_boost, I_nu_seed);
    compute_IC_spectrum(gamma, dN_e_boost, nu_seed, dnu_seed, nu_IC, I_nu_seed);

    generated = true;
}

template <typename Electrons, typename Photons>
typename ICPhoton<Electrons, Photons>::GridParams ICPhoton<Electrons, Photons>::compute_grid_params() const {
    GridParams params;

    constexpr Real eps_tail = 1e-2;     // target suppression level for single-cutoff tail
    constexpr Real safety_margin = 2.0; // retain headroom beyond estimated tail
    constexpr Real min_tail_factor = 5.0;

    const Real tail_factor = std::max(-std::log(eps_tail), min_tail_factor);

    params.gamma_min = std::min(electrons.gamma_m, electrons.gamma_c) / 3;
    params.gamma_max = std::max(electrons.gamma_M * tail_factor, params.gamma_min);
    params.gamma_breaks = {electrons.gamma_a, electrons.gamma_m, electrons.gamma_c, electrons.gamma_M};

    params.nu_min = std::min(photons.nu_a, photons.nu_m) / 10;
    params.nu_max = std::max(photons.nu_M * tail_factor, params.nu_min);

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

    std::array<Real, 4> dN_e_break{};
    for (size_t j = 0; j < gamma_breaks.size(); ++j) {
        dN_e_break[j] = electrons.compute_column_den(gamma_breaks[j]);
        params.gamma_break_weight[j] = dN_e_break[j] / (gamma_breaks[j] * gamma_breaks[j]);
    }

    size_t k = 0;
    for (size_t i = 0; i < nu_seed_breaks.size(); ++i) {
        const Real nu_seed = nu_seed_breaks[i];
        const Real I_seed = photons.compute_I_nu(nu_seed);
        params.nu_break_weight[i] = I_seed / (nu_seed * nu_seed);

        for (size_t j = 0; j < gamma_breaks.size(); ++j) {
            const Real gamma_b = gamma_breaks[j];
            const Real nu_IC = 2 * IC_x0 * nu_seed * gamma_b * gamma_b;
            const Real dN_e = dN_e_break[j];
            const Real KN_correction = KN ? compton_correction(gamma_b * nu_seed) : Real(1);

            params.nu_IC_break_weight[k] = dN_e * I_seed * KN_correction;
            params.nu_IC_breaks[k] = nu_IC;
            ++k;
        }
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::initialize_grids(GridParams const& params, Array& gamma, Array& dgamma,
                                                    Array& nu_seed, Array& dnu_seed, Array& nu_IC) {
    constexpr Real gamma_grid_per_order{5};
    constexpr Real nu_grid_per_order{3.5};
    constexpr Real ic_grid_per_order{2};

    adaptive_grid_with_breaks(std::log2(params.nu_IC_min), std::log2(params.nu_IC_max),
                              std::span<const Real>(params.nu_IC_breaks),
                              std::span<const Real>(params.nu_IC_break_weight), ic_grid_per_order, nu_IC, 1, 2, 3);
    log2_nu_IC = xt::log2(nu_IC);

    adaptive_grid_with_breaks(std::log2(params.nu_min), std::log2(params.nu_max),
                              std::span<const Real>(params.nu_breaks), std::span<const Real>(params.nu_break_weight),
                              nu_grid_per_order, nu_seed, 1, 0.5, 3);
    compute_bin_widths(nu_seed, dnu_seed);

    log2space(std::log2(params.gamma_min), std::log2(params.gamma_max), gamma_grid_per_order, gamma);
    compute_trapezoidal_weights(gamma, dgamma);
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::sample_distributions(Array const& gamma, Array const& dgamma, Array const& nu_seed,
                                                        Array& dN_e_boost, Array& I_nu_seed) {
    dN_e_boost = Array::from_shape({gamma.size()});
    const size_t g_size = gamma.size();
    for (size_t i = 0; i < g_size; ++i) {
        const Real gamma_i = gamma(i);
        dN_e_boost(i) = electrons.compute_column_den(gamma_i) * dgamma(i) / (gamma_i * gamma_i);
    }

    I_nu_seed = Array::from_shape({nu_seed.size()});
    const size_t nu_size = nu_seed.size();
    for (size_t j = 0; j < nu_size; ++j) {
        I_nu_seed(j) = photons.compute_I_nu(nu_seed(j));
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::build_cdf_thomson(Array const& I_nu_seed, Array const& nu_seed,
                                                     Array const& dnu_seed, Array& fv_buf, Array& cdf_buf) {
    const int nu_last = static_cast<int>(nu_seed.size()) - 1;

    fv_buf(nu_last) = I_nu_seed(nu_last) / (nu_seed(nu_last) * nu_seed(nu_last));
    cdf_buf(nu_last) = 0;
    for (int j = nu_last - 1; j >= 0; --j) {
        fv_buf(j) = I_nu_seed(j) / (nu_seed(j) * nu_seed(j));
        cdf_buf(j) = cdf_buf(j + 1) + 0.5 * (fv_buf(j) + fv_buf(j + 1)) * dnu_seed(j);
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::build_cdf_KN(Real gamma_i, Array const& I_nu_seed, Array const& nu_seed,
                                                Array const& dnu_seed, Array& fv_buf, Array& cdf_buf) {
    const int nu_last = static_cast<int>(nu_seed.size()) - 1;

    fv_buf(nu_last) =
        I_nu_seed(nu_last) / (nu_seed(nu_last) * nu_seed(nu_last)) * compton_correction(gamma_i * nu_seed(nu_last));
    cdf_buf(nu_last) = 0;
    for (int j = nu_last - 1; j >= 0; --j) {
        fv_buf(j) = I_nu_seed(j) / (nu_seed(j) * nu_seed(j)) * compton_correction(gamma_i * nu_seed(j));
        cdf_buf(j) = cdf_buf(j + 1) + 0.5 * (fv_buf(j) + fv_buf(j + 1)) * dnu_seed(j);
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::accumulate_IC(Real dN_e_boost, Real gamma_i2, Array const& nu_IC,
                                                 Array const& nu_seed, Array const& dnu_seed, Array const& fv_buf,
                                                 Array const& cdf_buf, Array& I_buf, size_t spec_size) {
    const int nu_last = static_cast<int>(nu_seed.size()) - 1;
    const Real inv_gi2 = inv_4x0 / gamma_i2;

    if (cdf_buf(0) <= 0)
        return;

    // Plateau: nu_seed < nu_seed_min, all seed photons overshoot, CDF = cdf_buf[0]
    const Real plateau = dN_e_boost * cdf_buf(0);
    size_t k = 0;
    while (k < spec_size && nu_IC(k) * inv_gi2 < nu_seed(0)) {
        I_buf(k++) += plateau;
    }

    for (int j = 0; j < nu_last && k < spec_size; ++j) {
        const Real f_lo = fv_buf(j);
        const Real f_hi = fv_buf(j + 1);
        const Real cdf_hi = cdf_buf(j + 1);
        const Real dnu = dnu_seed(j);
        const Real nu_lo = nu_seed(j);
        const Real inv_dnu = 1.0 / dnu;

        size_t k_hi = k;
        while (k_hi < spec_size && nu_IC(k_hi) * inv_gi2 < nu_seed(j + 1)) {
            ++k_hi;
        }

        for (size_t kk = k; kk < k_hi; ++kk) {
            const Real frac = (nu_IC(kk) * inv_gi2 - nu_lo) * inv_dnu;
            const Real rem = 1.0 - frac;
            const Real f_seed = f_lo * rem + f_hi * frac;
            I_buf(kk) += dN_e_boost * (cdf_hi + 0.5 * (f_seed + f_hi) * rem * dnu);
        }

        k = k_hi;
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::compute_IC_spectrum(Array const& gamma, Array const& dN_e_boost,
                                                       Array const& nu_seed, Array const& dnu_seed, Array const& nu_IC,
                                                       Array const& I_nu_seed) {
    const size_t spec_size = log2_nu_IC.size();
    const size_t nu_size = nu_seed.size();
    const int gamma_size = static_cast<int>(std::min(gamma.size(), dN_e_boost.size()));

    log2_I_nu_IC = Array::from_shape({spec_size});

    static thread_local Array I_buf, cdf_buf, fv_buf;
    I_buf.resize({spec_size});
    cdf_buf.resize({nu_size});
    fv_buf.resize({nu_size});
    I_buf.fill(0);

    if (KN) {
        for (int i = 0; i < gamma_size; ++i) {
            if (dN_e_boost(i) <= 0)
                continue;
            const Real gamma_i = gamma(i);
            build_cdf_KN(gamma_i, I_nu_seed, nu_seed, dnu_seed, fv_buf, cdf_buf);
            accumulate_IC(dN_e_boost(i), gamma_i * gamma_i, nu_IC, nu_seed, dnu_seed, fv_buf, cdf_buf, I_buf,
                          spec_size);
        }
    } else {
        build_cdf_thomson(I_nu_seed, nu_seed, dnu_seed, fv_buf, cdf_buf);
        for (int i = 0; i < gamma_size; ++i) {
            if (dN_e_boost(i) <= 0)
                continue;
            accumulate_IC(dN_e_boost(i), gamma(i) * gamma(i), nu_IC, nu_seed, dnu_seed, fv_buf, cdf_buf, I_buf,
                          spec_size);
        }
    }

    const Real log2_scale = fast_log2(0.25 * con::sigmaT);
    interp_slope = Array::from_shape({spec_size - 1});
    for (size_t i = 0; i < spec_size; ++i) {
        log2_I_nu_IC(i) = fast_log2(I_buf[i]) + log2_nu_IC(i) + log2_scale;
        if (i > 0) {
            const Real dl = log2_nu_IC(i) - log2_nu_IC(i - 1);
            interp_slope(i - 1) = (dl != 0) ? (log2_I_nu_IC(i) - log2_I_nu_IC(i - 1)) / dl : 0;
        }
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

    for (size_t i = 0; i < phi_compute; ++i) {
        for (size_t j : shock.theta_reps) {
            for (size_t k = 0; k < t_size; ++k) {
                IC_ph(i, j, k) = ICPhoton(electrons(i, j, k), photons(i, j, k), KN);
            }
        }
    }

    if (shock.symmetry != Symmetry::structured) {
        for (size_t i = 0; i < phi_compute; ++i) {
            for (size_t j : shock.theta_reps) {
                for (size_t k = 0; k < t_size; ++k)
                    IC_ph(i, j, k).compute_log2_I_nu(0);
            }
        }
    }

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
                    auto const& inj = electrons(i, j, k_inj - 1);
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
