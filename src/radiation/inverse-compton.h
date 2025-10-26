//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <array>

#include "../core/mesh.h"
#include "../core/physics.h"
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
     * @param B Magnetic field strength
     * @param Y_T Thomson Y parameter
     * <!-- ************************************************************************************** -->
     */
    InverseComptonY(Real gamma_m, Real gamma_c, Real B, Real Y_T, Real eps_e_on_eps_B) noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Simple constructor that initializes with only the Thomson Y parameter for special cases.
     * @param Y_T Thomson Y parameter
     * <!-- ************************************************************************************** -->
     */
    explicit InverseComptonY(Real Y_T) noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Default constructor that initializes all member variables to zero.
     * <!-- ************************************************************************************** -->
     */
    InverseComptonY() noexcept;

    // Member variables
    Real nu_m_hat{0};    ///< Frequency threshold for minimum electrons
    Real nu_c_hat{0};    ///< Frequency threshold for cooling electrons
    Real nu_self{0};     ///< Frequency threshold for nu_self = nu_self_hat
    Real nu0{0};         ///< Frequency threshold for Y(nu0) = 1
    Real gamma_m_hat{0}; ///< Lorentz factor threshold for minimum energy electrons
    Real gamma_c_hat{0}; ///< Lorentz factor threshold for cooling electrons
    Real gamma_self{0};  ///< Lorentz factor threshold for gamma_self = gamma_self_hat
    Real gamma0{0};      ///< Lorentz factor threshold for Y(gamma0) = 1
    Real Y_T{0};         ///< Thomson scattering Y parameter
    size_t regime{0};    ///< Indicator for the operating regime (1=fast IC cooling, 2=slow IC cooling, 3=special case)

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter for a given frequency and spectral index.
     * @details Different scaling relations apply depending on the cooling regime and frequency range.
     * @param nu Frequency at which to compute the Y parameter
     * @param p Spectral index of electron distribution
     * @return The effective Y parameter at the given frequency
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real evaluate_at_nu(Real nu, Real p) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter for a given Lorentz factor and spectral index.
     * @details Different scaling relations apply depending on the cooling regime and gamma value.
     * @param gamma Electron Lorentz factor
     * @param p Spectral index of electron distribution
     * @return The effective Y parameter at the given gamma
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real evaluate_at_gamma(Real gamma, Real p) const;
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the Compton scattering cross-section as a function of frequency (nu).
 * @param nu The frequency at which to compute the cross-section
 * @return Compton cross-section
 * <!-- ************************************************************************************** -->
 */
Real compton_cross_section(Real nu);

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
        size_t nu_size;
        Real nu_IC_min;
        Real nu_IC_max;
        size_t spectrum_resol;
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
     * @param dgamma Output array for gamma bin widths
     * @param nu_sync Output array for synchrotron frequencies
     * @param dnu_sync Output array for synchrotron frequency bin widths
     * @param dN_e Output array for electron column densities
     * @param I_nu_sync Output array for photon intensities
     * <!-- ************************************************************************************** -->
     */
    void preprocess_distributions(GridParams const& params, Array& gamma, Array& dgamma, Array& nu_sync,
                                  Array& dnu_sync, Array& dN_e, Array& I_nu_sync) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the IC spectrum through integration over electron and photon distributions
     * @param params Grid parameters
     * @param gamma Array of gamma values
     * @param nu_sync Array of synchrotron frequencies
     * @param dN_e Array of electron column densities
     * @param I_nu_sync Array of photon intensities
     * @param dnu_sync Array of synchrotron frequency bin widths
     * <!-- ************************************************************************************** -->
     */
    void compute_IC_spectrum(GridParams const& params, Array const& gamma, Array const& nu_sync, Array const& dN_e,
                             Array const& I_nu_sync, Array const& dnu_sync);

    Array log2_nu_IC;

    Array log2_I_nu_IC;

    Array interp_slope;

    Real inv_dlog2_nu{0};

    static constexpr size_t gamma_grid_per_order{7}; // Number of frequency bins

    static constexpr size_t nu_grid_per_order{5}; // Number of gamma bins

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
 * @return A 3D grid of IC photons
 * <!-- ************************************************************************************** -->
 */
template <typename Electrons, typename Photons>
ICPhotonGrid<Electrons, Photons> generate_IC_photons(ElectronGrid<Electrons> const& electrons,
                                                     PhotonGrid<Photons> const& photons, bool KN = true) noexcept;

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
    const auto params = compute_grid_params();
    initialize_grids(params);

    Array gamma, dgamma, nu_sync, dnu_sync, dN_e, I_nu_sync;
    preprocess_distributions(params, gamma, dgamma, nu_sync, dnu_sync, dN_e, I_nu_sync);

    compute_IC_spectrum(params, gamma, nu_sync, dN_e, I_nu_sync, dnu_sync);

    generated = true;
}

template <typename Electrons, typename Photons>
typename ICPhoton<Electrons, Photons>::GridParams ICPhoton<Electrons, Photons>::compute_grid_params() const {
    GridParams params;

    params.gamma_min = electrons.gamma_m;
    params.gamma_max = electrons.gamma_M * 30;
    params.gamma_size =
        static_cast<size_t>(std::max(std::log10(params.gamma_max / params.gamma_min), 3.) * gamma_grid_per_order);

    params.nu_min = std::min(photons.nu_a, photons.nu_m) / 10;
    params.nu_max = photons.nu_M * 30;
    params.nu_size = static_cast<size_t>(std::max(std::log10(params.nu_max / params.nu_min), 3.) * nu_grid_per_order);

    params.nu_IC_min = 4 * IC_x0 * params.nu_min * params.gamma_min * params.gamma_min;
    params.nu_IC_max = 4 * IC_x0 * params.nu_max * params.gamma_max * params.gamma_max;
    params.spectrum_resol = static_cast<size_t>(std::max(10 * std::log10(params.nu_IC_max / (params.nu_IC_min)), 30.));

    return params;
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::initialize_grids(GridParams const& params) {
    log2_nu_IC = xt::linspace(std::log2(params.nu_IC_min), std::log2(params.nu_IC_max), params.spectrum_resol);

    interp_slope = Array::from_shape({params.spectrum_resol - 1});

    inv_dlog2_nu = 1 / (log2_nu_IC(1) - log2_nu_IC(0));
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::preprocess_distributions(GridParams const& params, Array& gamma, Array& dgamma,
                                                            Array& nu_sync, Array& dnu_sync, Array& dN_e,
                                                            Array& I_nu_sync) const {

    logspace_boundary_center(std::log2(params.nu_min), std::log2(params.nu_max), params.nu_size, nu_sync, dnu_sync);
    logspace_boundary_center(std::log2(params.gamma_min), std::log2(params.gamma_max), params.gamma_size, gamma,
                             dgamma);

    dN_e = Array::from_shape({params.gamma_size});
    for (size_t i = 0; i < params.gamma_size; ++i) {
        dN_e(i) = electrons.compute_column_den(gamma(i)) * dgamma(i);
    }

    I_nu_sync = Array::from_shape({params.nu_size});
    for (size_t j = 0; j < params.nu_size; ++j) {
        I_nu_sync(j) = photons.compute_I_nu(nu_sync(j));
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::compute_IC_spectrum(GridParams const& params, Array const& gamma,
                                                       Array const& nu_sync, Array const& dN_e, Array const& I_nu_sync,
                                                       Array const& dnu_sync) {
    Array nu_IC = xt::exp2(log2_nu_IC);
    Array I_nu_IC = xt::zeros<Real>({params.spectrum_resol});

    //const auto sigma =
    //    KN ? +[](Real nu_comv) { return compton_cross_section(nu_comv); } : +[](Real) { return con::sigmaT; };

    if (KN) { //separate the KN and non-KN case to ensure inline
        for (size_t i = params.gamma_size; i-- > 0;) {
            const Real gamma_i = gamma(i);
            const Real Ndgamma = dN_e(i);
            const Real down_scatter = 1 / (4 * IC_x0 * gamma_i * gamma_i);

            const Real nu_comv_max = gamma_i * nu_sync.back();
            Real row_integral =
                I_nu_sync.back() * dnu_sync.back() * compton_cross_section(nu_comv_max) / (nu_comv_max * nu_comv_max);
            Real slope = 0;
            Real nu_sync_j = 0;

            int k = static_cast<int>(params.spectrum_resol) - 1;
            for (; k >= 0; --k) {
                if (nu_IC(k) * down_scatter < nu_sync.back()) {
                    break;
                }
            }

            for (size_t j = params.nu_size - 1; j-- > 0 && k >= 0;) {
                nu_sync_j = nu_sync(j);
                const Real nu_comv = gamma_i * nu_sync_j;
                const Real inv = 1 / nu_comv;

                slope = I_nu_sync(j) * compton_cross_section(nu_comv) * inv * inv;
                row_integral += slope * dnu_sync(j);

                for (; k >= 0; --k) {
                    if (const Real nu_seed = nu_IC(k) * down_scatter; nu_seed >= nu_sync_j) {
                        I_nu_IC(k) += Ndgamma * (row_integral + slope * (nu_sync_j - nu_seed));
                    } else {
                        break;
                    }
                }
            }

            // Handle the remaining low-frequency tail
            for (; k >= 0; --k) {
                const Real nu_seed = nu_IC(k) * down_scatter;
                I_nu_IC(k) += Ndgamma * (row_integral + slope * (nu_sync_j - nu_seed));
            }
        }
        log2_I_nu_IC = xt::log2(I_nu_IC * nu_IC * 0.25);
    } else {
        Array factor = I_nu_sync / (nu_sync * nu_sync);
        for (size_t i = params.gamma_size; i-- > 0;) {
            const Real gamma_i = gamma(i);
            const Real gamma_i2 = gamma_i * gamma_i;
            const Real Ndgamma = dN_e(i);
            const Real down_scatter = 1 / (4 * IC_x0 * gamma_i2);
            const Real nu_comv_max = gamma_i * nu_sync.back();

            Real row_integral = I_nu_sync.back() * dnu_sync.back() / (nu_comv_max * nu_comv_max);
            Real slope = 0;
            Real nu_sync_j = 0;

            int k = static_cast<int>(params.spectrum_resol) - 1;
            for (; k >= 0; --k) {
                if (nu_IC(k) * down_scatter < nu_sync.back()) {
                    break;
                }
            }

            for (size_t j = params.nu_size - 1; j-- > 0 && k >= 0;) {
                nu_sync_j = nu_sync(j);

                slope = factor(j) / gamma_i2;
                row_integral += slope * dnu_sync(j);

                for (; k >= 0; --k) {
                    if (const Real nu_seed = nu_IC(k) * down_scatter; nu_seed >= nu_sync_j) {
                        I_nu_IC(k) += Ndgamma * (row_integral + slope * (nu_sync_j - nu_seed));
                    } else {
                        break;
                    }
                }
            }

            // Handle the remaining low-frequency tail
            for (; k >= 0; --k) {
                const Real nu_seed = nu_IC(k) * down_scatter;
                I_nu_IC(k) += Ndgamma * (row_integral + slope * (nu_sync_j - nu_seed));
            }
        }
        log2_I_nu_IC = xt::log2(I_nu_IC * nu_IC * 0.25 * con::sigmaT);
    }

    for (size_t i = 0; i < params.spectrum_resol - 1; ++i) {
        interp_slope(i) = (log2_I_nu_IC(i + 1) - log2_I_nu_IC(i)) * inv_dlog2_nu;
    }
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

    size_t idx = 0;
    if (const Real off_set = (log2_nu - log2_nu_IC(0)); off_set >= 0) {
        idx = static_cast<size_t>(std::floor(off_set * inv_dlog2_nu));
    }

    if (idx >= interp_slope.size()) {
        return -con::inf;
    } else [[likely]] {
        return log2_I_nu_IC(idx) + (log2_nu - log2_nu_IC(idx)) * interp_slope(idx);
    }
}

template <typename Electrons, typename Photons>
ICPhotonGrid<Electrons, Photons> generate_IC_photons(ElectronGrid<Electrons> const& electrons,
                                                     PhotonGrid<Photons> const& photons, bool KN) noexcept {
    size_t phi_size = electrons.shape()[0];
    size_t theta_size = electrons.shape()[1];
    size_t t_size = electrons.shape()[2];
    ICPhotonGrid<Electrons, Photons> IC_ph({phi_size, theta_size, t_size});

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                IC_ph(i, j, k) = ICPhoton(electrons(i, j, k), photons(i, j, k), KN);
            }
        }
    }
    return IC_ph;
}

inline Real eta_rad(Real gamma_m, Real gamma_c, Real p) {
    return gamma_c < gamma_m ? 1 : fast_pow(gamma_c / gamma_m, (2 - p));
}

Real compute_gamma_c(Real t_comv, Real B, Real Y);

Real compute_syn_gamma_M(Real B, Real Y, Real p);

Real compute_syn_I_peak(Real B, Real p, Real column_den);

Real compute_syn_gamma_a(Real B, Real I_nu_peak, Real gamma_m, Real gamma_c, Real gamma_M, Real p);

size_t determine_regime(Real gamma_a, Real gamma_c, Real gamma_m);

void update_gamma_c_Thomson(Real& gamma_c, InverseComptonY& Ys, RadParams const& rad, Real B, Real t_com, Real gamma_m);

void update_gamma_c_KN(Real& gamma_c, InverseComptonY& Ys, RadParams const& rad, Real B, Real t_com, Real gamma_m);

void update_gamma_M(Real& gamma_M, InverseComptonY const& Ys, Real p, Real B);

template <typename Electrons, typename Photons, typename Updater>
void IC_cooling(ElectronGrid<Electrons>& electrons, PhotonGrid<Photons>& photons, Shock const& shock,
                Updater&& update_gamma_c) {
    const size_t phi_size = electrons.shape()[0];
    const size_t theta_size = electrons.shape()[1];
    const size_t t_size = electrons.shape()[2];

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            const size_t k_inj = shock.injection_idx(i, j);

            for (size_t k = 0; k < t_size; ++k) {
                if (shock.required(i, j, k) == 0)
                    continue;

                const Real t_com = shock.t_comv(i, j, k);
                const Real B = shock.B(i, j, k);

                auto& elec = electrons(i, j, k);
                auto& Ys = elec.Ys;
                const Real p = elec.p;

                update_gamma_c(elec.gamma_c, Ys, shock.rad, B, t_com, elec.gamma_m);

                update_gamma_M(elec.gamma_M, Ys, p, B);

                if (k >= k_inj) {
                    const auto& inj = electrons(i, j, k_inj);
                    elec.gamma_c = inj.gamma_c * elec.gamma_m / inj.gamma_m;
                    elec.gamma_M = elec.gamma_c;
                }

                const Real I_nu_peak = compute_syn_I_peak(B, p, elec.column_den);
                elec.gamma_a = compute_syn_gamma_a(B, I_nu_peak, elec.gamma_m, elec.gamma_c, elec.gamma_M, p);
                elec.regime = determine_regime(elec.gamma_a, elec.gamma_c, elec.gamma_m);
                elec.Y_c = Ys.evaluate_at_gamma(elec.gamma_c, p);
            }
        }
    }
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
