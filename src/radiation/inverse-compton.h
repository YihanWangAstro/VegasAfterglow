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
 * @param conical flag for conical solution
 * @return A 3D grid of IC photons
 * <!-- ************************************************************************************** -->
 */
template <typename Electrons, typename Photons>
ICPhotonGrid<Electrons, Photons> generate_IC_photons(ElectronGrid<Electrons> const& electrons,
                                                     PhotonGrid<Photons> const& photons, bool KN = true,
                                                     Symmetry symmetry = Symmetry::structured) noexcept;

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
    // Ensure minimum spectrum resolution to avoid division by zero and out-of-bounds access
    const size_t safe_resol = std::max(params.spectrum_resol, static_cast<size_t>(2));

    log2_nu_IC = xt::linspace(std::log2(params.nu_IC_min), std::log2(params.nu_IC_max), safe_resol);

    interp_slope = Array::from_shape({safe_resol - 1});

    const Real dlog2 = log2_nu_IC(1) - log2_nu_IC(0);
    inv_dlog2_nu = (dlog2 != 0) ? (1 / dlog2) : 0;
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
                                                     PhotonGrid<Photons> const& photons, bool KN,
                                                     Symmetry symmetry) noexcept {
    size_t phi_size = electrons.shape()[0];
    size_t theta_size = electrons.shape()[1];
    size_t t_size = electrons.shape()[2];
    ICPhotonGrid<Electrons, Photons> IC_ph({phi_size, theta_size, t_size});

    const size_t phi_compute = (symmetry != Symmetry::structured) ? 1 : phi_size;
    const size_t theta_compute = (symmetry == Symmetry::isotropic) ? 1 : theta_size;

    for (size_t i = 0; i < phi_compute; ++i) {
        for (size_t j = 0; j < theta_compute; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                IC_ph(i, j, k) = ICPhoton(electrons(i, j, k), photons(i, j, k), KN);
            }
        }
    }

    if (symmetry == Symmetry::isotropic) {
        for (size_t i = 0; i < phi_size; ++i)
            for (size_t j = 0; j < theta_size; ++j)
                if (i != 0 || j != 0)
                    xt::view(IC_ph, i, j, xt::all()) = xt::view(IC_ph, 0, 0, xt::all());
    } else if (symmetry >= Symmetry::phi_symmetric) {
        // piecewise: input data already group-identical from synchrotron broadcast, just broadcast phi
        for (size_t i = 1; i < phi_size; ++i)
            for (size_t j = 0; j < theta_size; ++j)
                xt::view(IC_ph, i, j, xt::all()) = xt::view(IC_ph, 0, j, xt::all());
    }

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
    const size_t theta_size = electrons.shape()[1];
    const size_t t_size = electrons.shape()[2];

    const size_t phi_compute = (shock.symmetry != Symmetry::structured) ? 1 : phi_size;

    for (size_t i = 0; i < phi_compute; ++i) {
        size_t rep_idx = 0;
        for (size_t j = 0; j < theta_size; ++j) {
            if (shock.symmetry >= Symmetry::piecewise) {
                if (rep_idx < shock.theta_reps.size() && shock.theta_reps[rep_idx] == j)
                    ++rep_idx;
                else
                    continue;
            }
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
                    elec.gamma_c = cool_after_crossing(inj.gamma_c, inj.gamma_m, elec.gamma_m, dt_comv, B, 0);
                    // gamma_M: step-by-step synchrotron cooling from previous step (continuous, uses local B)
                    const auto& prev = electrons(i, j, k - 1);
                    const Real dt_step = t_com - shock.t_comv(i, j, k - 1);
                    const Real gamma_c_step = compute_gamma_c(dt_step, B, 0);
                    elec.gamma_M = prev.gamma_M * gamma_c_step / (prev.gamma_M + gamma_c_step);
                }

                const Real I_nu_peak = compute_syn_I_peak(B, p, elec.column_den);
                elec.gamma_a = compute_syn_gamma_a(B, I_nu_peak, elec.gamma_m, elec.gamma_c, elec.gamma_M, p);
                elec.regime = determine_regime(elec.gamma_a, elec.gamma_c, elec.gamma_m);
                elec.Y_c = Ys.gamma_spectrum(elec.gamma_c);
            }
        }
    }
    if (shock.symmetry == Symmetry::isotropic) {
        for (size_t i = 0; i < phi_size; ++i)
            for (size_t j = 0; j < theta_size; ++j)
                if (i != 0 || j != 0)
                    xt::view(electrons, i, j, xt::all()) = xt::view(electrons, 0, 0, xt::all());
    } else if (shock.symmetry == Symmetry::piecewise) {
        for (size_t r = 0; r < shock.theta_reps.size(); ++r) {
            const size_t j_start = shock.theta_reps[r];
            const size_t j_end = (r + 1 < shock.theta_reps.size()) ? shock.theta_reps[r + 1] : theta_size;
            for (size_t j = j_start + 1; j < j_end; ++j)
                xt::view(electrons, 0, j, xt::all()) = xt::view(electrons, 0, j_start, xt::all());
        }
        for (size_t i = 1; i < phi_size; ++i)
            for (size_t j = 0; j < theta_size; ++j)
                xt::view(electrons, i, j, xt::all()) = xt::view(electrons, 0, j, xt::all());
    } else if (shock.symmetry == Symmetry::phi_symmetric) {
        for (size_t i = 1; i < phi_size; ++i)
            for (size_t j = 0; j < theta_size; ++j)
                xt::view(electrons, i, j, xt::all()) = xt::view(electrons, 0, j, xt::all());
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

template <typename Photon>
Real inverse_compton_correction(Photon const& ph, Real nu) {
    if (ph.Y_c < 1e-3) {
        return 1.;
    } else {
        return (1. + ph.Y_c) / (1 + ph.Ys.nu_spectrum(nu));
    }
}
