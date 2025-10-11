//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <array>

#include "macros.h"
#include "mesh.h"
#include "shock.h"
#include "utilities.h"
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
     * @param nu_m Characteristic frequency for the minimum Lorentz factor
     * @param nu_c Characteristic frequency for the cooling Lorentz factor
     * @param B Magnetic field strength
     * @param Y_T Thomson Y parameter
     * <!-- ************************************************************************************** -->
     */
    InverseComptonY(Real nu_m, Real nu_c, Real B, Real Y_T) noexcept;

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
    Real nu_hat_m{0};    ///< Frequency threshold for minimum electrons
    Real nu_hat_c{0};    ///< Frequency threshold for cooling electrons
    Real gamma_hat_m{0}; ///< Lorentz factor threshold for minimum energy electrons
    Real gamma_hat_c{0}; ///< Lorentz factor threshold for cooling electrons
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
    [[nodiscard]] Real compute_val_at_nu(Real nu, Real p) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter for a given Lorentz factor and spectral index.
     * @details Different scaling relations apply depending on the cooling regime and gamma value.
     * @param gamma Electron Lorentz factor
     * @param p Spectral index of electron distribution
     * @return The effective Y parameter at the given gamma
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real compute_val_at_gamma(Real gamma, Real p) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns the Thomson Y parameter from the provided InverseComptonY object.
     * @details Previously supported summing Y parameters from multiple objects.
     * @param Ys InverseComptonY object
     * @return The Thomson Y parameter
     * <!-- ************************************************************************************** -->
     */
    static Real compute_Y_Thompson(InverseComptonY const& Ys); ///< Returns Y_T parameter

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter at a specific Lorentz factor and spectral index.
     * @details Previously supported summing contributions from multiple InverseComptonY objects.
     * @param Ys InverseComptonY object
     * @param gamma Electron Lorentz factor
     * @param p Spectral index of electron distribution
     * @return The effective Y parameter at the given gamma
     * <!-- ************************************************************************************** -->
     */
    static Real compute_Y_tilt_at_gamma(InverseComptonY const& Ys, Real gamma, Real p);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter at a specific frequency and spectral index.
     * @details Previously supported summing contributions from multiple InverseComptonY objects.
     * @param Ys InverseComptonY object
     * @param nu Frequency at which to compute the Y parameter
     * @param p Spectral index of electron distribution
     * @return The effective Y parameter at the given frequency
     * <!-- ************************************************************************************** -->
     */
    static Real compute_Y_tilt_at_nu(InverseComptonY const& Ys, Real nu, Real p);
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
    void generate_grid();

    Array nu0; // input frequency nu0 grid value

    Array I_nu_dnu; // input I_nu

    Array gamma; // gamma grid boundary values

    Array dN_e; // electron column density

    MeshGrid IC_tab;

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
void ICPhoton<Electrons, Photons>::generate_grid() {
    const Real gamma_min = std::min(electrons.gamma_m, electrons.gamma_c);
    const Real gamma_max = electrons.gamma_M * 10;
    const size_t gamma_size = static_cast<size_t>(std::log10(gamma_max / gamma_min) * gamma_grid_per_order);

    const Real nu_min = std::min(photons.nu_a, photons.nu_m) / 10;
    const Real nu_max = photons.nu_M * 10;
    const size_t nu_size = static_cast<size_t>(std::log10(nu_max / nu_min) * nu_grid_per_order);

    Array dnu0, dgamma;

    logspace_boundary_center(std::log2(nu_min), std::log2(nu_max), nu_size, nu0, dnu0);

    logspace_boundary_center(std::log2(gamma_min), std::log2(gamma_max), gamma_size, gamma, dgamma);

    IC_tab = MeshGrid({gamma_size, nu_size}, -1);
    generated = true;

    I_nu_dnu = Array::from_shape({nu_size});
    dN_e = Array::from_shape({gamma_size});

    for (size_t i = 0; i < gamma_size; ++i) {
        dN_e(i) = electrons.compute_column_den(gamma(i)) * dgamma(i);
    }

    for (size_t j = 0; j < nu_size; ++j) {
        I_nu_dnu(j) = photons.compute_I_nu(nu0(j)) * dnu0(j);
    }
}

template <typename Electrons, typename Photons>
Real ICPhoton<Electrons, Photons>::compute_I_nu(Real nu) {
    if (generated == false)
        generate_grid();

    Real IC_I_nu = 0;

    const size_t gamma_size = gamma.size();
    const size_t nu_size = nu0.size();

    const auto sigma =
        KN ? +[](Real nu_comv) { return compton_cross_section(nu_comv); } : +[](Real) { return con::sigmaT; };

    for (size_t i = gamma_size; i-- > 0;) {
        const Real gamma_i = gamma(i);
        const Real upscatter = 4 * IC_x0 * gamma_i * gamma_i;
        const Real Ndgamma = dN_e(i);

        if (nu > upscatter * nu0.back())
            break;

        bool extrapolate = true;
        for (size_t j = nu_size; j-- > 0;) {
            const Real nu0_j = nu0(j);

            if (IC_tab(i, j) < 0) { // integral at (gamma(i), nu(j)) has not been evaluated
                Real nu_comv = gamma_i * nu0_j;
                Real inv = 1 / nu_comv;

                const Real grid_value = I_nu_dnu(j) * sigma(nu_comv) * inv * inv;
                IC_tab(i, j) = (j != nu_size - 1) ? (IC_tab(i, j + 1) + grid_value) : grid_value;
            }

            if (upscatter * nu0_j < nu) {
                IC_I_nu += Ndgamma * (IC_tab(i, j + 1) + (IC_tab(i, j) - IC_tab(i, j + 1)) / (nu0(j + 1) - nu0(j)) *
                                                             (nu0(j + 1) - nu / upscatter));

                extrapolate = false;
                break;
            }
        }

        if (extrapolate) {
            IC_I_nu += Ndgamma *
                       (IC_tab(i, 0) + (IC_tab(i, 0) - IC_tab(i, 1)) / (nu0(1) - nu0(0)) * (nu0(0) - nu / upscatter));
        }
    }

    return IC_I_nu * nu / 4;
}

template <typename Electrons, typename Photons>
Real ICPhoton<Electrons, Photons>::compute_log2_I_nu(Real log2_nu) {
    return std::log2(compute_I_nu(std::exp2(log2_nu)));
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

template <typename Electrons>
Real compute_Thomson_Y(Real B, Real t_com, Real eps_e, Real eps_B, Electrons const& e) {
    Real eta_e = eta_rad(e.gamma_m, e.gamma_c, e.p);
    Real b = eta_e * eps_e / eps_B;
    Real Y0 = (std::sqrt(1 + 4 * b) - 1) / 2;
    Real Y1 = 2 * Y0;
    while (std::fabs((Y1 - Y0) / Y0) > 1e-4) {
        Y1 = Y0;
        Real gamma_c = compute_gamma_c(t_com, B, e.Ys, e.p);
        eta_e = eta_rad(e.gamma_m, gamma_c, e.p);
        b = eta_e * eps_e / eps_B;
        Y0 = (std::sqrt(1 + 4 * b) - 1) / 2;
    }
    return Y0;
}

template <typename Electrons, typename Photons>
void Thomson_cooling(ElectronGrid<Electrons>& electrons, PhotonGrid<Photons>& photons, Shock const& shock) {
    const size_t phi_size = electrons.shape()[0];
    const size_t theta_size = electrons.shape()[1];
    const size_t t_size = electrons.shape()[2];

    for (size_t i = 0; i < phi_size; i++) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                const Real Y_T = compute_Thomson_Y(shock.B(i, j, k), shock.t_comv(i, j, k), shock.rad.eps_e,
                                                   shock.rad.eps_B, electrons(i, j, k));
                electrons(i, j, k).Ys = InverseComptonY(Y_T);
            }
        }
    }
    update_electrons_4Y(electrons, shock);
    generate_syn_photons(photons, shock, electrons);
}

template <typename Electrons, typename Photons>
void KN_cooling(ElectronGrid<Electrons>& electrons, PhotonGrid<Photons>& photons, Shock const& shock) {
    const size_t phi_size = electrons.shape()[0];
    const size_t theta_size = electrons.shape()[1];
    const size_t r_size = electrons.shape()[2];
    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < r_size; ++k) {
                const Real Y_T = compute_Thomson_Y(shock.B(i, j, k), shock.t_comv(i, j, k), shock.rad.eps_e,
                                                   shock.rad.eps_B, electrons(i, j, k));
                // Clear existing Ys and emplace a new InverseComptonY with additional synchrotron frequency parameters.
                electrons(i, j, k).Ys =
                    InverseComptonY(photons(i, j, k).nu_m, photons(i, j, k).nu_c, shock.B(i, j, k), Y_T);
            }
        }
    }
    update_electrons_4Y(electrons, shock);
    generate_syn_photons(photons, shock, electrons);
}
