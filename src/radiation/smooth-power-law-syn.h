//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once

#include "inverse-compton.h"

/**
 * <!-- ************************************************************************************** -->
 * @struct SmoothPowerLawSyn
 * @brief Represents synchrotron photons in the comoving frame and provides spectral functions.
 * <!-- ************************************************************************************** -->
 */
struct SmoothPowerLawSyn {
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
    [[nodiscard]] Real compute_I_nu(Real nu) const noexcept; ///< Linear power PER SOLID ANGLE

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the base-2 logarithm of comoving synchrotron specific intensity at a given frequency.
     * @details Optimized for numerical computation by using logarithmic arithmetic.
     * @param log2_nu Base-2 logarithm of the frequency
     * @return Base-2 logarithm of synchrotron specific intensity
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real compute_log2_I_nu(Real log2_nu) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Updates cached calculation constants used for efficiently computing synchrotron spectra.
     * @details Constants vary based on the electron regime (1-6) and involve different power laws.
     * <!-- ************************************************************************************** -->
     */
    void build() noexcept;

  private:
    Real log2_norm_{0};       ///< Cached spectral coefficient 0 in log2
    Real log2_thick_norm_{0}; ///< Cached spectral coefficient 1 in log2
    Real smooth_m_fast_{1};   ///<
    Real smooth_c_fast_{1};   ///<
    Real smooth_m_slow_{1};   ///<
    Real smooth_c_slow_{1};   ///<
    Real smooth_a_;           ///<

    // Precomputed s*(beta_lo - beta_hi) for broken power law transitions
    Real diff_slope_m_slow_{0}; ///< smooth_m_slow_ * (1/3 - (1-p)/2)
    Real diff_slope_c_slow_{0}; ///< smooth_c_slow_ * ((1-p)/2 - (-p/2))
    Real diff_slope_m_fast_{0}; ///< smooth_m_fast_ * (1/3 - (-1/2))
    Real diff_slope_c_fast_{0}; ///< smooth_c_fast_ * (-1/2 - (-p/2))
    Real inv_nu_M_{0};          ///< Cached 1/nu_M for division optimization

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the synchrotron spectrum at a given frequency based on the electron regime.
     * @details Implements the broken power-law with exponential cutoff formulae for different regimes.
     * @param nu The frequency at which to compute the spectrum
     * @return The normalized synchrotron spectrum value
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] inline Real compute_spectrum(Real nu) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the base-2 logarithm of synchrotron spectrum at a given frequency.
     * @details Uses logarithmic arithmetic for numerical stability in different spectral regimes.
     * @param log2_nu Base-2 logarithm of the frequency
     * @return Base-2 logarithm of the synchrotron spectrum
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] inline Real compute_log2_spectrum(Real log2_nu) const noexcept;

    [[nodiscard]] inline Real log2_optical_thin_sharp(Real log2_nu) const noexcept;
    [[nodiscard]] inline Real log2_optical_thick_sharp(Real log2_nu) const noexcept;

    [[nodiscard]] inline Real log2_optical_thin(Real log2_nu) const noexcept;
    [[nodiscard]] inline Real log2_optical_thick(Real log2_nu) const noexcept;
};
