//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "synchrotron.h"

#include "../core/physics.h"
#include "../util/macros.h"
#include "../util/utilities.h"
#include "inverse-compton.h"

//========================================================================================================
//                                  Helper Functions - Simple Utilities
//========================================================================================================

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Helper function that checks if three values are in non-decreasing order.
 * @param a First value
 * @param b Middle value
 * @param c Last value
 * @return True if a ≤ b ≤ c, false otherwise
 * <!-- ************************************************************************************** -->
 */
inline bool order(Real a, Real b, Real c) {
    return a <= b && b <= c;
};

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Determines the spectral regime (1-6) based on the ordering of characteristic Lorentz factors.
 * @details Classifies the regime based on the ordering of absorption (a), cooling (c),
 *          and minimum (m) Lorentz factors.
 * @param a Absorption Lorentz factor
 * @param c Cooling Lorentz factor
 * @param m Minimum Lorentz factor
 * @return Regime number (1-6) or 0 if no valid regime is found
 * <!-- ************************************************************************************** -->
 */
size_t determine_regime(Real a, Real c, Real m) {
    if (order(a, m, c)) {
        return 1;
    } else if (order(m, a, c)) {
        return 2;
    } else if (order(a, c, m)) {
        return 3;
    } else if (order(c, a, m)) {
        return 4;
    } else if (order(m, c, a)) {
        return 5;
    } else if (order(c, m, a)) {
        return 6;
    } else
        return 0;
}

//========================================================================================================
//                                  Helper Functions - Foundation Physics
//========================================================================================================

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the peak synchrotron power per electron in the comoving frame.
 * @details Based on magnetic field strength B and power-law index p of the electron distribution.
 * @param B Magnetic field strength
 * @param p Power-law index of electron distribution
 * @return Peak synchrotron power per electron
 * <!-- ************************************************************************************** -->
 */
Real compute_single_elec_P_nu_max(Real B, Real /*p*/) {
    constexpr Real sin_angle_ave = con::pi / 4;
    constexpr Real Fx_max = 0.92; // Bing's book 5.5
    return B * (sin_angle_ave * Fx_max * 1.73205080757 * con::e3 / (con::me * con::c2));
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the peak synchrotron intensity for a given column number density.
 * @details Uses the peak synchrotron power and column number density.
 * @param B Magnetic field strength
 * @param p Power-law index of electron distribution
 * @param column_den Column number density
 * @return Peak synchrotron intensity
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_I_peak(Real B, Real p, Real column_den) {
    return compute_single_elec_P_nu_max(B, p) * column_den / (4 * con::pi);
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the characteristic synchrotron frequency for electrons with a given Lorentz factor.
 * @details Uses the standard synchrotron formula.
 * @param gamma Electron Lorentz factor
 * @param B Magnetic field strength
 * @return Characteristic synchrotron frequency
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_freq(Real gamma, Real B) {
    return 3 * con::e / (4 * con::pi * con::me * con::c) * B * gamma * gamma;
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the electron Lorentz factor corresponding to a synchrotron frequency.
 * @details Inverse of the compute_syn_freq function.
 * @param nu Synchrotron frequency
 * @param B Magnetic field strength
 * @return Corresponding electron Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_gamma(Real nu, Real B) {
    return std::sqrt((4 * con::pi * con::me * con::c / (3 * con::e)) * (nu / B));
}

//========================================================================================================
//                                  Helper Functions - Gamma Computations
//========================================================================================================

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the maximum electron Lorentz factor for synchrotron emission.
 * @details Uses an iterative approach to account for inverse Compton cooling effects.
 * @param B Magnetic field strength
 * @param Ys InverseComptonY object
 * @param p Spectral index of electron distribution
 * @return Maximum electron Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_gamma_M(Real B, Real Y, Real /*p*/) {
    if (B == 0) {
        return std::numeric_limits<Real>::infinity();
    }
    return std::sqrt(6 * con::pi * con::e / con::sigmaT / (B * (1 + Y)));
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the minimum electron Lorentz factor for synchrotron emission.
 * @details Accounts for different power-law indices with special handling for the p=2 case.
 *          Uses the fraction of shock energy given to electrons (eps_e) and electron fraction (xi).
 * @param Gamma_th Downstream thermal Lorentz factor
 * @param gamma_M Maximum electron Lorentz factor
 * @param eps_e Fraction of shock energy given to electrons
 * @param p Power-law index of electron distribution
 * @param xi Fraction of electrons accelerated
 * @return Minimum electron Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_gamma_m(Real Gamma_th, Real gamma_M, Real eps_e, Real p, Real xi) {
    const Real gamma_ave_minus_1 = eps_e * (Gamma_th - 1) * (con::mp / con::me) / xi;
    Real gamma_m_minus_1 = 1;
    if (p > 2) {
        gamma_m_minus_1 = (p - 2) / (p - 1) * gamma_ave_minus_1;
    } else if (p < 2) {
        gamma_m_minus_1 = std::pow((2 - p) / (p - 1) * gamma_ave_minus_1 * std::pow(gamma_M, p - 2), 1 / (p - 1));
    } else {
        gamma_m_minus_1 = root_bisect(
            [=](Real x) -> Real {
                return (x * std::log(gamma_M) - (x + 1) * std::log(x) - gamma_ave_minus_1 - std::log(gamma_M));
            },
            0, gamma_M);
    }
    return gamma_m_minus_1 + 1;
}

Real compute_gamma_c(Real t_comv, Real B, Real Y) {
    constexpr Real ad_cooling = 1;
    //-sqrt(Gamma * Gamma - 1) * con::c* t_comv / r;  // adiabatic cooling

    const Real gamma_bar = (6 * con::pi * con::me * con::c / con::sigmaT) / (B * B * (1 + Y) * t_comv) * ad_cooling;
    const Real gamma_c = (gamma_bar + std::sqrt(gamma_bar * gamma_bar + 4)) / 2; // correction on newtonian regime

    return gamma_c;
}

Real cool_after_crossing(Real gamma_x, Real gamma_m_x, Real gamma_m, Real /*dt_comv*/, Real /*B*/, Real /*Y*/) {
    //Real gamma_c_dt = compute_gamma_c(dt_comv, B, Y);
    Real gamma_syn = gamma_x; //* gamma_c_dt / (gamma_x + gamma_c_dt);
    Real f_ad = (gamma_m - 1) / (gamma_m_x - 1);
    return (gamma_syn - 1) * f_ad + 1;
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the self-absorption Lorentz factor by equating synchrotron emission to blackbody.
 * @details Uses the peak intensity and shock parameters to determine where absorption becomes important.
 *          Handles both weak and strong absorption regimes.
 * @param B Magnetic field strength
 * @param I_syn_peak Peak synchrotron intensity
 * @param gamma_m Minimum electron Lorentz factor
 * @param gamma_c Cooling electron Lorentz factor
 * @param gamma_M Maximum electron Lorentz factor
 * @param p Power-law index of electron distribution
 * @return Self-absorption Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_gamma_a(Real B, Real I_syn_peak, Real gamma_m, Real gamma_c, Real /*gamma_M*/, Real p) {
    const Real gamma_peak = std::min(gamma_m, gamma_c);
    const Real nu_peak = compute_syn_freq(gamma_peak, B);

    const Real kT = (gamma_peak - 1) * (con::me * con::c2) / 3;
    // 2kT(nu_a/c)^2 = I_peak*(nu_a/nu_peak)^(1/3) // first assume nu_a is in the 1/3 segment
    Real nu_a = fast_pow(I_syn_peak * con::c2 / (std::cbrt(nu_peak) * 2 * kT), 0.6);

    if (nu_a > nu_peak) {        // nu_a is not in the 1/3 segment
        if (gamma_c > gamma_m) { // first assume nu_a is in the -(p-1)/2 segment
            const Real nu_m = compute_syn_freq(gamma_m, B);
            nu_a = fast_pow(I_syn_peak * con::c2 / (2 * kT) * fast_pow(nu_m, p / 2), 2 / (p + 4));
            const Real nu_c = compute_syn_freq(gamma_c, B);
            if (nu_a > nu_c) { //  nu_a is not in the -(p-1)/2 but -p/2 segment
                nu_a = fast_pow(I_syn_peak * con::c2 / (2 * kT) * std::sqrt(nu_c) * fast_pow(nu_m, p / 2), 2 / (p + 5));
            }
        } else { //first assume nu_a is in the -1/2 segment
            const Real nu_c = compute_syn_freq(gamma_c, B);
            nu_a = fast_pow(I_syn_peak * con::c2 / (2 * kT) * std::sqrt(nu_c), 0.4);
            const Real nu_m = compute_syn_freq(gamma_m, B);
            if (nu_a > nu_m) { // nu_a is not in the -1/2 segment but -p/2 segment
                nu_a = fast_pow(I_syn_peak * con::c2 / (2 * kT) * std::sqrt(nu_c) * fast_pow(nu_m, p / 2), 2 / (p + 5));
            }
        }
    }
    return compute_syn_gamma(nu_a, B) + 1;
}

Real compute_gamma_peak(Real gamma_a, Real gamma_m, Real gamma_c) {
    const Real gamma_peak = std::min(gamma_m, gamma_c);
    if (gamma_a > gamma_c) {
        return gamma_a;
    } else {
        return gamma_peak;
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Determines the peak Lorentz factor directly from a SynElectrons object.
 * @details Convenient wrapper around the three-parameter version.
 * @param e Synchrotron electron object
 * @return Peak Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_gamma_peak(SynElectrons const& e) {
    return compute_gamma_peak(e.gamma_a, e.gamma_m, e.gamma_c);
}

Real cyclotron_correction(Real gamma_m, Real p) {
    Real f = (gamma_m - 1) / gamma_m;
    if (p > 3) {
        f = fast_pow(f, (p - 1) / 2);
    }
    return f;
}

//========================================================================================================
//                                  SynElectrons Class Methods
//========================================================================================================

Real SynElectrons::compute_spectrum(Real gamma) const {
    // Smooth broken power law with sharpness s=1
    // General s: fast_pow(1 + fast_pow(gamma/gamma_break, s*delta), -1/s)
    // With s=1:  1 / (1 + fast_pow(gamma/gamma_break, delta))

    switch (regime) {
        case 1: // slow cooling: gamma_m < gamma_c
        case 2:
        case 5:
            // (p-1)/gamma_m * (gamma/gamma_m)^-p * gamma_c/(gamma + gamma_c) * exp cutoffs
            return (p - 1) / gamma_m * fast_exp(-gamma / gamma_M - gamma_m / gamma) * fast_pow(gamma / gamma_m, -p) *
                   gamma_c / (gamma + gamma_c);

        case 3: // fast cooling: gamma_c < gamma_m
        case 4:
        case 6:
            // gamma_c/gamma^2 / (1 + (gamma/gamma_m)^(p-1)) * exp cutoffs
            return fast_exp(-gamma / gamma_M - gamma_c / gamma) * gamma_c / (gamma * gamma) /
                   (1.0 + fast_pow(gamma / gamma_m, p - 1));

        default:
            return 0;
    }
}

Real SynElectrons::compute_N_gamma(Real gamma) const {
    if (gamma <= gamma_c) { // Below the cooling Lorentz factor: direct scaling
        return N_e * compute_spectrum(gamma);
    } else {
        return N_e * compute_spectrum(gamma) * (1 + Y_c) / (1 + Ys.gamma_spectrum(gamma));
    }
}

Real SynElectrons::compute_column_den(Real gamma) const {
    if (gamma <= gamma_c) { // Below the cooling Lorentz factor: direct scaling
        return column_den * compute_spectrum(gamma);
    } else {
        return column_den * compute_spectrum(gamma) * (1 + Y_c) / (1 + Ys.gamma_spectrum(gamma));
    }
}

//========================================================================================================
//                                  Factory Functions - Synchrotron Electrons
//========================================================================================================

SynElectronGrid generate_syn_electrons(Shock const& shock) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    SynElectronGrid electrons({phi_size, theta_size, t_size});

    generate_syn_electrons(electrons, shock);

    return electrons;
}

void generate_syn_electrons(SynElectronGrid& electrons, Shock const& shock) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    const RadParams rad = shock.rad;

    electrons.resize({phi_size, theta_size, t_size});

    const size_t phi_compute = (shock.symmetry != Symmetry::structured) ? 1 : phi_size;

    for (size_t i = 0; i < phi_compute; ++i) {
        for (size_t j : shock.theta_reps) {
            const size_t k_inj = shock.injection_idx(i, j);
            for (size_t k = 0; k < t_size; ++k) {
                const Real t_com = shock.t_comv(i, j, k);
                const Real B = shock.B(i, j, k);
                const Real r = shock.r(i, j, k);
                // Real Gamma = shock.Gamma(i, j, k);
                const Real Gamma_th = shock.Gamma_th(i, j, k);

                auto& elec = electrons(i, j, k);

                elec.gamma_M = compute_syn_gamma_M(B, 0., rad.p);
                elec.gamma_m = compute_syn_gamma_m(Gamma_th, elec.gamma_M, rad.eps_e, rad.p, rad.xi_e);

                // Fraction of synchrotron electrons; the rest are cyclotron
                const Real f_syn = cyclotron_correction(elec.gamma_m, rad.p);

                elec.N_e = shock.N_p(i, j, k) * rad.xi_e * f_syn;
                elec.column_den = elec.N_e / (r * r);
                const Real I_nu_peak = compute_syn_I_peak(B, rad.p, elec.column_den);

                // no new shocked electrons, cool the relic population from crossing
                if (k >= k_inj) {
                    auto const& inj = electrons(i, j, k_inj - 1);
                    const Real dt_comv = t_com - shock.t_comv(i, j, k_inj - 1);
                    elec.gamma_c = cool_after_crossing(inj.gamma_c, inj.gamma_m, elec.gamma_m, t_com, B, 0);

                    elec.gamma_M = cool_after_crossing(inj.gamma_M, inj.gamma_m, elec.gamma_m, dt_comv, B, 0);
                } else {
                    elec.gamma_c = compute_gamma_c(t_com, B, 0.);
                }

                elec.gamma_a = compute_syn_gamma_a(B, I_nu_peak, elec.gamma_m, elec.gamma_c, elec.gamma_M, rad.p);
                elec.regime = determine_regime(elec.gamma_a, elec.gamma_c, elec.gamma_m);
                elec.p = rad.p;
            }
        }
    }

    broadcast_symmetry(electrons, shock);
}

//========================================================================================================
//                                  Factory Functions - Synchrotron Photons
//========================================================================================================

SynPhotonGrid generate_syn_photons(Shock const& shock, SynElectronGrid const& electrons) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    SynPhotonGrid photons({phi_size, theta_size, t_size});

    generate_syn_photons(photons, shock, electrons);

    return photons;
}

void generate_syn_photons(SynPhotonGrid& photons, Shock const& shock, SynElectronGrid const& electrons) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    photons.resize({phi_size, theta_size, t_size});

    const size_t phi_compute = (shock.symmetry != Symmetry::structured) ? 1 : phi_size;

    for (size_t i = 0; i < phi_compute; ++i) {
        for (size_t j : shock.theta_reps) {
            for (size_t k = 0; k < t_size; ++k) {
                auto& ph = photons(i, j, k);
                auto& elec = electrons(i, j, k);
                ph.p = elec.p;
                ph.Ys = elec.Ys;
                ph.Y_c = elec.Y_c;
                ph.regime = elec.regime;

                const Real B = shock.B(i, j, k);

                ph.nu_M = compute_syn_freq(elec.gamma_M, B);
                ph.nu_m = compute_syn_freq(elec.gamma_m, B);
                ph.nu_c = compute_syn_freq(elec.gamma_c, B);
                ph.nu_a = compute_syn_freq(elec.gamma_a, B);
                ph.I_nu_max = compute_syn_I_peak(B, elec.p, elec.column_den);

                ph.build();
            }
        }
    }

    broadcast_symmetry(photons, shock);
}
