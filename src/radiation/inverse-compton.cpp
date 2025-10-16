//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "inverse-compton.h"

#include "IO.h"
#include "macros.h"
#include "physics.h"
#include "utilities.h"

InverseComptonY::InverseComptonY(Real gamma_m, Real gamma_c, Real B, Real Y_T) noexcept {
    const Real nu_m = compute_syn_freq(gamma_m, B);  // Compute minimum synchrotron frequency
    const Real nu_c = compute_syn_freq(gamma_c, B);  // Compute cooling synchrotron frequency
    gamma_m_hat = con::me * con::c2 / con::h / nu_m; // Compute minimum characteristic Lorentz factor
    gamma_c_hat = con::me * con::c2 / con::h / nu_c; // Compute cooling characteristic Lorentz factor
    this->Y_T = Y_T;                                 // Set the Thomson Y parameter
    nu_m_hat = compute_syn_freq(gamma_m_hat, B);     // Compute the corresponding synchrotron frequency for gamma_hat_m
    nu_c_hat = compute_syn_freq(gamma_c_hat, B);     // Compute the corresponding synchrotron frequency for gamma_hat_c

    if (nu_m_hat <= nu_c_hat) {
        regime = 1; // fast IC cooling regime
    } else {
        regime = 2; // slow IC cooling regime
    }
}

InverseComptonY::InverseComptonY(Real Y_T) noexcept {
    this->Y_T = Y_T; // Set the Thomson Y parameter
    regime = 3;      // Set regime to 3 (special case)
}

InverseComptonY::InverseComptonY() noexcept {
    nu_m_hat = 0;
    nu_c_hat = 0;
    gamma_m_hat = 0;
    gamma_c_hat = 0;
    Y_T = 0;
    regime = 0;
}

Real InverseComptonY::evaluate_at_gamma(Real gamma, Real p) const {
    switch (regime) {
        case 3:
            return Y_T; // In regime 3, simply return Y_T
            break;
        case 1:
            if (gamma <= gamma_m_hat) {
                return Y_T; // For gamma below gamma_hat_m, no modification
            } else if (gamma <= gamma_c_hat) {
                return Y_T / std::sqrt(gamma / gamma_m_hat); // Intermediate regime scaling
            } else {
                return Y_T * pow43(gamma_c_hat / gamma) * std::sqrt(gamma_m_hat / gamma_c_hat); // High gamma scaling
            }
            break;
        case 2:
            if (gamma <= gamma_c_hat) {
                return Y_T; // For gamma below gamma_hat_c, no modification
            } else if (gamma <= gamma_m_hat) {
                return Y_T * fast_pow(gamma / gamma_c_hat, (p - 3) / 2); // Scaling in intermediate regime
            } else {
                return Y_T * pow43(gamma_m_hat / gamma) *
                       fast_pow(gamma_m_hat / gamma_c_hat, (p - 3) / 2); // High gamma scaling
            }
            break;
        default:
            return 0;
            break;
    }
}

Real InverseComptonY::evaluate_at_nu(Real nu, Real p) const {
    switch (regime) {
        case 3:
            return Y_T; // In regime 3, simply return Y_T
            break;
        case 1:
            if (nu <= nu_m_hat) {
                return Y_T; // For frequencies below nu_hat_m, no modification
            } else if (nu <= nu_c_hat) {
                return Y_T * std::sqrt(std::sqrt(nu_m_hat / nu)); // Intermediate frequency scaling
            } else {
                return Y_T * pow23(nu_c_hat / nu) * std::sqrt(std::sqrt(nu_m_hat / nu_c_hat)); // High-frequency scaling
            }
            break;
        case 2:
            if (nu <= nu_c_hat) {
                return Y_T; // For frequencies below nu_hat_c, no modification
            } else if (nu <= nu_m_hat) {
                return Y_T * fast_pow(nu / nu_c_hat, (p - 3) / 4); // Intermediate frequency scaling
            } else {
                return Y_T * pow23(nu_m_hat / nu) *
                       fast_pow(nu_m_hat / nu_c_hat, (p - 3) / 4); // High-frequency scaling
            }
            break;
        default:
            return 0;
            break;
    }
}

Real compton_cross_section(Real nu) {
    const Real x = con::h / (con::me * con::c2) * nu;
    /*if (x <= 1) {
        return con::sigmaT;
    } else {
        return 0;
    }*/

    if (x < 1e-2) {
        return con::sigmaT * (1 - 2 * x);
    } else if (x > 1e2) {
        return 3. / 8 * con::sigmaT * (log(2 * x) + 0.5) / x;
    } else {
        const Real l = std::log1p(2.0 * x); // log(1+2x)
        const Real invx = 1.0 / x;
        const Real invx2 = invx * invx;
        const Real term1 = 1.0 + 2.0 * x;
        const Real invt1 = 1.0 / term1;
        const Real invt1_2 = invt1 * invt1;

        // ((1+x)/x^3) * (2x(1+x)/(1+2x) - log(1+2x)) + log(1+2x)/(2x) - (1+3x)/(1+2x)^2
        const Real a = (1.0 + x) * invx2 * invx;        // (1+x)/x^3
        const Real b = 2.0 * x * (1.0 + x) * invt1 - l; // bracket
        const Real c = 0.5 * l * invx;                  // log_term/(2x)
        const Real d = (1.0 + 3.0 * x) * invt1_2;       // (1+3x)/(1+2x)^2

        return 0.75 * con::sigmaT * (a * b + c - d);
    }
}
