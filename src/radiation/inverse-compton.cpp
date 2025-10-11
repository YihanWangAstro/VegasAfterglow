//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "inverse-compton.h"

#include "IO.h"
#include "macros.h"
#include "utilities.h"

InverseComptonY::InverseComptonY(Real nu_m, Real nu_c, Real B, Real Y_T) noexcept {
    gamma_hat_m = con::me * con::c2 / con::h / nu_m; // Compute minimum characteristic Lorentz factor
    gamma_hat_c = con::me * con::c2 / con::h / nu_c; // Compute cooling characteristic Lorentz factor
    this->Y_T = Y_T;                                 // Set the Thomson Y parameter
    nu_hat_m = compute_syn_freq(gamma_hat_m, B);     // Compute the corresponding synchrotron frequency for gamma_hat_m
    nu_hat_c = compute_syn_freq(gamma_hat_c, B);     // Compute the corresponding synchrotron frequency for gamma_hat_c

    if (nu_hat_m <= nu_hat_c) {
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
    nu_hat_m = 0;
    nu_hat_c = 0;
    gamma_hat_m = 0;
    gamma_hat_c = 0;
    Y_T = 0;
    regime = 0;
}

Real InverseComptonY::compute_val_at_gamma(Real gamma, Real p) const {
    switch (regime) {
        case 3:
            return Y_T; // In regime 3, simply return Y_T
            break;
        case 1:
            if (gamma <= gamma_hat_m) {
                return Y_T; // For gamma below gamma_hat_m, no modification
            } else if (gamma <= gamma_hat_c) {
                return Y_T / std::sqrt(gamma / gamma_hat_m); // Intermediate regime scaling
            } else
                return Y_T * pow43(gamma_hat_c / gamma) * std::sqrt(gamma_hat_m / gamma_hat_c); // High gamma scaling

            break;
        case 2:
            if (gamma <= gamma_hat_c) {
                return Y_T; // For gamma below gamma_hat_c, no modification
            } else if (gamma <= gamma_hat_m) {
                return Y_T * fast_pow(gamma / gamma_hat_c, (p - 3) / 2); // Scaling in intermediate regime
            } else
                return Y_T * pow43(gamma_hat_m / gamma) *
                       fast_pow(gamma_hat_m / gamma_hat_c, (p - 3) / 2); // High gamma scaling

            break;
        default:
            return 0;
            break;
    }
}

Real InverseComptonY::compute_val_at_nu(Real nu, Real p) const {
    switch (regime) {
        case 3:
            return Y_T; // In regime 3, simply return Y_T
            break;
        case 1:
            if (nu <= nu_hat_m) {
                return Y_T; // For frequencies below nu_hat_m, no modification
            } else if (nu <= nu_hat_c) {
                return Y_T * std::sqrt(std::sqrt(nu_hat_m / nu)); // Intermediate frequency scaling
            } else
                return Y_T * pow23(nu_hat_c / nu) * std::sqrt(std::sqrt(nu_hat_m / nu_hat_c)); // High-frequency scaling

            break;
        case 2:
            if (nu <= nu_hat_c) {
                return Y_T; // For frequencies below nu_hat_c, no modification
            } else if (nu <= nu_hat_m) {
                return Y_T * fast_pow(nu / nu_hat_c, (p - 3) / 4); // Intermediate frequency scaling
            } else
                return Y_T * pow23(nu_hat_m / nu) *
                       fast_pow(nu_hat_m / nu_hat_c, (p - 3) / 4); // High-frequency scaling

            break;
        default:
            return 0;
            break;
    }
}

Real InverseComptonY::compute_Y_Thompson(InverseComptonY const& Ys) {
    return Ys.Y_T;
}

Real InverseComptonY::compute_Y_tilt_at_gamma(InverseComptonY const& Ys, Real gamma, Real p) {
    return Ys.compute_val_at_gamma(gamma, p);
}

Real InverseComptonY::compute_Y_tilt_at_nu(InverseComptonY const& Ys, Real nu, Real p) {
    return Ys.compute_val_at_nu(nu, p);
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
