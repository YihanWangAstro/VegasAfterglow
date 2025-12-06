//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "inverse-compton.h"

#include "../core/physics.h"
#include "../util/IO.h"
#include "../util/macros.h"
#include "../util/utilities.h"

//========================================================================================================
//                                  InverseComptonY Constructors
//========================================================================================================

InverseComptonY::InverseComptonY(Real gamma_m, Real gamma_c, Real B, Real Y_T, Real eps_e_on_eps_B) noexcept {
    constexpr Real B_QED = 4 * con::pi * con::me * con::me * con::c2 * con::c / (3 * con::e * con::h);
    const Real nu_m = compute_syn_freq(gamma_m, B);
    const Real nu_c = compute_syn_freq(gamma_c, B);
    gamma_m_hat = con::me * con::c2 / con::h / nu_m;
    gamma_c_hat = con::me * con::c2 / con::h / nu_c;

    this->Y_T = Y_T; // Set the Thomson Y parameter
    nu_m_hat = compute_syn_freq(gamma_m_hat, B);
    nu_c_hat = compute_syn_freq(gamma_c_hat, B);
    if (gamma_m <= gamma_c) { // slow IC cooling regime
        regime = 1;
    } else {                          //fast IC cooling regime
        if (gamma_m <= gamma_m_hat) { //weak KN effects
            regime = 2;
        } else { //strong KN effects
            gamma_self = std::cbrt(B_QED / B);
            nu_self = compute_syn_freq(gamma_self, B);

            gamma0 = std::sqrt(eps_e_on_eps_B * gamma_m * gamma_m_hat);
            nu0 = compute_syn_freq(gamma0, B);

            if (gamma0 > gamma_m) {
                regime = 3;
                //.../
            } else if (gamma0 > gamma_self) {
                regime = 4;
                //../
            } else {
                regime = 5;
            }
        }
    }
}

InverseComptonY::InverseComptonY(Real Y_T) noexcept {
    this->Y_T = Y_T; // Set the Thomson Y parameter
    regime = 0;      // Set regime to 0 (Thomson only, no KN effects)
}

InverseComptonY::InverseComptonY() noexcept {
    nu_m_hat = 0;
    nu_c_hat = 0;
    gamma_m_hat = 0;
    gamma_c_hat = 0;
    Y_T = 0;
    regime = 0;
}

//========================================================================================================
//                                  InverseComptonY Public Methods
//========================================================================================================

Real InverseComptonY::evaluate_at_gamma(Real gamma, Real p) const {
    switch (regime) {
        case 0:
            return Y_T; // Thomas only, no KN effects
            break;
        case 1:
            if (gamma <= gamma_c_hat) {
                return Y_T; // For gamma below gamma_hat_c, no modification
            } else if (gamma <= gamma_m_hat) {
                return Y_T * fast_pow(gamma / gamma_c_hat, (p - 3) / 2); // Scaling in intermediate regime
            } else {
                return Y_T * pow43(gamma_m_hat / gamma) *
                       fast_pow(gamma_m_hat / gamma_c_hat, (p - 3) / 2); // High gamma scaling
            }
            break;
        case 2:
            if (gamma <= gamma_m_hat) {
                return Y_T; // For gamma below gamma_hat_m, no modification
            } else if (gamma <= gamma_c_hat) {
                return Y_T / std::sqrt(gamma / gamma_m_hat); // Intermediate regime scaling
            } else {
                return Y_T * pow43(gamma_c_hat / gamma) * std::sqrt(gamma_m_hat / gamma_c_hat); // High gamma scaling
            }
            break;

        default:
            return 0;
            break;
    }
}

Real InverseComptonY::evaluate_at_nu(Real nu, Real p) const {
    switch (regime) {
        case 0:
            return Y_T; // Thomas only, no KN effects
            break;
        case 1:
            if (nu <= nu_c_hat) {
                return Y_T; // For frequencies below nu_hat_c, no modification
            } else if (nu <= nu_m_hat) {
                return Y_T * fast_pow(nu / nu_c_hat, (p - 3) / 4); // Intermediate frequency scaling
            } else {
                return Y_T * pow23(nu_m_hat / nu) *
                       fast_pow(nu_m_hat / nu_c_hat, (p - 3) / 4); // High-frequency scaling
            }
            break;
        case 2:
            if (nu <= nu_m_hat) {
                return Y_T; // For frequencies below nu_hat_m, no modification
            } else if (nu <= nu_c_hat) {
                return Y_T * std::sqrt(std::sqrt(nu_m_hat / nu)); // Intermediate frequency scaling
            } else {
                return Y_T * pow23(nu_c_hat / nu) * std::sqrt(std::sqrt(nu_m_hat / nu_c_hat)); // High-frequency scaling
            }
            break;

        default:
            return 0;
            break;
    }
}

//========================================================================================================
//                                  Helper Functions for Update Functions
//========================================================================================================

void update_gamma_c_Thomson(Real& gamma_c, InverseComptonY& Ys, RadParams const& rad, Real B, Real t_com,
                            Real gamma_m) {
    Real eta_e = eta_rad(gamma_m, gamma_c, rad.p);
    Real b = eta_e * rad.eps_e / rad.eps_B;
    Real Y_T = (std::sqrt(1 + 4 * b) - 1) / 2;

    Real gamma_c_new = compute_gamma_c(t_com, B, Y_T);
    while (std::fabs((gamma_c_new - gamma_c) / gamma_c) > 1e-3) {
        gamma_c = gamma_c_new;
        eta_e = eta_rad(gamma_m, gamma_c, rad.p);
        b = eta_e * rad.eps_e / rad.eps_B;
        Y_T = (std::sqrt(1 + 4 * b) - 1) / 2;
        gamma_c_new = compute_gamma_c(t_com, B, Y_T);
    }
    gamma_c = gamma_c_new;
    Ys = InverseComptonY(Y_T);
}

void update_gamma_c_KN(Real& gamma_c, InverseComptonY& Ys, RadParams const& rad, Real B, Real t_com, Real gamma_m) {
    Real eta_e = eta_rad(gamma_m, gamma_c, rad.p);
    Real b = eta_e * rad.eps_e / rad.eps_B;
    Real Y_T = (std::sqrt(1 + 4 * b) - 1) / 2;
    Ys = InverseComptonY(gamma_m, gamma_c, B, Y_T, rad.eps_e / rad.eps_B);
    Real Y_c = 0; //Ys.evaluate_at_gamma(gamma_c, rad.p);
    Real gamma_c_new = compute_gamma_c(t_com, B, Y_c);

    while (std::fabs((gamma_c_new - gamma_c) / gamma_c) > 1e-3) {
        gamma_c = gamma_c_new;
        eta_e = eta_rad(gamma_m, gamma_c, rad.p);
        b = eta_e * rad.eps_e / rad.eps_B;
        Y_T = (std::sqrt(1 + 4 * b) - 1) / 2;
        Ys = InverseComptonY(gamma_m, gamma_c, B, Y_T, rad.eps_e / rad.eps_B);
        Y_c = Ys.evaluate_at_gamma(gamma_c, rad.p);
        gamma_c_new = compute_gamma_c(t_com, B, Y_c);
    }
    gamma_c = gamma_c_new;
}

void update_gamma_M(Real& gamma_M, InverseComptonY const& Ys, Real p, Real B) {
    if (B == 0) {
        gamma_M = std::numeric_limits<Real>::infinity();
        return;
    }

    Real Y_M = Ys.evaluate_at_gamma(gamma_M, p);
    Real gamma_M_new = compute_syn_gamma_M(B, Y_M, p);

    while (std::fabs((gamma_M - gamma_M_new) / gamma_M_new) > 1e-3) {
        gamma_M = gamma_M_new;
        Y_M = Ys.evaluate_at_gamma(gamma_M, p);
        gamma_M_new = compute_syn_gamma_M(B, Y_M, p);
    }
}

//========================================================================================================
//                                  Standalone Physics Functions
//========================================================================================================

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
