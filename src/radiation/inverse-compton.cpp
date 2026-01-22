//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "inverse-compton.h"

#include "../core/physics.h"
#include "../util/macros.h"
#include "../util/utilities.h"
#include "synchrotron.h"

//========================================================================================================
//                                  InverseComptonY Constructors
//========================================================================================================
InverseComptonY::InverseComptonY(Real gamma_m, Real gamma_c, Real p, Real B, Real Y_T, bool is_KN) noexcept {
    const Real nu_m = compute_syn_freq(gamma_m, B);
    this->gamma_m_hat = std::max(con::me * con::c2 / con::h / nu_m, 1.0);

    this->gamma_self = fast_pow(gamma_m_hat * gamma_m * gamma_m, 1.0 / 3.0);

    this->gamma_self3 = gamma_self * gamma_self * gamma_self;

    this->B_ = B;
    this->gamma_m_ = gamma_m;
    this->p_ = p;
    if (is_KN) {
        update_cooling_breaks(gamma_c, Y_T);
    } else {
        this->gamma_c_hat = compute_gamma_hat(gamma_c);
        this->Y_T = Y_T;
        regime = 0;
        build_segments();
    }
}

InverseComptonY::InverseComptonY() noexcept {
    gamma_m_hat = 1.0;
    gamma_c_hat = 1.0;
    Y_T = 0.0;
    regime = 0;
    active_segment_count = 0;
}

//========================================================================================================
//                                  InverseComptonY Methods
//========================================================================================================

void InverseComptonY::update_gamma0(Real gamma_c) noexcept {
    if (Y_T < 1) {
        gamma0 = 0.0;
        return;
    }

    if (gamma_m_ < gamma_c) { // slow cooling
        gamma0 = fast_pow(Y_T, 2.0 / (3.0 - p_)) * gamma_c_hat;
        if (gamma0 > gamma_m_hat) {
            gamma0 = gamma_m_hat * fast_pow(Y_T, 3.0 / 4.0) * fast_pow(gamma_c / gamma_m_, 0.75 * (p_ - 3.0));
        }

    } else {                          //fast cooling
        if (gamma_m_ < gamma_m_hat) { //weak KN regime (gamma_m < gamma_m_hat)
            gamma0 = Y_T * Y_T * gamma_m_hat;
            if (gamma0 > gamma_c_hat) {
                gamma0 = fast_pow(Y_T * gamma_c / gamma_m_, 3.0 / 4.0) * gamma_c_hat;
            }
        } else { //strong KN regime (gamma_m_hat < gamma_self < gamma_m)
            gamma0 = Y_T * Y_T * gamma_m_hat;

            if (gamma0 > gamma_self) {
                gamma0 = std::sqrt(Y_T * gamma_m_ * gamma_m_hat);
            }

            /*if (gamma0 > gamma_m_) {
                gamma0 = gamma_m_hat * fast_pow(Y_T, 1.0 / (4.0 - p_)) *
                         fast_pow(gamma_m_hat / gamma_m_, (p_ - 3.0) / (4.0 - p_));
            }
            Real gamma_m_hat_hat = compute_gamma_hat(gamma_m_hat);
            if (gamma0 > gamma_m_hat_hat) {

                Real Y_B3 = Y_T * fast_pow(gamma_m_hat_hat / gamma_m_, p_ - 3.0) * gamma_m_hat / gamma_m_hat_hat;
                gamma0 = gamma_m_hat_hat * fast_pow(Y_B3, 4.0 / (5.0 - p_));
            }*/
        }
    }
}

void InverseComptonY::update_cooling_breaks(Real gamma_c, Real Y_T) noexcept {
    gamma_c_hat = compute_gamma_hat(gamma_c);

    this->Y_T = Y_T; // Set the normalization

    update_gamma0(gamma_c);

    if (gamma_m_ < gamma_c) { // slow cooling
        regime = 1;
    } else {
        if (gamma_m_ < gamma_m_hat) { //weak KN regime (gamma_m < gamma_m_hat)
            regime = 2;
        } else {                 //strong KN regime (gamma_m_hat < gamma_self < gamma_m)
            if (gamma0 == 0.0) { // No gamma_0 because Y(gamma)<1 everywhere
                regime = 2;
            } else if (gamma0 < gamma_self) { // (gamma_m_hat <= gamma0 < gamma_self)
                regime = 3;
            } else if (gamma0 < gamma_m_) { // (gamma_self <= gamma0 < gamma_m)
                regime = 4;
            } else { // (gamma_m < gamma0)
                regime = 5;
            }
        }
    }
    build_segments();
}

void InverseComptonY::build_segments() noexcept {
    active_segment_count = 0;

    Real current_val = Y_T;
    Real current_gamma = 1.0;

    auto chain = [&](Real target_gamma, Real slope) {
        segments[active_segment_count++] = {target_gamma, current_val, slope, current_gamma};

        if (target_gamma > current_gamma) {
            current_val *= fast_pow(target_gamma / current_gamma, slope);
            current_gamma = target_gamma;
        }
    };

    constexpr Real INF_GAMMA = 1.0e100;

    switch (regime) {
        case 0: // Thomson
            chain(INF_GAMMA, 0.0);
            break;

        case 1: // Slow Cooling
            chain(gamma_c_hat, 0.0);
            chain(gamma_m_hat, 0.5 * (p_ - 3.0));
            chain(INF_GAMMA, -4.0 / 3.0);
            break;

        case 2: // Fast Cooling, Weak KN
            chain(gamma_m_hat, 0.0);
            chain(gamma_c_hat, -0.5);
            chain(INF_GAMMA, -4.0 / 3.0);
            break;

        case 3: // Fast Cooling, Strong KN (gamma_m_hat <= gamma0 < gamma_self)
        {
            Real gamma0_hat = compute_gamma_hat(gamma0);
            Real gamma_m_hat_hat = compute_gamma_hat(gamma_m_hat);

            chain(gamma_m_hat, 0.0);
            chain(gamma0_hat, -0.5);
            chain(gamma_m_hat_hat, -0.75);
            chain(INF_GAMMA, -0.5);
            break;
        }
        case 4: // Fast Cooling, Strong KN (gamma_self <= gamma0 < gamma_m)
        {
            Real gamma0_hat = compute_gamma_hat(gamma0);
            Real gamma0_hat_hat = compute_gamma_hat(gamma0_hat);
            Real gamma_m_hat_hat = compute_gamma_hat(gamma_m_hat);

            chain(gamma_m_hat, 0.0);
            chain(gamma0_hat, -0.5);
            chain(gamma0_hat_hat, -1.0);
            chain(gamma_m_hat_hat, -0.75);
            chain(INF_GAMMA, -0.5);
            break;
        }
        case 5: // Fast Cooling, Strong KN (gamma_m < gamma0)
        {
            Real gamma0_hat = compute_gamma_hat(gamma0);
            Real gamma_m_hat_hat = compute_gamma_hat(gamma_m_hat);
            Real gamma0_hat_hat = compute_gamma_hat(gamma0_hat);

            chain(gamma0_hat, 0.0);
            chain(gamma_m_hat, 0.5 * (p_ - 3.0));
            chain(gamma_m_hat_hat, -1.0);
            chain(gamma0_hat_hat, 0.25 * (p_ - 5.0));
            chain(INF_GAMMA, -0.5);
            break;
        }
    }
}

Real InverseComptonY::gamma_spectrum(Real gamma) const {
    for (int i = 0; i < active_segment_count; ++i) {
        if (gamma <= segments[i].gamma_max) {
            return segments[i].eval(gamma);
        }
    }
    return 0.0;
}

Real InverseComptonY::nu_spectrum(Real nu) const {
    const Real gamma = compute_syn_gamma(nu, B_);
    return gamma_spectrum(gamma);
}

Real InverseComptonY::compute_gamma_hat(Real gamma) const noexcept {
    return std::max(gamma_self3 / (gamma * gamma), 1.0);
}
//========================================================================================================
//                                  Helper Functions for Update Functions
//========================================================================================================
void update_gamma_c_Thomson(Real& gamma_c, InverseComptonY& Ys, RadParams const& rad, Real B, Real t_com, Real gamma_m,
                            Real gamma_c_last) {
    Real Y_T = compute_Thomson_Y(rad, gamma_m, gamma_c);
    Real gamma_c_new = gamma_c_last; //compute_gamma_c(t_com, B, Y_T);

    while (std::fabs((gamma_c_new - gamma_c) / gamma_c) > 1e-3) {
        gamma_c = gamma_c_new;
        Y_T = compute_Thomson_Y(rad, gamma_m, gamma_c);
        gamma_c_new = compute_gamma_c(t_com, B, Y_T);
    }
    gamma_c = gamma_c_new;
    Ys = InverseComptonY(gamma_m, gamma_c, rad.p, B, Y_T, false);
}

void update_gamma_c_KN(Real& gamma_c, InverseComptonY& Ys, RadParams const& rad, Real B, Real t_com, Real gamma_m,
                       Real gamma_c_last) {
    // the iteration may converge to just one solution or three solutions, when there are three solutions:
    // 1. gamma_c_min (IC cooling dominant, Y(gamma_c) > 1). solution for gamma_m < gamma_c_trans
    // 2. gamma_c_trans (unstable state),
    // 3. gamma_c_max (synchrotron cooling dominant, Y(gamma_c) < 1). solution for gamma_m > gamma_c_trans

    Real gamma_c_sync = compute_gamma_c(t_com, B, 0);
    Real gamma_c_trans = std::max(gamma_c_sync * 0.5, 1.0); // where Y(gamma_c) = 1
    Real gamma_c_new = 2 * gamma_c_sync;

    if (gamma_m < gamma_c_trans) {
        //if (gamma_m < gamma_c_sync) {
        gamma_c_new = 1;
    }
    gamma_c_new = gamma_c_last; //initial guess

    Real Y_T = compute_Thomson_Y(rad, gamma_m, gamma_c_new);
    Ys = InverseComptonY(gamma_m, gamma_c_new, rad.p, B, Y_T, true);
    size_t max_iter = 100;
    size_t iter = 0;
    do {
        gamma_c = gamma_c_new;
        Y_T = compute_Thomson_Y(rad, gamma_m, gamma_c);
        Ys.update_cooling_breaks(gamma_c, Y_T);
        Real Y_c = Ys.gamma_spectrum(gamma_c);
        gamma_c_new = compute_gamma_c(t_com, B, Y_c);
        iter++;
    } while (std::fabs((gamma_c_new - gamma_c) / gamma_c) > 1e-3 && iter < max_iter);
    gamma_c = gamma_c_new;
}

void update_gamma_M(Real& gamma_M, InverseComptonY const& Ys, Real p, Real B) {
    if (B == 0) {
        gamma_M = std::numeric_limits<Real>::infinity();
        return;
    }

    Real Y_M = Ys.gamma_spectrum(gamma_M);
    Real gamma_M_new = compute_syn_gamma_M(B, Y_M, p);

    while (std::fabs((gamma_M - gamma_M_new) / gamma_M_new) > 1e-3) {
        gamma_M = gamma_M_new;
        Y_M = Ys.gamma_spectrum(gamma_M);
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
