//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "power-law-syn.h"

#include "radiation/inverse-compton.h"

//========================================================================================================
//                                  PowerLawSyn Class Methods
//========================================================================================================

Real PowerLawSyn::compute_spectrum(Real nu) const {
    switch (regime) {
        case 1:
            if (nu <= nu_a) {
                return C1_ * (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_m) {
                return std::cbrt(nu / nu_m);
            }
            if (nu <= nu_c) {
                return fast_pow(nu / nu_m, -(p - 1) / 2);
            }

            return C2_ * fast_pow(nu / nu_c, -p / 2);

            break;
        case 2:
            if (nu <= nu_m) {
                return C1_ * (nu / nu_m) * (nu / nu_m);
            }
            if (nu <= nu_a) {
                return C2_ * pow52(nu / nu_a); // Using pow52 for (nu / nu_a)^(5/2)
            }
            if (nu <= nu_c) {
                return fast_pow(nu / nu_m, -(p - 1) / 2);
            }

            return C3_ * fast_pow(nu / nu_c, -p / 2);

            break;
        case 3:
            if (nu <= nu_a) {
                return C1_ * (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_c) {
                return std::cbrt(nu / nu_c);
            }
            if (nu <= nu_m) {
                return std::sqrt(nu_c / nu);
            }
            return C2_ * fast_pow(nu / nu_m, -p / 2);

            break;
        case 4:
            if (nu <= nu_a) {
                return 3 * C2_ * (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_m) {
                return 3 * C2_ * std::sqrt(nu_a / nu);
            }
            return 3 * C2_ * C1_ * fast_pow(nu / nu_m, -p / 2);

            break;
        case 5:
        case 6:
            if (nu <= nu_m) {
                return C3_ * C2_ * C1_ * (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_a) {
                return C3_ * C2_ * pow52(nu / nu_a);
            }
            return C3_ * C2_ * fast_pow(nu / nu_a, -p / 2);

            break;

        default:
            return 0;
            break;
    }
}

Real PowerLawSyn::compute_log2_spectrum(Real log2_nu) const {
    constexpr Real log2_3 = 1.5849625007; // log2(3)
    switch (regime) {
        case 1:
            if (log2_nu <= log2_nu_a) {
                return log2_C1_ + 2. * log2_nu;
            }
            if (log2_nu <= log2_nu_m) {
                return log2_C2_ + log2_nu / 3.;
            }
            if (log2_nu <= log2_nu_c) {
                return log2_C3_ - (p - 1.) / 2. * log2_nu;
            }
            return log2_C4_ - p / 2. * log2_nu;

            break;
        case 2:
            if (log2_nu <= log2_nu_m) {
                return log2_C1_ + 2. * log2_nu;
            }
            if (log2_nu <= log2_nu_a) {
                return log2_C2_ + 2.5 * log2_nu;
            }
            if (log2_nu <= log2_nu_c) {
                return log2_C3_ - (p - 1.) / 2. * log2_nu;
            }

            return log2_C4_ - p / 2. * log2_nu;

            break;
        case 3:
            if (log2_nu <= log2_nu_a) {
                return log2_C1_ + 2. * log2_nu;
            }
            if (log2_nu <= log2_nu_c) {
                return log2_C2_ + log2_nu / 3.;
            }
            if (log2_nu <= log2_nu_m) {
                return log2_C3_ - log2_nu / 2.;
            }

            return log2_C4_ - p / 2. * log2_nu;

            break;
        case 4:

            if (log2_nu <= log2_nu_a) {
                return log2_3 + log2_C4_ + log2_C1_ + 2. * log2_nu;
            }
            if (log2_nu <= log2_nu_m) {
                return log2_3 + log2_C2_ - log2_nu / 2.;
            }

            return log2_3 + log2_C3_ - p / 2. * log2_nu;

            break;
        case 5:
        case 6:
            if (log2_nu <= log2_nu_m) {
                return log2_C3_ + log2_C4_ + 0.5 * log2_nu_m + log2_C1_ + 2. * log2_nu;
            }
            if (log2_nu <= log2_nu_a) {
                return log2_C3_ + log2_C4_ + log2_C1_ + 2.5 * log2_nu;
            }

            return log2_C3_ + log2_C2_ - p / 2. * log2_nu;

            break;
        default:
            return -con::inf;
            break;
    }
}

void PowerLawSyn::update_constant() {
    // Update constants based on spectral parameters
    if (regime == 1) {
        // a_m_1_3 = std::cbrt(nu_a / nu_m);  // (nu_a / nu_m)^(1/3)
        // c_m_mpa1_2 = fastPow(nu_c / nu_m, (-p + 1) / 2);  // (nu_c / nu_m)^((-p+1)/2)
        C1_ = std::cbrt(nu_a / nu_m);
        C2_ = fast_pow(nu_c / nu_m, (-p + 1) / 2);

        log2_C1_ = (log2_nu_a - log2_nu_m) / 3 - 2 * log2_nu_a;
        log2_C2_ = -log2_nu_m / 3;
        log2_C3_ = (p - 1) / 2 * log2_nu_m;
        log2_C4_ = (p - 1) / 2 * (log2_nu_m - log2_nu_c) + p / 2 * log2_nu_c;
    } else if (regime == 2) {
        // m_a_pa4_2 = fastPow(nu_m / nu_a, (p + 4) / 2);    // (nu_m / nu_a)^((p+4)/2)
        // a_m_mpa1_2 = fastPow(nu_a / nu_m, (-p + 1) / 2);  // (nu_a / nu_m)^((-p+1)/2)
        // c_m_mpa1_2 = fastPow(nu_c / nu_m, (-p + 1) / 2);  // (nu_c / nu_m)^((-p+1)/2)
        C1_ = fast_pow(nu_m / nu_a, (p + 4) / 2);
        C2_ = fast_pow(nu_a / nu_m, (-p + 1) / 2);
        C3_ = fast_pow(nu_c / nu_m, (-p + 1) / 2);

        log2_C1_ = (p + 4) / 2 * (log2_nu_m - log2_nu_a) - 2 * log2_nu_m;
        log2_C2_ = (p - 1) / 2 * (log2_nu_m - log2_nu_a) - 2.5 * log2_nu_a;
        log2_C3_ = (p - 1) / 2 * log2_nu_m;
        log2_C4_ = (p - 1) / 2 * (log2_nu_m - log2_nu_c) + p / 2 * log2_nu_c;
    } else if (regime == 3) {
        // a_c_1_3 = std::cbrt(nu_a / nu_c);  // (nu_a / nu_c)^(1/3)
        // c_m_1_2 = std::sqrt(nu_c / nu_m);  // (nu_c / nu_m)^(1/2)
        C1_ = std::cbrt(nu_a / nu_c);
        C2_ = std::sqrt(nu_c / nu_m);

        log2_C1_ = (log2_nu_a - log2_nu_c) / 3 - 2 * log2_nu_a;
        log2_C2_ = -log2_nu_c / 3;
        log2_C3_ = log2_nu_c / 2;
        log2_C4_ = (log2_nu_c - log2_nu_m) / 2 + p / 2 * log2_nu_m;
    } else if (regime == 4) {
        C1_ = std::sqrt(nu_a / nu_m);
        C3_ = 3;
        C2_ = std::sqrt(nu_c / nu_a) / C3_;

        log2_C4_ = fast_log2(C2_);

        log2_C1_ = -2 * log2_nu_a;
        log2_C2_ = log2_C4_ + log2_nu_a / 2;
        log2_C3_ = log2_C4_ + (log2_nu_a - log2_nu_m) / 2 + p / 2 * log2_nu_m;

    } else if (regime == 5) {
        C1_ = std::sqrt(nu_m / nu_a);
        C3_ = 3 / (p - 1);
        C2_ = std::sqrt(nu_c / nu_a) * fast_pow(nu_m / nu_a, (p - 1) / 2) / C3_;

        log2_C4_ = fast_log2(C2_);

        log2_C1_ = -2.5 * log2_nu_a;
        log2_C2_ = log2_C4_ + p / 2 * log2_nu_a;

        log2_C3_ = fast_log2(C3_);
    } else if (regime == 6) {
        C1_ = std::sqrt(nu_m / nu_a);
        C3_ = 3;
        C2_ = std::sqrt(nu_c / nu_a) * fast_pow(nu_m / nu_a, (p - 1) / 2) / C3_;

        log2_C4_ = fast_log2(C2_);

        log2_C1_ = -2.5 * log2_nu_a;
        log2_C2_ = log2_C4_ + p / 2 * log2_nu_a;

        log2_C3_ = 1.5849625007; // log2(3)
    }
}

Real PowerLawSyn::compute_I_nu(Real nu) const {
    if (nu <= nu_c) { // Below cooling frequency, simple scaling
        return fast_exp(-nu / nu_M) * I_nu_max * compute_spectrum(nu);
    } else {
        return fast_exp(-nu / nu_M) * I_nu_max * compute_spectrum(nu) * inverse_compton_correction(*this, nu);
    }
}

Real PowerLawSyn::compute_log2_I_nu(Real log2_nu) const {
    if (log2_nu <= log2_nu_c) { // Below cooling frequency, simple scaling
        return log2_I_nu_max + compute_log2_spectrum(log2_nu) - 1.442695 * fast_exp2(log2_nu) / nu_M;
    } else {
        const Real nu = fast_exp2(log2_nu);
        const Real cooling_factor = inverse_compton_correction(*this, nu);
        return log2_I_nu_max + compute_log2_spectrum(log2_nu) + fast_log2(cooling_factor) - 1.442695 * nu / nu_M;
    }
}
