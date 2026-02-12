//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "power-law-syn.h"

#include "inverse-compton.h"

//========================================================================================================
//                                  PowerLawSyn Class Methods
//========================================================================================================

Real PowerLawSyn::compute_spectrum(Real nu) const {
    return segments_.eval(nu);
}

Real PowerLawSyn::compute_log2_spectrum(Real log2_nu) const {
    return segments_.log2_eval(log2_nu);
}

void PowerLawSyn::build() {
    log2_I_nu_max = fast_log2(I_nu_max);

    inv_nu_M_ = 1.0 / nu_M;

    switch (regime) {
        case 1: // a ≤ m ≤ c
            segments_.first_segment(std::cbrt(nu_a / nu_m), nu_a, 2.0);
            segments_.add_segment(nu_a, 1.0 / 3);
            segments_.add_segment(nu_m, -(p - 1) / 2);
            segments_.add_segment(nu_c, -p / 2);
            break;

        case 2: // m ≤ a ≤ c
            segments_.first_segment(fast_pow(nu_m / nu_a, (p + 4) / 2), nu_m, 2.0);
            segments_.add_segment(nu_m, 2.5);
            segments_.add_segment(nu_a, -(p - 1) / 2);
            segments_.add_segment(nu_c, -p / 2);
            break;

        case 3: // a ≤ c ≤ m
            segments_.first_segment(std::cbrt(nu_a / nu_c), nu_a, 2.0);
            segments_.add_segment(nu_a, 1.0 / 3);
            segments_.add_segment(nu_c, -0.5);
            segments_.add_segment(nu_m, -p / 2);
            break;

        case 4: // c ≤ a ≤ m
            segments_.first_segment(std::sqrt(nu_c / nu_a), nu_a, 2.0);
            segments_.add_segment(nu_a, -0.5);
            segments_.add_segment(nu_m, -p / 2);
            break;

        case 5: // m ≤ c ≤ a
        case 6: // c ≤ m ≤ a
            segments_.first_segment(std::sqrt(nu_c / nu_a) * fast_pow(nu_m / nu_a, (p + 4) / 2), nu_m, 2.0);
            segments_.add_segment(nu_m, 2.5);
            segments_.add_segment(nu_a, -p / 2);
            break;

        default:
            segments_.clear();
            break;
    }
}

Real PowerLawSyn::compute_I_nu(Real nu) const {
    if (nu <= nu_c) { // Below cooling frequency, simple scaling
        return fast_exp(-nu * inv_nu_M_) * I_nu_max * compute_spectrum(nu);
    } else {
        return fast_exp(-nu * inv_nu_M_) * I_nu_max * compute_spectrum(nu) * inverse_compton_correction(*this, nu);
    }
}

Real PowerLawSyn::compute_log2_I_nu(Real log2_nu) const {
    constexpr Real ln2 = 1.442695040888963407359924681001892137;
    if (log2_nu <= log2_nu_c) { // Below cooling frequency, simple scaling
        return log2_I_nu_max + compute_log2_spectrum(log2_nu) - ln2 * fast_exp2(log2_nu) * inv_nu_M_;
    } else {
        const Real nu = fast_exp2(log2_nu);
        return log2_I_nu_max + compute_log2_spectrum(log2_nu) + fast_log2(inverse_compton_correction(*this, nu)) -
               ln2 * nu * inv_nu_M_;
    }
}
