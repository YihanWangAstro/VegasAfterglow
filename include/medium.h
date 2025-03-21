//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _MEDIUM_
#define _MEDIUM_

#include <iostream>

#include "macros.h"
#include "mesh.h"
#include "utilities.h"

/********************************************************************************************************************
 * CLASS: Medium
 * DESCRIPTION: Represents the interstellar medium (ISM) or any surrounding medium.
 *              The medium is characterized by a density profile with index k, a characteristic number density n_c,
 *              and a characteristic radius r_c. The class provides methods to compute the mass density
 *              (rho) and the enclosed mass up to radius r.
 ********************************************************************************************************************/
class Medium {
   public:
    TernaryFunc rho{func::zero};  // density(phi, theta, r)
};

namespace evn {
    inline auto ISM(Real n_ism) {
        return [n_ism](Real phi, Real theta, Real r) { return n_ism * con::mp; };
    }

    inline auto wind(Real A_star) {
        Real A = A_star * 5e11 * con::g / con::cm;
        return [A](Real phi, Real theta, Real r) { return A / (r * r); };
    }
}  // namespace evn

#endif