//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "shock.h"

#include "shock-physics.h"

Shock::Shock(size_t phi_size, size_t theta_size, size_t t_size, RadParams const& rad_params)
    : t_comv({phi_size, theta_size, t_size}, 0),     // Initialize comoving time grid with 0
      r({phi_size, theta_size, t_size}, 0),          // Initialize radius grid with 0
      theta({phi_size, theta_size, t_size}, 0),      // Initialize theta grid with 0
      Gamma({phi_size, theta_size, t_size}, 1),      // Initialize Gamma grid with 1
      Gamma_th({phi_size, theta_size, t_size}, 1),   // Initialize Gamma_th grid with 1
      B({phi_size, theta_size, t_size}, 0),          // Initialize magnetic field grid with 0
      N_p({phi_size, theta_size, t_size}, 0),        // Initialize column density grid with 0
      injection_idx({phi_size, theta_size}, t_size), // Initialize a cross-index grid with t_size
      // required({phi_size, theta_size, t_size}, 1),    // Initialize the required grid with all-true
      rad(rad_params),        // Set radiation parameters
      phi_size(phi_size),     // Store phi grid size
      theta_size(theta_size), // Store theta grid size
      t_size(t_size) {}

void Shock::resize(size_t phi_size, size_t theta_size, size_t t_size) {
    this->phi_size = phi_size;
    this->theta_size = theta_size;
    this->t_size = t_size;
    t_comv.resize({phi_size, theta_size, t_size});
    r.resize({phi_size, theta_size, t_size});
    theta.resize({phi_size, theta_size, t_size});
    Gamma.resize({phi_size, theta_size, t_size});
    Gamma_th.resize({phi_size, theta_size, t_size});
    B.resize({phi_size, theta_size, t_size});
    N_p.resize({phi_size, theta_size, t_size});
    injection_idx.resize({phi_size, theta_size});
    injection_idx.fill(t_size);
    // required.resize({phi_size, theta_size, t_size});
    // required.fill(1);
}

void Shock::broadcast_groups(Coord const& coord) {
    if (coord.symmetry == Symmetry::isotropic) {
        for (size_t j = 1; j < theta_size; ++j) {
            injection_idx(0, j) = injection_idx(0, 0);
            for (size_t k = 0; k < t_size; ++k) {
                t_comv(0, j, k) = t_comv(0, 0, k);
                this->r(0, j, k) = this->r(0, 0, k);
                theta(0, j, k) = coord.theta(j);
                Gamma(0, j, k) = Gamma(0, 0, k);
                Gamma_th(0, j, k) = Gamma_th(0, 0, k);
                B(0, j, k) = B(0, 0, k);
                N_p(0, j, k) = N_p(0, 0, k);
            }
        }
    } else if (coord.symmetry == Symmetry::piecewise) {
        for (size_t r = 0; r < coord.theta_reps.size(); ++r) {
            const size_t j_start = coord.theta_reps[r];
            const size_t j_end = (r + 1 < coord.theta_reps.size()) ? coord.theta_reps[r + 1] : theta_size;
            for (size_t j = j_start + 1; j < j_end; ++j) {
                injection_idx(0, j) = injection_idx(0, j_start);
                for (size_t k = 0; k < t_size; ++k) {
                    t_comv(0, j, k) = t_comv(0, j_start, k);
                    this->r(0, j, k) = this->r(0, j_start, k);
                    theta(0, j, k) = coord.theta(j);
                    Gamma(0, j, k) = Gamma(0, j_start, k);
                    Gamma_th(0, j, k) = Gamma_th(0, j_start, k);
                    B(0, j, k) = B(0, j_start, k);
                    N_p(0, j, k) = N_p(0, j_start, k);
                }
            }
        }
    }

    for (size_t i = 1; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            injection_idx(i, j) = injection_idx(0, j);
            for (size_t k = 0; k < t_size; ++k) {
                t_comv(i, j, k) = t_comv(0, j, k);
                this->r(i, j, k) = this->r(0, j, k);
                theta(i, j, k) = coord.theta(j);
                Gamma(i, j, k) = Gamma(0, j, k);
                Gamma_th(i, j, k) = Gamma_th(0, j, k);
                B(i, j, k) = B(0, j, k);
                N_p(i, j, k) = N_p(0, j, k);
            }
        }
    }
}

Real compute_downstr_4vel(Real gamma_rel, Real sigma) noexcept {
    const Real ad_idx = physics::thermo::adiabatic_idx(gamma_rel);
    const Real gamma_m_1 = gamma_rel - 1; // (gamma_rel - 1)
    const Real ad_idx_m_2 = ad_idx - 2;   // (ad_idx - 2)
    const Real ad_idx_m_1 = ad_idx - 1;   // (ad_idx - 1)
    if (sigma <= con::sigma_cut) {
        return std::sqrt(std::fabs(gamma_m_1 * ad_idx_m_1 * ad_idx_m_1 / (-ad_idx * ad_idx_m_2 * gamma_m_1 + 2)));
    } else {
        const Real gamma_sq = gamma_rel * gamma_rel; // gamma_rel^2
        const Real gamma_p_1 = gamma_rel + 1;        // (gamma_rel + 1)

        // Precompute common terms
        const Real term1 = -ad_idx * ad_idx_m_2;
        const Real term2 = gamma_sq - 1;

        // Compute coefficients
        const Real A = term1 * gamma_m_1 + 2;
        const Real B = -gamma_p_1 * (-ad_idx_m_2 * (ad_idx * gamma_sq + 1) + ad_idx * ad_idx_m_1 * gamma_rel) * sigma -
                       gamma_m_1 * (term1 * (gamma_sq - 2) + 2 * gamma_rel + 3);
        const Real C = gamma_p_1 * (ad_idx * (1 - ad_idx / 4) * term2 + 1) * sigma * sigma +
                       term2 * (2 * gamma_rel + ad_idx_m_2 * (ad_idx * gamma_rel - 1)) * sigma +
                       gamma_p_1 * gamma_m_1 * gamma_m_1 * ad_idx_m_1 * ad_idx_m_1;
        const Real D = -gamma_m_1 * gamma_p_1 * gamma_p_1 * ad_idx_m_2 * ad_idx_m_2 * sigma * sigma / 4;

        const Real b = B / A;
        const Real c = C / A;
        const Real d = D / A;
        const Real P = c - b * b / 3;
        const Real Q = 2 * b * b * b / 27 - b * c / 3 + d;
        const Real u = std::sqrt(-P / 3);
        const Real v = std::clamp(3 * Q / (2 * P * u), -1.0, 1.0);
        const Real uds = 2 * u * std::cos((std::acos(v) - 2 * con::pi) / 3) - b / 3;
        return std::sqrt(uds);
    }
}
