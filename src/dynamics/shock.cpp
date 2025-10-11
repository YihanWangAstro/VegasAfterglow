//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "shock.h"

Shock::Shock(size_t phi_size, size_t theta_size, size_t t_size, RadParams const& rad_params)
    : t_comv({phi_size, theta_size, t_size}, 0),     // Initialize comoving time grid with 0
      r({phi_size, theta_size, t_size}, 0),          // Initialize radius grid with 0
      theta({phi_size, theta_size, t_size}, 0),      // Initialize theta grid with 0
      Gamma({phi_size, theta_size, t_size}, 1),      // Initialize Gamma grid with 1
      Gamma_th({phi_size, theta_size, t_size}, 1),   // Initialize Gamma_th grid with 1
      B({phi_size, theta_size, t_size}, 0),          // Initialize magnetic field grid with 0
      N_p({phi_size, theta_size, t_size}, 0),        // Initialize column density grid with 0
      injection_idx({phi_size, theta_size}, t_size), // Initialize a cross-index grid with t_size
      required({phi_size, theta_size, t_size}, 1),   // Initialize the required grid with true
      rad(rad_params),                               // Set radiation parameters
      phi_size(phi_size),                            // Store phi grid size
      theta_size(theta_size),                        // Store theta grid size
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
    required.resize({phi_size, theta_size, t_size});
    required.fill(1);
}

Real compute_downstr_4vel(Real gamma_rel, Real sigma) {
    const Real ad_idx = adiabatic_idx(gamma_rel);
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
