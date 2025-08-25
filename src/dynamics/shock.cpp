//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "shock.h"

Shock::Shock(size_t phi_size, size_t theta_size, size_t t_size, RadParams const& rad_params)
    : t_comv({phi_size, theta_size, t_size}, 0),      // Initialize comoving time grid with 0
      r({phi_size, theta_size, t_size}, 0),           // Initialize radius grid with 0
      theta({phi_size, theta_size, t_size}, 0),       // Initialize theta grid with 0
      Gamma({phi_size, theta_size, t_size}, 1),       // Initialize Gamma grid with 1
      Gamma_th({phi_size, theta_size, t_size}, 1),    // Initialize Gamma_th grid with 1
      B({phi_size, theta_size, t_size}, 0),           // Initialize magnetic field grid with 0
      N_p({phi_size, theta_size, t_size}, 0),         // Initialize column density grid with 0
      injection_idx({phi_size, theta_size}, t_size),  // Initialize cross index grid with t_size
      required({phi_size, theta_size, t_size}, 1),    // Initialize required grid with true
      rad(rad_params),                                // Set radiation parameters
      phi_size(phi_size),                             // Store phi grid size
      theta_size(theta_size),                         // Store theta grid size
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
    Real ad_idx = adiabatic_idx(gamma_rel);
    Real gamma_m_1 = gamma_rel - 1;  // (gamma_rel - 1)
    Real ad_idx_m_2 = ad_idx - 2;    // (ad_idx - 2)
    Real ad_idx_m_1 = ad_idx - 1;    // (ad_idx - 1)
    if (std::abs(sigma) <= 1e-6) {
        return std::sqrt(std::fabs(gamma_m_1 * ad_idx_m_1 * ad_idx_m_1 / (-ad_idx * ad_idx_m_2 * gamma_m_1 + 2)));
    } else {
        Real gamma_sq = gamma_rel * gamma_rel;  // gamma_rel^2
        Real gamma_p_1 = gamma_rel + 1;         // (gamma_rel + 1)

        // Precompute common terms
        Real term1 = -ad_idx * ad_idx_m_2;
        Real term2 = gamma_sq - 1;

        // Compute coefficients
        Real A = term1 * gamma_m_1 + 2;
        Real B = -gamma_p_1 * (-ad_idx_m_2 * (ad_idx * gamma_sq + 1) + ad_idx * ad_idx_m_1 * gamma_rel) * sigma -
                 gamma_m_1 * (term1 * (gamma_sq - 2) + 2 * gamma_rel + 3);
        Real C = gamma_p_1 * (ad_idx * (1 - ad_idx / 4) * term2 + 1) * sigma * sigma +
                 term2 * (2 * gamma_rel + ad_idx_m_2 * (ad_idx * gamma_rel - 1)) * sigma +
                 gamma_p_1 * gamma_m_1 * gamma_m_1 * ad_idx_m_1 * ad_idx_m_1;
        Real D = -gamma_m_1 * gamma_p_1 * gamma_p_1 * ad_idx_m_2 * ad_idx_m_2 * sigma * sigma / 4;

        Real b = B / A;
        Real c = C / A;
        Real d = D / A;
        Real P = c - b * b / 3;
        Real Q = 2 * b * b * b / 27 - b * c / 3 + d;
        Real u = std::sqrt(-P / 3);
        Real v = std::clamp(3 * Q / (2 * P * u), -1.0, 1.0);
        Real uds = 2 * u * std::cos((std::acos(v) - 2 * con::pi) / 3) - b / 3;
        return std::sqrt(uds);
    }
}

void save_shock_state(Shock& shock, size_t i, size_t j, size_t k, Real t_comv, Real r, Real theta, Real Gamma,
                      Real Gamma_th, Real B, Real mass) {
    shock.t_comv(i, j, k) = t_comv;
    shock.r(i, j, k) = r;
    shock.theta(i, j, k) = theta;
    shock.Gamma(i, j, k) = Gamma;
    shock.Gamma_th(i, j, k) = Gamma_th;
    shock.B(i, j, k) = B;
    shock.N_p(i, j, k) = mass / con::mp;
}

Real compute_4vel_jump(Real gamma_rel, Real sigma_upstr) {
    Real u_down_s_ = compute_downstr_4vel(gamma_rel, sigma_upstr);
    Real u_up_s_ = compute_upstr_4vel(u_down_s_, gamma_rel);
    Real ratio_u = u_up_s_ / u_down_s_;
    if (u_down_s_ == 0.) {
        ratio_u = 4 * gamma_rel;  // (g_hat*gamma_rel+1)/(g_hat-1)
    }
    return ratio_u;
}

Real compute_compression(Real Gamma_upstr, Real Gamma_downstr, Real sigma_upstr) {
    Real Gamma_rel = compute_rel_Gamma(Gamma_upstr, Gamma_downstr);
    return compute_4vel_jump(Gamma_rel, sigma_upstr);
}

Real compute_downstr_B(Real eps_B, Real rho_upstr, Real B_upstr, Real Gamma_th, Real comp_ratio) {
    Real rho_downstr = rho_upstr * comp_ratio;

    Real e_th = (Gamma_th - 1) * rho_downstr * con::c2;

    return compute_comv_weibel_B(eps_B, e_th) + B_upstr * comp_ratio;
}
