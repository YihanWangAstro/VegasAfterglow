//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "prompt.h"

CoastingShock::CoastingShock(size_t phi_size, size_t theta_size, size_t t_size)
    : r({phi_size, theta_size, t_size}, 0),       // Initialize engine time grid with 0
      theta({phi_size, theta_size, t_size}, 0),   // Initialize theta grid with 0
      Gamma({phi_size, theta_size, t_size}, 1),   // Initialize Gamma_rel grid with 1
      epsilon({phi_size, theta_size, t_size}, 0), // Initialize column density grid with 0
      phi_size(phi_size),                         // Store phi grid size
      theta_size(theta_size),                     // Store theta grid size
      t_size(t_size) {}

PromptPhotonsGrid gen_prompt_photons(CoastingShock const& shock, Real R0, Real nu_0, Real alpha) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    PromptPhotonsGrid ph({phi_size, theta_size, t_size});

    const Real Gamma_c = shock.Gamma(0, 0, 0);
    const Real beta_c = gamma_to_beta(Gamma_c);

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            const Real Gamma = shock.Gamma(i, j, 0);
            const Real beta = gamma_to_beta(Gamma);
            const Real R = R0 * beta / (1 - beta) * (1 - beta_c) / beta_c;
            const Real Rmin = R * 0;
            const Real Rmax = R * 1.;
            for (size_t k = 0; k < t_size; ++k) {
                if (shock.r(i, j, k) >= Rmin && shock.r(i, j, k) <= Rmax) {
                    ph(i, j, k).E_nu_peak = shock.epsilon(i, j, k);
                    ph(i, j, k).nu_0 = nu_0;
                    ph(i, j, k).alpha = alpha;
                }
            }
        }
    }
    return ph;
}

Real PromptPhotons::I_nu(Real nu) const {
    return E_nu_peak * std::pow(nu / nu_0, -alpha);
}
