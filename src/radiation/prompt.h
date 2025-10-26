//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include "../core/mesh.h"
#include "../core/physics.h"

struct PromptPhotons {
    Real E_nu_peak{0};
    Real nu_0{0};
    Real alpha{0};

    [[nodiscard]] Real I_nu(Real nu) const;
};

using PromptPhotonsGrid = xt::xtensor<PromptPhotons, 3>;

class CoastingShock {
  public:
    CoastingShock(size_t phi_size, size_t theta_size, size_t t_size);
    CoastingShock() = delete;

    MeshGrid3d r;       // radius
    MeshGrid3d theta;   // theta for jet spreading
    MeshGrid3d Gamma;   // relative Lorentz factor between down stream and upstream
    MeshGrid3d epsilon; // relative energy per solid angle

    [[nodiscard]] auto shape() const {
        return std::make_tuple(phi_size, theta_size, t_size);
    } // Returns grid dimensions

  private:
    size_t const phi_size{0};   // Number of grid points in the phi direction
    size_t const theta_size{0}; // Number of grid points in the theta direction
    size_t const t_size{0};     // Number of grid points in the time direction
};

template <typename Ejecta>
CoastingShock gen_coasting_shock(Coord const& coord, Ejecta const& jet) {
    auto [phi_size, theta_size, t_size] = coord.shape();
    CoastingShock shock(1, theta_size, t_size);

    for (size_t j = 0; j < theta_size; ++j) {
        const Real Gamma = jet.Gamma0(coord.phi(0), coord.theta(j));
        const Real beta = gamma_to_beta(Gamma);
        const Real epsilon = jet.eps_k(coord.phi(0), coord.theta(j));
        for (size_t k = 0; k < t_size; ++k) {
            shock.Gamma(0, j, k) = Gamma;
            shock.epsilon(0, j, k) = epsilon;
            shock.r(0, j, k) = (beta * con::c) / std::fabs(1 - beta) * coord.t(k);
            shock.theta(0, j, k) = coord.theta(j);
        }
    }

    return shock;
}

PromptPhotonsGrid gen_prompt_photons(CoastingShock const& shock, Real R0, Real nu_0, Real alpha);
