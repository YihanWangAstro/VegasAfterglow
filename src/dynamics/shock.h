//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <numeric>
#include <tuple>
#include <vector>

#include "../core/mesh.h"
#include "../core/physics.h"
#include "../environment/jet.h"
#include "../environment/medium.h"

/**
 * <!-- ************************************************************************************** -->
 * @class Shock
 * @brief Represents a shock wave in an astrophysical environment.
 * @details The class stores physical properties of the shock across a 3D grid defined by azimuthal angle (phi),
 *          polar angle (theta), and time bins. Provides methods for shock calculations, including relativistic
 *          jump conditions, magnetic field calculations, and energy density computations.
 * <!-- ************************************************************************************** -->
 */
class Shock {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Constructs a Shock object with the given grid dimensions and energy fractions.
     * @details Initializes various 3D grids for storing physical properties of the shock, including comoving time,
     *          radius, Lorentz factors, magnetic fields, and downstream densities.
     * @param phi_size Number of grid points in the phi direction
     * @param theta_size Number of grid points in theta direction
     * @param t_size Number of grid points in time direction
     * @param rad_params Radiation parameters
     * <!-- ************************************************************************************** -->
     */
    Shock(size_t phi_size, size_t theta_size, size_t t_size, RadParams const& rad_params);

    Shock() noexcept = default;

    MeshGrid3d t_comv;       ///< Comoving time
    MeshGrid3d r;            ///< Radius
    MeshGrid3d theta;        ///< Theta for jet spreading
    MeshGrid3d Gamma;        ///< Bulk Lorentz factor
    MeshGrid3d Gamma_th;     ///< Downstream internal Lorentz factor
    MeshGrid3d B;            ///< Comoving magnetic field
    MeshGrid3d N_p;          ///< Downstream proton number per solid angle
    IndexGrid injection_idx; ///< Beyond which grid index there is no electron injection
    // MaskGrid required;        ///< Grid points actually required for final flux calculation
    RadParams rad;                           ///< Radiation parameters
    Symmetry symmetry{Symmetry::structured}; ///< Auto-detected symmetry level
    std::vector<size_t> theta_reps;          ///< Representative theta indices (contiguous groups)

    /// Returns grid dimensions as a tuple
    [[nodiscard]] auto shape() const { return std::make_tuple(phi_size, theta_size, t_size); }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Resizes all grid components of the Shock object to new dimensions.
     * @details Also resets injection_idx to t_size for all grid points.
     * @param phi_size New number of grid points in the phi direction
     * @param theta_size New number of grid points in theta direction
     * @param t_size New number of grid points in time direction
     * <!-- ************************************************************************************** -->
     */
    void resize(size_t phi_size, size_t theta_size, size_t t_size);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Auto-detect symmetry by comparing jet initial conditions at consecutive theta grid points.
     * @details Sets symmetry and theta_reps. Falls back to structured if jet is spreading or medium is anisotropic.
     * <!-- ************************************************************************************** -->
     */
    template <typename Ejecta, typename Medium>
    void detect_symmetry(Coord const& coord, Ejecta const& jet, Medium const& medium);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Broadcasts computed dynamics to all grid points based on detected symmetry.
     * @details For isotropic, copies (0,0) to all (0,j), then phi=0 to all phi.
     *          For piecewise, copies each representative theta to its group, then phi=0 to all phi.
     *          For phi_symmetric, only broadcasts phi=0 to all phi (all thetas already computed).
     * @param theta_coords 1D array of theta coordinate values.
     * <!-- ************************************************************************************** -->
     */
    void broadcast_groups(Array const& theta_coords);

  private:
    size_t phi_size{0};   ///< Number of grid points in phi direction
    size_t theta_size{0}; ///< Number of grid points in theta direction
    size_t t_size{0};     ///< Number of grid points in time direction
};

template <typename Ejecta, typename Medium>
void Shock::detect_symmetry(Coord const& coord, Ejecta const& jet, Medium const& medium) {
    auto [phi_size_, theta_size_, t_size_] = coord.shape();

    if (jet.spreading || !medium.isotropic) {
        symmetry = Symmetry::structured;
        theta_reps.resize(theta_size_);
        std::iota(theta_reps.begin(), theta_reps.end(), size_t(0));
        return;
    }

    const Real phi0 = coord.phi(0);
    theta_reps.clear();
    theta_reps.reserve(theta_size_);
    theta_reps.push_back(0);

    auto jet_ic_differs = [&](size_t ja, size_t jb) {
        const Real theta_a = coord.theta(ja), theta_b = coord.theta(jb);
        if (jet.eps_k(phi0, theta_a) != jet.eps_k(phi0, theta_b))
            return true;
        if (jet.Gamma0(phi0, theta_a) != jet.Gamma0(phi0, theta_b))
            return true;
        if constexpr (HasSigma<Ejecta>) {
            if (jet.sigma0(phi0, theta_a) != jet.sigma0(phi0, theta_b))
                return true;
        }
        if constexpr (HasDedt<Ejecta>) {
            for (size_t k = 0; k < t_size_; ++k) {
                const Real ta = coord.t(0, ja, k), tb = coord.t(0, jb, k);
                if (jet.deps_dt(phi0, theta_a, ta) != jet.deps_dt(phi0, theta_b, ta))
                    return true;
                if (jet.deps_dt(phi0, theta_a, tb) != jet.deps_dt(phi0, theta_b, tb))
                    return true;
            }
        }
        if constexpr (HasDmdt<Ejecta>) {
            for (size_t k = 0; k < t_size_; ++k) {
                const Real ta = coord.t(0, ja, k), tb = coord.t(0, jb, k);
                if (jet.dm_dt(phi0, theta_a, ta) != jet.dm_dt(phi0, theta_b, ta))
                    return true;
                if (jet.dm_dt(phi0, theta_a, tb) != jet.dm_dt(phi0, theta_b, tb))
                    return true;
            }
        }
        return false;
    };

    for (size_t j = 1; j < theta_size_; ++j) {
        if (jet_ic_differs(j - 1, j)) {
            theta_reps.push_back(j);
        }
    }

    if (theta_reps.size() == 1)
        symmetry = Symmetry::isotropic;
    else if (theta_reps.size() < theta_size_)
        symmetry = Symmetry::piecewise;
    else
        symmetry = Symmetry::phi_symmetric;
}
