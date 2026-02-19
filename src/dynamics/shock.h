//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <tuple>

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
    RadParams rad; ///< Radiation parameters

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
     * @brief Broadcasts computed dynamics to all grid points based on coord symmetry.
     * @details For isotropic, copies (0,0) to all (0,j), then phi=0 to all phi.
     *          For piecewise, copies each representative theta to its group, then phi=0 to all phi.
     *          For phi_symmetric, only broadcasts phi=0 to all phi (all thetas already computed).
     * @param coord Coordinate object containing symmetry and theta representatives.
     * <!-- ************************************************************************************** -->
     */
    void broadcast_groups(Coord const& coord);

  private:
    size_t phi_size{0};   ///< Number of grid points in phi direction
    size_t theta_size{0}; ///< Number of grid points in theta direction
    size_t t_size{0};     ///< Number of grid points in time direction
};
