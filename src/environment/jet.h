//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <cmath>
#include <utility>

#include "../util/utilities.h"
/**
 * <!-- ************************************************************************************** -->
 * @class Ejecta
 * @brief Represents generic ejecta properties for a simulation.
 * @details Uses ternary functions (TernaryFunc) to accept user-defined ejecta that describes various quantities
 *          as functions of phi, theta, and time. This class encapsulates all the properties of the material
 *          ejected in a gamma-ray burst, including its energy, magnetization, and Lorentz factor profiles.
 * <!-- ************************************************************************************** -->
 */
class Ejecta {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Constructor: Initialize with core angle, isotropic energy, and initial Lorentz factor
     * @param eps_k Energy per unit solid angle as a function of (phi, theta)
     * @param Gamma0 Lorentz factor as a function of (phi, theta)
     * @param sigma0 Magnetization parameter as a function of (phi, theta)
     * @param deps_dt Energy injection rate per solid angle as a function of (phi, theta, t)
     * @param dm_dt Mass injection rate per unit solid angle as a function of (phi, theta, t)
     * @param T0 Duration of the ejecta
     * @param spreading Flag indicating if the ejecta spreads laterally during evolution
     * <!-- ************************************************************************************** -->
     */
    Ejecta(BinaryFunc eps_k, BinaryFunc Gamma0, BinaryFunc sigma0 = func::zero_2d, TernaryFunc deps_dt = func::zero_3d,
           TernaryFunc dm_dt = func::zero_3d, bool spreading = false, Real T0 = 1 * unit::sec) noexcept
        : eps_k(std::move(eps_k)),
          Gamma0(std::move(Gamma0)),
          sigma0(std::move(sigma0)),
          deps_dt(std::move(deps_dt)),
          dm_dt(std::move(dm_dt)),
          T0(T0),
          spreading(spreading) {}

    Ejecta() = default;

    /// Initial energy per unit solid angle as a function of (phi, theta)
    BinaryFunc eps_k{func::zero_2d};

    /// Lorentz factor profile in the ejecta as a function of (phi, theta)
    /// Default is uniform (one) across all angles
    BinaryFunc Gamma0{func::one_2d};

    /// Initial magnetization parameter as a function of (phi, theta)
    BinaryFunc sigma0{func::zero_2d};

    /// Energy injection rate per solid angle as a function of (phi, theta, t)
    /// Default is no energy injection (zero)
    TernaryFunc deps_dt{func::zero_3d};

    /// Mass injection rate per unit solid angle as a function of (phi, theta, t)
    /// Default is no mass injection (zero)
    TernaryFunc dm_dt{func::zero_3d};

    /// Duration of the ejecta in seconds
    Real T0{1 * unit::sec};

    /// Flag indicating if the ejecta spreads laterally during evolution
    bool spreading{false};
};

/**
 * <!-- ************************************************************************************** -->
 * @class TophatJet
 * @brief Implements a tophat jet profile where properties are uniform within a core angle theta_c and zero outside.
 * @details This class provides a simple model for GRB jets with sharp edges, characterized by
 *          isotropic equivalent energy E_iso and initial Lorentz factor Gamma0 within the core angle.
 * <!-- ************************************************************************************** -->
 */
class TophatJet {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Constructor: Initialize with core angle, isotropic energy, and initial Lorentz factor
     * @param theta_c Core angle of the jet
     * @param E_iso Isotropic equivalent energy
     * @param Gamma0 Initial Lorentz factor
     * @param T0 Duration of the ejecta
     * @param spreading Flag indicating if the ejecta spreads laterally during evolution
     * <!-- ************************************************************************************** -->
     */
    TophatJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading = false, Real T0 = 1 * unit::sec) noexcept
        : T0(T0), spreading(spreading), theta_c_(theta_c), eps_k_(E_iso / (4 * con::pi)), Gamma0_(Gamma0) {}

    /**
     * <!-- ************************************************************************************** -->
     * @brief Energy per solid angle as a function of phi and theta
     * @param phi Azimuthal angle (unused)
     * @param theta Polar angle
     * @return Energy per solid angle (eps_k_ if theta < theta_c_, 0 otherwise)
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] inline Real eps_k(Real phi, Real theta) const noexcept { return theta < theta_c_ ? eps_k_ : 0; }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Initial Lorentz factor as a function of phi and theta
     * @param phi Azimuthal angle (unused)
     * @param theta Polar angle
     * @return Lorentz factor (Gamma0_ if theta < theta_c_, 1 otherwise)
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] inline Real Gamma0(Real phi, Real theta) const noexcept { return theta < theta_c_ ? Gamma0_ : 1; }

    Real T0{1 * unit::sec}; ///< Duration of the ejecta in seconds
    bool spreading{false};  ///< Flag indicating if the ejecta spreads laterally during evolution

  private:
    Real const theta_c_{0}; ///< Core angle of the jet
    Real const eps_k_{0};   ///< Energy per solid angle within the core
    Real const Gamma0_{1};  ///< Initial Lorentz factor within the core
};

/**
 * <!-- ************************************************************************************** -->
 * @class GaussianJet
 * @brief Implements a Gaussian jet profile where properties follow a Gaussian distribution with angles.
 * @details This class provides a smooth model for GRB jets, characterized by core angle theta_c,
 *          isotropic equivalent energy E_iso, and initial Lorentz factor Gamma0 at the center.
 * <!-- ************************************************************************************** -->
 */
class GaussianJet {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Constructor: Initialize with core angle, isotropic energy, and initial Lorentz factor
     * @param theta_c Core angle of the jet
     * @param E_iso Isotropic equivalent energy
     * @param Gamma0 Initial Lorentz factor
     * @param T0 Duration of the ejecta
     * @param spreading Flag indicating if the ejecta spreads laterally during evolution
     * <!-- ************************************************************************************** -->
     */
    GaussianJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading = false, Real T0 = 1 * unit::sec) noexcept
        : T0(T0),
          spreading(spreading),
          norm_(-1 / (2 * theta_c * theta_c)),
          eps_k_(E_iso / (4 * con::pi)),
          Gamma0_(Gamma0) {}

    /**
     * <!-- ************************************************************************************** -->
     * @brief Energy per solid angle as a function of phi and theta, with Gaussian falloff
     * @param phi Azimuthal angle (unused)
     * @param theta Polar angle
     * @return Energy per solid angle with Gaussian angular dependence
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] inline Real eps_k(Real phi, Real theta) const noexcept {
        return eps_k_ * fast_exp(theta * theta * norm_);
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Initial Lorentz factor as a function of phi and theta, with Gaussian falloff
     * @param phi Azimuthal angle (unused)
     * @param theta Polar angle
     * @return Lorentz factor with Gaussian angular dependence
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] inline Real Gamma0(Real phi, Real theta) const noexcept {
        return (Gamma0_ - 1) * fast_exp(theta * theta * norm_) + 1;
    }

    /// Duration of the ejecta in seconds
    Real T0{1 * unit::sec};
    /// Flag indicating if the ejecta spreads laterally during evolution
    bool spreading{false};

  private:
    Real const norm_{0};   ///< Normalization factor for Gaussian distribution
    Real const eps_k_{0};  ///< Peak energy per solid angle at center
    Real const Gamma0_{1}; ///< Peak Lorentz factor at the center
};

/**
 * <!-- ************************************************************************************** -->
 * @class PowerLawJet
 * @brief Implements a power-law jet profile where properties follow a power-law distribution with angles.
 * @details This class provides a model for GRB jets with a power-law decay, characterized by core angle theta_c,
 *          isotropic equivalent energy E_iso, initial Lorentz factor Gamma0, and power-law index k.
 * <!-- ************************************************************************************** -->
 */
class PowerLawJet {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Constructor: Initialize with core angle, isotropic energy, initial Lorentz factor, and power-law index
     * @param theta_c Core angle of the jet
     * @param E_iso Isotropic equivalent energy
     * @param Gamma0 Initial Lorentz factor
     * @param k_e Power-law index
     * @param k_g Power-law index
     * @param T0 Duration of the ejecta
     * @param spreading Flag indicating if the ejecta spreads laterally during evolution
     * <!-- ************************************************************************************** -->
     */
    PowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k_e, Real k_g, bool spreading = false,
                Real T0 = 1 * unit::sec) noexcept
        : T0(T0), spreading(spreading), theta_c_(theta_c), eps_k_(E_iso / (4 * con::pi)), Gamma0_(Gamma0), k_e_(k_e) {}

    /**
     * <!-- ************************************************************************************** -->
     * @brief Energy per solid angle as a function of phi and theta, with power-law falloff
     * @param phi Azimuthal angle (unused)
     * @param theta Polar angle
     * @return Energy per solid angle with power-law angular dependence
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] inline Real eps_k(Real phi, Real theta) const noexcept {
        return eps_k_ / (1 + fast_pow(theta / theta_c_, k_e_));
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Initial Lorentz factor as a function of phi and theta, with power-law falloff
     * @param phi Azimuthal angle (unused)
     * @param theta Polar angle
     * @return Lorentz factor with power-law angular dependence
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] inline Real Gamma0(Real phi, Real theta) const noexcept {
        return (Gamma0_ - 1) / (1 + fast_pow(theta / theta_c_, k_g_)) + 1;
    }

    /// Duration of the ejecta in seconds
    Real T0{1 * unit::sec};
    /// Flag indicating if the ejecta spreads laterally during evolution
    bool spreading{false};

  private:
    Real const theta_c_{0}; ///< Core angle of the jet
    Real const eps_k_{0};   ///< Energy per solid angle at the core
    Real const Gamma0_{1};  ///< Initial Lorentz factor at the core
    Real const k_e_{2};     ///< Power-law index for energy angular dependence
    Real const k_g_{2};     ///< Power-law index for Lorentz factor angular dependence
};

/**
 * <!-- ************************************************************************************** -->
 * @namespace math
 * @brief Mathematical helper functions for constructing jet profiles.
 * @details Provides a collection of functions for combining functions, constructing various injection profiles
 *          (e.g., tophat, Gaussian, power-law), and computing integrals. These functions are used to create
 *          complex jet profiles by combining spatial and temporal dependencies.
 * <!-- ************************************************************************************** -->
 */
namespace math {
    /**
     * <!-- ************************************************************************************** -->
     * @brief Combines a spatial function and a temporal function into one function of (phi, theta, t)
     * @tparam F1 Type of the spatial function
     * @tparam F2 Type of the temporal function
     * @param f_spatial Function of (phi, theta)
     * @param f_temporal Function of time
     * @return Combined function where a result is the product of both input functions
     * <!-- ************************************************************************************** -->
     */
    template <typename F1, typename F2>
    inline auto combine(F1 f_spatial, F2 f_temporal) noexcept {
        return [=](Real phi, Real theta, Real t) noexcept { return f_spatial(phi, theta) * f_temporal(t); };
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Creates a time-independent function from a spatial function by ignoring the time parameter t
     * @tparam F1 Type of the spatial function
     * @param f_spatial Function of (phi, theta)
     * @return Function of (phi, theta, t) that ignores t and returns the result of f_spatial
     * <!-- ************************************************************************************** -->
     */
    template <typename F1>
    inline auto t_indep(F1 f_spatial) noexcept {
        return [=](Real phi, Real theta, Real t) noexcept { return f_spatial(phi, theta); };
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns a constant (isotropic) function with a fixed height
     * @param height The constant value to return
     * @return Function that returns the same value regardless of angles
     * <!-- ************************************************************************************** -->
     */
    inline auto isotropic(Real height) noexcept {
        return [=](Real phi, Real theta) noexcept { return height; };
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns a tophat function: constant within the core angle, zero outside
     * @param theta_c Core angle
     * @param height Height (value inside the core)
     * @return Function implementing a tophat profile
     * <!-- ************************************************************************************** -->
     */
    inline auto tophat(Real theta_c, Real height) noexcept {
        return [=](Real phi, Real theta) noexcept { return theta < theta_c ? height : 0; };
    }

    inline auto tophat_plus_one(Real theta_c, Real height) noexcept {
        return [=](Real phi, Real theta) noexcept { return theta < theta_c ? height + 1 : 1; };
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns a Gaussian profile function for jet properties
     * @param theta_c Core angle (standard deviation of the Gaussian)
     * @param height Peak height at the center
     * @return Function implementing a Gaussian profile
     * <!-- ************************************************************************************** -->
     */
    inline auto gaussian(Real theta_c, Real height) noexcept {
        const Real spread = -2 * theta_c * theta_c;
        return [=](Real phi, Real theta) noexcept { return height * fast_exp(theta * theta / spread); };
    }

    inline auto gaussian_plus_one(Real theta_c, Real height) noexcept {
        const Real spread = -2 * theta_c * theta_c;
        return [=](Real phi, Real theta) noexcept { return height * fast_exp(theta * theta / spread) + 1; };
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns a power-law profile function for jet properties
     * @param theta_c Core angle
     * @param height Height at the center
     * @param k Power-law index
     * @return Function implementing a power-law profile
     * <!-- ************************************************************************************** -->
     */
    inline auto powerlaw(Real theta_c, Real height, Real k) noexcept {
        return [=](Real phi, Real theta) noexcept { return height / (1 + fast_pow(theta / theta_c, k)); };
    }

    inline auto powerlaw_plus_one(Real theta_c, Real height, Real k) noexcept {
        return [=](Real phi, Real theta) noexcept { return height / (1 + fast_pow(theta / theta_c, k)) + 1; };
    }

    inline auto powerlaw_wing(Real theta_c, Real height, Real k) noexcept {
        return [=](Real phi, Real theta) noexcept {
            if (theta <= theta_c) {
                return 0.;
            } else {
                return height * fast_pow(theta / theta_c, -k);
            }
        };
    }

    inline auto powerlaw_wing_plus_one(Real theta_c, Real height, Real k) noexcept {
        return [=](Real phi, Real theta) noexcept {
            if (theta <= theta_c) {
                return 1.;
            } else {
                return height * fast_pow(theta / theta_c, -k) + 1;
            }
        };
    }

    inline auto step_powerlaw(Real theta_c, Real height_c, Real height_w, Real k) noexcept {
        return [=](Real phi, Real theta) noexcept {
            if (theta <= theta_c) {
                return height_c;
            } else {
                return height_w * fast_pow(theta / theta_c, -k);
            }
        };
    }

    inline auto step_powerlaw_plus_one(Real theta_c, Real height_c, Real height_w, Real k) noexcept {
        return [=](Real phi, Real theta) noexcept {
            if (theta <= theta_c) {
                return height_c + 1;
            } else {
                return height_w * fast_pow(theta / theta_c, -k) + 1;
            }
        };
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns a two-component profile function for jet properties
     * @param theta_c Core angle
     * @param theta_w Wing angle
     * @param height_c Height at core
     * @param height_w Height at wing
     * @return Function implementing a two-component profile
     * <!-- ************************************************************************************** -->
     */
    inline auto two_component(Real theta_c, Real theta_w, Real height_c, Real height_w) noexcept {
        return [=](Real phi, Real theta) noexcept {
            if (theta <= theta_c) {
                return height_c;
            } else if (theta <= theta_w) {
                return height_w;
            } else {
                return 0.;
            }
        };
    }

    inline auto two_component_plus_one(Real theta_c, Real theta_w, Real height_c, Real height_w) noexcept {
        return [=](Real phi, Real theta) noexcept {
            if (theta <= theta_c) {
                return height_c + 1;
            } else if (theta <= theta_w) {
                return height_w + 1;
            } else {
                return 1.;
            }
        };
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Creates a constant injection profile: returns 1 regardless of time
     * @return Function that always returns 1
     * <!-- ************************************************************************************** -->
     */
    inline auto const_injection(Real height) noexcept {
        return [=](Real phi, Real theta, Real t) noexcept { return height; };
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Creates a step injection profile: returns 1 if t > t0, else 0
     * @param t0 Step time
     * @param height_low height of low
     * @param height_high height of high
     * @return Function implementing a step function
     * <!-- ************************************************************************************** -->
     */
    inline auto step_injection(Real t0, Real height_low, Real height_high) noexcept {
        return [=](Real t) noexcept { return t > t0 ? height_high : height_low; };
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Creates a square injection profile: returns 1 if t is between t0 and t1, else 0
     * @param t0 Start time
     * @param t1 End time
     * @param height_high height of high
     * @param height_low height of low
     * @return Function implementing a square wave
     * <!-- ************************************************************************************** -->
     */
    inline auto square_injection(Real t0, Real t1, Real height_high, Real height_low) noexcept {
        return [=](Real phi, Real theta, Real t) noexcept { return t > t0 && t < t1 ? height_high : height_low; };
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Creates a power-law injection profile: decaying with power-law index q
     * @param t0 Reference time
     * @param q Power-law decay index
     * @param height height
     * @return Function implementing power-law decay
     * <!-- ************************************************************************************** -->
     */
    inline auto powerlaw_injection(Real t0, Real q, Real height) noexcept {
        return [=](Real phi, Real theta, Real t) noexcept { return height * fast_pow(1 + t / t0, -q); };
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Creates a magnetar injection profile: returns power-law decay for theta < theta_c
     * @param t0 Reference time
     * @param q Power-law decay index
     * @param L0 Normalization factor
     * @param theta_c Critical angle
     * @return Function implementing a magnetar injection profile
     * <!-- ************************************************************************************** -->
     */
    inline auto magnetar_injection(Real t0, Real q, Real L0, Real theta_c) {
        return [=](Real phi, Real theta, Real t) noexcept {
            if (theta <= theta_c) {
                const Real tt = 1 + t / t0;
                return L0 * fast_pow(tt, -q);
            } else {
                return 0.;
            }
        };
    }
} // namespace math

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates a Lorentz factor distribution using the Liang & Ghirlanda (2010) prescription
 * @tparam F Type of the energy function
 * @param energy_func Function that computes energy as a function of (phi, theta, t)
 * @param e_max Maximum energy
 * @param gamma_max Maximum Lorentz factor
 * @param idx Power-law index relating energy to Lorentz factor
 * @return Function that calculates the Lorentz factor based on energy
 * @details Returns a lambda function that computes a Lorentz factor using the Liang & Ghirlanda (2010)
 *          prescription. This is a commonly used model for the Lorentz factor distribution in GRB jets.
 *          The returned function calculates:
 *              e = energy_func(phi, theta, 0)
 *              u = (e / e_max)^idx * gamma_max
 *              return sqrt(1 + u^2)
 *          This creates a power-law relationship between the energy and Lorentz factor of the jet.
 * <!-- ************************************************************************************** -->
 */
template <typename F>
auto LiangGhirlanda2010(F energy_func, Real e_max, Real gamma_max, Real idx) {
    return [=](Real phi, Real theta) noexcept {
        // Get the energy at the given angle
        const Real e = energy_func(phi, theta);

        // Calculate the velocity parameter u using power-law scaling
        const Real u = fast_pow(e / e_max, idx) * gamma_max;

        // Convert to the Lorentz factor
        return std::sqrt(1 + u * u);
    };
}
