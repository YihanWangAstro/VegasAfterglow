//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <concepts>
#include <cstddef>
#include <type_traits>

#include "../util/macros.h"

struct InverseComptonY;

/**
 * <!-- ************************************************************************************** -->
 * @concept SynElectronModel
 * @brief Contract a synchrotron electron-distribution type must satisfy to drive ICPhoton and
 *        the IC cooling pipeline.
 * @details This is the interface consumed by ICPhoton, IC_cooling/Thomson_cooling/KN_cooling/
 *          CMB_cooling, and the photon factory (generate_syn_photons). Mutability matters:
 *          IC_cooling writes gamma_a/gamma_c/gamma_M/Y_c/regime and mutates Ys in place, so
 *          these must be plain public data members (genuine lvalues), not accessors.
 *          n_breaks must be a constant expression equal to 4 — ICPhoton sizes its break arrays
 *          with it and fills them with {a, m, c, M}.
 * <!-- ************************************************************************************** -->
 */
template <typename T>
concept SynElectronModel =
    std::default_initializable<T> && std::copyable<T> && requires(T elec, T const c_elec, Real gamma) {
        typename std::integral_constant<size_t, T::n_breaks>;
        { elec.gamma_a } -> std::same_as<Real&>;
        { elec.gamma_m } -> std::same_as<Real&>;
        { elec.gamma_c } -> std::same_as<Real&>;
        { elec.gamma_M } -> std::same_as<Real&>;
        { elec.p } -> std::same_as<Real&>;
        { elec.N_e } -> std::same_as<Real&>;
        { elec.column_den } -> std::same_as<Real&>;
        { elec.Y_c } -> std::same_as<Real&>;
        { elec.regime } -> std::same_as<size_t&>;
        { elec.Ys } -> std::same_as<InverseComptonY&>;
        { c_elec.compute_column_den(gamma) } -> std::same_as<Real>;
    };

/**
 * <!-- ************************************************************************************** -->
 * @concept SynPhotonModel
 * @brief Contract a synchrotron photon-spectrum type must satisfy to be usable as SynPhotons.
 * @details This is the interface consumed by ICPhoton (seed spectrum), Observer::specific_flux
 *          (compute_log2_I_nu), the pybind details/spectrum paths, and the factory
 *          generate_syn_photons — which default-constructs elements, assigns the public fields,
 *          then calls build(). build() must be re-runnable: IC_cooling repopulates already-built
 *          grids in place. compute_I_nu/compute_log2_I_nu must be const (spectrum accessors hold
 *          const grids) and are only valid after build(). noexcept is not part of the contract.
 *          Swapping the active model is a one-line change of the SynPhotons alias in
 *          synchrotron.h; this concept is the compiler-checked checklist for any new model.
 * <!-- ************************************************************************************** -->
 */
template <typename T>
concept SynPhotonModel = std::default_initializable<T> && std::copyable<T> && requires(T ph, T const c_ph, Real nu) {
    typename std::integral_constant<size_t, T::n_breaks>;
    { ph.nu_a } -> std::same_as<Real&>;
    { ph.nu_m } -> std::same_as<Real&>;
    { ph.nu_c } -> std::same_as<Real&>;
    { ph.nu_M } -> std::same_as<Real&>;
    { ph.I_nu_max } -> std::same_as<Real&>;
    { ph.p } -> std::same_as<Real&>;
    { ph.Y_c } -> std::same_as<Real&>;
    { ph.regime } -> std::same_as<size_t&>;
    { ph.Ys } -> std::same_as<InverseComptonY&>;
    { ph.build() };
    { c_ph.compute_I_nu(nu) } -> std::same_as<Real>;
    { c_ph.compute_log2_I_nu(nu) } -> std::same_as<Real>;
};
