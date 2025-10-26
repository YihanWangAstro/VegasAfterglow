//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <string>

#include "../core/mesh.h"
#include "../dynamics/shock.h"
#include "../radiation/prompt.h"
#include "../radiation/synchrotron.h"

#ifndef NO_XTENSOR_IO
    #include "xtensor-io/xnpz.hpp"
#endif

/**
 * <!-- ************************************************************************************** -->
 * @defgroup IO_Functions Output and Printing Functions
 * @brief Functions for printing and outputting simulation data to files.
 * @details These functions handle output of various data types to files, including 1D arrays, grids, and
 *          different model components like SynPhotonGrid, SynElectronGrid, Shock, etc.
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Write Array to a CSV file with an optional unit for value scaling
 * @param filename The output filename
 * @param array The array to write
 * @param unit Optional scaling factor (default: 1.0)
 * <!-- ************************************************************************************** -->
 */
void write_csv(std::string const& filename, Array const& array, Real unit = 1.0);

/**
 * <!-- ************************************************************************************** -->
 * @brief Write MeshGrid to a CSV file with an optional unit for value scaling
 * @param filename The output filename
 * @param grid The grid to write
 * @param unit Optional scaling factor (default: 1.0)
 * <!-- ************************************************************************************** -->
 */
void write_csv(std::string const& filename, MeshGrid const& grid, Real unit = 1.0);

/**
 * <!-- ************************************************************************************** -->
 * @brief Write MeshGrid3d to a CSV file with an optional unit for value scaling
 * @param filename The output filename
 * @param grid3d The 3D grid to write
 * @param unit Optional scaling factor (default: 1.0)
 * <!-- ************************************************************************************** -->
 */
void write_csv(std::string const& filename, MeshGrid3d const& grid3d, Real unit = 1.0);

/**
 * <!-- ************************************************************************************** -->
 * @brief Write an array to an NPY file with an optional unit for value scaling
 * @tparam T Type of the array
 * @param filename The output filename
 * @param array The array to write
 * @param unit Optional scaling factor (default: 1.0)
 * <!-- ************************************************************************************** -->
 */
template <typename T>
void write_npy(std::string const& filename, const T& array, Real unit = 1.0);

#ifndef NO_XTENSOR_IO
/**
 * <!-- ************************************************************************************** -->
 * @brief Write SynPhotonGrid to an NPZ file
 * @param filename The output filename
 * @param syn_ph The synchrotron photon grid to write
 * <!-- ************************************************************************************** -->
 */
void write_npz(std::string const& filename, SynPhotonGrid const& syn_ph);

/**
 * <!-- ************************************************************************************** -->
 * @brief Write SynElectronGrid to an NPZ file
 * @param filename The output filename
 * @param syn_e The synchrotron electron grid to write
 * <!-- ************************************************************************************** -->
 */
void write_npz(std::string const& filename, SynElectronGrid const& syn_e);

/**
 * <!-- ************************************************************************************** -->
 * @brief Write Shock to an NPZ file
 * @param filename The output filename
 * @param shock The shock object to write
 * <!-- ************************************************************************************** -->
 */
void write_npz(std::string const& filename, Shock const& shock);

/**
 * <!-- ************************************************************************************** -->
 * @brief Write Coord to an NPZ file
 * @param filename The output filename
 * @param coord The coordinate object to write
 * <!-- ************************************************************************************** -->
 */
void write_npz(std::string const& filename, Coord const& coord);

/**
 * <!-- ************************************************************************************** -->
 * @brief Helper recursive function for writing multiple arrays to a single NPZ file
 * @tparam T Type of the current array
 * @tparam Rest Types of the rest of the arguments
 * @param filename The output filename
 * @param first Boolean indicating if this is the first array being written
 * @param name Name for the current array
 * @param array The current array to write
 * @param rest Rest of the arguments (name-array pairs)
 * <!-- ************************************************************************************** -->
 */
template <typename T, typename... Rest>
void write_npz_recursive(std::string const& filename, bool first, std::string const& name, const T& array,
                         const Rest&... rest);

/**
 * <!-- ************************************************************************************** -->
 * @brief Write an array to an NPZ file with an optional unit for value scaling
 * @tparam T Type of the array
 * @param filename The output filename
 * @param array The array to write
 * @param unit Optional scaling factor (default: 1.0)
 * <!-- ************************************************************************************** -->
 */
template <typename T>
void write_npz(std::string const& filename, const T& array, Real unit = 1.0) {
    xt::dump_npz(filename + ".npz", "array", xt::eval(array / unit), false, false);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Write multiple named arrays to a single NPZ file
 * @tparam Args Types of the arguments (alternating string names and arrays)
 * @param filename The output filename
 * @param args Variable arguments consisting of name-array pairs
 * @note Usage: write_npz("filename", "name1", array1, "name2", array2, ...)
 * <!-- ************************************************************************************** -->
 */
template <typename... Args>
void write_npz(std::string const& filename, Args const&... args);

//========================================================================================================
//                                  template function implementation
//========================================================================================================

template <typename T, typename... Rest>
void write_npz_recursive(std::string const& filename, bool first, std::string const& name, const T& array,
                         const Rest&... rest) {
    auto arr = xt::eval(array);                                // ensure evaluated
    xt::dump_npz(filename + ".npz", name, arr, false, !first); // append after first write

    if constexpr (sizeof...(rest) > 0) {
        write_npz_recursive(filename, false, rest...); // continue with rest
    }
}

template <typename... Args>
void write_npz(std::string const& filename, Args const&... args) {
    static_assert(sizeof...(args) % 2 == 0, "Arguments must be pairs: name1, array1, name2, array2, ...");

    write_npz_recursive(filename, true, args...);
}
#endif

template <typename T>
void write_npy(std::string const& filename, const T& array, Real unit) {
    xt::dump_npy(filename + ".npy", xt::eval(array / unit));
}
