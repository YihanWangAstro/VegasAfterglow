//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "mesh.h"

#include "../util/macros.h"

bool is_linear_scale(Array const& arr, Real tolerance) {
    if (arr.size() < 2)
        return false; // At least two elements are needed.

    const Real diff = arr[1] - arr[0];
    for (size_t i = 2; i < arr.size(); ++i) {
        if (std::fabs((arr[i] - arr[i - 1] - diff) / diff) > tolerance) {
            return false;
        }
    }
    return true;
}

bool is_log_scale(Array const& arr, Real tolerance) {
    if (arr.size() < 2)
        return false; // At least two elements are needed.

    const Real ratio = arr[1] / arr[0];
    for (size_t i = 2; i < arr.size(); ++i) {
        if (std::fabs((arr[i] / arr[i - 1] - ratio) / ratio) > tolerance) {
            return false;
        }
    }
    return true;
}

Array boundary_to_center(Array const& boundary) {
    Array center({boundary.size() - 1}, 0);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
    return center;
}

Array boundary_to_center_log(Array const& boundary) {
    Array center({boundary.size() - 1}, 0);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
    return center;
}
