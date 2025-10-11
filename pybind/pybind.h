//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <pybind11/stl.h>

#include "xtensor-python/pytensor.hpp"

namespace py = pybind11;
using PyArray = xt::pytensor<double, 1>;
using PyGrid = xt::pytensor<double, 2>;

template <typename Array>
bool is_ascending(Array const& arr) {
    if (arr.size() <= 1)
        return true;
    for (size_t i = 1; i < arr.size(); ++i) {
        if (arr(i) < arr(i - 1)) {
            return false;
        }
    }
    return true;
}
