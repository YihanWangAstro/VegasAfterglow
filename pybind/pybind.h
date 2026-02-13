//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "xtensor-python/pytensor.hpp"

namespace py = pybind11;
using PyArray = xt::pytensor<double, 1>;
using PyGrid = xt::pytensor<double, 2>;
using PyGrid3d = xt::pytensor<double, 3>;

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

inline std::vector<size_t> logscale_screen(PyArray const& data, size_t num_order) {
    const size_t total_size = data.size();

    // Handle edge cases: empty or single-element arrays
    if (total_size == 0) {
        return {};
    }
    if (total_size == 1) {
        return {0};
    }

    if (num_order == 0) {
        std::vector<size_t> indices(total_size);
        std::iota(indices.begin(), indices.end(), 0);
        return indices;
    }

    const double log_start = std::log10(static_cast<double>(data(0)));
    const double log_end = std::log10(static_cast<double>(data(total_size - 1)));
    const double log_range = log_end - log_start;
    const size_t total_points = static_cast<size_t>(std::ceil(log_range * static_cast<double>(num_order))) + 1;

    std::vector<size_t> indices;
    indices.reserve(total_points);

    // Always include the first point
    indices.push_back(0);

    // Handle case where total_points <= 1 (avoid division by zero)
    if (total_points <= 1) {
        if (total_size > 1) {
            indices.push_back(total_size - 1);
        }
        return indices;
    }

    const double step = log_range / static_cast<double>(total_points - 1);

    for (size_t i = 1; i < total_points - 1; ++i) {
        const double log_target = log_start + static_cast<double>(i) * step;
        const double target_value = std::pow(10.0, log_target);

        size_t best_idx = 1;
        double min_diff = std::abs(static_cast<double>(data(1)) - target_value);

        for (size_t j = 2; j < total_size - 1; ++j) {
            const double diff = std::abs(static_cast<double>(data(j)) - target_value);
            if (diff < min_diff) {
                min_diff = diff;
                best_idx = j;
            }
        }

        if (std::ranges::find(indices, best_idx) == indices.end()) {
            indices.push_back(best_idx);
        }
    }

    // Always include the last point
    if (total_size > 1) {
        indices.push_back(total_size - 1);
    }

    std::ranges::sort(indices);
    indices.erase(std::ranges::unique(indices).begin(), indices.end());

    return indices;
}
