//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <array>

#include "macros.h"

// c++20 concept
template <typename T>
concept HasDmdt = requires(T t) {
    { t.dm_dt(0.0, 0.0, 0.0) };
};

template <typename T>
concept HasDedt = requires(T t) {
    { t.deps_dt(0.0, 0.0, 0.0) };
};

template <typename T>
concept HasSigma = requires(T t) {
    { t.sigma0(0.0, 0.0) };
};

template <typename T>
concept HasU = requires(T t) {
    { t.U2_th };
};

template <typename T>
concept HasMass = requires(T t) {
    { t.mass(0.0, 0.0, 0.0) };
};

#define MAKE_THIS_ODEINT_STATE(classname, data, array_size)                                                            \
    using array_type = std::array<Real, array_size>;                                                                   \
    using value_type = typename array_type::value_type;                                                                \
    using iterator = typename array_type::iterator;                                                                    \
    using const_iterator = typename array_type::const_iterator;                                                        \
    classname() : data{} {};                                                                                           \
    constexpr size_t size() const noexcept {                                                                           \
        return array_size;                                                                                             \
    }                                                                                                                  \
    constexpr iterator begin() noexcept {                                                                              \
        return data.begin();                                                                                           \
    }                                                                                                                  \
    constexpr iterator end() noexcept {                                                                                \
        return data.end();                                                                                             \
    }                                                                                                                  \
    constexpr const_iterator begin() const noexcept {                                                                  \
        return data.begin();                                                                                           \
    }                                                                                                                  \
    constexpr const_iterator end() const noexcept {                                                                    \
        return data.end();                                                                                             \
    }                                                                                                                  \
    constexpr value_type& operator[](size_t i) noexcept {                                                              \
        return data[i];                                                                                                \
    }                                                                                                                  \
    constexpr const value_type& operator[](size_t i) const noexcept {                                                  \
        return data[i];                                                                                                \
    }
