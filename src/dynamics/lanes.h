//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <algorithm>
#include <array>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <cmath>
#include <cstddef>

#include "../util/macros.h"

/// Fixed-width value pack for the batched shock-dynamics solver: one shock per lane,
/// integrated in lockstep. Native vector storage on clang/gcc (guaranteed SIMD),
/// elementwise loops on MSVC; both compute identical IEEE per-lane values.
#if defined(_MSC_VER) && !defined(__clang__)
    #define LANES_LOOP_OPS 1
#endif

/// Per-lane boolean result of comparing Lanes packs. Combine with &&/||/! (eager,
/// no short-circuit); consume with where() / all_of() / any_of(). For Value = Real
/// the same generic code produces a plain bool instead.
template <std::size_t K>
struct LaneMask {
#ifdef LANES_LOOP_OPS
    std::array<bool, K> m;
    bool operator[](std::size_t w) const noexcept { return m[w]; }

    #define LANES_MASK_OP(name, op)                                                                                    \
        friend LaneMask name(LaneMask const& a, LaneMask const& b) noexcept {                                          \
            LaneMask r;                                                                                                \
            for (std::size_t w = 0; w < K; ++w) {                                                                      \
                r.m[w] = a.m[w] op b.m[w];                                                                             \
            }                                                                                                          \
            return r;                                                                                                  \
        }
    LANES_MASK_OP(operator&&, &&)
    LANES_MASK_OP(operator||, ||)
    #undef LANES_MASK_OP
#else
    typedef long long storage __attribute__((vector_size(K * sizeof(long long))));
    storage m;
    bool operator[](std::size_t w) const noexcept { return m[w] != 0; }

    friend LaneMask operator&&(LaneMask const& a, LaneMask const& b) noexcept { return LaneMask{a.m & b.m}; }
    friend LaneMask operator||(LaneMask const& a, LaneMask const& b) noexcept { return LaneMask{a.m | b.m}; }
#endif
};

template <std::size_t K>
inline bool all_of(LaneMask<K> const& m) noexcept {
    bool all = true;
    for (std::size_t w = 0; w < K; ++w) {
        all = all && m[w];
    }
    return all;
}

template <std::size_t K>
inline bool any_of(LaneMask<K> const& m) noexcept {
    bool any = false;
    for (std::size_t w = 0; w < K; ++w) {
        any = any || m[w];
    }
    return any;
}

/// Scalar counterparts so Value-generic code compiles unchanged for Value = Real.
inline constexpr bool all_of(bool m) noexcept {
    return m;
}
inline constexpr bool any_of(bool m) noexcept {
    return m;
}
inline constexpr Real where(bool c, Real a, Real b) noexcept {
    return c ? a : b;
}

template <std::size_t K>
struct Lanes {
#ifdef LANES_LOOP_OPS
    std::array<Real, K> v;
#else
    typedef Real storage __attribute__((vector_size(K * sizeof(Real))));
    storage v;
#endif

    Lanes() = default;
    Lanes(Real s) noexcept { ///< broadcast (implicit: mirrors Real -> Value in scalar code)
        for (std::size_t w = 0; w < K; ++w) {
            v[w] = s;
        }
    }

    Real operator[](std::size_t w) const noexcept { return v[w]; }

    // Mixed Real operands broadcast through the implicit constructor.
#ifdef LANES_LOOP_OPS
    #define LANES_BINARY_OP(op)                                                                                        \
        friend Lanes operator op(Lanes const& a, Lanes const& b) noexcept {                                            \
            Lanes r;                                                                                                   \
            for (std::size_t w = 0; w < K; ++w) {                                                                      \
                r.v[w] = a.v[w] op b.v[w];                                                                             \
            }                                                                                                          \
            return r;                                                                                                  \
        }
    #define LANES_CMP_OP(op)                                                                                           \
        friend LaneMask<K> operator op(Lanes const& a, Lanes const& b) noexcept {                                      \
            LaneMask<K> r;                                                                                             \
            for (std::size_t w = 0; w < K; ++w) {                                                                      \
                r.m[w] = a.v[w] op b.v[w];                                                                             \
            }                                                                                                          \
            return r;                                                                                                  \
        }
#else
    #define LANES_BINARY_OP(op)                                                                                        \
        friend Lanes operator op(Lanes const& a, Lanes const& b) noexcept {                                            \
            Lanes r;                                                                                                   \
            r.v = a.v op b.v;                                                                                          \
            return r;                                                                                                  \
        }
    #define LANES_CMP_OP(op)                                                                                           \
        friend LaneMask<K> operator op(Lanes const& a, Lanes const& b) noexcept {                                      \
            return LaneMask<K>{a.v op b.v};                                                                            \
        }
#endif

    LANES_BINARY_OP(+)
    LANES_BINARY_OP(-)
    LANES_BINARY_OP(*)
    LANES_BINARY_OP(/)
    LANES_CMP_OP(<)
    LANES_CMP_OP(<=)
    LANES_CMP_OP(>)
    LANES_CMP_OP(>=)
    LANES_CMP_OP(==)
    LANES_CMP_OP(!=)
#undef LANES_BINARY_OP
#undef LANES_CMP_OP

    friend Lanes operator-(Lanes const& a) noexcept { return Lanes(0) - a; }
    Lanes& operator+=(Lanes const& b) noexcept { return *this = *this + b; }
    Lanes& operator-=(Lanes const& b) noexcept { return *this = *this - b; }
    Lanes& operator*=(Lanes const& b) noexcept { return *this = *this * b; }
    Lanes& operator/=(Lanes const& b) noexcept { return *this = *this / b; }

    /// Per-lane conditional: lane w of the result is a's where c[w] holds, else b's.
    friend Lanes where(LaneMask<K> const& c, Lanes const& a, Lanes const& b) noexcept {
        Lanes r;
#ifdef LANES_LOOP_OPS
        for (std::size_t w = 0; w < K; ++w) {
            r.v[w] = c.m[w] ? a.v[w] : b.v[w];
        }
#else
        r.v = c.m ? a.v : b.v;
#endif
        return r;
    }

    /// Elementwise min/max/clamp (per-lane std semantics; mixed Real operands broadcast).
    friend Lanes min(Lanes const& a, Lanes const& b) noexcept {
        Lanes r;
        for (std::size_t w = 0; w < K; ++w) {
            r.v[w] = std::min(a.v[w], b.v[w]);
        }
        return r;
    }
    friend Lanes max(Lanes const& a, Lanes const& b) noexcept {
        Lanes r;
        for (std::size_t w = 0; w < K; ++w) {
            r.v[w] = std::max(a.v[w], b.v[w]);
        }
        return r;
    }
    friend Lanes clamp(Lanes const& x, Lanes const& lo, Lanes const& hi) noexcept {
        Lanes r;
        for (std::size_t w = 0; w < K; ++w) {
            r.v[w] = std::clamp(x.v[w], lo.v[w], hi.v[w]);
        }
        return r;
    }
};

namespace lane {

    /// Lanes per batch: a full 512-bit register on AVX-512, otherwise 4 (one AVX register,
    /// or two NEON/SSE registers). Define AFTERGLOW_LANE_WIDTH to override.
#if defined(AFTERGLOW_LANE_WIDTH)
    inline constexpr std::size_t width = AFTERGLOW_LANE_WIDTH;
#elif defined(__AVX512F__)
    inline constexpr std::size_t width = 8;
#else
    inline constexpr std::size_t width = 4;
#endif

    template <typename Value>
    inline constexpr std::size_t width_of = 1;
    template <std::size_t K>
    inline constexpr std::size_t width_of<Lanes<K>> = K;

    using Batch = Lanes<width>;

    /// Scalar lane `w` of a value (identity for Real).
    inline Real get(Real v, std::size_t) noexcept {
        return v;
    }
    template <std::size_t K>
    inline Real get(Lanes<K> const& v, std::size_t w) noexcept {
        return v[w];
    }

    /// Writes scalar lane `w` of a value (the whole value for Real).
    inline void set(Real& v, std::size_t, Real x) noexcept {
        v = x;
    }
    template <std::size_t K>
    inline void set(Lanes<K>& v, std::size_t w, Real x) noexcept {
        v.v[w] = x;
    }

    /// Applies a scalar function elementwise across one or more packs.
    template <std::size_t K, typename F, typename... L>
    inline Lanes<K> map(F&& f, Lanes<K> const& first, L const&... rest) noexcept {
        Lanes<K> r;
        for (std::size_t w = 0; w < K; ++w) {
            r.v[w] = f(first[w], rest[w]...);
        }
        return r;
    }

    /// Builds a value by evaluating a scalar function once per lane: the bridge for
    /// per-lane constants and for scalar-only callables (profiles, branchy helpers).
    template <typename Value, typename F>
    inline Value from_lanes(F&& f) noexcept {
        if constexpr (width_of<Value> == 1) {
            return f(std::size_t{0});
        } else {
            Value r;
            for (std::size_t w = 0; w < width_of<Value>; ++w) {
                r.v[w] = f(w);
            }
            return r;
        }
    }

    /// Calls a scalar function once per lane with that lane of each argument:
    /// r[w] = f(get(args, w)...). For scalar-only helpers (e.g. the sigma > 0 cubic).
    template <typename Value, typename F, typename... Args>
    inline Value apply(F&& f, Args const&... args) noexcept {
        return from_lanes<Value>([&](std::size_t w) { return f(get(args, w)...); });
    }

    /// all/any over the drivers' per-lane bookkeeping flags.
    template <std::size_t K>
    inline bool any_of(std::array<bool, K> const& m) noexcept {
        bool any = false;
        for (std::size_t w = 0; w < K; ++w) {
            any = any || m[w];
        }
        return any;
    }

    /// Largest |lane| (the scalar reduction for the lockstep ODE error norm).
    inline Real max_abs(Real v) noexcept {
        return std::fabs(v);
    }
    template <std::size_t K>
    inline Real max_abs(Lanes<K> const& v) noexcept {
        Real m = 0;
        for (std::size_t w = 0; w < K; ++w) {
            m = std::max(m, std::fabs(v[w]));
        }
        return m;
    }

    /// Copies lane `from` onto lane `to` (parks dead/padded lanes on live values).
    inline void copy_lane(Real&, std::size_t, std::size_t) noexcept {}
    template <std::size_t K>
    inline void copy_lane(Lanes<K>& v, std::size_t from, std::size_t to) noexcept {
        v.v[to] = v.v[from];
    }

} // namespace lane

template <std::size_t K>
inline Lanes<K> sqrt(Lanes<K> const& a) noexcept {
    return lane::map([](Real x) { return std::sqrt(x); }, a);
}
template <std::size_t K>
inline Lanes<K> sin(Lanes<K> const& a) noexcept {
    return lane::map([](Real x) { return std::sin(x); }, a);
}
template <std::size_t K>
inline Lanes<K> cos(Lanes<K> const& a) noexcept {
    return lane::map([](Real x) { return std::cos(x); }, a);
}
template <std::size_t K>
inline Lanes<K> abs(Lanes<K> const& a) noexcept {
    return lane::map([](Real x) { return std::fabs(x); }, a);
}
template <std::size_t K>
inline LaneMask<K> isfinite(Lanes<K> const& a) noexcept {
    LaneMask<K> r{};
    for (std::size_t w = 0; w < K; ++w) {
#ifdef LANES_LOOP_OPS
        r.m[w] = std::isfinite(a[w]);
#else
        r.m[w] = std::isfinite(a[w]) ? -1 : 0;
#endif
    }
    return r;
}

// Value-generic bodies call these unqualified with a function-local `using std::sqrt;`
// (etc.), so the Real instantiation resolves to libm while ADL finds the Lanes
// overloads above.

/// Scalar view of one lane of a batched odeint state.
template <template <typename, typename> class State, typename E, typename Value>
inline State<E, Real> extract_lane(State<E, Value> const& state, std::size_t w) noexcept {
    State<E, Real> out;
    for (std::size_t e = 0; e < out.size(); ++e) {
        out[e] = lane::get(state[e], w);
    }
    return out;
}

/// Writes a scalar state into lane `w` of a batched odeint state (inverse of extract_lane).
template <template <typename, typename> class State, typename E, typename Value>
inline void insert_lane(State<E, Value>& state, std::size_t w, State<E, Real> const& scalar) noexcept {
    for (std::size_t e = 0; e < scalar.size(); ++e) {
        lane::set(state[e], w, scalar[e]);
    }
}

/// Copies lane `src`'s values onto every non-live lane: dead rows and padding must carry
/// finite values so they don't poison the batch's shared step-size control.
template <typename State, std::size_t K>
inline void park_dead_lanes(State& state, std::array<bool, K> const& live, std::size_t src) noexcept {
    for (std::size_t w = 0; w < K; ++w) {
        if (!live[w]) {
            for (std::size_t e = 0; e < state.size(); ++e) {
                lane::copy_lane(state[e], src, w);
            }
        }
    }
}

/// odeint range algebra whose error norm reduces across lanes, so a batch steps at the
/// most demanding lane's step size (lockstep).
struct LaneAlgebra : boost::numeric::odeint::range_algebra {
    template <typename S>
    Real norm_inf(S const& s) const noexcept {
        Real m = 0;
        for (auto const& v : s) {
            m = std::max(m, lane::max_abs(v));
        }
        return m;
    }
};
