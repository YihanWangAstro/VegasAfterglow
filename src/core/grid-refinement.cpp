#include "grid-refinement.h"

size_t adaptive_grid_with_breaks(Real lg2_min, Real lg2_max, std::span<const Real> breaks,
                                 std::span<const Real> break_weights, Real pts_per_decade, Array& grid,
                                 size_t max_refined_breaks, Real refine_radius_decades, Real refine_factor) {
    constexpr size_t MAX_PTS = 256;
    const Real lg2_per_decade = log2_10;
    const Real refine_radius = refine_radius_decades * lg2_per_decade;
    const Real coarse_step = lg2_per_decade / std::max(pts_per_decade, Real(1e-6));
    const Real fine_step = coarse_step / refine_factor;
    const Real merge_eps = 0.5 * fine_step;

    // Collect valid breaks within range with weights
    struct Break {
        Real lg2;
        Real weight;
    };
    std::vector<Break> valid_breaks;
    valid_breaks.reserve(breaks.size());
    for (size_t i = 0; i < breaks.size(); ++i) {
        const Real b = breaks[i];
        if (!(b > 0) || !std::isfinite(b)) {
            continue;
        }
        const Real lg2_b = std::log2(b);
        if (lg2_b > lg2_min && lg2_b < lg2_max) {
            Real w = 1;
            if (i < break_weights.size()) {
                const Real wi = break_weights[i];
                if (std::isfinite(wi) && wi >= 0) {
                    w = wi;
                }
            }
            valid_breaks.push_back({lg2_b, w});
        }
    }

    // Sort by weight descending; top N get refinement
    std::ranges::sort(valid_breaks, [](const auto& a, const auto& b) { return a.weight > b.weight; });
    const size_t n_refine = std::min(max_refined_breaks, valid_breaks.size());

    // Build grid points with break tracking for merge protection
    struct Pt {
        Real lg2;
        bool is_break;
    };
    std::vector<Pt> pts;
    pts.reserve(MAX_PTS * 2);

    // Coarse uniform grid
    const Real span = std::max(lg2_max - lg2_min, Real(0));
    const size_t n_coarse = (span > 0 && coarse_step > 0)
                                ? std::max<size_t>(1, static_cast<size_t>(std::ceil(span / coarse_step)))
                                : size_t(1);
    for (size_t i = 0; i <= n_coarse; ++i) {
        pts.push_back({lg2_min + span * static_cast<Real>(i) / static_cast<Real>(n_coarse), false});
    }

    // Break anchors
    for (const auto& br : valid_breaks) {
        pts.push_back({br.lg2, true});
    }

    // Refinement around top-weighted breaks
    for (size_t i = 0; i < n_refine; ++i) {
        const Real lo = std::max(lg2_min, valid_breaks[i].lg2 - refine_radius);
        const Real hi = std::min(lg2_max, valid_breaks[i].lg2 + refine_radius);
        const size_t n = std::max<size_t>(1, static_cast<size_t>((hi - lo) / fine_step));
        for (size_t k = 0; k <= n; ++k) {
            pts.push_back({lo + (hi - lo) * static_cast<Real>(k) / static_cast<Real>(n), false});
        }
    }

    // Sort and merge close points, protecting break positions
    std::ranges::sort(pts, [](const auto& a, const auto& b) { return a.lg2 < b.lg2; });

    std::vector<Pt> merged;
    merged.reserve(pts.size());
    for (const auto& p : pts) {
        if (!merged.empty() && p.lg2 - merged.back().lg2 < merge_eps) {
            if (p.is_break && !merged.back().is_break) {
                merged.back() = p;
            } else {
                merged.back().is_break = merged.back().is_break || p.is_break;
            }
        } else {
            merged.push_back(p);
        }
    }

    // Ensure endpoints
    if (merged.empty()) {
        merged.push_back({lg2_min, true});
        merged.push_back({lg2_max, true});
    } else {
        merged.front().lg2 = lg2_min;
        if (merged.size() == 1) {
            merged.push_back({lg2_max, true});
        } else {
            merged.back().lg2 = lg2_max;
        }
    }

    // Hard cap with uniform subsampling
    if (merged.size() > MAX_PTS) {
        std::vector<Pt> capped;
        capped.reserve(MAX_PTS);
        for (size_t i = 0; i < MAX_PTS; ++i) {
            capped.push_back(merged[i * (merged.size() - 1) / (MAX_PTS - 1)]);
        }
        merged.swap(capped);
    }

    grid = Array::from_shape({merged.size()});
    for (size_t i = 0; i < merged.size(); ++i) {
        grid(i) = std::exp2(merged[i].lg2);
    }
    return merged.size();
}

Array jump_refinement_grid(Array const& jumps, Real theta_min, Real theta_max, Real avg_spacing) {
    std::vector<Real> points;
    points.reserve(jumps.size() * 3);
    const Real tight = avg_spacing / 8;
    for (size_t idx = 0; idx < jumps.size(); ++idx) {
        const Real jump_theta = jumps(idx);
        if (jump_theta >= con::pi / 2 - 0.01) {
            continue;
        }
        // Tight neighbors to minimize quadrature error at discontinuity
        if (jump_theta - tight >= theta_min) {
            points.push_back(jump_theta - tight);
        }
        if (jump_theta + tight <= theta_max) {
            points.push_back(jump_theta + tight);
        }
        // Add the jump position itself
        if (jump_theta >= theta_min && jump_theta <= theta_max) {
            points.push_back(jump_theta);
        }
    }
    std::ranges::sort(points);
    points.erase(std::ranges::unique(points).begin(), points.end());
    return xt::adapt(points, {points.size()});
}

Array shock_crossing_refinement_logspace(Real t_start, Real t_end, Real t_cross, size_t t_num, size_t base_t_num) {
    t_cross = std::clamp(t_cross, t_start, t_end);
    if (t_cross <= t_start || t_cross >= t_end) {
        return xt::logspace(std::log10(t_start), std::log10(t_end), t_num);
    }

    const Real log_total = std::log10(t_end / t_start);
    const Real log_after = std::log10(t_end / t_cross);

    // Post-crossing: same density as base (fwd-equivalent) logspace
    size_t n_post = std::max<size_t>(static_cast<size_t>(base_t_num * log_after / log_total), 2);
    if (n_post >= t_num) {
        n_post = t_num / 2;
    }
    // Pre-crossing: gets all remaining points (base + extra)
    size_t n_pre = t_num + 1 - n_post;

    Array result = xt::zeros<Real>({t_num});
    size_t idx = 0;

    auto seg1 = xt::logspace(std::log10(t_start), std::log10(t_cross), n_pre);
    for (size_t k = 0; k < n_pre; ++k) {
        result(idx++) = seg1(k);
    }

    auto seg2 = xt::logspace(std::log10(t_cross), std::log10(t_end), n_post);
    for (size_t k = 1; k < n_post && idx < t_num; ++k) {
        result(idx++) = seg2(k);
    }

    return result;
}
