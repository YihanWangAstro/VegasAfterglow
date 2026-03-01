#include <boost/test/unit_test.hpp>
#include <cmath>
#include <span>

#include "core/mesh.h"
#include "util/macros.h"

BOOST_AUTO_TEST_SUITE(Mesh)

// ---------------------------------------------------------------------------
// boundary_to_center (non-template, returns Array)
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(boundary_to_center_linear) {
    Array boundary = {1.0, 3.0, 5.0};
    Array center = boundary_to_center(boundary);
    BOOST_REQUIRE_EQUAL(center.size(), 2u);
    BOOST_CHECK_CLOSE(center(0), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(center(1), 4.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(boundary_to_center_single) {
    Array boundary = {1.0, 3.0};
    Array center = boundary_to_center(boundary);
    BOOST_REQUIRE_EQUAL(center.size(), 1u);
    BOOST_CHECK_CLOSE(center(0), 2.0, 1e-10);
}

// ---------------------------------------------------------------------------
// boundary_to_center_log (non-template, returns Array)
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(boundary_to_center_log_values) {
    Array boundary = {1.0, 4.0, 16.0};
    Array center = boundary_to_center_log(boundary);
    BOOST_REQUIRE_EQUAL(center.size(), 2u);
    BOOST_CHECK_CLOSE(center(0), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(center(1), 8.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(boundary_to_center_log_single) {
    Array boundary = {1.0, 9.0};
    Array center = boundary_to_center_log(boundary);
    BOOST_REQUIRE_EQUAL(center.size(), 1u);
    BOOST_CHECK_CLOSE(center(0), 3.0, 1e-10);
}

// ---------------------------------------------------------------------------
// boundary_to_center (template version) edge cases
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(boundary_to_center_empty) {
    // Empty boundary should not crash
    Array boundary = Array::from_shape({0});
    Array center = boundary_to_center(boundary);
    BOOST_CHECK_EQUAL(center.size(), 0u);
}

BOOST_AUTO_TEST_CASE(boundary_to_center_size_mismatch) {
    // Template version: if center.size()+1 != boundary.size(), early return
    Array boundary = {1.0, 3.0, 5.0};
    Array center = Array::from_shape({5}); // wrong size: need 2, got 5
    xt::xarray<Real> old_center = center;
    boundary_to_center(boundary, center);
    // center should be untouched (early return)
    for (size_t i = 0; i < center.size(); ++i) {
        BOOST_CHECK_EQUAL(center(i), old_center(i));
    }
}

// ---------------------------------------------------------------------------
// logspace_center
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(logspace_center_count) {
    Array center;
    logspace_center(0.0, 10.0, 20, center);
    BOOST_CHECK_EQUAL(center.size(), 20u);
}

BOOST_AUTO_TEST_CASE(logspace_center_bounds) {
    Array center;
    const Real lg2_min = 2.0;
    const Real lg2_max = 8.0;
    logspace_center(lg2_min, lg2_max, 50, center);
    for (size_t i = 0; i < center.size(); ++i) {
        BOOST_CHECK_GE(center(i), std::exp2(lg2_min));
        BOOST_CHECK_LE(center(i), std::exp2(lg2_max));
    }
}

BOOST_AUTO_TEST_CASE(logspace_center_monotonic) {
    Array center;
    logspace_center(0.0, 20.0, 100, center);
    for (size_t i = 1; i < center.size(); ++i) {
        BOOST_CHECK_GT(center(i), center(i - 1));
    }
}

BOOST_AUTO_TEST_CASE(logspace_center_ratio) {
    Array center;
    logspace_center(0.0, 10.0, 40, center);
    // In logspace, consecutive centers should have a uniform ratio
    const Real expected_ratio = center(1) / center(0);
    for (size_t i = 2; i < center.size(); ++i) {
        const Real ratio = center(i) / center(i - 1);
        BOOST_CHECK_CLOSE(ratio, expected_ratio, 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(logspace_center_zero_size) {
    Array center;
    logspace_center(0.0, 10.0, 0, center);
    BOOST_CHECK_EQUAL(center.size(), 0u);
}

BOOST_AUTO_TEST_CASE(logspace_center_size_one) {
    Array center;
    logspace_center(0.0, 10.0, 1, center);
    BOOST_REQUIRE_EQUAL(center.size(), 1u);
    // Single center should be the geometric mean of endpoints
    const Real lo = std::exp2(0.0);
    const Real hi = std::exp2(10.0);
    BOOST_CHECK_CLOSE(center(0), std::sqrt(lo * hi), 1e-8);
}

// ---------------------------------------------------------------------------
// log2space
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(log2space_count) {
    Array arr;
    log2space(0.0, 10.0 * log2_10, 4, arr);
    // 10 decades * 4 points per decade = 40, plus 1 -> 41
    const Real decades = 10.0;
    const size_t expected = static_cast<size_t>(std::ceil(decades * 4)) + 1;
    BOOST_CHECK_EQUAL(arr.size(), expected);
}

BOOST_AUTO_TEST_CASE(log2space_bounds) {
    Array arr;
    const Real lg2_min = 5.0;
    const Real lg2_max = 15.0;
    log2space(lg2_min, lg2_max, 6, arr);
    BOOST_CHECK_CLOSE(arr(0), std::exp2(lg2_min), 1e-8);
    BOOST_CHECK_CLOSE(arr(arr.size() - 1), std::exp2(lg2_max), 1e-8);
}

// ---------------------------------------------------------------------------
// logspace_boundary_center
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(logspace_boundary_center_consistency) {
    Array center;
    Array bin_width;
    logspace_boundary_center(2.0, 10.0, 30, center, bin_width);
    BOOST_REQUIRE_EQUAL(center.size(), 30u);
    BOOST_REQUIRE_EQUAL(bin_width.size(), 30u);

    // Reconstruct boundaries from center and bin_width:
    // center(i) = left * sqrt(r), bin_width(i) = left * (r - 1)
    // so left = center(i) / sqrt(r), right = left * r
    // We check that center is approximately the geometric mean of adjacent boundaries.
    // A simpler check: build boundaries from the logspace and compare.
    const Real dlg2 = (10.0 - 2.0) / 30.0;
    const Real r = std::exp2(dlg2);
    Real left = std::exp2(2.0);
    for (size_t i = 0; i < 30; ++i) {
        const Real right = left * r;
        const Real geo_center = std::sqrt(left * right);
        BOOST_CHECK_CLOSE(center(i), geo_center, 1e-8);
        left = right;
    }
}

BOOST_AUTO_TEST_CASE(logspace_boundary_center_width) {
    Array center;
    Array bin_width;
    logspace_boundary_center(0.0, 20.0, 50, center, bin_width);
    for (size_t i = 0; i < bin_width.size(); ++i) {
        BOOST_CHECK_GT(bin_width(i), 0.0);
    }
}

BOOST_AUTO_TEST_CASE(logspace_boundary_center_zero_size) {
    Array center;
    Array bin_width;
    logspace_boundary_center(0.0, 10.0, 0, center, bin_width);
    BOOST_CHECK_EQUAL(center.size(), 0u);
    BOOST_CHECK_EQUAL(bin_width.size(), 0u);
}

// ---------------------------------------------------------------------------
// adaptive_grid_with_breaks
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(adaptive_grid_endpoints) {
    Array grid;
    Real breaks[] = {1e5, 1e8};
    Real weights[] = {1.0, 1.0};
    size_t n = adaptive_grid_with_breaks(0.0, 40.0, std::span(breaks), std::span(weights), 4, grid);
    BOOST_REQUIRE_GT(n, 1u);
    BOOST_CHECK_CLOSE(grid(0), std::exp2(0.0), 1e-8);
    BOOST_CHECK_CLOSE(grid(n - 1), std::exp2(40.0), 1e-4);
}

BOOST_AUTO_TEST_CASE(adaptive_grid_contains_breaks) {
    Array grid;
    Real breaks[] = {1e5, 1e8};
    Real weights[] = {1.0, 1.0};
    size_t n = adaptive_grid_with_breaks(0.0, 40.0, std::span(breaks), std::span(weights), 4, grid);

    // Each break should appear in the grid approximately (within a factor of 2)
    for (Real brk : breaks) {
        bool found = false;
        for (size_t i = 0; i < n; ++i) {
            if (std::fabs(std::log2(grid(i)) - std::log2(brk)) < 1.0) {
                found = true;
                break;
            }
        }
        BOOST_CHECK_MESSAGE(found, "Break " << brk << " not found (approximately) in grid");
    }
}

BOOST_AUTO_TEST_CASE(adaptive_grid_max_256) {
    Array grid;
    Real breaks[] = {1e5};
    Real weights[] = {1.0};
    size_t n = adaptive_grid_with_breaks(0.0, 40.0, std::span(breaks), std::span(weights), 100, grid);
    BOOST_CHECK_LE(n, 256u);
    BOOST_CHECK_LE(grid.size(), 256u);
}

BOOST_AUTO_TEST_CASE(adaptive_grid_refined_denser) {
    Array grid;
    Real breaks[] = {1e10};
    Real weights[] = {1.0};
    size_t n = adaptive_grid_with_breaks(0.0, 40.0, std::span(breaks), std::span(weights), 4, grid);

    // Compute average spacing near break vs far from break
    const Real lg2_break = std::log2(1e10);
    Real near_sum = 0;
    size_t near_count = 0;
    Real far_sum = 0;
    size_t far_count = 0;

    for (size_t i = 1; i < n; ++i) {
        const Real lg2_mid = 0.5 * (std::log2(grid(i)) + std::log2(grid(i - 1)));
        const Real step = std::log2(grid(i)) - std::log2(grid(i - 1));
        if (std::fabs(lg2_mid - lg2_break) < 2.0) {
            near_sum += step;
            near_count++;
        } else if (std::fabs(lg2_mid - lg2_break) > 10.0) {
            far_sum += step;
            far_count++;
        }
    }

    if (near_count > 0 && far_count > 0) {
        const Real near_avg = near_sum / near_count;
        const Real far_avg = far_sum / far_count;
        BOOST_CHECK_LT(near_avg, far_avg);
    }
}

BOOST_AUTO_TEST_CASE(adaptive_grid_no_breaks) {
    Array grid;
    std::span<const Real> no_breaks;
    std::span<const Real> no_weights;
    size_t n = adaptive_grid_with_breaks(0.0, 20.0, no_breaks, no_weights, 4, grid);
    BOOST_CHECK_GT(n, 1u);
    BOOST_CHECK_CLOSE(grid(0), std::exp2(0.0), 1e-8);
    BOOST_CHECK_CLOSE(grid(n - 1), std::exp2(20.0), 1e-4);
}

BOOST_AUTO_TEST_CASE(adaptive_grid_min_equals_max) {
    Array grid;
    std::span<const Real> no_breaks;
    std::span<const Real> no_weights;
    size_t n = adaptive_grid_with_breaks(10.0, 10.0, no_breaks, no_weights, 4, grid);
    // Degenerate range: grid should still have valid endpoints
    BOOST_CHECK_GE(n, 1u);
    BOOST_CHECK_CLOSE(grid(0), std::exp2(10.0), 1e-8);
}

BOOST_AUTO_TEST_CASE(adaptive_grid_huge_pts_per_decade) {
    Array grid;
    Real breaks[] = {1e5};
    Real weights[] = {1.0};
    size_t n = adaptive_grid_with_breaks(0.0, 10.0, std::span(breaks), std::span(weights), 1000, grid);
    BOOST_CHECK_LE(n, 256u);
}

BOOST_AUTO_TEST_CASE(adaptive_grid_invalid_breaks) {
    Array grid;
    Real breaks[] = {std::numeric_limits<Real>::quiet_NaN(), -5.0, 1e5};
    Real weights[] = {1.0, 1.0, 1.0};
    // NaN and negative break values should be skipped without crashing
    size_t n = adaptive_grid_with_breaks(0.0, 40.0, std::span(breaks), std::span(weights), 4, grid);
    BOOST_CHECK_GT(n, 1u);
    BOOST_CHECK_CLOSE(grid(0), std::exp2(0.0), 1e-8);
    BOOST_CHECK_CLOSE(grid(n - 1), std::exp2(40.0), 1e-4);

    // Only the valid break (1e5) should appear; NaN and negative should not
    // Verify the grid is monotonically increasing and finite
    for (size_t i = 0; i < n; ++i) {
        BOOST_CHECK(std::isfinite(grid(i)));
        if (i > 0) {
            BOOST_CHECK_GE(grid(i), grid(i - 1));
        }
    }
}

// ---------------------------------------------------------------------------
// merge_grids
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(merge_grids_sorted) {
    Array a = {5.0, 1.0, 3.0};
    // merge_grids assumes sorted input, so sort first
    Array a_sorted = {1.0, 3.0, 5.0};
    Array b_sorted = {2.0, 4.0, 6.0};
    Array merged = merge_grids(a_sorted, b_sorted);
    for (size_t i = 1; i < merged.size(); ++i) {
        BOOST_CHECK_LE(merged(i - 1), merged(i));
    }
}

BOOST_AUTO_TEST_CASE(merge_grids_unique) {
    Array a = {1.0, 2.0, 3.0};
    Array b = {2.0, 3.0, 4.0};
    Array merged = merge_grids(a, b);
    for (size_t i = 1; i < merged.size(); ++i) {
        BOOST_CHECK_NE(merged(i - 1), merged(i));
    }
}

BOOST_AUTO_TEST_CASE(merge_grids_disjoint) {
    Array a = {1.0, 3.0, 5.0};
    Array b = {2.0, 4.0, 6.0};
    Array merged = merge_grids(a, b);
    BOOST_REQUIRE_EQUAL(merged.size(), 6u);
    BOOST_CHECK_CLOSE(merged(0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(merged(1), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(merged(2), 3.0, 1e-10);
    BOOST_CHECK_CLOSE(merged(3), 4.0, 1e-10);
    BOOST_CHECK_CLOSE(merged(4), 5.0, 1e-10);
    BOOST_CHECK_CLOSE(merged(5), 6.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(merge_grids_overlapping) {
    Array a = {1.0, 2.0, 3.0};
    Array b = {2.0, 3.0, 4.0};
    Array merged = merge_grids(a, b);
    BOOST_REQUIRE_EQUAL(merged.size(), 4u);
    BOOST_CHECK_CLOSE(merged(0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(merged(1), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(merged(2), 3.0, 1e-10);
    BOOST_CHECK_CLOSE(merged(3), 4.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(merge_grids_empty_first) {
    Array a = Array::from_shape({0});
    Array b = {2.0, 4.0, 6.0};
    Array merged = merge_grids(a, b);
    BOOST_REQUIRE_EQUAL(merged.size(), 3u);
    BOOST_CHECK_CLOSE(merged(0), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(merged(1), 4.0, 1e-10);
    BOOST_CHECK_CLOSE(merged(2), 6.0, 1e-10);
}

// ---------------------------------------------------------------------------
// structure_weight
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(structure_weight_one) {
    // Gamma=1: structure_weight = 1 * sqrt(|0*1|) = 0
    BOOST_CHECK_SMALL(structure_weight(1.0), 1e-15);
}

BOOST_AUTO_TEST_CASE(structure_weight_positive) {
    // Gamma=2: 2 * sqrt(|1*2|) = 2*sqrt(2) > 0
    BOOST_CHECK_GT(structure_weight(2.0), 0.0);
    BOOST_CHECK_CLOSE(structure_weight(2.0), 2.0 * std::sqrt(2.0), 1e-10);
}

BOOST_AUTO_TEST_CASE(structure_weight_monotonic) {
    Real prev = structure_weight(1.5);
    for (Real Gamma = 2.0; Gamma <= 100.0; Gamma += 1.0) {
        Real val = structure_weight(Gamma);
        BOOST_CHECK_GT(val, prev);
        prev = val;
    }
}

// ---------------------------------------------------------------------------
// refined_time_grid
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(refined_time_grid_count) {
    Array grid = refined_time_grid(100.0, 1e6, 1e4, 100);
    BOOST_CHECK_EQUAL(grid.size(), 100u);
}

BOOST_AUTO_TEST_CASE(refined_time_grid_bounds) {
    const Real t_start = 100.0;
    const Real t_end = 1e6;
    Array grid = refined_time_grid(t_start, t_end, 1e4, 100);
    // First element should be close to t_start, last close to t_end
    BOOST_CHECK_CLOSE(grid(0), t_start, 1.0);
    BOOST_CHECK_CLOSE(grid(99), t_end, 1.0);
}

BOOST_AUTO_TEST_CASE(refined_time_grid_denser_around_refine) {
    const Real t_start = 100.0;
    const Real t_end = 1e8;
    const Real t_refine = 1e5;
    Array grid = refined_time_grid(t_start, t_end, t_refine, 200);

    // Measure spacing in log space around the refinement region vs far away
    Real near_spacing = 0;
    size_t near_count = 0;
    Real far_spacing = 0;
    size_t far_count = 0;

    for (size_t i = 1; i < grid.size(); ++i) {
        if (grid(i) <= 0 || grid(i - 1) <= 0) {
            continue;
        }
        const Real log_mid = 0.5 * (std::log10(grid(i)) + std::log10(grid(i - 1)));
        const Real log_step = std::log10(grid(i)) - std::log10(grid(i - 1));
        if (log_step <= 0) {
            continue;
        }
        const Real log_refine = std::log10(t_refine);
        if (std::fabs(log_mid - log_refine) < 0.5) {
            near_spacing += log_step;
            near_count++;
        } else if (std::fabs(log_mid - log_refine) > 2.0) {
            far_spacing += log_step;
            far_count++;
        }
    }

    if (near_count > 0 && far_count > 0) {
        const Real near_avg = near_spacing / near_count;
        const Real far_avg = far_spacing / far_count;
        BOOST_CHECK_LT(near_avg, far_avg);
    }
}

BOOST_AUTO_TEST_CASE(refined_time_grid_refine_outside_range) {
    // t_refine far above t_end -> falls back to simple logspace
    const Real t_start = 100.0;
    const Real t_end = 1e6;
    const Real t_refine = 1e20; // far outside range
    Array grid = refined_time_grid(t_start, t_end, t_refine, 50);
    BOOST_CHECK_EQUAL(grid.size(), 50u);
    // Should still span the correct range
    BOOST_CHECK_CLOSE(grid(0), t_start, 1.0);
    BOOST_CHECK_CLOSE(grid(49), t_end, 1.0);
    // Should be a simple logspace: uniform ratio
    if (grid.size() >= 3) {
        const Real ratio_0 = grid(1) / grid(0);
        const Real ratio_1 = grid(2) / grid(1);
        BOOST_CHECK_CLOSE(ratio_0, ratio_1, 1.0);
    }
}

// ---------------------------------------------------------------------------
// jump_refinement_grid
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(jump_refinement_grid_near_jump) {
    Array jumps = {0.1};
    Array result = jump_refinement_grid(jumps, 0.01, 1.0, 0.01, 1.0);
    BOOST_CHECK_GT(result.size(), 0u);

    // Check that points cluster near the jump location (0.1)
    size_t near_count = 0;
    for (size_t i = 0; i < result.size(); ++i) {
        if (std::fabs(result(i) - 0.1) < 0.03) {
            near_count++;
        }
    }
    // Most points should be close to the jump
    BOOST_CHECK_GT(near_count, result.size() / 3);

    // Check all points are in range
    for (size_t i = 0; i < result.size(); ++i) {
        BOOST_CHECK_GE(result(i), 0.01);
        BOOST_CHECK_LE(result(i), 1.0);
    }

    // Check sorted
    for (size_t i = 1; i < result.size(); ++i) {
        BOOST_CHECK_LE(result(i - 1), result(i));
    }
}

BOOST_AUTO_TEST_SUITE_END()
