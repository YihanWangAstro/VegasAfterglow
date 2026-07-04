#include <boost/test/unit_test.hpp>
#include <cmath>
#include <span>

#include "core/grid-refinement.h"
#include "util/macros.h"

BOOST_AUTO_TEST_SUITE(Mesh)

// ---------------------------------------------------------------------------
// boundary_to_center (non-template, returns Array)
// ---------------------------------------------------------------------------

// Centers are the arithmetic midpoints of adjacent boundaries: {1,3,5} -> {2,4}.
BOOST_AUTO_TEST_CASE(boundary_to_center_linear) {
    Array boundary = {1.0, 3.0, 5.0};
    Array center = boundary_to_center(boundary);
    BOOST_REQUIRE_EQUAL(center.size(), 2u);
    BOOST_CHECK_CLOSE(center(0), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(center(1), 4.0, 1e-10);
}

// A two-point boundary yields a single center at the arithmetic midpoint.
BOOST_AUTO_TEST_CASE(boundary_to_center_single) {
    Array boundary = {1.0, 3.0};
    Array center = boundary_to_center(boundary);
    BOOST_REQUIRE_EQUAL(center.size(), 1u);
    BOOST_CHECK_CLOSE(center(0), 2.0, 1e-10);
}

// ---------------------------------------------------------------------------
// boundary_to_center_log (non-template, returns Array)
// ---------------------------------------------------------------------------

// Log-space centers are the geometric means of adjacent boundaries: {1,4,16} -> {2,8}.
BOOST_AUTO_TEST_CASE(boundary_to_center_log_values) {
    Array boundary = {1.0, 4.0, 16.0};
    Array center = boundary_to_center_log(boundary);
    BOOST_REQUIRE_EQUAL(center.size(), 2u);
    BOOST_CHECK_CLOSE(center(0), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(center(1), 8.0, 1e-10);
}

// A two-point boundary yields a single geometric-mean center: {1,9} -> {3}.
BOOST_AUTO_TEST_CASE(boundary_to_center_log_single) {
    Array boundary = {1.0, 9.0};
    Array center = boundary_to_center_log(boundary);
    BOOST_REQUIRE_EQUAL(center.size(), 1u);
    BOOST_CHECK_CLOSE(center(0), 3.0, 1e-10);
}

// ---------------------------------------------------------------------------
// boundary_to_center (template version) edge cases
// ---------------------------------------------------------------------------

// An empty boundary array yields an empty center array without crashing.
BOOST_AUTO_TEST_CASE(boundary_to_center_empty) {
    // Empty boundary should not crash
    Array boundary = Array::from_shape({0});
    Array center = boundary_to_center(boundary);
    BOOST_CHECK_EQUAL(center.size(), 0u);
}

// When center.size() + 1 != boundary.size(), the template overload returns early and leaves center untouched.
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
// log2space
// ---------------------------------------------------------------------------

// Point count is ceil(decades * pts_per_decade) + 1: 10 decades at 4 per decade gives 41 points.
BOOST_AUTO_TEST_CASE(log2space_count) {
    Array arr;
    log2space(0.0, 10.0 * log2_10, 4, arr);
    // 10 decades * 4 points per decade = 40, plus 1 -> 41
    const Real decades = 10.0;
    const size_t expected = static_cast<size_t>(std::ceil(decades * 4)) + 1;
    BOOST_CHECK_EQUAL(arr.size(), expected);
}

// First and last grid values equal exp2 of the requested log2 bounds.
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

// Each returned center equals the geometric mean of its log-spaced bin boundaries.
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
    // A simpler check: build boundaries from the logspace and compare. The mesh
    // builds its grid through fast_exp2, so under AFTERGLOW_FAST_MATH the
    // comparison against this std::exp2 reference carries the kernel's ~1e-7
    // relative error COMPOUNDED across the bins (the boundaries are built
    // multiplicatively, so the drift grows linearly with bin index: ~30 bins
    // x 1e-7 = 3e-6 relative). BOOST_CHECK_CLOSE tolerances are in percent.
#ifdef AFTERGLOW_FAST_MATH
    constexpr Real MESH_TOL = 1e-3;
#else
    constexpr Real MESH_TOL = 1e-8;
#endif
    const Real dlg2 = (10.0 - 2.0) / 30.0;
    const Real r = std::exp2(dlg2);
    Real left = std::exp2(2.0);
    for (size_t i = 0; i < 30; ++i) {
        const Real right = left * r;
        const Real geo_center = std::sqrt(left * right);
        BOOST_CHECK_CLOSE(center(i), geo_center, MESH_TOL);
        left = right;
    }
}

// All bin widths are strictly positive.
BOOST_AUTO_TEST_CASE(logspace_boundary_center_width) {
    Array center;
    Array bin_width;
    logspace_boundary_center(0.0, 20.0, 50, center, bin_width);
    for (size_t i = 0; i < bin_width.size(); ++i) {
        BOOST_CHECK_GT(bin_width(i), 0.0);
    }
}

// Requesting zero bins yields empty center and bin_width arrays.
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

// First and last grid points land on exp2(lg2_min) and exp2(lg2_max) even with breaks present.
BOOST_AUTO_TEST_CASE(adaptive_grid_endpoints) {
    Array grid;
    Real breaks[] = {1e5, 1e8};
    Real weights[] = {1.0, 1.0};
    size_t n = adaptive_grid_with_breaks(0.0, 40.0, std::span(breaks), std::span(weights), 4, grid);
    BOOST_REQUIRE_GT(n, 1u);
    BOOST_CHECK_CLOSE(grid(0), std::exp2(0.0), 1e-8);
    BOOST_CHECK_CLOSE(grid(n - 1), std::exp2(40.0), 1e-4);
}

// Every requested break frequency appears in the grid to within a factor of 2 (1 in log2).
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

// Grid size is capped at 256 points even when pts_per_decade requests more.
BOOST_AUTO_TEST_CASE(adaptive_grid_max_256) {
    Array grid;
    Real breaks[] = {1e5};
    Real weights[] = {1.0};
    size_t n = adaptive_grid_with_breaks(0.0, 40.0, std::span(breaks), std::span(weights), 100, grid);
    BOOST_CHECK_LE(n, 256u);
    BOOST_CHECK_LE(grid.size(), 256u);
}

// Average log2 spacing within 2 log2 units of the break is finer than more than 10 log2 units
// away, compared only when both regions contain grid intervals.
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

// With empty break spans the grid is still built and spans exp2(lg2_min) to exp2(lg2_max).
BOOST_AUTO_TEST_CASE(adaptive_grid_no_breaks) {
    Array grid;
    std::span<const Real> no_breaks;
    std::span<const Real> no_weights;
    size_t n = adaptive_grid_with_breaks(0.0, 20.0, no_breaks, no_weights, 4, grid);
    BOOST_CHECK_GT(n, 1u);
    BOOST_CHECK_CLOSE(grid(0), std::exp2(0.0), 1e-8);
    BOOST_CHECK_CLOSE(grid(n - 1), std::exp2(20.0), 1e-4);
}

// Degenerate range (lg2_min == lg2_max) still yields at least one point, located at exp2(lg2_min).
BOOST_AUTO_TEST_CASE(adaptive_grid_min_equals_max) {
    Array grid;
    std::span<const Real> no_breaks;
    std::span<const Real> no_weights;
    size_t n = adaptive_grid_with_breaks(10.0, 10.0, no_breaks, no_weights, 4, grid);
    // Degenerate range: grid should still have valid endpoints
    BOOST_CHECK_GE(n, 1u);
    BOOST_CHECK_CLOSE(grid(0), std::exp2(10.0), 1e-8);
}

// The 256-point cap holds even at an extreme density of 1000 points per decade.
BOOST_AUTO_TEST_CASE(adaptive_grid_huge_pts_per_decade) {
    Array grid;
    Real breaks[] = {1e5};
    Real weights[] = {1.0};
    size_t n = adaptive_grid_with_breaks(0.0, 10.0, std::span(breaks), std::span(weights), 1000, grid);
    BOOST_CHECK_LE(n, 256u);
}

// NaN and negative breaks do not crash grid construction: the result stays finite,
// non-decreasing, and spans exp2(lg2_min) to exp2(lg2_max).
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

// Merging two sorted grids produces a non-decreasing result.
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

// Values shared between the inputs appear only once: no adjacent duplicates in the merged grid.
BOOST_AUTO_TEST_CASE(merge_grids_unique) {
    Array a = {1.0, 2.0, 3.0};
    Array b = {2.0, 3.0, 4.0};
    Array merged = merge_grids(a, b);
    for (size_t i = 1; i < merged.size(); ++i) {
        BOOST_CHECK_NE(merged(i - 1), merged(i));
    }
}

// Disjoint inputs {1,3,5} and {2,4,6} interleave into the sorted 6-element union.
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

// Overlapping inputs {1,2,3} and {2,3,4} merge to the 4-element union {1,2,3,4}.
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

// Merging an empty first grid with a non-empty second returns the second grid unchanged.
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

// Gamma = 1 (no bulk motion) gives zero structure weight: 1 * sqrt(0 * 1) = 0.
BOOST_AUTO_TEST_CASE(structure_weight_one) {
    // Gamma=1: structure_weight = 1 * sqrt(|0*1|) = 0
    BOOST_CHECK_SMALL(structure_weight(1.0), 1e-15);
}

// Gamma = 2 gives a positive weight matching the closed form Gamma * sqrt((Gamma-1) * Gamma) = 2*sqrt(2).
BOOST_AUTO_TEST_CASE(structure_weight_positive) {
    // Gamma=2: 2 * sqrt(|1*2|) = 2*sqrt(2) > 0
    BOOST_CHECK_GT(structure_weight(2.0), 0.0);
    BOOST_CHECK_CLOSE(structure_weight(2.0), 2.0 * std::sqrt(2.0), 1e-10);
}

// structure_weight increases strictly monotonically with Gamma over [1.5, 100].
BOOST_AUTO_TEST_CASE(structure_weight_monotonic) {
    Real prev = structure_weight(1.5);
    for (Real Gamma = 2.0; Gamma <= 100.0; Gamma += 1.0) {
        Real val = structure_weight(Gamma);
        BOOST_CHECK_GT(val, prev);
        prev = val;
    }
}

// ---------------------------------------------------------------------------
// jump_refinement_grid
// ---------------------------------------------------------------------------

// Refinement points cluster near the jump (over 1/3 of them within 0.03 of theta = 0.1),
// stay within [theta_min, theta_max], and come out sorted.
BOOST_AUTO_TEST_CASE(jump_refinement_grid_near_jump) {
    Array jumps = {0.1};
    Array result = jump_refinement_grid(jumps, 0.01, 1.0, 0.01);
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
