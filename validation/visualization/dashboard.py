#!/usr/bin/env python3
"""VegasAfterglow Visualization Dashboard - generates PDF reports from test results."""

import argparse
import os
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib
matplotlib.use('Agg')

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

DEFAULT_WORKERS = os.cpu_count() or 4

from validation.visualization.common import (ReportBuilder, extract_session_metadata, find_benchmark_file,
                                        find_regression_file, get_runtime_build_info, load_json,
                                        setup_plot_style)
from validation.visualization.benchmark import generate_convergence_pages, plot_benchmark_overview, plot_convergence_summary, plot_error_distribution
from validation.visualization.regression import (plot_dynamics_summary_grid, plot_spectrum_summary_grid,
                                            plot_frequencies_summary_grid,
                                            plot_shock_quantities_combined, plot_frequency_quantities_combined,
                                            plot_spectrum_shapes,
                                            plot_rvs_shock_quantities_combined,
                                            plot_rvs_frequency_quantities_combined)


class ComprehensiveDashboard:
    def __init__(self, output_dir: str = "output"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        setup_plot_style()

    def generate_full_report(self, benchmark_file: Optional[str] = None, regression_file: Optional[str] = None,
                             n_workers: int = 0):
        output_path = self.output_dir / "comprehensive_report.pdf"

        if benchmark_file is None:
            benchmark_file = find_benchmark_file()
        if regression_file is None:
            regression_file = find_regression_file()

        benchmark_data = None
        regression_data = None

        if benchmark_file and Path(benchmark_file).exists():
            benchmark_data = load_json(benchmark_file)
            print(f"Loaded benchmark data from: {benchmark_file}")
        else:
            print("No benchmark data found")

        if regression_file and Path(regression_file).exists():
            regression_data = load_json(regression_file)
            print(f"Loaded regression data from: {regression_file}")
        else:
            print("No regression data found")

        if benchmark_data and benchmark_data.get("configs"):
            metadata = extract_session_metadata(benchmark_data)
        else:
            metadata = get_runtime_build_info()

        builder = ReportBuilder(output_path)
        section_num = 0
        toc_nodes: List[Dict] = []

        # Title page
        builder.add_title_page("Comprehensive Test Report",
                               "Regression Tests \u2022 Benchmarks \u2022 Physics Diagnostics", metadata)

        # TOC placeholder
        builder.set_toc_position()

        # =============================================================
        # Section 1: Regression Tests
        # =============================================================
        if regression_data:
            section_num += 1
            section_node = {"title": f"{section_num}. Regression Tests", "page": builder.current_page, "children": []}

            builder.add_section_header(section_num, "Regression Tests",
                                       "Power-law scaling verification against theoretical predictions")

            print("\nGenerating regression section...")
            print("  - Adding regression guide pages")
            section_node["children"].append({"title": "How to Read", "page": builder.current_page, "children": []})
            builder.add_guide_pages("regression_guide")

            has_rvs = any(k.startswith("rvs_") for k in regression_data)

            # Track detail page numbers for internal links: {viz_key: {medium: page_num}}
            detail_pages = {}
            # Back-link region: top-left title area (clicking title goes back to summary)
            BACK_LINK_RECT = (0.0, 0.92, 1.0, 0.08)  # full width at top

            def _add_medium_pages(viz_key, label_fmt, plot_fn, node, summary_page=None, **plot_kw):
                """Add detail pages and track page for each medium."""
                if viz_key in regression_data:
                    medium_pages = {}
                    pages_added = []
                    for med in ["ISM", "Wind"]:
                        if med in regression_data[viz_key]:
                            label = label_fmt.format(medium=med)
                            print(f"  - {label}")
                            node["children"].append({"title": label, "page": builder.current_page, "children": []})
                            medium_pages[med] = builder.current_page
                            pages_added.append(builder.current_page)
                            builder.add_fig(plot_fn(regression_data[viz_key], med, **plot_kw))
                    if medium_pages:
                        detail_pages[viz_key] = medium_pages
                    # Register back-links from detail pages to summary
                    if summary_page and pages_added:
                        for page in pages_added:
                            builder.add_internal_links(page, [(BACK_LINK_RECT, summary_page)])

            # --- 1.1 Dynamics ---
            dynamics_node = {"title": "1.1 Dynamics", "page": builder.current_page, "children": []}
            dynamics_summary_page = None
            dynamics_section_rects = {}

            if "shock_grid" in regression_data or (has_rvs and any(
                    f"rvs_shock_grid_{r}" in regression_data for r in ["thin", "thick"])):
                print("  - Dynamics summary")
                dynamics_node["children"].append({"title": "Dynamics Summary", "page": builder.current_page, "children": []})
                dynamics_summary_page = builder.current_page
                fig, dynamics_section_rects = plot_dynamics_summary_grid(regression_data)
                builder.add_fig(fig)

            _add_medium_pages("viz_shock_dynamics", "Forward Shock ({medium})",
                              plot_shock_quantities_combined, dynamics_node,
                              summary_page=dynamics_summary_page)

            if has_rvs:
                for regime in ["thin", "thick"]:
                    _add_medium_pages(f"viz_rvs_shock_dynamics_{regime}",
                                      f"Reverse Shock {regime.title()} ({{medium}})",
                                      plot_rvs_shock_quantities_combined, dynamics_node,
                                      summary_page=dynamics_summary_page, regime=regime)

            # Register dynamics summary links (per medium)
            if dynamics_summary_page and dynamics_section_rects:
                dyn_links = []
                grid_to_viz = {
                    "shock_grid": "viz_shock_dynamics",
                    "rvs_shock_grid_thin": "viz_rvs_shock_dynamics_thin",
                    "rvs_shock_grid_thick": "viz_rvs_shock_dynamics_thick",
                }
                for grid_key, medium_rects in dynamics_section_rects.items():
                    viz_key = grid_to_viz.get(grid_key)
                    if viz_key and viz_key in detail_pages:
                        viz_medium_pages = detail_pages[viz_key]
                        for medium, rect in medium_rects.items():
                            if medium in viz_medium_pages:
                                dyn_links.append((rect, viz_medium_pages[medium]))
                if dyn_links:
                    builder.add_internal_links(dynamics_summary_page, dyn_links)

            section_node["children"].append(dynamics_node)

            # --- 1.2 Radiation ---
            radiation_node = {"title": "1.2 Radiation", "page": builder.current_page, "children": []}

            # --- 1.2.1 Synchrotron Spectrum Shapes ---
            spectrum_node = {"title": "1.2.1 Synchrotron Spectrum Shapes", "page": builder.current_page, "children": []}
            spectrum_summary_page = None
            spectrum_section_rects = {}

            if "spectrum_grid" in regression_data:
                print("  - Spectrum summary")
                spectrum_node["children"].append({"title": "Spectrum Summary", "page": builder.current_page, "children": []})
                spectrum_summary_page = builder.current_page
                fig, spectrum_section_rects = plot_spectrum_summary_grid(regression_data)
                builder.add_fig(fig)

            spectrum_detail_page = None
            if "viz_spectrum_shapes" in regression_data:
                print("  - Spectrum detail")
                spectrum_node["children"].append({"title": "Spectrum Detail", "page": builder.current_page, "children": []})
                spectrum_detail_page = builder.current_page
                builder.add_fig(plot_spectrum_shapes(regression_data["viz_spectrum_shapes"]))

            # Register spectrum summary link (forward and back)
            if spectrum_summary_page and spectrum_detail_page and "spectrum_grid" in spectrum_section_rects:
                # Spectrum uses "_all" key since it doesn't have ISM/Wind split
                spectrum_rect = spectrum_section_rects["spectrum_grid"].get("_all")
                if spectrum_rect:
                    builder.add_internal_links(spectrum_summary_page, [(spectrum_rect, spectrum_detail_page)])
                # Back-link from detail to summary
                builder.add_internal_links(spectrum_detail_page, [(BACK_LINK_RECT, spectrum_summary_page)])

            radiation_node["children"].append(spectrum_node)

            # --- 1.2.2 Characteristic Frequencies ---
            freq_node = {"title": "1.2.2 Characteristic Frequencies", "page": builder.current_page, "children": []}
            freq_summary_page = None
            freq_section_rects = {}

            if "freq_grid" in regression_data or (has_rvs and any(
                    f"rvs_freq_grid_{r}" in regression_data for r in ["thin", "thick"])):
                print("  - Frequencies summary")
                freq_node["children"].append({"title": "Frequencies Summary", "page": builder.current_page, "children": []})
                freq_summary_page = builder.current_page
                fig, freq_section_rects = plot_frequencies_summary_grid(regression_data)
                builder.add_fig(fig)

            _add_medium_pages("viz_frequencies", "Forward Shock ({medium})",
                              plot_frequency_quantities_combined, freq_node,
                              summary_page=freq_summary_page)

            if has_rvs:
                for regime in ["thin", "thick"]:
                    _add_medium_pages(f"viz_rvs_frequencies_{regime}",
                                      f"Reverse Shock {regime.title()} ({{medium}})",
                                      plot_rvs_frequency_quantities_combined, freq_node,
                                      summary_page=freq_summary_page, regime=regime)

            # Register frequencies summary links (per medium)
            if freq_summary_page and freq_section_rects:
                freq_links = []
                grid_to_viz = {
                    "freq_grid": "viz_frequencies",
                    "rvs_freq_grid_thin": "viz_rvs_frequencies_thin",
                    "rvs_freq_grid_thick": "viz_rvs_frequencies_thick",
                }
                for grid_key, medium_rects in freq_section_rects.items():
                    viz_key = grid_to_viz.get(grid_key)
                    if viz_key and viz_key in detail_pages:
                        viz_medium_pages = detail_pages[viz_key]
                        for medium, rect in medium_rects.items():
                            if medium in viz_medium_pages:
                                freq_links.append((rect, viz_medium_pages[medium]))
                if freq_links:
                    builder.add_internal_links(freq_summary_page, freq_links)

            radiation_node["children"].append(freq_node)
            section_node["children"].append(radiation_node)
            toc_nodes.append(section_node)

        # =============================================================
        # Section 2: Benchmark Results
        # =============================================================
        if benchmark_data and benchmark_data.get("configs"):
            section_num += 1
            bench_node = {"title": f"{section_num}. Benchmark Results", "page": builder.current_page, "children": []}

            builder.add_section_header(section_num, "Benchmark Results",
                                       "Performance timing and resolution convergence analysis")

            print("  - Adding benchmark guide pages")
            bench_node["children"].append({"title": "How to Read", "page": builder.current_page, "children": []})
            builder.add_guide_pages("benchmark_guide")

            configs = benchmark_data.get("configs", [])

            print(f"\nGenerating benchmark section ({len(configs)} configs)...")
            print("  - Overview plots")

            on_axis_configs = [c for c in configs if c.get("theta_obs_ratio", 0) <= 1]
            off_axis_configs = [c for c in configs if c.get("theta_obs_ratio", 0) > 1]

            if on_axis_configs:
                print("    - On-axis overview (θ_v/θ_c≤1)")
                bench_node["children"].append({"title": "On-axis Overview", "page": builder.current_page, "children": []})
                builder.add_fig(plot_benchmark_overview(benchmark_data, angle_filter="on-axis"))

            if off_axis_configs:
                print("    - Off-axis overview (θ_v/θ_c>1)")
                bench_node["children"].append({"title": "Off-axis Overview", "page": builder.current_page, "children": []})
                builder.add_fig(plot_benchmark_overview(benchmark_data, angle_filter="off-axis"))

            has_convergence = any(c.get("phi_convergence") for c in configs)
            if has_convergence:
                print("  - Convergence summary")
                bench_node["children"].append({"title": "Convergence Summary", "page": builder.current_page, "children": []})
                summary_page = builder.current_page
                fig, model_id_map, models_list, cell_rects = plot_convergence_summary(benchmark_data)
                builder.add_fig(fig)

                print("  - Error distribution")
                bench_node["children"].append({"title": "Error Distribution", "page": builder.current_page, "children": []})
                builder.add_fig(plot_error_distribution(benchmark_data))

                # Generate convergence pages (parallel or sequential)
                detail_start_page = builder.current_page
                convergence_pdf_files = generate_convergence_pages(models_list, n_workers)
                builder.add_fig_files(convergence_pdf_files)

                # Register jump links from summary grid cells to detail pages
                links = []
                for i, model in enumerate(models_list):
                    if i in cell_rects and i < len(convergence_pdf_files):
                        links.append((cell_rects[i], detail_start_page + i))
                if links:
                    builder.add_internal_links(summary_page, links)

                # Register back-links from each detail page to summary
                # Back-link region: top title area
                BACK_LINK_RECT = (0.0, 0.92, 1.0, 0.08)
                for i in range(len(convergence_pdf_files)):
                    detail_page = detail_start_page + i
                    builder.add_internal_links(detail_page, [(BACK_LINK_RECT, summary_page)])

            toc_nodes.append(bench_node)

        # Save final document
        print(f"\n{'='*70}")
        builder.save(toc_nodes=toc_nodes if toc_nodes else None)
        print(f"{'='*70}")
        return output_path

    def generate_benchmark_report(self, benchmark_file: str, n_workers: int = 0):
        if not Path(benchmark_file).exists():
            print(f"Benchmark file not found: {benchmark_file}")
            return

        data = load_json(benchmark_file)
        output_path = self.output_dir / "benchmark_report.pdf"

        if data.get("configs"):
            metadata = extract_session_metadata(data)
        else:
            metadata = get_runtime_build_info()

        builder = ReportBuilder(output_path)
        current_page = 1

        builder.add_title_page("Benchmark Report", metadata=metadata)
        current_page += 1

        print("  Adding benchmark guide pages...")
        builder.add_guide_pages("benchmark_guide")
        current_page += 3  # Approximate guide page count

        if data.get("configs"):
            configs = data.get("configs", [])

            on_axis_configs = [c for c in configs if c.get("theta_obs_ratio", 0) <= 1]
            off_axis_configs = [c for c in configs if c.get("theta_obs_ratio", 0) > 1]

            if on_axis_configs:
                builder.add_fig(plot_benchmark_overview(data, angle_filter="on-axis"))
                current_page += 1
            if off_axis_configs:
                builder.add_fig(plot_benchmark_overview(data, angle_filter="off-axis"))
                current_page += 1

            has_convergence = any(c.get("phi_convergence") for c in configs)
            if has_convergence:
                summary_page = current_page
                fig, model_id_map, models_list, cell_rects = plot_convergence_summary(data)
                builder.add_fig(fig)
                current_page += 1

                builder.add_fig(plot_error_distribution(data))
                current_page += 1

                detail_start_page = current_page
                convergence_pdf_files = generate_convergence_pages(models_list, n_workers)
                builder.add_fig_files(convergence_pdf_files)
                current_page += len(convergence_pdf_files)

                links = []
                for i in range(len(models_list)):
                    if i in cell_rects and i < len(convergence_pdf_files):
                        links.append((cell_rects[i], detail_start_page + i))
                if links:
                    builder.add_internal_links(summary_page, links)

        builder.save()

    def generate_regression_report(self, regression_file: str):
        if not Path(regression_file).exists():
            print(f"Regression file not found: {regression_file}")
            return

        data = load_json(regression_file)
        output_path = self.output_dir / "regression_report.pdf"

        metadata = data.get("metadata", {})
        if not metadata.get("Version"):
            metadata = get_runtime_build_info()

        summary = data.get("summary", {})
        by_model = summary.get("by_model", {})
        n_total = sum(by_model.get(m, {}).get("pass", 0) + by_model.get(m, {}).get("fail", 0) for m in by_model)
        n_pass = sum(by_model.get(m, {}).get("pass", 0) for m in by_model)

        has_rvs = any(k.startswith("rvs_") for k in data)

        builder = ReportBuilder(output_path)
        subtitle = f"{n_pass}/{n_total} tests passed" if n_total > 0 else ""
        builder.add_title_page("Regression Test Report", subtitle=subtitle, metadata=metadata)

        print("  Adding regression guide pages...")
        builder.add_guide_pages("regression_guide")

        # --- Dynamics ---
        if "shock_grid" in data or (has_rvs and any(
                f"rvs_shock_grid_{r}" in data for r in ["thin", "thick"])):
            print("  Generating dynamics summary...")
            builder.add_fig(plot_dynamics_summary_grid(data))

        def _gen_medium(viz_key, label, plot_fn, **kw):
            if viz_key in data:
                for med in ["ISM", "Wind"]:
                    if med in data[viz_key]:
                        print(f"  Generating {label} ({med})...")
                        builder.add_fig(plot_fn(data[viz_key], med, **kw))

        _gen_medium("viz_shock_dynamics", "shock dynamics", plot_shock_quantities_combined)
        if has_rvs:
            for regime in ["thin", "thick"]:
                _gen_medium(f"viz_rvs_shock_dynamics_{regime}", f"reverse shock dynamics {regime}",
                            plot_rvs_shock_quantities_combined, regime=regime)

        # --- Spectrum Shapes ---
        if "spectrum_grid" in data:
            print("  Generating spectrum summary...")
            builder.add_fig(plot_spectrum_summary_grid(data))

        if "viz_spectrum_shapes" in data:
            print("  Generating spectrum detail...")
            builder.add_fig(plot_spectrum_shapes(data["viz_spectrum_shapes"]))

        # --- Frequencies ---
        if "freq_grid" in data or (has_rvs and any(
                f"rvs_freq_grid_{r}" in data for r in ["thin", "thick"])):
            print("  Generating frequencies summary...")
            builder.add_fig(plot_frequencies_summary_grid(data))

        _gen_medium("viz_frequencies", "frequencies", plot_frequency_quantities_combined)
        if has_rvs:
            for regime in ["thin", "thick"]:
                _gen_medium(f"viz_rvs_frequencies_{regime}", f"reverse shock frequencies {regime}",
                            plot_rvs_frequency_quantities_combined, regime=regime)

        builder.save()


def main():
    parser = argparse.ArgumentParser(description="VegasAfterglow Visualization Dashboard")
    parser.add_argument("--full", action="store_true", help="Full report")
    parser.add_argument("--benchmark", action="store_true", help="Benchmark-only")
    parser.add_argument("--regression", action="store_true", help="Regression-only")
    parser.add_argument("-j", "--parallel", type=int, default=DEFAULT_WORKERS, metavar="N",
                        help=f"Parallel workers (default: {DEFAULT_WORKERS})")
    parser.add_argument("--benchmark-file", type=str, help="Benchmark JSON path")
    parser.add_argument("--regression-file", type=str, help="Regression JSON path")
    parser.add_argument("--output", type=str, default=str(Path(__file__).parent.parent), help="Output dir")
    args = parser.parse_args()
    if not any([args.full, args.benchmark, args.regression]):
        args.full = True
    print("=" * 70 + "\nVegasAfterglow Visualization Dashboard\n" + "=" * 70)

    dashboard = ComprehensiveDashboard(output_dir=args.output)

    benchmark_file = args.benchmark_file or find_benchmark_file()
    regression_file = args.regression_file or find_regression_file()

    if benchmark_file:
        print(f"Benchmark data: {benchmark_file}")
    if regression_file:
        print(f"Regression data: {regression_file}")
    if args.parallel > 0:
        print(f"Parallel workers: {args.parallel}")

    if args.full:
        dashboard.generate_full_report(benchmark_file=benchmark_file, regression_file=regression_file, n_workers=args.parallel)

    if args.benchmark and benchmark_file:
        dashboard.generate_benchmark_report(benchmark_file, n_workers=args.parallel)

    if args.regression and regression_file:
        dashboard.generate_regression_report(regression_file)


if __name__ == "__main__":
    main()
