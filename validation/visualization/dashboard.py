#!/usr/bin/env python3
"""VegasAfterglow Visualization Dashboard - generates PDF reports from test results."""

import argparse
import os
import shutil
from io import BytesIO
from pathlib import Path
from typing import Dict, List, Optional

# Default to number of CPUs for parallel execution
DEFAULT_WORKERS = os.cpu_count() or 4

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend before importing pyplot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader, PdfWriter

import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from validation.visualization.common import (create_benchmark_guide_pages, create_regression_guide_pages,
                                        create_section_header, create_toc_page_reportlab,
                                        create_title_page, extract_session_metadata, find_benchmark_file,
                                        find_regression_file, get_runtime_build_info, load_json, setup_plot_style)
from validation.visualization.benchmark import generate_convergence_pages, plot_benchmark_overview, plot_convergence_summary
from validation.visualization.regression import (plot_dynamics_summary_grid, plot_spectrum_summary_grid,
                                            plot_frequencies_summary_grid,
                                            plot_shock_quantities_combined, plot_frequency_quantities_combined,
                                            plot_spectrum_shapes,
                                            plot_rvs_shock_quantities_combined,
                                            plot_rvs_frequency_quantities_combined)


def _create_toc_page_matplotlib(toc_nodes: List[Dict]) -> plt.Figure:
    """Fallback: create TOC page using matplotlib when reportlab is unavailable."""
    fig = plt.figure(figsize=(8.5, 11))
    fig.patch.set_facecolor("white")
    fig.text(0.5, 0.92, "Table of Contents", ha="center", va="top", fontsize=18, fontweight="bold")

    y = 0.85
    line_height = 0.024

    def _render(node, level):
        nonlocal y
        if y < 0.05:
            return
        indent = 0.1 + level * 0.05
        fs = max(8, 12 - level)
        weight = "bold" if level <= 1 else "normal"
        color = "#000000" if level == 0 else ("#333333" if level == 1 else "#666666")
        fig.text(indent, y, node["title"], ha="left", va="top", fontsize=fs, fontweight=weight, color=color)
        fig.text(0.9, y, str(node.get("page", "")), ha="right", va="top", fontsize=fs, color=color)
        y -= line_height * (1.3 if level == 0 else 1.0)
        for child in node.get("children", []):
            _render(child, level + 1)
        if level == 0:
            y -= line_height * 0.3

    for node in toc_nodes:
        _render(node, 0)
    return fig


def _add_hierarchical_bookmarks(writer: PdfWriter, toc_nodes: List[Dict]):
    """Add hierarchical PDF bookmarks from TOC node structure."""
    def _add(node, parent=None):
        page = node.get("page", 1) - 1  # 0-indexed
        bookmark = writer.add_outline_item(node["title"], max(0, page), parent=parent)
        for child in node.get("children", []):
            _add(child, parent=bookmark)
    for node in toc_nodes:
        _add(node)


def _merge_pending_pdfs(output_path: Path, pending_pdf_merges: List[Dict], insert_at: int = 5):
    """Merge pending PDFs into the main output.

    Each merge_info dict can have an 'insert_at' key to specify its insertion position.
    If not specified, uses the default insert_at parameter.
    Merges are processed from highest to lowest position to maintain correct indices.
    """
    if not pending_pdf_merges:
        return

    print("Merging parallel-generated pages...")

    # Group merges by their insert position
    merges_by_position = {}
    for merge_info in pending_pdf_merges:
        pos = merge_info.get("insert_at", insert_at)
        if pos not in merges_by_position:
            merges_by_position[pos] = []
        merges_by_position[pos].append(merge_info)

    # Process from highest to lowest position (to maintain correct page indices)
    sorted_positions = sorted(merges_by_position.keys(), reverse=True)

    reader = PdfReader(str(output_path))
    main_pages = list(reader.pages)
    total_merged = 0

    for pos in sorted_positions:
        merge_infos = merges_by_position[pos]

        # Collect all pages to insert at this position
        insert_pages = []
        for merge_info in merge_infos:
            for pdf_file in merge_info["files"]:
                merge_reader = PdfReader(pdf_file)
                insert_pages.extend(merge_reader.pages)

        actual_pos = min(pos, len(main_pages))

        # Insert pages at this position
        new_pages = main_pages[:actual_pos] + insert_pages + main_pages[actual_pos:]
        main_pages = new_pages
        total_merged += len(insert_pages)
        print(f"  Merged {len(insert_pages)} pages at position {actual_pos}")

    # Write final result
    writer = PdfWriter()
    for page in main_pages:
        writer.add_page(page)

    temp_output = str(output_path) + ".tmp"
    with open(temp_output, "wb") as f:
        writer.write(f)
    shutil.move(temp_output, str(output_path))

    # Cleanup temp directories
    for merge_info in pending_pdf_merges:
        if merge_info.get("temp_dir"):
            shutil.rmtree(merge_info["temp_dir"], ignore_errors=True)

    print(f"  Total: {total_merged} pages merged")


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

        section_num = 0
        pending_pdf_merges = []
        convergence_insert_at = None
        # Hierarchical TOC structure: [{"title", "page", "children": [...]}]
        toc_nodes: List[Dict] = []
        current_page = 1
        actual_page = 0  # 0-indexed position in actual PDF

        def _add_page(fig):
            nonlocal current_page, actual_page
            pdf.savefig(fig)
            plt.close(fig)
            current_page += 1
            actual_page += 1

        with PdfPages(output_path) as pdf:
            create_title_page(pdf, "Comprehensive Test Report", "Regression Tests \u2022 Benchmarks \u2022 Physics Diagnostics", metadata)
            current_page += 1
            actual_page += 1

            # Placeholder for TOC page (will be replaced)
            fig_toc_placeholder = plt.figure(figsize=(8.5, 11))
            fig_toc_placeholder.patch.set_facecolor("white")
            pdf.savefig(fig_toc_placeholder)
            plt.close(fig_toc_placeholder)
            current_page += 1
            actual_page += 1

            # =============================================================
            # Section 1: Regression Tests
            # =============================================================
            if regression_data:
                section_num += 1
                section_node = {"title": f"{section_num}. Regression Tests", "page": current_page, "children": []}

                create_section_header(pdf, section_num, "Regression Tests", "Power-law scaling verification against theoretical predictions")
                current_page += 1
                actual_page += 1

                print("\nGenerating regression section...")
                print("  - Adding regression guide pages")
                section_node["children"].append({"title": "How to Read", "page": current_page, "children": []})
                guide_pages, guide_merge_info = create_regression_guide_pages(pdf)
                if guide_merge_info:
                    guide_merge_info["insert_at"] = actual_page
                    pending_pdf_merges.append(guide_merge_info)
                current_page += guide_pages

                has_rvs = any(k.startswith("rvs_") for k in regression_data)

                # --- 1.1 Dynamics ---
                dynamics_node = {"title": "1.1 Dynamics", "page": current_page, "children": []}

                if "shock_grid" in regression_data or (has_rvs and any(
                        f"rvs_shock_grid_{r}" in regression_data for r in ["thin", "thick"])):
                    print("  - Dynamics summary")
                    dynamics_node["children"].append({"title": "Dynamics Summary", "page": current_page, "children": []})
                    _add_page(plot_dynamics_summary_grid(regression_data))

                # Forward shock dynamics detail
                if "viz_shock_dynamics" in regression_data:
                    for medium in ["ISM", "Wind"]:
                        if medium in regression_data["viz_shock_dynamics"]:
                            label = f"Forward Shock ({medium})"
                            print(f"  - {label}")
                            dynamics_node["children"].append({"title": label, "page": current_page, "children": []})
                            _add_page(plot_shock_quantities_combined(regression_data["viz_shock_dynamics"], medium))

                # Reverse shock dynamics detail (thin then thick)
                if has_rvs:
                    for regime in ["thin", "thick"]:
                        for medium in ["ISM", "Wind"]:
                            viz_key = f"viz_rvs_shock_dynamics_{regime}"
                            if viz_key in regression_data and medium in regression_data[viz_key]:
                                label = f"Reverse Shock {regime.title()} ({medium})"
                                print(f"  - {label}")
                                dynamics_node["children"].append({"title": label, "page": current_page, "children": []})
                                _add_page(plot_rvs_shock_quantities_combined(regression_data[viz_key], medium, regime))

                section_node["children"].append(dynamics_node)

                # --- 1.2 Radiation ---
                radiation_node = {"title": "1.2 Radiation", "page": current_page, "children": []}

                # --- 1.2.1 Synchrotron Spectrum Shapes ---
                spectrum_node = {"title": "1.2.1 Synchrotron Spectrum Shapes", "page": current_page, "children": []}

                if "spectrum_grid" in regression_data:
                    print("  - Spectrum summary")
                    spectrum_node["children"].append({"title": "Spectrum Summary", "page": current_page, "children": []})
                    _add_page(plot_spectrum_summary_grid(regression_data))

                if "viz_spectrum_shapes" in regression_data:
                    print("  - Spectrum detail")
                    spectrum_node["children"].append({"title": "Spectrum Detail", "page": current_page, "children": []})
                    _add_page(plot_spectrum_shapes(regression_data["viz_spectrum_shapes"]))

                radiation_node["children"].append(spectrum_node)

                # --- 1.2.2 Characteristic Frequencies ---
                freq_node = {"title": "1.2.2 Characteristic Frequencies", "page": current_page, "children": []}

                if "freq_grid" in regression_data or (has_rvs and any(
                        f"rvs_freq_grid_{r}" in regression_data for r in ["thin", "thick"])):
                    print("  - Frequencies summary")
                    freq_node["children"].append({"title": "Frequencies Summary", "page": current_page, "children": []})
                    _add_page(plot_frequencies_summary_grid(regression_data))

                # Forward shock frequencies detail
                if "viz_frequencies" in regression_data:
                    for medium in ["ISM", "Wind"]:
                        if medium in regression_data["viz_frequencies"]:
                            label = f"Forward Shock ({medium})"
                            print(f"  - Freq {label}")
                            freq_node["children"].append({"title": label, "page": current_page, "children": []})
                            _add_page(plot_frequency_quantities_combined(regression_data["viz_frequencies"], medium))

                # Reverse shock frequencies detail (thin then thick)
                if has_rvs:
                    for regime in ["thin", "thick"]:
                        for medium in ["ISM", "Wind"]:
                            freq_key = f"viz_rvs_frequencies_{regime}"
                            if freq_key in regression_data and medium in regression_data[freq_key]:
                                label = f"Reverse Shock {regime.title()} ({medium})"
                                print(f"  - Freq {label}")
                                freq_node["children"].append({"title": label, "page": current_page, "children": []})
                                _add_page(plot_rvs_frequency_quantities_combined(regression_data[freq_key], medium, regime))

                radiation_node["children"].append(freq_node)
                section_node["children"].append(radiation_node)
                toc_nodes.append(section_node)

            # =============================================================
            # Section 2: Benchmark Results
            # =============================================================
            if benchmark_data and benchmark_data.get("configs"):
                section_num += 1
                bench_node = {"title": f"{section_num}. Benchmark Results", "page": current_page, "children": []}

                create_section_header(pdf, section_num, "Benchmark Results", "Performance timing and resolution convergence analysis")
                current_page += 1
                actual_page += 1

                print("  - Adding benchmark guide pages")
                bench_node["children"].append({"title": "How to Read", "page": current_page, "children": []})
                guide_pages, guide_merge_info = create_benchmark_guide_pages(pdf)
                if guide_merge_info:
                    guide_merge_info["insert_at"] = actual_page
                    pending_pdf_merges.insert(0, guide_merge_info)
                current_page += guide_pages

                configs = benchmark_data.get("configs", [])

                print(f"\nGenerating benchmark section ({len(configs)} configs)...")
                print("  - Overview plots")

                on_axis_configs = [c for c in configs if c.get("theta_obs_ratio", 0) <= 1]
                off_axis_configs = [c for c in configs if c.get("theta_obs_ratio", 0) > 1]

                if on_axis_configs:
                    print("    - On-axis overview (\u03b8_v/\u03b8_c\u22641)")
                    bench_node["children"].append({"title": "On-axis Overview", "page": current_page, "children": []})
                    _add_page(plot_benchmark_overview(benchmark_data, angle_filter="on-axis"))

                if off_axis_configs:
                    print("    - Off-axis overview (\u03b8_v/\u03b8_c>1)")
                    bench_node["children"].append({"title": "Off-axis Overview", "page": current_page, "children": []})
                    _add_page(plot_benchmark_overview(benchmark_data, angle_filter="off-axis"))

                has_convergence = any(c.get("phi_convergence") for c in configs)
                if has_convergence:
                    print("  - Convergence summary")
                    bench_node["children"].append({"title": "Convergence Summary", "page": current_page, "children": []})
                    fig, model_id_map, models_list = plot_convergence_summary(benchmark_data)
                    _add_page(fig)
                    convergence_insert_at = actual_page
                    pending_pdf_merges.extend(generate_convergence_pages(pdf, models_list, n_workers))

                toc_nodes.append(bench_node)

        # Merge convergence detail pages at the correct position
        if pending_pdf_merges and convergence_insert_at is not None:
            _merge_pending_pdfs(output_path, pending_pdf_merges, insert_at=convergence_insert_at)

        # Replace TOC placeholder page with actual TOC
        if toc_nodes:
            print("  Adding table of contents...")
            import tempfile
            toc_temp = Path(tempfile.mkdtemp()) / "toc.pdf"
            toc_pages = create_toc_page_reportlab(toc_nodes, toc_temp)

            reader = PdfReader(str(output_path))
            writer = PdfWriter()

            # Page 1: Title
            writer.add_page(reader.pages[0])

            # Page 2+: TOC (replace placeholder)
            if toc_pages > 0 and toc_temp.exists():
                toc_reader = PdfReader(str(toc_temp))
                for tp in toc_reader.pages:
                    writer.add_page(tp)
            else:
                # Fallback to matplotlib TOC
                toc_fig = _create_toc_page_matplotlib(toc_nodes)
                toc_buffer = BytesIO()
                toc_fig.savefig(toc_buffer, format="pdf")
                plt.close(toc_fig)
                toc_buffer.seek(0)
                toc_reader = PdfReader(toc_buffer)
                writer.add_page(toc_reader.pages[0])

            # Remaining pages (skip placeholder at index 1)
            for i in range(2, len(reader.pages)):
                writer.add_page(reader.pages[i])

            # Add hierarchical PDF bookmarks
            _add_hierarchical_bookmarks(writer, toc_nodes)

            temp_output = str(output_path) + ".tmp"
            with open(temp_output, "wb") as f:
                writer.write(f)
            shutil.move(temp_output, str(output_path))

            # Cleanup
            if toc_temp.exists():
                shutil.rmtree(toc_temp.parent, ignore_errors=True)

        print(f"\n{'='*70}")
        print(f"Comprehensive report saved to: {output_path}")
        print(f"{'='*70}")
        return output_path

    def generate_benchmark_report(self, benchmark_file: str, n_workers: int = 0):
        if not Path(benchmark_file).exists():
            print(f"Benchmark file not found: {benchmark_file}")
            return

        data = load_json(benchmark_file)
        output_path = self.output_dir / "benchmark_report.pdf"
        pending_pdf_merges = []

        with PdfPages(output_path) as pdf:
            if data.get("configs"):
                metadata = extract_session_metadata(data)
            else:
                metadata = get_runtime_build_info()

            create_title_page(pdf, "Benchmark Report", metadata=metadata)

            # Add guide pages
            print("  Adding benchmark guide pages...")
            guide_pages, guide_merge_info = create_benchmark_guide_pages(pdf)
            if guide_merge_info:
                guide_merge_info["insert_at"] = 1  # After title page
                pending_pdf_merges.insert(0, guide_merge_info)

            if data.get("configs"):
                configs = data.get("configs", [])

                on_axis_configs = [c for c in configs if c.get("theta_obs_ratio", 0) <= 1]
                off_axis_configs = [c for c in configs if c.get("theta_obs_ratio", 0) > 1]

                if on_axis_configs:
                    fig = plot_benchmark_overview(data, angle_filter="on-axis")
                    pdf.savefig(fig)
                    plt.close(fig)

                if off_axis_configs:
                    fig = plot_benchmark_overview(data, angle_filter="off-axis")
                    pdf.savefig(fig)
                    plt.close(fig)

                has_convergence = any(c.get("phi_convergence") for c in configs)
                if has_convergence:
                    fig, model_id_map, models_list = plot_convergence_summary(data)
                    pdf.savefig(fig)
                    plt.close(fig)

                    pending_pdf_merges.extend(generate_convergence_pages(pdf, models_list, n_workers))

        _merge_pending_pdfs(output_path, pending_pdf_merges)
        print(f"Benchmark report saved to: {output_path}")

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

        pending_pdf_merges = []
        has_rvs = any(k.startswith("rvs_") for k in data)

        with PdfPages(output_path) as pdf:
            subtitle = f"{n_pass}/{n_total} tests passed" if n_total > 0 else ""
            create_title_page(pdf, "Regression Test Report", subtitle=subtitle, metadata=metadata)

            # Add guide pages
            print("  Adding regression guide pages...")
            guide_pages, guide_merge_info = create_regression_guide_pages(pdf)
            if guide_merge_info:
                guide_merge_info["insert_at"] = 1
                pending_pdf_merges.append(guide_merge_info)

            # --- Dynamics ---
            if "shock_grid" in data or (has_rvs and any(
                    f"rvs_shock_grid_{r}" in data for r in ["thin", "thick"])):
                print("  Generating dynamics summary...")
                fig = plot_dynamics_summary_grid(data)
                pdf.savefig(fig)
                plt.close(fig)

            if "viz_shock_dynamics" in data:
                for medium in ["ISM", "Wind"]:
                    if medium in data["viz_shock_dynamics"]:
                        print(f"  Generating shock dynamics ({medium})...")
                        fig = plot_shock_quantities_combined(data["viz_shock_dynamics"], medium)
                        pdf.savefig(fig)
                        plt.close(fig)

            if has_rvs:
                for regime in ["thin", "thick"]:
                    for medium in ["ISM", "Wind"]:
                        viz_key = f"viz_rvs_shock_dynamics_{regime}"
                        if viz_key in data and medium in data[viz_key]:
                            print(f"  Generating reverse shock dynamics {regime} ({medium})...")
                            fig = plot_rvs_shock_quantities_combined(data[viz_key], medium, regime)
                            pdf.savefig(fig)
                            plt.close(fig)

            # --- Spectrum Shapes ---
            if "spectrum_grid" in data:
                print("  Generating spectrum summary...")
                fig = plot_spectrum_summary_grid(data)
                pdf.savefig(fig)
                plt.close(fig)

            if "viz_spectrum_shapes" in data:
                print("  Generating spectrum detail...")
                fig = plot_spectrum_shapes(data["viz_spectrum_shapes"])
                pdf.savefig(fig)
                plt.close(fig)

            # --- Frequencies ---
            if "freq_grid" in data or (has_rvs and any(
                    f"rvs_freq_grid_{r}" in data for r in ["thin", "thick"])):
                print("  Generating frequencies summary...")
                fig = plot_frequencies_summary_grid(data)
                pdf.savefig(fig)
                plt.close(fig)

            if "viz_frequencies" in data:
                for medium in ["ISM", "Wind"]:
                    if medium in data["viz_frequencies"]:
                        print(f"  Generating frequencies ({medium})...")
                        fig = plot_frequency_quantities_combined(data["viz_frequencies"], medium)
                        pdf.savefig(fig)
                        plt.close(fig)

            if has_rvs:
                for regime in ["thin", "thick"]:
                    for medium in ["ISM", "Wind"]:
                        freq_key = f"viz_rvs_frequencies_{regime}"
                        if freq_key in data and medium in data[freq_key]:
                            print(f"  Generating reverse shock frequencies {regime} ({medium})...")
                            fig = plot_rvs_frequency_quantities_combined(data[freq_key], medium, regime)
                            pdf.savefig(fig)
                            plt.close(fig)

        if pending_pdf_merges:
            _merge_pending_pdfs(output_path, pending_pdf_merges, insert_at=pending_pdf_merges[0].get("insert_at", 1))

        print(f"Regression report saved to: {output_path}")


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
