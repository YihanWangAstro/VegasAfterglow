#!/usr/bin/env python3
"""VegasAfterglow Visualization Dashboard - generates PDF reports from test results."""

import argparse
import os
import shutil
from io import BytesIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Default to number of CPUs for parallel execution
DEFAULT_WORKERS = os.cpu_count() or 4

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend before importing pyplot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader, PdfWriter

import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from validation.visualization.common import (create_benchmark_guide_pages, create_regression_guide_pages, create_section_header,
                                        create_title_page, extract_session_metadata, find_benchmark_file,
                                        find_regression_file, get_runtime_build_info, load_json, setup_plot_style)
from validation.visualization.benchmark import generate_convergence_pages, plot_benchmark_overview, plot_convergence_summary
from validation.visualization.regression import (plot_combined_summary_grid, plot_frequency_quantities_combined,
                                            plot_shock_quantities_combined, plot_spectrum_shapes)


def _create_toc_page(sections: List[Tuple[str, int, List[Tuple[str, int]]]]) -> plt.Figure:
    """Create a table of contents page.

    Args:
        sections: List of (section_title, page_num, subsections) where subsections is [(name, page), ...]
    """
    fig = plt.figure(figsize=(8.5, 11))
    fig.patch.set_facecolor("white")

    # Title
    fig.text(0.5, 0.92, "Table of Contents", ha="center", va="top", fontsize=18, fontweight="bold")

    y = 0.85
    line_height = 0.028

    for section_title, section_page, subsections in sections:
        # Section header
        fig.text(0.1, y, section_title, ha="left", va="top", fontsize=12, fontweight="bold")
        fig.text(0.9, y, str(section_page), ha="right", va="top", fontsize=12)
        # Dotted line
        fig.text(0.5, y - 0.005, "." * 80, ha="center", va="top", fontsize=8, color="#CCCCCC")
        y -= line_height * 1.2

        # Subsections
        for sub_name, sub_page in subsections:
            fig.text(0.15, y, sub_name, ha="left", va="top", fontsize=10, color="#444444")
            fig.text(0.9, y, str(sub_page), ha="right", va="top", fontsize=10, color="#444444")
            y -= line_height

        y -= line_height * 0.5  # Extra space between sections

    return fig


def _add_pdf_outline(output_path: Path, sections: List[Tuple[str, int, List[Tuple[str, int]]]]):
    """Add PDF bookmarks/outline to an existing PDF.

    Args:
        output_path: Path to the PDF file
        sections: List of (section_title, page_num, subsections)
    """
    reader = PdfReader(str(output_path))
    writer = PdfWriter()

    # Copy all pages
    for page in reader.pages:
        writer.add_page(page)

    # Add bookmarks (page numbers are 0-indexed for pypdf)
    for section_title, section_page, subsections in sections:
        parent = writer.add_outline_item(section_title, section_page - 1)
        for sub_name, sub_page in subsections:
            writer.add_outline_item(sub_name, sub_page - 1, parent=parent)

    # Write back
    temp_output = str(output_path) + ".tmp"
    with open(temp_output, "wb") as f:
        writer.write(f)
    shutil.move(temp_output, str(output_path))


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

        if benchmark_data and benchmark_data.get("sessions"):
            metadata = extract_session_metadata(benchmark_data["sessions"][-1])
        else:
            metadata = get_runtime_build_info()

        section_num = 0
        pending_pdf_merges = []
        convergence_insert_at = None  # Track where to insert convergence detail pages
        # Track sections for TOC: [(title, page_num, [(subsection_name, page_num), ...])]
        sections: List[Tuple[str, int, List[Tuple[str, int]]]] = []
        # current_page: logical page number for TOC (includes pages that will be merged later)
        # actual_page: physical page count in PDF (for calculating insert positions)
        current_page = 1
        actual_page = 0  # 0-indexed position in actual PDF

        with PdfPages(output_path) as pdf:
            create_title_page(pdf, "Comprehensive Test Report", "Benchmarks \u2022 Regression Tests \u2022 Physics Diagnostics", metadata)
            current_page += 1
            actual_page += 1

            # Placeholder for TOC page (will be page 2)
            toc_page_num = current_page
            fig_toc_placeholder = plt.figure(figsize=(8.5, 11))
            fig_toc_placeholder.patch.set_facecolor("white")
            pdf.savefig(fig_toc_placeholder)
            plt.close(fig_toc_placeholder)
            current_page += 1
            actual_page += 1

            # Section 1: Benchmark Results
            if benchmark_data and benchmark_data.get("sessions"):
                section_num += 1
                section_start = current_page
                subsections = []

                create_section_header(pdf, section_num, "Benchmark Results", "Performance timing and resolution convergence analysis")
                current_page += 1
                actual_page += 1

                # Add guide pages explaining how to read benchmark results
                print("  - Adding benchmark guide pages")
                subsections.append(("How to Read", current_page))
                guide_pages, guide_merge_info = create_benchmark_guide_pages(pdf)
                if guide_merge_info:
                    guide_merge_info["insert_at"] = actual_page  # Use actual position
                    pending_pdf_merges.insert(0, guide_merge_info)
                current_page += guide_pages  # Logical page count includes guide pages

                latest = benchmark_data["sessions"][-1]
                configs = latest.get("configs", latest.get("results", []))

                print(f"\nGenerating benchmark section ({len(configs)} configs)...")
                print("  - Overview plots")

                on_axis_configs = [c for c in configs if c.get("theta_obs_ratio", 0) <= 1]
                off_axis_configs = [c for c in configs if c.get("theta_obs_ratio", 0) > 1]

                if on_axis_configs:
                    print("    - On-axis overview (\u03b8_v/\u03b8_c\u22641)")
                    subsections.append(("On-axis Overview", current_page))
                    fig = plot_benchmark_overview(latest, angle_filter="on-axis")
                    pdf.savefig(fig)
                    plt.close(fig)
                    current_page += 1
                    actual_page += 1

                if off_axis_configs:
                    print("    - Off-axis overview (\u03b8_v/\u03b8_c>1)")
                    subsections.append(("Off-axis Overview", current_page))
                    fig = plot_benchmark_overview(latest, angle_filter="off-axis")
                    pdf.savefig(fig)
                    plt.close(fig)
                    current_page += 1
                    actual_page += 1

                has_convergence = any(c.get("phi_convergence") for c in configs)
                if has_convergence:
                    print("  - Convergence summary")
                    subsections.append(("Convergence Summary", current_page))
                    fig, model_id_map, models_list = plot_convergence_summary(latest)
                    pdf.savefig(fig)
                    plt.close(fig)
                    current_page += 1
                    actual_page += 1
                    # Remember where to insert convergence detail pages
                    convergence_insert_at = actual_page

                    pending_pdf_merges.extend(generate_convergence_pages(pdf, models_list, n_workers))

                sections.append(("1. Benchmark Results", section_start, subsections))

            # Section 2: Regression Tests
            if regression_data:
                section_num += 1
                section_start = current_page
                subsections = []

                create_section_header(pdf, section_num, "Regression Tests", "Power-law scaling verification against theoretical predictions")
                current_page += 1
                actual_page += 1

                # Add guide pages explaining how to read regression results
                print("\nGenerating regression section...")
                print("  - Adding regression guide pages")
                subsections.append(("How to Read", current_page))
                guide_pages, guide_merge_info = create_regression_guide_pages(pdf)
                if guide_merge_info:
                    guide_merge_info["insert_at"] = actual_page  # Use actual position
                    pending_pdf_merges.append(guide_merge_info)
                current_page += guide_pages  # Logical page count includes guide pages

                if "shock_grid" in regression_data or "freq_grid" in regression_data:
                    print("  - Summary grid")
                    subsections.append(("Summary Grid", current_page))
                    fig = plot_combined_summary_grid(regression_data)
                    pdf.savefig(fig)
                    plt.close(fig)
                    current_page += 1

                for medium in ["ISM", "Wind"]:
                    if "viz_shock_dynamics" in regression_data:
                        viz_data = regression_data["viz_shock_dynamics"]
                        if medium in viz_data:
                            print(f"  - Shock dynamics ({medium})")
                            subsections.append((f"Shock Dynamics ({medium})", current_page))
                            fig = plot_shock_quantities_combined(viz_data, medium)
                            pdf.savefig(fig)
                            plt.close(fig)
                            current_page += 1

                    if "viz_frequencies" in regression_data:
                        viz_data = regression_data["viz_frequencies"]
                        if medium in viz_data:
                            print(f"  - Frequencies ({medium})")
                            subsections.append((f"Frequencies ({medium})", current_page))
                            fig = plot_frequency_quantities_combined(viz_data, medium)
                            pdf.savefig(fig)
                            plt.close(fig)
                            current_page += 1

                if "viz_spectrum_shapes" in regression_data:
                    print("  - Spectrum shapes")
                    subsections.append(("Spectrum Shapes", current_page))
                    fig = plot_spectrum_shapes(regression_data["viz_spectrum_shapes"])
                    pdf.savefig(fig)
                    plt.close(fig)
                    current_page += 1

                sections.append((f"{section_num}. Regression Tests", section_start, subsections))

        # Merge convergence detail pages at the correct position (after summary grid)
        if pending_pdf_merges and convergence_insert_at is not None:
            _merge_pending_pdfs(output_path, pending_pdf_merges, insert_at=convergence_insert_at)

        # Replace TOC placeholder page with actual TOC
        if sections:
            print("  Adding table of contents...")
            reader = PdfReader(str(output_path))
            writer = PdfWriter()

            # Page 1: Title
            writer.add_page(reader.pages[0])

            # Page 2: TOC (replace placeholder)
            toc_fig = _create_toc_page(sections)
            toc_buffer = BytesIO()
            toc_fig.savefig(toc_buffer, format="pdf")
            plt.close(toc_fig)
            toc_buffer.seek(0)
            toc_reader = PdfReader(toc_buffer)
            writer.add_page(toc_reader.pages[0])

            # Remaining pages (skip placeholder at index 1)
            for i in range(2, len(reader.pages)):
                writer.add_page(reader.pages[i])

            # Add PDF bookmarks
            for section_title, section_page, subsections in sections:
                parent = writer.add_outline_item(section_title, section_page - 1)
                for sub_name, sub_page in subsections:
                    writer.add_outline_item(sub_name, sub_page - 1, parent=parent)

            temp_output = str(output_path) + ".tmp"
            with open(temp_output, "wb") as f:
                writer.write(f)
            shutil.move(temp_output, str(output_path))

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
            if data.get("sessions"):
                metadata = extract_session_metadata(data["sessions"][-1])
            else:
                metadata = get_runtime_build_info()

            create_title_page(pdf, "Benchmark Report", metadata=metadata)

            # Add guide pages
            print("  Adding benchmark guide pages...")
            guide_pages, guide_merge_info = create_benchmark_guide_pages(pdf)
            if guide_merge_info:
                guide_merge_info["insert_at"] = 1  # After title page
                pending_pdf_merges.insert(0, guide_merge_info)

            if data.get("sessions"):
                latest = data["sessions"][-1]
                configs = latest.get("configs", latest.get("results", []))

                on_axis_configs = [c for c in configs if c.get("theta_obs_ratio", 0) <= 1]
                off_axis_configs = [c for c in configs if c.get("theta_obs_ratio", 0) > 1]

                if on_axis_configs:
                    fig = plot_benchmark_overview(latest, angle_filter="on-axis")
                    pdf.savefig(fig)
                    plt.close(fig)

                if off_axis_configs:
                    fig = plot_benchmark_overview(latest, angle_filter="off-axis")
                    pdf.savefig(fig)
                    plt.close(fig)

                has_convergence = any(c.get("phi_convergence") for c in configs)
                if has_convergence:
                    fig, model_id_map, models_list = plot_convergence_summary(latest)
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

        with PdfPages(output_path) as pdf:
            subtitle = f"{n_pass}/{n_total} tests passed" if n_total > 0 else ""
            create_title_page(pdf, "Regression Test Report", subtitle=subtitle, metadata=metadata)

            # Add guide pages
            print("  Adding regression guide pages...")
            guide_pages, guide_merge_info = create_regression_guide_pages(pdf)
            if guide_merge_info:
                guide_merge_info["insert_at"] = 1  # After title page
                pending_pdf_merges.append(guide_merge_info)

            if "shock_grid" in data or "freq_grid" in data:
                print("  Generating summary grid...")
                fig = plot_combined_summary_grid(data)
                pdf.savefig(fig)
                plt.close(fig)

            if "viz_shock_dynamics" in data:
                viz_data = data["viz_shock_dynamics"]
                for medium in ["ISM", "Wind"]:
                    if medium in viz_data:
                        print(f"  Generating shock dynamics ({medium})...")
                        fig = plot_shock_quantities_combined(viz_data, medium)
                        pdf.savefig(fig)
                        plt.close(fig)

            if "viz_frequencies" in data:
                viz_data = data["viz_frequencies"]
                for medium in ["ISM", "Wind"]:
                    if medium in viz_data:
                        print(f"  Generating frequencies ({medium})...")
                        fig = plot_frequency_quantities_combined(viz_data, medium)
                        pdf.savefig(fig)
                        plt.close(fig)

            if "viz_spectrum_shapes" in data:
                print("  Generating spectrum shapes...")
                fig = plot_spectrum_shapes(data["viz_spectrum_shapes"])
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
