"""VegasAfterglow visualization package for generating PDF reports from test results."""

from .dashboard import ComprehensiveDashboard, main
from .common import setup_plot_style, find_benchmark_file, find_regression_file
from .benchmark import plot_benchmark_overview, plot_convergence_summary, plot_single_model_convergence_page
from .regression import (plot_dynamics_summary_grid, plot_spectrum_summary_grid, plot_frequencies_summary_grid,
                         plot_shock_quantities_combined, plot_frequency_quantities_combined,
                         plot_spectrum_shapes,
                         plot_rvs_shock_quantities_combined, plot_rvs_frequency_quantities_combined)

__all__ = [
    "ComprehensiveDashboard",
    "main",
    "setup_plot_style",
    "find_benchmark_file",
    "find_regression_file",
    "plot_benchmark_overview",
    "plot_convergence_summary",
    "plot_single_model_convergence_page",
    "plot_dynamics_summary_grid",
    "plot_spectrum_summary_grid",
    "plot_frequencies_summary_grid",
    "plot_shock_quantities_combined",
    "plot_frequency_quantities_combined",
    "plot_spectrum_shapes",
    "plot_rvs_shock_quantities_combined",
    "plot_rvs_frequency_quantities_combined",
]
