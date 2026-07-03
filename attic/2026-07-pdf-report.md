# PDF validation report generator removed 2026-07-02

The matplotlib/reportlab PDF pipeline (`comprehensive_report.pdf`) was replaced
by the unified self-contained HTML report (`tests/report.py`), which renders the
same validation figures (evolution panels, spectral snapshots, convergence and
timing) natively as SVG. Too large to archive verbatim (~4,500 lines); recover
from git history at the paths below (last present at the commit before this
note was added).

Removed files, original locations:

- `validation/visualization/dashboard.py` — `ComprehensiveDashboard`, page
  composition, TOC/link plumbing
- `validation/visualization/benchmark.py` — benchmark overview, convergence
  summary/detail pages, error-distribution plots
- `validation/visualization/regression.py` — dynamics/frequency/spectrum
  summary grids and diagnostic panels
- `validation/visualization/common.py` — plot style, `ReportBuilder`,
  reportlab markdown rendering; the non-plotting helpers (fiducial constants,
  `load_json_safe`, `get_runtime_build_info`, …) live on in
  `tests/validation/common.py`
- `validation/visualization/__init__.py`
- `validation/comprehensive_report.pdf` — the committed artifact

Kept and relocated:

- `validation/visualization/benchmark_svg.py` → `tests/validation/benchmark_svg.py`
  (README performance SVGs)
- `validation/visualization/guides/` → `tests/validation/guides/`
  (methodology docs)

Also removed with the pipeline: the `[test]` extras `pypdf`, `reportlab`,
`pdfrw`, `Pillow`, `pylatexenc`; the `--no-report`/`--output`/`--version`
flags of `run_validation.py`; the PDF steps in `PyPI-build.yml` and
`deploy-report.yml` (both now publish the HTML report to GitHub Pages
`reports/latest/`).
