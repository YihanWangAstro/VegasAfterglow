"""Common utilities, constants, and helper functions for visualization."""

import json
import math
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

# Try to import VegasAfterglow
try:
    import VegasAfterglow as va
    HAS_VA = True
except ImportError:
    HAS_VA = False

# Page sizes
PAGE_PORTRAIT, PAGE_LANDSCAPE = (8.5, 11), (11, 8.5)

# Color schemes
COLORS = {"pass": "#2ECC71", "fail": "#E74C3C", "primary": "#3498DB", "secondary": "#9B59B6", "neutral": "#7F8C8D", "background": "#2C3E50"}
PHASE_COLORS = {"coasting": "#E74C3C", "crossing": "#F39C12", "BM": "#3498DB", "post_crossing": "#3498DB", "deep_newtonian": "#2ECC71"}
PHASE_NAMES = {"coasting": "Coasting", "crossing": "Crossing", "BM": "Blandford-McKee", "post_crossing": "Post-crossing", "deep_newtonian": "Sedov-Taylor"}
BAND_COLORS = {"Radio": "firebrick", "Optical": "yellowgreen", "X-ray": "royalblue", "TeV": "purple"}
MEDIUM_STYLES, MEDIUM_MARKERS = {"ISM": "-", "wind": "--"}, {"ISM": "o", "wind": "s"}

# Convergence thresholds
MAX_ERROR_THRESHOLD, MEAN_ERROR_THRESHOLD = 0.15, 0.05

# Fiducial resolution settings (single source of truth, used by benchmark and visualization)
FIDUCIAL_RESOLUTION = (0.1, 0.5, 5)
DIM_INDEX = {"phi": 0, "theta": 1, "t": 2}
FIDUCIAL_VALUES = {dim: FIDUCIAL_RESOLUTION[idx] for dim, idx in DIM_INDEX.items()}

QTY_SYMBOLS = {"u": r"$\Gamma\beta$", "Gamma": r"$\Gamma$", "r": r"$r$", "B": r"$B$", "N_p": r"$N_p$",
               "nu_m": r"$\nu_m$", "nu_c": r"$\nu_c$", "nu_a": r"$\nu_a$", "nu_M": r"$\nu_M$"}


def to_float(value):
    """Convert a value to float, handling Fraction and string representations."""
    from fractions import Fraction
    if value is None:
        return None
    if isinstance(value, str):
        try:
            return float(Fraction(value))
        except ValueError:
            return float(value)
    return float(value)


def format_slope(value):
    from fractions import Fraction
    if value is None:
        return ""
    if isinstance(value, (str, Fraction)):
        return str(value)
    return str(Fraction(value).limit_denominator(100))


def get_flux_time(config: Dict) -> float:
    timing = config.get("timing", {})
    return timing.get("flux_single_ms", 0) or timing.get("flux_density_ms", 0)


def find_fiducial_index(values: List, fiducial: float) -> int:
    if not values:
        return -1
    for i, v in enumerate(values):
        if abs(v - fiducial) < 1e-9:
            return i
    diffs = [abs(v - fiducial) for v in values]
    return diffs.index(min(diffs))


def get_unique_short_names(names: List[str], max_len: int = 6) -> Dict[str, str]:
    short = {}
    for name in names:
        for length in range(3, len(name) + 1):
            prefix = name[:length]
            if prefix not in [short.get(n, "")[:length] for n in names if n != name]:
                short[name] = prefix[:max_len]
                break
        else:
            short[name] = name[:max_len]
    return short


def worst_fiducial_error(errors_by_band: Dict, fid_idx: int) -> float:
    """Return the worst (max absolute) error across all bands at the fiducial index."""
    worst = 0.0
    for errs in errors_by_band.values():
        if errs and 0 <= fid_idx < len(errs):
            v = errs[fid_idx]
            if v is not None and not math.isnan(v):
                worst = max(worst, abs(v))
    return worst


def get_git_commit() -> str:
    """Get the short git commit hash of HEAD."""
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def load_json_safe(path: Path) -> Tuple[Optional[Dict], Optional[str]]:
    """Load JSON with error handling. Returns (data, None) on success or (None, error_msg) on failure."""
    if not path.exists():
        return None, f"Results not found: {path}"
    try:
        return json.loads(path.read_text()), None
    except json.JSONDecodeError as e:
        return None, f"Invalid JSON: {e}"


def setup_plot_style():
    plt.rcParams.update({
        "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 12, "legend.fontsize": 9,
        "xtick.labelsize": 9, "ytick.labelsize": 9, "figure.dpi": 72, "savefig.dpi": 72,
        "axes.grid": True, "grid.alpha": 0.3, "lines.linewidth": 1.5,
        "pdf.compression": 9, "mathtext.fontset": "dejavusans",
    })


def _parse_cmake_flags(cmake_path: Path) -> str:
    """Parse compiler flags from CMakeLists.txt based on current platform."""
    import platform
    import re

    if not cmake_path.exists():
        return "unknown"

    try:
        content = cmake_path.read_text()
        system = platform.system()

        if system == "Windows":
            # Look for MSVC flags
            match = re.search(r'elseif\s*\(\s*MSVC\s*\)\s*\n\s*add_compile_options\s*\(([^)]+)\)', content)
            if match:
                flags = match.group(1).strip()
                return flags.replace('\n', ' ').strip()
        else:
            # GNU/Clang flags (macOS/Linux)
            flags = []
            in_gnu_section = False
            for line in content.split('\n'):
                if 'CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang"' in line:
                    in_gnu_section = True
                elif in_gnu_section:
                    if 'elseif' in line or 'endif' in line:
                        break
                    match = re.search(r'add_compile_options\s*\(([^)]+)\)', line)
                    if match:
                        flags.append(match.group(1).strip())

            if flags:
                return ' '.join(flags)
    except Exception:
        pass

    return "unknown"


def get_runtime_build_info() -> Dict[str, str]:
    """Get build/system information at runtime."""
    import platform
    import subprocess
    info = {"Version": "unknown", "Python": sys.version.split()[0], "Platform": f"{platform.system()} {platform.machine()}",
            "Commit": "unknown", "Compiler": "unknown", "Flags": "unknown", "CPU": "unknown"}

    # Get version from VegasAfterglow package (trim commit suffix like +gb28e78a26.d20260204)
    if HAS_VA:
        full_version = getattr(va, "__version__", "unknown")
        # Remove the +commit.date suffix for cleaner display
        info["Version"] = full_version.split("+")[0] if "+" in full_version else full_version

    # Get CPU info
    try:
        cpu_brand = None
        # On macOS, use sysctl for detailed CPU info (platform.processor() only returns "arm" or "i386")
        if platform.system() == "Darwin":
            result = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True)
            if result.returncode == 0 and result.stdout.strip():
                cpu_brand = result.stdout.strip()
        # On Linux, parse /proc/cpuinfo
        elif platform.system() == "Linux":
            try:
                with open("/proc/cpuinfo") as f:
                    for line in f:
                        if line.startswith("model name"):
                            cpu_brand = line.split(":")[1].strip()
                            break
            except (IOError, OSError):
                pass
        # Fallback to platform.processor()
        if not cpu_brand:
            cpu_brand = platform.processor()
        if cpu_brand and cpu_brand not in ("", "arm", "i386", "x86_64", "aarch64"):
            info["CPU"] = cpu_brand
        elif cpu_brand:
            # Use machine type as fallback with arch
            info["CPU"] = f"{platform.machine()} ({cpu_brand})"
    except Exception:
        pass

    # Get git commit
    try:
        result = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True)
        info["Commit"] = result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    # Get compiler info from compile_commands.json (if available)
    compile_commands = Path(__file__).parent.parent.parent / "compile_commands.json"
    if compile_commands.exists():
        try:
            data = json.loads(compile_commands.read_text())
            if data:
                parts = data[0].get("command", "").split()
                if parts:
                    compiler_name = Path(parts[0]).name
                    try:
                        result = subprocess.run([compiler_name, "--version"], capture_output=True, text=True, check=True)
                        info["Compiler"] = result.stdout.split("\n")[0].strip()
                    except (subprocess.CalledProcessError, FileNotFoundError):
                        info["Compiler"] = compiler_name
        except Exception:
            pass

    # Get compiler flags from CMakeLists.txt (more reliable than compile_commands.json)
    cmake_path = Path(__file__).parent.parent.parent / "CMakeLists.txt"
    cmake_flags = _parse_cmake_flags(cmake_path)
    if cmake_flags != "unknown":
        info["Flags"] = cmake_flags

    return info


def extract_session_metadata(session: Dict) -> Dict[str, str]:
    metadata = get_runtime_build_info()
    for key, field in [("vegasafterglow_version", "Version"), ("python_version", "Python"), ("commit", "Commit"),
                       ("platform", "Platform"), ("compiler", "Compiler"), ("compile_flags", "Flags"), ("cpu", "CPU")]:
        if (val := session.get(key)) and val != "unknown":
            # Trim version string if it has commit suffix
            if field == "Version" and "+" in val:
                val = val.split("+")[0]
            metadata[field] = val
    return metadata


def _find_result_file(subdir: str, filename: str) -> Optional[str]:
    """Search common locations for a validation result file."""
    candidates = [
        Path(__file__).parent.parent / subdir / "results" / filename,
        Path(f"tests/{subdir}/results/{filename}"),
        Path(f"{subdir}/results/{filename}"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate.resolve())
    return None


def find_benchmark_file() -> Optional[str]:
    return _find_result_file("benchmark", "benchmark_history.json")


def find_regression_file() -> Optional[str]:
    return _find_result_file("regression", "regression_results.json")


def load_json(filepath: str) -> Dict:
    return json.loads(Path(filepath).read_text())


def _parse_markdown_content(markdown_path: Path) -> List[Dict]:
    """Parse markdown file into structured elements for rendering.

    Returns a list of dicts with 'type' and 'content' keys.
    Types: 'h1', 'h2', 'h3', 'paragraph', 'bullet', 'table', 'hr', 'code'
    """
    import re

    if not markdown_path.exists():
        return []

    content = markdown_path.read_text()

    # Skip YAML frontmatter if present
    if content.startswith('---'):
        end_idx = content.find('---', 3)
        if end_idx != -1:
            content = content[end_idx + 3:].strip()

    elements = []
    lines = content.split('\n')
    i = 0

    while i < len(lines):
        line = lines[i]

        # Horizontal rule
        if re.match(r'^-{3,}$', line.strip()):
            elements.append({'type': 'hr', 'content': ''})
            i += 1
            continue

        # Headers
        _header_match = None
        for prefix, htype in [('### ', 'h3'), ('## ', 'h2'), ('# ', 'h1')]:
            if line.startswith(prefix):
                elements.append({'type': htype, 'content': line[len(prefix):].strip()})
                i += 1
                _header_match = True
                break
        if _header_match:
            continue

        # Table (starts with |)
        if line.strip().startswith('|'):
            table_lines = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                table_lines.append(lines[i])
                i += 1
            # Parse table
            rows = []
            for tl in table_lines:
                # Skip separator row (|---|---|)
                if re.match(r'^\s*\|[-:\s|]+\|\s*$', tl):
                    continue
                cells = [c.strip() for c in tl.strip().strip('|').split('|')]
                if cells and any(c for c in cells):
                    rows.append(cells)
            if rows:
                elements.append({'type': 'table', 'content': rows})
            continue

        # Bullet points
        if line.strip().startswith('- '):
            bullet_text = line.strip()[2:]
            # Handle multi-line bullets
            while i + 1 < len(lines) and lines[i + 1].strip() and not lines[i + 1].strip().startswith('-') and not lines[i + 1].startswith('#'):
                if lines[i + 1].startswith('  ') or lines[i + 1].startswith('\t'):
                    i += 1
                    bullet_text += ' ' + lines[i].strip()
                else:
                    break
            elements.append({'type': 'bullet', 'content': bullet_text})
            i += 1
            continue

        # Empty line
        if not line.strip():
            i += 1
            continue

        # Regular paragraph (collect consecutive non-empty lines)
        para_lines = [line]
        while i + 1 < len(lines):
            next_line = lines[i + 1]
            if not next_line.strip() or next_line.startswith('#') or next_line.strip().startswith('|') or next_line.strip().startswith('-'):
                break
            para_lines.append(next_line)
            i += 1
        paragraph = ' '.join(l.strip() for l in para_lines)
        if paragraph.strip():
            elements.append({'type': 'paragraph', 'content': paragraph})
        i += 1

    return elements


def _render_markdown_with_reportlab(markdown_path: Path, output_path: Path) -> bool:
    """Render markdown to PDF using reportlab. Returns True on success."""
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import inch
        from reportlab.platypus import (
            HRFlowable,
            Paragraph,
            SimpleDocTemplate,
            Spacer,
            Table,
            TableStyle,
        )

        elements_data = _parse_markdown_content(markdown_path)
        if not elements_data:
            return False

        # Create document
        doc = SimpleDocTemplate(
            str(output_path),
            pagesize=letter,
            rightMargin=0.75*inch,
            leftMargin=0.75*inch,
            topMargin=0.75*inch,
            bottomMargin=0.75*inch
        )

        # Define styles
        styles = getSampleStyleSheet()
        _ps = ParagraphStyle
        title_style = _ps('CustomTitle', parent=styles['Heading1'], fontSize=18, spaceAfter=12,
                          textColor=colors.HexColor('#2C3E50'), fontName='Helvetica-Bold')
        h2_style = _ps('CustomH2', parent=styles['Heading2'], fontSize=14, spaceBefore=16, spaceAfter=8,
                       textColor=colors.HexColor('#34495E'), fontName='Helvetica-Bold')
        h3_style = _ps('CustomH3', parent=styles['Heading3'], fontSize=12, spaceBefore=12, spaceAfter=6,
                       textColor=colors.HexColor('#4A5568'), fontName='Helvetica-Bold')
        body_style = _ps('CustomBody', parent=styles['Normal'], fontSize=10, spaceBefore=4, spaceAfter=4,
                         leading=14, fontName='Helvetica')
        bullet_style = _ps('CustomBullet', parent=body_style, leftIndent=20, bulletIndent=10,
                           spaceBefore=2, spaceAfter=2)

        # Build document content
        story = []

        for elem in elements_data:
            elem_type = elem['type']
            content = elem['content']

            # Convert markdown formatting to reportlab tags
            def convert_formatting(text):
                import re

                # LaTeX math: $...$ - convert to Unicode using pylatexenc if available
                def latex_to_unicode(match):
                    latex = match.group(1)
                    try:
                        from pylatexenc.latex2text import LatexNodes2Text
                        # Convert LaTeX to Unicode text
                        unicode_text = LatexNodes2Text().latex_to_text(latex)
                        # Handle subscripts and superscripts for reportlab
                        unicode_text = re.sub(r'_\{([^}]+)\}', r'<sub>\1</sub>', unicode_text)
                        unicode_text = re.sub(r'_([a-zA-Z0-9])', r'<sub>\1</sub>', unicode_text)
                        unicode_text = re.sub(r'\^\{([^}]+)\}', r'<super>\1</super>', unicode_text)
                        unicode_text = re.sub(r'\^([a-zA-Z0-9\-]+)', r'<super>\1</super>', unicode_text)
                        return unicode_text
                    except ImportError:
                        # Fallback: basic substitutions if pylatexenc not available
                        greek = {
                            r'\\nu': 'ν', r'\\alpha': 'α', r'\\beta': 'β', r'\\gamma': 'γ',
                            r'\\Gamma': 'Γ', r'\\theta': 'θ', r'\\phi': 'φ', r'\\tau': 'τ',
                            r'\\sigma': 'σ', r'\\pi': 'π', r'\\rho': 'ρ', r'\\lambda': 'λ',
                            r'\\approx': '≈', r'\\propto': '∝', r'\\sim': '~',
                            r'\\leq': '≤', r'\\geq': '≥', r'\\times': '×', r'\\pm': '±',
                        }
                        for pattern, replacement in greek.items():
                            latex = re.sub(pattern, replacement, latex)
                        latex = re.sub(r'_\{([^}]+)\}', r'<sub>\1</sub>', latex)
                        latex = re.sub(r'_([a-zA-Z0-9])', r'<sub>\1</sub>', latex)
                        latex = re.sub(r'\^\{([^}]+)\}', r'<super>\1</super>', latex)
                        latex = re.sub(r'\^([a-zA-Z0-9\-]+)', r'<super>\1</super>', latex)
                        latex = re.sub(r'\\([a-zA-Z]+)', r'\1', latex)
                        return latex

                text = re.sub(r'\$([^$]+)\$', latex_to_unicode, text)

                # Superscript: ^(...) or ^X for single char - process BEFORE other formatting
                text = re.sub(r'\^[\(\[]([^\)\]]+)[\)\]]', r'<super>\1</super>', text)
                text = re.sub(r'\^(\d+)', r'<super>\1</super>', text)
                text = re.sub(r'\^(-\d+)', r'<super>\1</super>', text)
                # Bold: **text** or __text__
                text = re.sub(r'\*\*([^*]+)\*\*', r'<b>\1</b>', text)
                text = re.sub(r'__([^_]+)__', r'<b>\1</b>', text)
                # Italic: *text* or _text_
                text = re.sub(r'\*([^*]+)\*', r'<i>\1</i>', text)
                text = re.sub(r'(?<![_\w])_([^_]+)_(?![_\w])', r'<i>\1</i>', text)
                # Code: `text`
                text = re.sub(r'`([^`]+)`', r'<font face="Courier" size="9">\1</font>', text)
                return text

            _style_map = {'h1': title_style, 'h2': h2_style, 'h3': h3_style, 'paragraph': body_style}
            if elem_type in _style_map:
                story.append(Paragraph(convert_formatting(content), _style_map[elem_type]))
            elif elem_type == 'bullet':
                story.append(Paragraph(f"• {convert_formatting(content)}", bullet_style))
            elif elem_type == 'hr':
                story.append(Spacer(1, 6))
                story.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor('#BDC3C7')))
                story.append(Spacer(1, 6))
            elif elem_type == 'table' and content:
                table_data = [[Paragraph(convert_formatting(cell), body_style) for cell in row] for row in content]
                t = Table(table_data, repeatRows=1)
                t.setStyle(TableStyle([
                    ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#ECF0F1')),
                    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                    ('FONTSIZE', (0, 0), (-1, -1), 9), ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                    ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                    ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#BDC3C7')),
                    ('TOPPADDING', (0, 0), (-1, -1), 4), ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
                    ('LEFTPADDING', (0, 0), (-1, -1), 6), ('RIGHTPADDING', (0, 0), (-1, -1), 6),
                ]))
                story.extend([Spacer(1, 6), t, Spacer(1, 6)])

        doc.build(story)
        return True

    except ImportError as e:
        print(f"    Warning: reportlab not available ({e})")
        return False
    except Exception as e:
        print(f"    Warning: reportlab rendering failed ({e})")
        return False


def render_markdown_to_pdf_file(markdown_path: Path) -> Optional[Dict]:
    """Render markdown to a standalone PDF file using reportlab.

    Returns dict with 'pdf_path' and 'num_pages' on success, None on failure.
    The caller is responsible for cleaning up the temp file.
    """
    if not markdown_path.exists():
        return None

    import tempfile
    from pypdf import PdfReader

    # Create temp file for the PDF
    temp_dir = tempfile.mkdtemp(prefix="guide_pdf_")
    output_path = Path(temp_dir) / f"{markdown_path.stem}.pdf"

    if _render_markdown_with_reportlab(markdown_path, output_path):
        # Count pages in generated PDF
        try:
            reader = PdfReader(str(output_path))
            num_pages = len(reader.pages)
            return {"pdf_path": str(output_path), "num_pages": num_pages, "temp_dir": temp_dir}
        except Exception:
            return None

    return None


def create_toc_page_reportlab(toc_nodes: List[Dict], output_path: Path) -> int:
    """Create hierarchical table of contents using reportlab.

    Args:
        toc_nodes: List of {"title": str, "page": int, "children": [...]}
        output_path: Path to save standalone TOC PDF

    Returns:
        Number of pages generated (0 if reportlab unavailable).
    """
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import ParagraphStyle
        from reportlab.lib.units import inch
        from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer
    except ImportError:
        return 0

    try:
        styles_by_level = {
            0: ParagraphStyle("TOC0", fontSize=12, fontName="Helvetica-Bold",
                              textColor=colors.HexColor("#2C3E50"), spaceBefore=10, spaceAfter=2, leading=16),
            1: ParagraphStyle("TOC1", fontSize=10, fontName="Helvetica-Bold",
                              textColor=colors.HexColor("#34495E"), leftIndent=20, spaceBefore=6, spaceAfter=1, leading=14),
            2: ParagraphStyle("TOC2", fontSize=9, fontName="Helvetica-Bold",
                              textColor=colors.HexColor("#4A5568"), leftIndent=40, spaceBefore=4, spaceAfter=1, leading=13),
            3: ParagraphStyle("TOC3", fontSize=9, fontName="Helvetica",
                              textColor=colors.HexColor("#666666"), leftIndent=60, spaceBefore=1, spaceAfter=1, leading=12),
        }

        story = []
        title_style = ParagraphStyle("TOCTitle", fontSize=18, fontName="Helvetica-Bold",
                                     textColor=colors.HexColor("#2C3E50"), spaceAfter=20, alignment=1)
        story.append(Paragraph("Table of Contents", title_style))
        story.append(Spacer(1, 12))

        def _render_node(node, level):
            style = styles_by_level.get(min(level, 3))
            title = node["title"]
            story.append(Paragraph(title, style))
            for child in node.get("children", []):
                _render_node(child, level + 1)

        for node in toc_nodes:
            _render_node(node, 0)

        doc = SimpleDocTemplate(str(output_path), pagesize=letter,
                                rightMargin=0.75 * inch, leftMargin=0.75 * inch,
                                topMargin=0.75 * inch, bottomMargin=0.75 * inch)
        doc.build(story)

        from pypdf import PdfReader
        reader = PdfReader(str(output_path))
        return len(reader.pages)
    except Exception as e:
        print(f"  Warning: reportlab TOC failed ({e})")
        return 0


# ---------------------------------------------------------------------------
# ReportBuilder: hybrid PDF assembly (reportlab text + matplotlib figures)
# ---------------------------------------------------------------------------

def _get_reportlab():
    """Import reportlab components. Returns dict of classes."""
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.platypus import (Image, PageBreak, Paragraph,
                                     SimpleDocTemplate, Spacer)
    return {
        "colors": colors, "letter": letter, "inch": inch,
        "ParagraphStyle": ParagraphStyle,
        "Image": Image, "PageBreak": PageBreak,
        "SimpleDocTemplate": SimpleDocTemplate,
        "Paragraph": Paragraph, "Spacer": Spacer,
    }


class ReportBuilder:
    """Assembles a PDF report from matplotlib figures and reportlab flowables.

    Uses a hybrid approach for optimal file size:
    - reportlab for text-heavy pages (title, TOC, section headers, guides)
    - matplotlib PdfPages for figure pages (shared font resources)
    - pypdf for final assembly in correct page order

    Usage:
        builder = ReportBuilder("output.pdf")
        builder.add_title_page("Report Title", metadata={...})
        builder.add_section_header(1, "Section 1")
        builder.add_fig(matplotlib_figure)
        builder.save()
    """

    def __init__(self, output_path: Union[str, Path]):
        self.output_path = Path(output_path)
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        # Each entry: (type, data, page_count)
        # Types: "reportlab", "fig", "pdf_file", "pdf_files", "ext_pdf", "toc_placeholder"
        self._entries = []
        self._toc_insert_idx = None
        self._current_page = 1  # Track current page number (1-based)
        # Internal links: list of (source_page_1based, fig_rect, target_page_1based)
        # fig_rect = (x, y, w, h) in figure coordinates (0-1)
        self._internal_links = []

    @property
    def current_page(self) -> int:
        """Return the current page number (1-based) for the next entry."""
        return self._current_page

    def add_fig(self, fig: plt.Figure):
        """Add a matplotlib figure as a full page."""
        self._entries.append(("fig", fig))
        self._current_page += 1

    def add_fig_file(self, pdf_path: Union[str, Path]):
        """Add a pre-rendered single-page PDF file."""
        self._entries.append(("pdf_file", str(pdf_path)))
        self._current_page += 1

    def add_fig_files(self, pdf_paths: List[Union[str, Path]]):
        """Add multiple pre-rendered PDF files as consecutive pages."""
        self._entries.append(("pdf_files", [str(p) for p in pdf_paths]))
        self._current_page += len(pdf_paths)

    def add_pdf_pages(self, pdf_path: Union[str, Path], num_pages: int = 1):
        """Add all pages from an external PDF file."""
        self._entries.append(("ext_pdf", str(pdf_path)))
        self._current_page += num_pages

    def add_flowables(self, flowables: list):
        """Add reportlab flowables as a page (with automatic page break)."""
        self._entries.append(("reportlab", list(flowables)))
        self._current_page += 1

    def set_toc_position(self, estimated_pages: int = 1):
        """Mark where TOC should be inserted during save()."""
        self._toc_insert_idx = len(self._entries)
        self._entries.append(("toc_placeholder", None))
        self._current_page += estimated_pages

    def add_internal_links(self, source_page: int, links: List[tuple]):
        """Register clickable jump links on a page.

        Args:
            source_page: 1-based page number of the source page.
            links: List of (fig_rect, target_page) where fig_rect is
                   (x, y, w, h) in matplotlib figure coordinates (0-1 range)
                   and target_page is the 1-based destination page number.
        """
        for fig_rect, target_page in links:
            self._internal_links.append((source_page, fig_rect, target_page))

    def add_title_page(self, title: str, subtitle: str = "", metadata: Optional[Dict] = None):
        """Add a title page using reportlab flowables."""
        rl = _get_reportlab()
        inch = rl["inch"]
        Paragraph = rl["Paragraph"]
        Spacer = rl["Spacer"]
        ParagraphStyle = rl["ParagraphStyle"]
        colors = rl["colors"]

        story = []
        story.append(Spacer(1, 2 * inch))

        # Logo or title text
        logo_path = None
        assets_dir = Path(__file__).parent.parent.parent / "assets"
        if (assets_dir / "logo.png").exists():
            logo_path = str(assets_dir / "logo.png")

        if logo_path:
            try:
                # Preserve aspect ratio: read actual size, scale to target height
                from reportlab.lib.utils import ImageReader
                img_reader = ImageReader(logo_path)
                iw, ih = img_reader.getSize()
                target_h = 1.8 * inch
                scale = target_h / ih
                img = rl["Image"](logo_path, width=iw * scale, height=target_h)
                img.hAlign = "CENTER"
                story.append(img)
                story.append(Spacer(1, 0.4 * inch))
            except Exception:
                logo_path = None

        if not logo_path:
            style = ParagraphStyle("TitleLogo", fontSize=36, fontName="Helvetica-Bold",
                                   textColor=colors.HexColor("#2C3E50"), alignment=1)
            story.append(Paragraph("VegasAfterglow", style))
            story.append(Spacer(1, 0.3 * inch))

        title_style = ParagraphStyle("ReportTitle", fontSize=24, fontName="Helvetica",
                                     textColor=colors.HexColor("#34495E"), alignment=1, spaceAfter=16)
        story.append(Paragraph(title, title_style))

        if subtitle:
            sub_style = ParagraphStyle("Subtitle", fontSize=14, fontName="Helvetica",
                                       textColor=colors.HexColor("#7F8C8D"), alignment=1, spaceAfter=20)
            story.append(Paragraph(subtitle, sub_style))

        story.append(Spacer(1, 0.4 * inch))
        ts_style = ParagraphStyle("Timestamp", fontSize=11, fontName="Helvetica",
                                  textColor=colors.HexColor("#95A5A6"), alignment=1, spaceAfter=8)
        story.append(Paragraph(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", ts_style))

        if metadata:
            story.append(Spacer(1, 0.4 * inch))
            version_info = " | ".join(f"{k} {metadata[k]}" for k in ("Version", "Python") if k in metadata)
            if version_info:
                vs = ParagraphStyle("VersionInfo", fontSize=10, fontName="Helvetica-Bold",
                                    textColor=colors.HexColor("#34495E"), alignment=1, spaceAfter=6)
                story.append(Paragraph(version_info.replace("Version", "VegasAfterglow"), vs))
            env_info = " | ".join(f"{k}: {metadata[k]}" for k in ("Commit", "Platform") if k in metadata)
            if env_info:
                es = ParagraphStyle("EnvInfo", fontSize=9, fontName="Helvetica",
                                    textColor=colors.HexColor("#7F8C8D"), alignment=1, spaceAfter=6)
                story.append(Paragraph(env_info, es))
            if "CPU" in metadata and metadata["CPU"] != "unknown":
                cpu_style = ParagraphStyle("CPUInfo", fontSize=8, fontName="Helvetica",
                                           textColor=colors.HexColor("#95A5A6"), alignment=1, spaceAfter=6)
                story.append(Paragraph(f"CPU: {metadata['CPU']}", cpu_style))

        story.append(Spacer(1, 1 * inch))
        ps = ParagraphStyle("Powered", fontSize=9, fontName="Helvetica-Oblique",
                            textColor=colors.HexColor("#95A5A6"), alignment=1)
        story.append(Paragraph("Powered by Claude Code", ps))

        self.add_flowables(story)

    def add_section_header(self, section_num: int, title: str, description: str = ""):
        """Add a section header page using reportlab flowables."""
        rl = _get_reportlab()
        inch = rl["inch"]
        Paragraph = rl["Paragraph"]
        Spacer = rl["Spacer"]
        ParagraphStyle = rl["ParagraphStyle"]
        colors = rl["colors"]

        story = []
        story.append(Spacer(1, 3.5 * inch))

        sec_style = ParagraphStyle("SectionNum", fontSize=14, fontName="Helvetica",
                                   textColor=colors.HexColor("#7F8C8D"), alignment=1, spaceAfter=20)
        story.append(Paragraph(f"Section {section_num}", sec_style))

        title_style = ParagraphStyle("SectionTitle", fontSize=28, fontName="Helvetica-Bold",
                                     textColor=colors.HexColor("#2C3E50"), alignment=1, spaceAfter=24)
        story.append(Paragraph(title, title_style))

        if description:
            desc_style = ParagraphStyle("SectionDesc", fontSize=12, fontName="Helvetica",
                                        textColor=colors.HexColor("#7F8C8D"), alignment=1, leading=18)
            story.append(Paragraph(description, desc_style))

        self.add_flowables(story)

    def add_guide_pages(self, guide_name: str):
        """Add guide pages from a markdown file, rendered via reportlab."""
        guide_path = Path(__file__).parent / "guides" / f"{guide_name}.md"
        result = render_markdown_to_pdf_file(guide_path)
        if result:
            num_pages = result.get("num_pages", 1)
            self.add_pdf_pages(result["pdf_path"], num_pages=num_pages)
            # Don't clean up yet - save() will read the file
        else:
            rl = _get_reportlab()
            label = guide_name.replace("_", " ").title()
            style = rl["ParagraphStyle"]("GuideNA", fontSize=14, fontName="Helvetica",
                                         alignment=1, textColor=rl["colors"].HexColor("#666666"))
            self.add_flowables([
                rl["Spacer"](1, 4 * rl["inch"]),
                rl["Paragraph"](f"{label} - See {guide_name}.md for details", style),
            ])

    def _build_toc_pages(self, toc_nodes: List[Dict], output_path: Path) -> int:
        """Build TOC as a standalone PDF. Returns number of pages."""
        return create_toc_page_reportlab(toc_nodes, output_path)

    def save(self, toc_nodes: Optional[List[Dict]] = None):
        """Build the final PDF by assembling all components.

        Strategy:
        1. Group consecutive matplotlib figures into shared PdfPages (font dedup)
        2. Build reportlab text pages into standalone PDFs
        3. Merge everything in order using pypdf
        """
        from pypdf import PdfReader, PdfWriter

        temp_dir = tempfile.mkdtemp(prefix="report_build_")
        temp_pdfs = []  # List of (pdf_path, page_count) in order

        try:
            # Resolve TOC if needed
            entries = list(self._entries)
            if self._toc_insert_idx is not None and toc_nodes:
                toc_path = Path(temp_dir) / "toc.pdf"
                toc_pages = self._build_toc_pages(toc_nodes, toc_path)
                if toc_pages > 0:
                    entries[self._toc_insert_idx] = ("ext_pdf", str(toc_path))

            # Process entries into temp PDFs
            fig_batch = []  # Collect consecutive figures

            def _flush_fig_batch():
                if not fig_batch:
                    return
                pdf_path = Path(temp_dir) / f"figs_{len(temp_pdfs):04d}.pdf"
                with PdfPages(str(pdf_path)) as pdf:
                    for fig in fig_batch:
                        pdf.savefig(fig)
                        plt.close(fig)
                reader = PdfReader(str(pdf_path))
                temp_pdfs.append((str(pdf_path), len(reader.pages)))
                fig_batch.clear()

            for entry_type, entry_data in entries:
                if entry_type == "fig":
                    fig_batch.append(entry_data)

                elif entry_type == "pdf_file":
                    _flush_fig_batch()
                    temp_pdfs.append((entry_data, 1))

                elif entry_type == "pdf_files":
                    _flush_fig_batch()
                    # Merge batch into single PDF for font dedup
                    if entry_data:
                        merged_path = Path(temp_dir) / f"merged_{len(temp_pdfs):04d}.pdf"
                        writer = PdfWriter()
                        for path in entry_data:
                            reader = PdfReader(str(path))
                            for page in reader.pages:
                                writer.add_page(page)
                        with open(merged_path, "wb") as f:
                            writer.write(f)
                        reader = PdfReader(str(merged_path))
                        temp_pdfs.append((str(merged_path), len(reader.pages)))

                elif entry_type == "ext_pdf":
                    _flush_fig_batch()
                    reader = PdfReader(str(entry_data))
                    temp_pdfs.append((str(entry_data), len(reader.pages)))

                elif entry_type == "reportlab":
                    _flush_fig_batch()
                    rl = _get_reportlab()
                    pdf_path = Path(temp_dir) / f"rl_{len(temp_pdfs):04d}.pdf"
                    doc = rl["SimpleDocTemplate"](
                        str(pdf_path),
                        pagesize=rl["letter"],
                        rightMargin=0.75 * rl["inch"],
                        leftMargin=0.75 * rl["inch"],
                        topMargin=0.75 * rl["inch"],
                        bottomMargin=0.75 * rl["inch"],
                    )
                    doc.build(entry_data + [rl["PageBreak"]()])
                    reader = PdfReader(str(pdf_path))
                    temp_pdfs.append((str(pdf_path), len(reader.pages)))

                elif entry_type == "toc_placeholder":
                    _flush_fig_batch()
                    # No TOC provided, skip

            _flush_fig_batch()

            # Assemble final PDF
            writer = PdfWriter()
            for pdf_path, _ in temp_pdfs:
                reader = PdfReader(pdf_path)
                for page in reader.pages:
                    writer.add_page(page)

            # Deduplicate identical objects (fonts, images) across merged pages
            print("  Compressing and deduplicating PDF objects...")
            writer.compress_identical_objects(remove_identicals=True, remove_orphans=True)

            with open(self.output_path, "wb") as f:
                writer.write(f)

            # Post-process: add page numbers, bookmarks, and internal links
            self._post_process(toc_nodes)

        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

        print(f"Report saved to: {self.output_path}")

    def _post_process(self, toc_nodes: Optional[List[Dict]] = None):
        """Add bookmarks and internal links.

        Uses append() rather than clone_from to avoid shared object references
        which cause annotations to appear on all pages.
        """
        from pypdf import PdfReader, PdfWriter
        from pypdf.generic import (
            ArrayObject,
            DictionaryObject,
            NameObject,
            NumberObject,
        )

        reader = PdfReader(str(self.output_path))
        writer = PdfWriter()
        n_pages = len(reader.pages)

        # Add pages without modification
        for page in reader.pages:
            writer.add_page(page)

        # --- Bookmarks ---
        if toc_nodes:
            def _add_bookmark(node, parent=None):
                page = node.get("page", 1) - 1
                bookmark = writer.add_outline_item(node["title"], max(0, page), parent=parent)
                for child in node.get("children", []):
                    _add_bookmark(child, parent=bookmark)

            for node in toc_nodes:
                _add_bookmark(node)

        # Deduplicate fonts from page number overlays BEFORE adding links
        # (doing it after would invalidate link references)
        writer.compress_identical_objects(remove_identicals=True, remove_orphans=False)

        # --- Internal links ---
        # Build annotations manually using DictionaryObject to ensure
        # each annotation is a unique object not shared between pages

        # First, group links by source page to minimize page modifications
        from collections import defaultdict
        links_by_page = defaultdict(list)
        for source_page_1, fig_rect, target_page_1 in self._internal_links:
            src_idx = source_page_1 - 1
            tgt_idx = target_page_1 - 1
            if 0 <= src_idx < n_pages and 0 <= tgt_idx < n_pages:
                links_by_page[src_idx].append((fig_rect, tgt_idx))

        n_links = 0
        for src_idx, page_links in links_by_page.items():
            src_page = writer.pages[src_idx]

            # CRITICAL: Create a new /Annots array for this page to break
            # shared references with other pages from the same source PDF
            new_annots = ArrayObject()
            src_page[NameObject("/Annots")] = new_annots

            # Get page dimensions from MediaBox
            media_box = src_page.mediabox
            page_w = float(media_box.width)
            page_h = float(media_box.height)

            for fig_rect, tgt_idx in page_links:
                tgt_page = writer.pages[tgt_idx]

                # Convert figure coords (0-1) to PDF points
                fx, fy, fw, fh = fig_rect
                x1 = fx * page_w
                y1 = fy * page_h
                x2 = (fx + fw) * page_w
                y2 = (fy + fh) * page_h

                # Build link annotation manually
                link_annot = DictionaryObject()
                link_annot[NameObject("/Type")] = NameObject("/Annot")
                link_annot[NameObject("/Subtype")] = NameObject("/Link")
                link_annot[NameObject("/Rect")] = ArrayObject([
                    NumberObject(x1), NumberObject(y1),
                    NumberObject(x2), NumberObject(y2)
                ])
                # Invisible border (set to [0,0,1] for debug)
                link_annot[NameObject("/Border")] = ArrayObject([
                    NumberObject(0), NumberObject(0), NumberObject(0)
                ])
                # Destination: fit page in window
                link_annot[NameObject("/Dest")] = ArrayObject([
                    tgt_page.indirect_reference,
                    NameObject("/Fit")
                ])

                # Add as indirect object for proper PDF structure
                annot_ref = writer._add_object(link_annot)
                new_annots.append(annot_ref)
                n_links += 1

        if n_links:
            print(f"  Added {n_links} internal jump links")

        temp_output = str(self.output_path) + ".tmp"
        with open(temp_output, "wb") as f:
            writer.write(f)

        shutil.move(temp_output, str(self.output_path))
