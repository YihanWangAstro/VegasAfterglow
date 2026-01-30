"""Common utilities, constants, and helper functions for visualization."""

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, str(Path(__file__).parent.parent / "benchmark"))
from benchmark_suite import FIDUCIAL_VALUES

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
PHASE_COLORS = {"coasting": "#E74C3C", "crossing": "#F39C12", "BM": "#3498DB", "deep_newtonian": "#2ECC71"}
PHASE_NAMES = {"coasting": "Coasting", "crossing": "Crossing", "BM": "Blandford-McKee", "deep_newtonian": "Sedov-Taylor"}
BAND_COLORS = {
    "Radio": "firebrick", "Optical": "yellowgreen", "X-ray": "royalblue",
    "Radio (fwd)": "firebrick", "Optical (fwd)": "yellowgreen", "X-ray (fwd)": "royalblue",
    "Radio (rvs)": "salmon", "Optical (rvs)": "darkseagreen", "X-ray (rvs)": "cornflowerblue",
}
MEDIUM_STYLES, MEDIUM_MARKERS = {"ISM": "-", "wind": "--"}, {"ISM": "o", "wind": "s"}

# Convergence thresholds
MAX_ERROR_THRESHOLD, MEAN_ERROR_THRESHOLD = 0.15, 0.05

QTY_SYMBOLS = {"u": r"$\Gamma\beta$", "Gamma": r"$\Gamma$", "r": r"$r$", "B": r"$B$", "N_p": r"$N_p$",
               "nu_m": r"$\nu_m$", "nu_c": r"$\nu_c$", "nu_a": r"$\nu_a$", "nu_M": r"$\nu_M$"}


def to_float(value):
    from fractions import Fraction
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, Fraction):
        return float(value)
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
    if isinstance(value, str):
        return value
    if isinstance(value, Fraction):
        return str(value)
    frac = Fraction(value).limit_denominator(100)
    return str(frac)


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


def setup_plot_style():
    plt.rcParams.update({
        "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 12, "legend.fontsize": 9,
        "xtick.labelsize": 9, "ytick.labelsize": 9, "figure.dpi": 100, "savefig.dpi": 100,
        "axes.grid": True, "grid.alpha": 0.3, "lines.linewidth": 1.5,
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
    info = {"Version": "unknown", "Python": sys.version.split()[0], "Platform": f"{platform.system()} {platform.machine()}",
            "Commit": "unknown", "Compiler": "unknown", "Flags": "unknown"}

    # Get version from VegasAfterglow package
    if HAS_VA:
        info["Version"] = getattr(va, "__version__", "unknown")

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
            with open(compile_commands, "r") as f:
                data = json.load(f)
            if data:
                cmd = data[0].get("command", "")
                parts = cmd.split()
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
                       ("platform", "Platform"), ("compiler", "Compiler"), ("compile_flags", "Flags")]:
        session_value = session.get(key)
        # Only override if session has a valid (non-empty, non-unknown) value
        if session_value and session_value != "unknown":
            metadata[field] = session_value
    return metadata


def find_benchmark_file() -> Optional[str]:
    candidates = [Path(__file__).parent.parent / "benchmark" / "results" / "benchmark_history.json",
                  Path("tests/benchmark/results/benchmark_history.json"), Path("benchmark/results/benchmark_history.json")]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate.resolve())
    return None


def find_regression_file() -> Optional[str]:
    candidates = [Path(__file__).parent.parent / "regression" / "results" / "regression_results.json",
                  Path("tests/regression/results/regression_results.json"), Path("regression/results/regression_results.json")]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate.resolve())
    return None


def load_json(filepath: str) -> Dict:
    with open(filepath) as f:
        return json.load(f)


def create_title_page(pdf: PdfPages, title: str, subtitle: str = "", metadata: Optional[Dict] = None):
    fig = plt.figure(figsize=(8.5, 11))

    # Try to load and display logo (prefer PNG, fall back to SVG conversion)
    logo_displayed = False
    assets_dir = Path(__file__).parent.parent.parent / "assets"
    logo_png = assets_dir / "logo.png"
    logo_svg = assets_dir / "logo.svg"

    img = None
    # Try PNG first (no external dependencies)
    if logo_png.exists():
        try:
            from PIL import Image
            img = Image.open(logo_png)
        except ImportError:
            pass

    # Try SVG with cairosvg if PNG not available
    if img is None and logo_svg.exists():
        try:
            import cairosvg
            from io import BytesIO
            from PIL import Image
            png_data = cairosvg.svg2png(url=str(logo_svg), output_width=400)
            img = Image.open(BytesIO(png_data))
        except (ImportError, OSError):
            pass

    if img is not None:
        try:
            # Add axes for logo (centered, upper portion of page)
            ax = fig.add_axes([0.25, 0.62, 0.5, 0.25])  # [left, bottom, width, height]
            ax.imshow(img)
            ax.axis("off")
            logo_displayed = True
        except Exception:
            pass

    if not logo_displayed:
        fig.text(0.5, 0.70, "VegasAfterglow", fontsize=36, ha="center", fontweight="bold", color="#2C3E50")

    # Adjust vertical positions based on whether logo is displayed
    title_y = 0.55 if logo_displayed else 0.62
    subtitle_y = 0.49 if logo_displayed else 0.56
    timestamp_y = 0.43 if logo_displayed else 0.48
    metadata_y = 0.36 if logo_displayed else 0.40

    fig.text(0.5, title_y, title, fontsize=24, ha="center", color="#34495E")
    if subtitle:
        fig.text(0.5, subtitle_y, subtitle, fontsize=14, ha="center", color="#7F8C8D")
    fig.text(0.5, timestamp_y, f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", fontsize=11, ha="center", color="#95A5A6")
    if metadata:
        y_start = metadata_y
        version_info = []
        if "Version" in metadata:
            version_info.append(f"VegasAfterglow {metadata['Version']}")
        if "Python" in metadata:
            version_info.append(f"Python {metadata['Python']}")
        if version_info:
            fig.text(0.5, y_start, " | ".join(version_info), fontsize=10, ha="center", color="#34495E", fontweight="bold")
            y_start -= 0.04
        env_info = []
        if "Commit" in metadata:
            env_info.append(f"Commit: {metadata['Commit']}")
        if "Platform" in metadata:
            env_info.append(f"Platform: {metadata['Platform']}")
        if env_info:
            fig.text(0.5, y_start, " | ".join(env_info), fontsize=9, ha="center", color="#7F8C8D")
    fig.text(0.5, 0.08, "Powered by Claude Code", fontsize=9, ha="center", color="#95A5A6", style="italic")
    pdf.savefig(fig)
    plt.close(fig)


def create_section_header(pdf: PdfPages, section_num: int, title: str, description: str = ""):
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.5, 0.6, f"Section {section_num}", fontsize=14, ha="center", color="#7F8C8D")
    fig.text(0.5, 0.5, title, fontsize=28, ha="center", fontweight="bold", color="#2C3E50")
    if description:
        fig.text(0.5, 0.4, description, fontsize=12, ha="center", color="#7F8C8D", wrap=True)
    pdf.savefig(fig)
    plt.close(fig)


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
        if line.startswith('# '):
            elements.append({'type': 'h1', 'content': line[2:].strip()})
            i += 1
            continue
        if line.startswith('## '):
            elements.append({'type': 'h2', 'content': line[3:].strip()})
            i += 1
            continue
        if line.startswith('### '):
            elements.append({'type': 'h3', 'content': line[4:].strip()})
            i += 1
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

        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=18,
            spaceAfter=12,
            textColor=colors.HexColor('#2C3E50'),
            fontName='Helvetica-Bold'
        )

        h2_style = ParagraphStyle(
            'CustomH2',
            parent=styles['Heading2'],
            fontSize=14,
            spaceBefore=16,
            spaceAfter=8,
            textColor=colors.HexColor('#34495E'),
            fontName='Helvetica-Bold'
        )

        h3_style = ParagraphStyle(
            'CustomH3',
            parent=styles['Heading3'],
            fontSize=12,
            spaceBefore=12,
            spaceAfter=6,
            textColor=colors.HexColor('#4A5568'),
            fontName='Helvetica-Bold'
        )

        body_style = ParagraphStyle(
            'CustomBody',
            parent=styles['Normal'],
            fontSize=10,
            spaceBefore=4,
            spaceAfter=4,
            leading=14,
            fontName='Helvetica'
        )

        bullet_style = ParagraphStyle(
            'CustomBullet',
            parent=body_style,
            leftIndent=20,
            bulletIndent=10,
            spaceBefore=2,
            spaceAfter=2
        )

        # Build document content
        story = []

        for elem in elements_data:
            elem_type = elem['type']
            content = elem['content']

            # Convert markdown formatting to reportlab tags
            def convert_formatting(text):
                import re
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

            if elem_type == 'h1':
                story.append(Paragraph(convert_formatting(content), title_style))
            elif elem_type == 'h2':
                story.append(Paragraph(convert_formatting(content), h2_style))
            elif elem_type == 'h3':
                story.append(Paragraph(convert_formatting(content), h3_style))
            elif elem_type == 'paragraph':
                story.append(Paragraph(convert_formatting(content), body_style))
            elif elem_type == 'bullet':
                story.append(Paragraph(f"• {convert_formatting(content)}", bullet_style))
            elif elem_type == 'hr':
                story.append(Spacer(1, 6))
                story.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor('#BDC3C7')))
                story.append(Spacer(1, 6))
            elif elem_type == 'table':
                rows = content
                if rows:
                    # Convert cell content
                    table_data = []
                    for row in rows:
                        table_data.append([Paragraph(convert_formatting(cell), body_style) for cell in row])

                    # Create table with styling
                    t = Table(table_data, repeatRows=1)
                    t.setStyle(TableStyle([
                        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#ECF0F1')),
                        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                        ('FONTSIZE', (0, 0), (-1, -1), 9),
                        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                        ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#BDC3C7')),
                        ('TOPPADDING', (0, 0), (-1, -1), 4),
                        ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
                        ('LEFTPADDING', (0, 0), (-1, -1), 6),
                        ('RIGHTPADDING', (0, 0), (-1, -1), 6),
                    ]))
                    story.append(Spacer(1, 6))
                    story.append(t)
                    story.append(Spacer(1, 6))

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
            page = node.get("page", "")
            # Use a table-like approach with dots leader
            dots = "." * max(1, 80 - len(title) - len(str(page)) - level * 4)
            text = f'{title} <font color="#CCCCCC">{dots}</font> {page}'
            story.append(Paragraph(text, style))
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


def _create_guide_pages(pdf: PdfPages, guide_name: str) -> tuple:
    """Create guide pages from a markdown file.

    Args:
        pdf: PdfPages object to write fallback page into.
        guide_name: Base name of the guide file (e.g., "benchmark_guide" or "regression_guide").

    Returns (num_pages, pending_merge_info) where pending_merge_info is either
    None (pages added inline) or a dict with PDF to merge later.
    """
    guide_path = Path(__file__).parent / "guides" / f"{guide_name}.md"

    result = render_markdown_to_pdf_file(guide_path)
    if result:
        return (result["num_pages"], {"files": [result["pdf_path"]], "temp_dir": result["temp_dir"]})

    # Fallback: create a simple placeholder page
    label = guide_name.replace("_", " ").title()
    fig = plt.figure(figsize=PAGE_PORTRAIT)
    fig.patch.set_facecolor("white")
    fig.text(0.5, 0.5, f"{label}\n\nSee {guide_name}.md for details",
             ha="center", va="center", fontsize=14)
    pdf.savefig(fig)
    plt.close(fig)
    return (1, None)


def create_benchmark_guide_pages(pdf: PdfPages) -> tuple:
    """Create guide pages for benchmark results."""
    return _create_guide_pages(pdf, "benchmark_guide")


def create_regression_guide_pages(pdf: PdfPages) -> tuple:
    """Create guide pages for regression test results."""
    return _create_guide_pages(pdf, "regression_guide")
