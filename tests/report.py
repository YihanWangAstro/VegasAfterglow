#!/usr/bin/env python3
"""Unified VegasAfterglow test & validation report.

One self-contained HTML page, organized by what each check establishes:

  code correctness     C++ unit tests (Boost.Test) + Python API tests (pytest)
  physics tests        pytest closure relations / invariants / golden baselines
  physics validation   tests/validation/regression power-law scalings vs analytic theory
  performance          tests/validation/benchmark timing, stage breakdown, convergence

Stdlib only. Chart palette follows the README performance SVGs.

Usage:
    python tests/report.py cpp-junit.xml py-junit.xml \
        --label cpp-junit.xml="C++ (Boost.Test)" --label py-junit.xml="Python (pytest)" \
        --validation tests/validation/regression/results/regression_results.json \
        --benchmark tests/validation/benchmark/results/benchmark_history.json \
        -o test-report.html
"""
import argparse
import ast
import html
import itertools
import json
import math
import re
import platform
import statistics
import subprocess
import sys
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

# pytest modules asserting quantitative physics (everything else = code/API correctness)
PHYSICS_MODULES = ("test_closure_relations", "test_shock_scalings",
                   "test_physics_invariants", "test_golden", "test_fit_smoke")

STATUS = {  # icon, label, css class — icon+label always paired, never color alone
    "passed": ("✓", "pass", "st-pass"),
    "failed": ("✕", "fail", "st-fail"),
    "skipped": ("○", "skip", "st-skip"),
}

# README benchmark-SVG palette — fixed assignment per radiation type, never cycled
RAD_COLORS = {
    "synchrotron": "#5B8ADB",
    "full_ssc": "#45AB8A",
    "ssc_kn": "#4DA5BB",
    "rvs_sync_thin": "#D49B40",
    "rvs_sync_thick": "#D46565",
}
RAD_LABELS = {
    "synchrotron": "synchrotron",
    "full_ssc": "SSC",
    "ssc_kn": "SSC + KN",
    "rvs_sync_thin": "reverse shock (thin)",
    "rvs_sync_thick": "reverse shock (thick)",
}
# Stage order and colors match the README benchmark SVGs (tests/validation/benchmark_svg.py)
STAGES = [
    ("dynamics", "#D49B40", "shock dynamics"),
    ("EAT_grid", "#4DA5BB", "EAT grid"),
    ("syn_electrons", "#45AB8A", "syn. electrons"),
    ("syn_photons", "#D46565", "syn. photons"),
    ("sync_flux", "#5B8ADB", "flux integration"),
    ("cooling", "#9B72CF", "IC cooling"),
    ("ic_photons", "#D4693A", "IC photons"),
    ("ssc_flux", "#52A87C", "SSC flux"),
]

# Benchmark config axes, in README display order
BENCH_JETS = [("tophat", "Tophat"), ("two_component", "Two-comp."),
              ("gaussian", "Gaussian"), ("powerlaw", "Power-law")]
BENCH_MEDIA = ["ISM", "wind"]
BENCH_ANGLES = [0.0, 1.0, 2.0, 4.0]
BAND_COLORS = {"Radio": "#D49B40", "Optical": "#45AB8A", "X-ray": "#5B8ADB",
               "TeV": "#9B72CF"}
CONV_ERR_FLOOR = 1e-4  # log-scale floor: errors at/below 0.01% sit on the axis

# Convergence classification at the fiducial resolution (mirrors run_validation.py)
FIDUCIAL = {"phi": 0.1, "theta": 0.25, "t": 10}
MEAN_ERR_THRESHOLD, MAX_ERR_THRESHOLD = 0.05, 0.15


def esc(s):
    return html.escape(str(s), quote=True)


def fmt_t(seconds):
    if seconds < 1e-3:
        return f"{seconds * 1e6:.0f} µs"
    if seconds < 1:
        return f"{seconds * 1e3:.1f} ms"
    return f"{seconds:.2f} s"


def fmt_ms(ms):
    return f"{ms:.2f} ms" if ms < 100 else f"{ms / 1e3:.2f} s" if ms >= 1e3 else f"{ms:.0f} ms"


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                       text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "?"


def parse_junit(path):
    """Inputs are self-generated (Boost.Test / pytest); still, reject DTDs and
    entity declarations outright so hostile XML can't XXE/billion-laughs us."""
    text = Path(path).read_text(errors="replace")
    if "<!DOCTYPE" in text or "<!ENTITY" in text:
        raise ValueError(f"{path}: DTD/entity declarations are not allowed")
    root = ET.fromstring(text)
    suites = [root] if root.tag == "testsuite" else root.findall("testsuite")
    cases, total_time = [], 0.0
    for suite in suites:
        total_time += float(suite.get("time", 0) or 0)
        for tc in suite.findall("testcase"):
            status, detail = "passed", ""
            for child in tc:
                if child.tag in ("failure", "error"):
                    status = "failed"
                    detail = (child.get("message", "") or "") + "\n" + (child.text or "")
                elif child.tag == "skipped":
                    status = "skipped"
                    detail = child.get("message", "") or (child.text or "")
            cases.append({"classname": tc.get("classname", "") or "(no group)",
                          "name": tc.get("name", "?"),
                          "time": float(tc.get("time", 0) or 0),
                          "status": status, "detail": detail.strip()})
    return {"time": total_time}, cases


def counts(cases):
    return {s: sum(1 for c in cases if c["status"] == s) for s in STATUS}


# --------------------------------------------------------------------- svg helpers

def svg_open(width, height, aria):
    return (f'<svg viewBox="0 0 {width} {height}" role="img" '
            f'aria-label="{esc(aria)}">')


def stacked_status_bar(rows_data, aria):
    width, bar_h, gap, label_w = 780, 26, 14, 220
    out, y = [], 0
    for label, segs in rows_data:  # segs: [(count, css_class, icon, word)]
        total = max(1, sum(c for c, *_ in segs))
        x = label_w
        for count, cls, icon, word in segs:
            if count == 0:
                continue
            w = (width - label_w) * count / total
            text = (f'<text x="{x + w / 2:.1f}" y="{y + bar_h / 2:.1f}" class="seg-label" '
                    f'text-anchor="middle" dominant-baseline="central">{icon} {count}</text>'
                    if w > 46 else "")
            out.append(f'<rect x="{x:.1f}" y="{y}" width="{max(w - 2, 2):.1f}" '
                       f'height="{bar_h}" rx="4" class="{cls}">'
                       f'<title>{esc(label)}: {count} {word}</title></rect>{text}')
            x += w
        out.insert(0, f'<text x="{label_w - 12}" y="{y + bar_h / 2:.1f}" class="axis-label" '
                      f'text-anchor="end" dominant-baseline="central">{esc(label)}</text>')
        y += bar_h + gap
    return svg_open(width, y - gap, aria) + "".join(out) + "</svg>"


def outcome_rows(rows):
    return [(label,
             [(counts(cases)["passed"], "st-pass", "✓", "pass"),
              (counts(cases)["failed"], "st-fail", "✕", "fail"),
              (counts(cases)["skipped"], "st-skip", "○", "skip")])
            for label, cases in rows]


def hbar_chart(rows, aria):
    """Simple horizontal bars: rows = [(label, value, color, tooltip)]."""
    width, bar_h, gap, label_w = 780, 22, 9, 220
    v_max = max(v for _, v, *_ in rows) or 1
    out, y = [], 0
    for label, value, *rest in rows:
        col = rest[0] if rest and rest[0] else "#5B8ADB"
        tip = rest[1] if len(rest) > 1 else f"{label}: {value:.2f} ms"
        w = max((width - label_w - 90) * value / v_max, 2)
        out.append(
            f'<text x="{label_w - 12}" y="{y + bar_h / 2:.1f}" class="axis-label" '
            f'text-anchor="end" dominant-baseline="central">{esc(label)}</text>'
            f'<rect x="{label_w}" y="{y}" width="{w:.1f}" height="{bar_h}" rx="4" '
            f'fill="{col}"><title>{esc(tip)}</title></rect>'
            f'<text x="{label_w + w + 8:.1f}" y="{y + bar_h / 2:.1f}" class="val-label" '
            f'dominant-baseline="central">{value:.2f} ms</text>')
        y += bar_h + gap
    return svg_open(width, y - gap, aria) + "".join(out) + "</svg>"


def slowest_chart(all_cases, top=10):
    slow = [c for c in sorted(all_cases, key=lambda c: -c["time"])[:top] if c["time"] > 0]
    if not slow:
        return ""
    rows = [(c["name"][:44] + ("…" if len(c["name"]) > 44 else ""),
             c["time"] * 1e3, None,
             f'{c["classname"]}::{c["name"]} — {fmt_t(c["time"])}') for c in slow]
    return hbar_chart(rows, "Slowest tests")


def interval_chart(tests):
    """Measured-vs-expected rows: tolerance band, expected tick, measured dot."""
    width, row_h, gap, label_w, val_w = 780, 24, 8, 250, 160
    plot_w = width - label_w - val_w
    out, y = [], 0
    for t in tests:
        exp, meas, tol = t["expected"], t.get("measured"), t["tolerance"]
        ok = t.get("passed", False)
        span = max(tol * 2.6, abs((meas if meas is not None else exp) - exp) * 1.3
                   + tol * 0.4)
        lo = exp - span
        scale = plot_w / (2 * span)

        def X(v):
            return label_w + (v - lo) * scale

        if meas is None:
            cls, icon = "st-skip", "○"
        else:
            cls = "st-pass" if ok else "st-fail"
            icon = "✓" if ok else "✕"
        meas_txt = "—" if meas is None else f"{meas:+.3f}"
        dot = ("" if meas is None else
               f'<circle cx="{X(max(min(meas, exp + span * 0.96), lo + span * 0.04)):.1f}" '
               f'cy="{y + row_h / 2:.1f}" r="5" class="{cls}">'
               f'<title>{esc(t["name"])}: measured {meas_txt}, expected {exp:+.3f} ± {tol}'
               f'</title></circle>')
        out.append(
            f'<text x="{label_w - 12}" y="{y + row_h / 2:.1f}" class="axis-label" '
            f'text-anchor="end" dominant-baseline="central">{esc(t["name"])}</text>'
            f'<rect x="{X(exp - tol):.1f}" y="{y + row_h / 2 - 5:.1f}" '
            f'width="{2 * tol * scale:.1f}" height="10" rx="5" class="tol-band"/>'
            f'<line x1="{X(exp):.1f}" x2="{X(exp):.1f}" y1="{y + 2}" y2="{y + row_h - 2}" '
            f'class="exp-tick"/>{dot}'
            f'<text x="{width - val_w + 12}" y="{y + row_h / 2:.1f}" class="val-label" '
            f'dominant-baseline="central">{icon} {meas_txt} <tspan class="mut">vs '
            f'{exp:+.3f}±{tol:g}</tspan></text>')
        y += row_h + gap
    return svg_open(width, y - gap, "Measured vs expected") + "".join(out) + "</svg>"


# ------------------------------------------------------------------- figures
#
# Native SVG renderings of the validation suite's viz_* payloads: evolution
# panels and spectral snapshots. Shaded bands mark the fitted phase windows;
# dashed guides carry the expected slope.

_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
_FIG_ID = itertools.count()

FIG_W, FIG_H = 386, 258
FIG_ML, FIG_MR, FIG_MT, FIG_MB = 54, 12, 26, 34


def p10(n):
    return "10" + str(n).translate(_SUP)


def sci(v):
    """5000.0 -> 5×10³ (unicode, for SVG text)."""
    if v is None or v <= 0 or not math.isfinite(v):
        return f"{v}"
    e = math.floor(math.log10(v))
    m = f"{v / 10 ** e:.1f}".rstrip("0").rstrip(".")
    return p10(e) if m == "1" else f"{m}×{p10(e)}"


def loglog_panel(title, series, phases=(), guides=(), vlines=(),
                 xlabel="t [s]", ylabel=""):
    """One small-multiple log-log panel.

    series: [(color, xs, ys)] — drawn as 2px lines, README palette
    phases: [(name, t1, t2, tooltip)] — shaded fit windows
    guides: [(x1, x2, slope, tooltip)] — dashed expected-slope segments,
            anchored to the first series inside the window
    vlines: [(x, label, color)] — dashed verticals (spectral breaks)
    """
    clean = []
    all_lx, all_ly = [], []
    for color, xs, ys in series:
        pts = [(math.log10(x), math.log10(y)) for x, y in zip(xs, ys)
               if x > 0 and y > 0 and math.isfinite(x) and math.isfinite(y)]
        if len(pts) >= 2:
            clean.append((color, pts))
            all_lx += [p[0] for p in pts]
            all_ly += [p[1] for p in pts]
    if not clean:
        return ""
    lx0, lx1 = min(all_lx), max(all_lx)
    ly0, ly1 = min(all_ly), max(all_ly)
    lx1, ly1 = max(lx1, lx0 + 1e-9), max(ly1, ly0 + 1e-9)
    pad = 0.05 * (ly1 - ly0)
    ly0, ly1 = ly0 - pad, ly1 + pad
    pw, ph = FIG_W - FIG_ML - FIG_MR, FIG_H - FIG_MT - FIG_MB

    def X(lx):
        return FIG_ML + (lx - lx0) / (lx1 - lx0) * pw

    def Y(ly):
        return FIG_MT + (ly1 - ly) / (ly1 - ly0) * ph

    cid = f"figclip{next(_FIG_ID)}"
    out = [f'<defs><clipPath id="{cid}"><rect x="{FIG_ML}" y="{FIG_MT}" '
           f'width="{pw}" height="{ph}"/></clipPath></defs>',
           f'<rect x="{FIG_ML}" y="{FIG_MT}" width="{pw}" height="{ph}" '
           f'class="fig-frame"/>']

    for name, t1, t2, tip in phases:
        if not (t1 and t2 and t1 > 0 and t2 > 0):
            continue
        a, b = max(math.log10(t1), lx0), min(math.log10(t2), lx1)
        if b <= a:
            continue
        out.append(f'<rect x="{X(a):.1f}" y="{FIG_MT}" width="{X(b) - X(a):.1f}" '
                   f'height="{ph}" class="phase-band"><title>{esc(tip)}</title></rect>')
        if X(b) - X(a) > 40:
            out.append(f'<text x="{(X(a) + X(b)) / 2:.1f}" y="{FIG_MT + 12}" '
                       f'class="phase-label" text-anchor="middle">{esc(name)}</text>')

    xstep = max(1, math.ceil((lx1 - lx0) / 5))
    for n in range(math.ceil(lx0), math.floor(lx1) + 1):
        if n % xstep:
            continue
        out.append(f'<line x1="{X(n):.1f}" x2="{X(n):.1f}" y1="{FIG_MT}" '
                   f'y2="{FIG_MT + ph}" class="grid-line"/>'
                   f'<text x="{X(n):.1f}" y="{FIG_MT + ph + 15}" class="tick-label" '
                   f'text-anchor="middle">{p10(n)}</text>')
    ystep = max(1, math.ceil((ly1 - ly0) / 4))
    for n in range(math.ceil(ly0), math.floor(ly1) + 1):
        if n % ystep:
            continue
        out.append(f'<line x1="{FIG_ML}" x2="{FIG_ML + pw}" y1="{Y(n):.1f}" '
                   f'y2="{Y(n):.1f}" class="grid-line"/>'
                   f'<text x="{FIG_ML - 6}" y="{Y(n):.1f}" class="tick-label" '
                   f'text-anchor="end" dominant-baseline="central">{p10(n)}</text>')

    marks = []
    ref = clean[0][1]

    def y_at(lx):
        if lx <= ref[0][0]:
            return ref[0][1]
        for (xa, ya), (xb, yb) in zip(ref, ref[1:]):
            if xa <= lx <= xb:
                return ya if xb <= xa else ya + (yb - ya) * (lx - xa) / (xb - xa)
        return ref[-1][1]

    for x1, x2, slope, tip in guides:
        if slope is None or not (x1 and x2 and x1 > 0 and x2 > 0):
            continue
        a, b = max(math.log10(x1), lx0), min(math.log10(x2), lx1)
        if b <= a:
            continue
        mid = (a + b) / 2
        off = 0.07 * (ly1 - ly0)
        ya = y_at(mid) + slope * (a - mid) + off
        yb = y_at(mid) + slope * (b - mid) + off
        marks.append(f'<line x1="{X(a):.1f}" y1="{Y(ya):.1f}" x2="{X(b):.1f}" '
                     f'y2="{Y(yb):.1f}" class="guide"><title>{esc(tip)}</title></line>')

    for color, pts in clean:
        stride = max(1, len(pts) // 140)
        pp = pts[::stride]
        if pp[-1] != pts[-1]:
            pp.append(pts[-1])
        path = "M" + " L".join(f"{X(a):.1f} {Y(b):.1f}" for a, b in pp)
        marks.append(f'<path d="{path}" fill="none" stroke="{color}" '
                     f'stroke-width="2" stroke-linejoin="round"/>')
    out.append(f'<g clip-path="url(#{cid})">{"".join(marks)}</g>')

    for vi, (xv, lab, color) in enumerate(vlines):
        if not xv or xv <= 0 or not math.isfinite(xv):
            continue
        lx = math.log10(xv)
        if not (lx0 <= lx <= lx1):
            continue
        out.append(f'<line x1="{X(lx):.1f}" x2="{X(lx):.1f}" y1="{FIG_MT}" '
                   f'y2="{FIG_MT + ph}" stroke="{color}" stroke-width="1.4" '
                   f'stroke-dasharray="4 4"><title>{esc(lab)} = {sci(xv)} Hz</title>'
                   f'</line><text x="{X(lx) + 4:.1f}" y="{FIG_MT + 13 + (vi % 2) * 11}" '
                   f'class="vline-label" fill="{color}">{esc(lab)}</text>')

    out.append(f'<text x="{FIG_ML}" y="16" class="fig-title">{esc(title)}</text>'
               f'<text x="{FIG_ML + pw / 2:.1f}" y="{FIG_H - 4}" class="ax-title" '
               f'text-anchor="middle">{esc(xlabel)}</text>')
    if ylabel:
        out.append(f'<text transform="translate(12,{FIG_MT + ph / 2:.1f}) rotate(-90)" '
                   f'class="ax-title" text-anchor="middle">{esc(ylabel)}</text>')
    return svg_open(FIG_W, FIG_H, title) + "".join(out) + "</svg>"


def phase_windows(mdata, qty=None):
    """Shaded windows (all phases) + expected-slope guides for one quantity."""
    phases, guides = [], []
    for pname, pdat in (mdata.get("phases") or {}).items():
        t1, t2 = pdat.get("t_range", (None, None))
        fits = pdat.get("fits", {})
        if qty is None:
            tip = pname + " — " + ", ".join(
                f"{q}: {f['measured']:+.3f}"
                + (f" (exp {f['expected']:+.3f})" if f.get("expected") is not None else "")
                for q, f in fits.items() if f.get("measured") is not None)
        else:
            f = fits.get(qty) or {}
            tip = f"{pname}: measured {f.get('measured', float('nan')):+.3f}"
            if f.get("expected") is not None:
                tip += f", expected {f['expected']:+.3f}"
                guides.append((t1, t2, f["expected"], tip))
        phases.append((pname.replace("_", " "), t1, t2, tip))
    return phases, guides


DYN_QTYS = [("u", "u = Γβ", ""), ("r", "r", "cm"), ("B", "B′", "G"), ("N_p", "N_p", "")]
FREQ_QTYS = [("nu_m", "ν_m", "#5B8ADB"), ("nu_c", "ν_c", "#45AB8A"),
             ("nu_a", "ν_a", "#D49B40"), ("nu_M", "ν_M", "#D46565")]


def dynamics_figures(viz):
    panels = []
    for medium, md in viz.items():
        for q, qlab, unit in DYN_QTYS:
            if q not in md:
                continue
            phases, guides = phase_windows(md, q)
            panels.append(loglog_panel(
                f"{medium} — {qlab}(t)", [("#5B8ADB", md["t"], md[q])],
                phases, guides, ylabel=f"{qlab} [{unit}]" if unit else qlab))
    return panels, ('<span class="lg"><span class="sw guide-sw"></span>'
                    'expected slope (hover for fit)</span>')


def frequency_figures(viz):
    panels = []
    for medium, md in viz.items():
        series = [(col, md["t"], md[q]) for q, _, col in FREQ_QTYS if q in md]
        phases, _ = phase_windows(md)
        panels.append(loglog_panel(f"{medium} — characteristic frequencies",
                                   series, phases, ylabel="ν [Hz]"))
    legend = "".join(legend_swatch(col, lab) for _, lab, col in FREQ_QTYS)
    return panels, legend


def spectrum_figures(viz, grid):
    panels = []
    for rkey, rd in sorted(viz.items()):
        regime = rd.get("actual_regime", rkey.replace("regime_", ""))
        med = rd.get("medium", "")
        guides = []
        for s in (grid.get(med, {}).get(regime, {}) or {}).get("segments", []):
            if s.get("expected") is not None and s.get("nu_low") and s.get("nu_high"):
                guides.append((s["nu_low"], s["nu_high"], s["expected"],
                               f"{s['segment']}: measured {s['measured']:+.3f}, "
                               f"expected {s.get('expected_expr', s['expected'])}"))
        vlines = [(rd.get("nu_a"), "ν_a", "#45AB8A"), (rd.get("nu_m"), "ν_m", "#4DA5BB"),
                  (rd.get("nu_c"), "ν_c", "#D46565")]
        panels.append(loglog_panel(
            f"Regime {regime} — {med}, t = {sci(rd.get('t'))} s",
            [("#5B8ADB", rd["nu"], rd["flux"])], guides=guides, vlines=vlines,
            xlabel="ν [Hz]", ylabel="F_ν"))
    return panels, ('<span class="lg"><span class="sw guide-sw"></span>'
                    'expected segment slope (hover for fit)</span>')


def fig_cells(panels):
    return "".join(f'<div class="fig-cell">{p}</div>' for p in panels if p)


def fig_grid(panels, legend=""):
    cells = fig_cells(panels)
    if not cells:
        return ""
    lg = f'<div class="legend">{legend}</div>' if legend else ""
    return f'<div class="fig-grid">{cells}</div>{lg}'


def legend_swatch(color, label):
    return (f'<span class="lg"><span class="sw" style="background:{color}"></span>'
            f'{label}</span>')


def stage_legend(present):
    return "".join(legend_swatch(col, esc(label))
                   for key, col, label in STAGES if key in present)


def grouped_timing_chart(cfgs, rad):
    """README-style chart for one radiation type: jet x medium groups, one
    stacked bar per viewing angle (theta_v/theta_c = 0|1|2|4), colored by stage."""
    index = {(c["jet_type"], c["medium"], float(c.get("theta_obs_ratio", 0))): c
             for c in cfgs if c["radiation"] == rad}
    if not index:
        return ""
    groups = [(f"{jlabel}", med, [index.get((jet, med, a)) for a in BENCH_ANGLES])
              for jet, jlabel in BENCH_JETS for med in BENCH_MEDIA]

    W, H = 950, 300
    XL, XR, YB, YT = 56, 936, 236, 30
    n_g, n_b, intra, inter = len(groups), len(BENCH_ANGLES), 3, 16
    bar_w = (XR - XL - (n_g - 1) * inter - n_g * (n_b - 1) * intra) / (n_g * n_b)

    def totals(c):
        sb = (c.get("timing") or {}).get("stage_breakdown") or {}
        return sum(sb.get(k, 0) for k, *_ in STAGES)

    t_max = max((totals(c) for row in groups for c in row[2] if c), default=1) or 1
    # snap the axis to 4 round gridline intervals (same rule as the README SVGs)
    ideal = max(t_max * 1.08 / 4, 1e-6)
    mag = 10 ** math.floor(math.log10(ideal))
    interval = next(m * mag for m in (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 7.5, 8, 10)
                    if m >= ideal / mag)
    y_max = interval * 4

    def gx(gi):
        return XL + gi * (n_b * bar_w + (n_b - 1) * intra + inter)

    def fy(ms):
        return YB - ms / y_max * (YB - YT)

    out, present = [], set()
    # y gridlines at round intervals (values in ms; unit lives on the axis title)
    for n in range(1, 5):
        v = interval * n
        out.append(f'<line x1="{XL}" x2="{XR}" y1="{fy(v):.1f}" y2="{fy(v):.1f}" '
                   f'class="grid-line"/><text x="{XL - 6}" y="{fy(v):.1f}" '
                   f'class="tick-label" text-anchor="end" '
                   f'dominant-baseline="central">{v:g}</text>')
    out.append(f'<line x1="{XL}" x2="{XR}" y1="{YB}" y2="{YB}" class="axis-line"/>'
               f'<text x="{XL - 6}" y="{YB}" class="tick-label" text-anchor="end" '
               f'dominant-baseline="central">0</text>')

    for gi, (jlabel, med, row) in enumerate(groups):
        for bi, c in enumerate(row):
            if c is None:
                continue
            sb = (c.get("timing") or {}).get("stage_breakdown") or {}
            x = gx(gi) + bi * (bar_w + intra)
            name = (f'{jlabel} / {med} / {RAD_LABELS.get(rad, rad)}, '
                    f'θ_v/θ_c = {BENCH_ANGLES[bi]:g}')
            bot = 0.0
            for key, col, slabel in STAGES:
                v = sb.get(key, 0)
                if v <= 0:
                    continue
                present.add(key)
                out.append(f'<rect x="{x:.1f}" y="{fy(bot + v):.1f}" width="{bar_w:.1f}" '
                           f'height="{max(v / y_max * (YB - YT), 0.5):.1f}" fill="{col}">'
                           f'<title>{esc(name)} · {esc(slabel)}: {v:.2f} ms</title></rect>')
                bot += v
            # invisible full-column overlay: hover anywhere on the bar for the
            # complete stage breakdown (sub-pixel segments are unhoverable alone)
            full_tip = esc(name) + " — " + ", ".join(
                f"{slabel} {sb.get(key, 0):.2f}" for key, _, slabel in STAGES
                if sb.get(key, 0) > 0) + f" ms (total {bot:.2f} ms)"
            out.append(f'<rect x="{x:.1f}" y="{fy(bot):.1f}" width="{bar_w:.1f}" '
                       f'height="{max(YB - fy(bot), 1):.1f}" fill="transparent">'
                       f'<title>{full_tip}</title></rect>')
            total_lab = f"{bot:.0f}" if bot >= 100 else f"{bot:.1f}"
            out.append(f'<text x="{x + bar_w / 2:.1f}" y="{fy(bot) - 4:.1f}" '
                       f'class="bar-total" text-anchor="middle">{total_lab}</text>'
                       f'<text x="{x + bar_w / 2:.1f}" y="{YB + 13}" class="tick-label" '
                       f'text-anchor="middle">{BENCH_ANGLES[bi]:g}</text>')
        cx = gx(gi) + (n_b * bar_w + (n_b - 1) * intra) / 2
        out.append(f'<text x="{cx:.1f}" y="{YB + 30}" class="group-label" '
                   f'text-anchor="middle">{esc(jlabel)}</text>'
                   f'<text x="{cx:.1f}" y="{YB + 44}" class="axis-label" '
                   f'text-anchor="middle">{esc(med)}</text>')
    out.append(f'<text transform="translate(14,{(YT + YB) / 2:.0f}) rotate(-90)" '
               f'class="ax-title" text-anchor="middle">wall time [ms]</text>'
               f'<text x="{XL}" y="{YB + 58}" class="ax-title">bar labels: '
               f'θ_v/θ_c ratio · single core · default resolution</text>')
    svg = (svg_open(W, H, f"Timing by configuration — {RAD_LABELS.get(rad, rad)}")
           + "".join(out) + "</svg>")
    return svg + f'<div class="legend">{stage_legend(present)}</div>'


def conv_panel(title, xvalues, band_series, worst_max, fid_x, y_max, xlabel):
    """Log-y convergence panel: per-band mean error lines, worst max-error
    dashed, 5%/15% threshold hairlines, fiducial-resolution marker.
    Errors at/below CONV_ERR_FLOOR (0.01%) are clamped onto the baseline."""
    W, H = 386, 232
    XL, XR, YT, YB = 50, 374, 30, 196
    x0, x1 = min(xvalues), max(xvalues)
    span = (x1 - x0) or 1
    ly0 = math.log10(CONV_ERR_FLOOR)
    ly1 = math.log10(max(y_max, CONV_ERR_FLOOR * 10))

    def X(v):
        return XL + (v - x0) / span * (XR - XL)

    def Y(e):
        le = math.log10(min(max(e, CONV_ERR_FLOOR), y_max))
        return YB - (le - ly0) / (ly1 - ly0) * (YB - YT)

    out = [f'<rect x="{XL}" y="{YT}" width="{XR - XL}" height="{YB - YT}" '
           f'class="fig-frame"/>']
    for n in range(math.ceil(ly0), math.floor(ly1) + 1):
        e = 10 ** n
        lab = f"≤{e * 100:g}%" if e <= CONV_ERR_FLOOR else f"{e * 100:g}%"
        out.append(f'<text x="{XL - 5}" y="{Y(e):.1f}" class="tick-label" '
                   f'text-anchor="end" dominant-baseline="central">{lab}</text>')
        if e > CONV_ERR_FLOOR:
            out.append(f'<line x1="{XL}" x2="{XR}" y1="{Y(e):.1f}" y2="{Y(e):.1f}" '
                       f'class="grid-line"/>')
    for pct, cls, lab in ((MEAN_ERR_THRESHOLD, "thresh-mean", "5% mean"),
                          (MAX_ERR_THRESHOLD, "thresh-max", "15% max")):
        if pct < y_max:
            out.append(f'<line x1="{XL}" x2="{XR}" y1="{Y(pct):.1f}" y2="{Y(pct):.1f}" '
                       f'class="{cls}"><title>{lab} threshold</title></line>'
                       f'<text x="{XR - 2}" y="{Y(pct) - 3:.1f}" class="thresh-label" '
                       f'text-anchor="end">{lab}</text>')
    for v in xvalues:
        out.append(f'<text x="{X(v):.1f}" y="{YB + 14}" class="tick-label" '
                   f'text-anchor="middle">{v:g}</text>')
    if fid_x is not None:
        fx = X(fid_x)
        anchor = ("start" if fx < XL + 30 else
                  "end" if fx > XR - 30 else "middle")
        out.append(f'<line x1="{fx:.1f}" x2="{fx:.1f}" y1="{YT}" '
                   f'y2="{YB}" class="fid-tick"><title>fiducial resolution'
                   f'</title></line><text x="{fx + (3 if anchor == "start" else 0):.1f}" '
                   f'y="{YT - 4}" class="phase-label" text-anchor="{anchor}">'
                   f'fiducial</text>')
    for band, ys in band_series:
        pairs = [(x, e) for x, e in zip(xvalues, ys) if e is not None]
        if len(pairs) < 2:
            continue
        col = BAND_COLORS.get(band, "#5B8ADB")
        path = "M" + " L".join(f"{X(x):.1f} {Y(e):.1f}" for x, e in pairs)
        out.append(f'<path d="{path}" fill="none" stroke="{col}" stroke-width="2" '
                   f'stroke-linejoin="round"/>')
        out.extend(f'<circle cx="{X(x):.1f}" cy="{Y(e):.1f}" r="3.4" fill="{col}">'
                   f'<title>{esc(band)} mean error: {e * 100:.3g}% at '
                   f'{x:g}</title></circle>'
                   for x, e in pairs)
    if worst_max:
        pairs = [(x, e) for x, e in zip(xvalues, worst_max) if e is not None]
        if len(pairs) >= 2:
            path = "M" + " L".join(f"{X(x):.1f} {Y(e):.1f}" for x, e in pairs)
            out.append(f'<path d="{path}" fill="none" class="worst-line">'
                       f'<title>worst max error across bands</title></path>')
    out.append(f'<text x="{XL}" y="16" class="fig-title">{esc(title)}</text>'
               f'<text x="{(XL + XR) / 2:.0f}" y="{H - 6}" class="ax-title" '
               f'text-anchor="middle">{esc(xlabel)}</text>')
    return svg_open(W, H, title) + "".join(out) + "</svg>"


# --------------------------------------------------------------------- tables

def collect_test_docs():
    """Map test name -> one-line description, from Python docstrings and the
    // comment block directly above each BOOST_AUTO_TEST_CASE."""
    root = Path(__file__).resolve().parents[1]
    docs = {}
    for p in (root / "tests").rglob("test_*.py"):
        try:
            tree = ast.parse(p.read_text(encoding="utf-8"))
        except (SyntaxError, OSError):
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name.startswith("test"):
                doc = ast.get_docstring(node)
                if doc:
                    para = doc.strip().split("\n\n")[0]
                    docs[node.name] = " ".join(line.strip()
                                               for line in para.splitlines())
    cpp_pat = re.compile(r"((?://[^\n]*\n)+)\s*BOOST_AUTO_TEST_CASE\((\w+)\)")
    for p in (root / "tests" / "cpp").glob("*.cpp"):
        try:
            src = p.read_text(encoding="utf-8")
        except OSError:
            continue
        for m in cpp_pat.finditer(src):
            lines = [line.lstrip("/ ").strip()
                     for line in m.group(1).strip().splitlines()]
            # drop decorative banner lines (----, ====) around section headers
            lines = [ln for ln in lines if ln.strip("-= ")]
            text = re.sub(r"^\d+\.\s*", "", " ".join(lines))
            if text:
                docs[m.group(2)] = text
    return docs


def suite_table(cases, docs=None):
    t_max = max((c["time"] for c in cases), default=1) or 1
    groups = {}
    for c in cases:
        groups.setdefault(c["classname"].split(".")[-1], []).append(c)
    out = []
    for group, gcases in sorted(groups.items()):
        n = counts(gcases)
        g_time = sum(c["time"] for c in gcases)
        badge = (f'<span class="chip st-fail-text">✕ {n["failed"]} failing</span>'
                 if n["failed"] else '<span class="chip st-pass-text">✓ all pass</span>')
        rows = []
        for c in sorted(gcases, key=lambda c: (c["status"] == "passed", c["name"])):
            icon, word, cls = STATUS[c["status"]]
            detail = (f'<details><summary>details</summary><pre>{esc(c["detail"][:4000])}'
                      f'</pre></details>') if c["detail"] else ""
            doc = (docs or {}).get(c["name"].split("[")[0], "")
            desc = f'<div class="tdesc">{esc(doc)}</div>' if doc else ""
            rows.append(
                f'<tr class="case-row" data-name="{esc(group.lower())} '
                f'{esc(c["name"].lower())} {esc(doc.lower())}">'
                f'<td><span class="pill {cls}"><span aria-hidden="true">{icon}</span> '
                f'{word}</span></td>'
                f'<td class="tname">{esc(c["name"])}{desc}{detail}</td>'
                f'<td class="tdur"><span class="mini-track"><span class="mini-bar" '
                f'style="width:{100 * c["time"] / t_max:.1f}%"></span></span>'
                f'{fmt_t(c["time"])}</td></tr>')
        out.append(
            f'<details class="group"{" open" if n["failed"] else ""}>'
            f'<summary><code>{esc(group)}</code>'
            f'<span class="group-meta">{len(gcases)} tests · {fmt_t(g_time)} · {badge}'
            f'</span></summary>'
            f'<table><thead><tr><th>status</th><th>test</th><th>duration</th></tr>'
            f'</thead><tbody>{"".join(rows)}</tbody></table></details>')
    return "".join(out)


# --------------------------------------------------------------------- sections

def how_to_read(items):
    lis = "".join(f"<li>{i}</li>" for i in items)
    return (f'<details class="howto"><summary>How to read this section</summary>'
            f'<ul>{lis}</ul></details>')


VALIDATION_HOWTO = [
    "Each check fits the power-law exponent of one quantity over a phase window "
    "and compares it with the analytic expectation from synchrotron afterglow "
    "theory: dot = measured, tick = expected, band = tolerance (±0.1 for "
    "dynamics and frequencies, ±0.15 for spectral slopes). The expected "
    "exponents live in <code>tests/validation/regression/run_regression.py</code>.",
    "Phases — forward shock: coasting, Blandford-McKee, deep-Newtonian "
    "(Sedov-Taylor); reverse shock: crossing, post-crossing, deep-Newtonian. "
    "Shaded bands in the evolution figures mark the fitted windows; dashed "
    "guides carry the expected slope (hover for the fit numbers).",
    "Spectral regimes I–V are the orderings of ν_a, ν_m, ν_c (I: ν_a &lt; ν_m "
    "&lt; ν_c, slow cooling … V: ν_c &lt; ν_m &lt; ν_a, absorbed fast cooling). "
    "Breaks appear as dashed verticals; each segment's expected slope is a "
    "dashed guide.",
]

PERFORMANCE_HOWTO = [
    "Timing: each bar is one 30-point broadband light curve on a single core, "
    "stacked by profiled pipeline stage; groups are jet × medium and the bars "
    "within a group are viewing angles θ_v/θ_c = 0, 1, 2, 4.",
    "Convergence: for each configuration and grid dimension (φ, θ in points "
    "per degree; t in points per decade), light curves at four resolutions are "
    "compared against a reference run 20% finer than the finest scanned value; "
    "errors are relative flux deviations per band (Radio, Optical, X-ray, plus "
    "TeV for SSC configurations).",
    "Status at the fiducial resolution (0.1, 0.25, 10): pass = mean error "
    "&lt; 5% and max &lt; 15%; acceptable = only the max exceeds 15%; fail = "
    "mean ≥ 5%. The curve panels show the per-band median mean error across "
    "configurations; the dashed line is the worst max error of any "
    "configuration or band.",
]


def validation_section(vpath):
    d = json.loads(Path(vpath).read_text(encoding="utf-8"))
    tests = d.get("tests", [])
    cats = d.get("categories", {})
    n_pass = sum(1 for t in tests if t.get("passed"))
    n_fail = len(tests) - n_pass
    ts = d.get("timestamp", "?")[:16].replace("T", " ")
    verdict = (f'<span class="chip st-pass-text">✓ {n_pass}/{len(tests)} within '
               f'tolerance</span>' if n_fail == 0 else
               f'<span class="chip st-fail-text">✕ {n_fail} of {len(tests)} outside '
               f'tolerance</span>')
    titles = {
        "shock_dynamics": "Forward shock dynamics (u, r, B, N_p × coasting / BM / deep Newtonian)",
        "frequencies": "Characteristic frequencies (ν_m, ν_c, ν_M)",
        "spectrum_shapes": "Spectral segment slopes (regimes I–V)",
        "rvs_shock_dynamics_thin": "Reverse shock dynamics — thin shell",
        "rvs_shock_dynamics_thick": "Reverse shock dynamics — thick shell",
        "rvs_frequencies_thin": "Reverse shock frequencies — thin shell",
        "rvs_frequencies_thick": "Reverse shock frequencies — thick shell",
    }
    viz_map = {
        "shock_dynamics": ("viz_shock_dynamics", dynamics_figures),
        "frequencies": ("viz_frequencies", frequency_figures),
        "rvs_shock_dynamics_thin": ("viz_rvs_shock_dynamics_thin", dynamics_figures),
        "rvs_shock_dynamics_thick": ("viz_rvs_shock_dynamics_thick", dynamics_figures),
        "rvs_frequencies_thin": ("viz_rvs_frequencies_thin", frequency_figures),
        "rvs_frequencies_thick": ("viz_rvs_frequencies_thick", frequency_figures),
    }
    blocks = []
    for cat, items in cats.items():
        valid = [t for t in items if t.get("measured") is not None]
        skipped = len(items) - len(valid)
        if not items:
            continue
        c_pass = sum(1 for t in valid if t.get("passed"))
        extra = f" · {skipped} not measurable" if skipped else ""
        figs = ""
        if cat == "spectrum_shapes" and "viz_spectrum_shapes" in d:
            figs = fig_grid(*spectrum_figures(d["viz_spectrum_shapes"],
                                              d.get("spectrum_grid", {})))
        elif cat in viz_map and viz_map[cat][0] in d:
            figs = fig_grid(*viz_map[cat][1](d[viz_map[cat][0]]))
        blocks.append(
            f'<details class="group" open><summary>'
            f'<code>{esc(titles.get(cat, cat))}</code>'
            f'<span class="group-meta">{c_pass}/{len(valid)} within tolerance{extra}'
            f'</span></summary><div class="pad">{figs}{interval_chart(items)}'
            f'</div></details>')
    return (f'<section id="validation"><div class="sec-head"><h2>Physics validation'
            f'<span class="sec-sub">power-law scalings vs analytic theory</span></h2>'
            f'{verdict}</div>'
            f'<p class="note">Full suite from <code>tests/validation/run_validation.py'
            f'</code>, run {esc(ts)}.</p>{how_to_read(VALIDATION_HOWTO)}'
            f'{"".join(blocks)}</section>', n_pass, n_fail)


def classify_convergence(conv, dim_name):
    values = conv.get("values") or []
    if not values:
        return None
    fid = FIDUCIAL.get(dim_name, values[-1])
    idx = min(range(len(values)), key=lambda i: abs(values[i] - fid))
    mean_by = conv.get("mean_errors_by_band") or {}
    err_by = conv.get("errors_by_band") or {}
    mean_err = max((v[idx] for v in mean_by.values() if len(v) > idx), default=0)
    max_err = max((v[idx] for v in err_by.values() if len(v) > idx), default=0)
    if mean_err >= MEAN_ERR_THRESHOLD:
        return "fail", mean_err, max_err
    if max_err >= MAX_ERR_THRESHOLD:
        return "acceptable", mean_err, max_err
    return "pass", mean_err, max_err


def worst_max_errors(convs, n_points):
    """Worst max-error across all bands of the given convergence dicts,
    one value per resolution step (None where no band has data)."""
    out = []
    for i in range(n_points):
        pts = [arr[i] for cv in convs
               for arr in (cv.get("errors_by_band") or {}).values()
               if len(arr) > i and arr[i] is not None]
        out.append(max(pts) if pts else None)
    return out


def benchmark_section(bpath):
    d = json.loads(Path(bpath).read_text(encoding="utf-8"))
    cfgs = d.get("configs", [])
    ts = d.get("timestamp", "?")[:16].replace("T", " ")
    cpu = d.get("cpu", "?")

    # --- timing by radiation type, split on-axis vs off-axis ---
    def is_on_axis(c):
        return float(c.get("theta_obs_ratio", c.get("theta_obs", 0))) == 0.0

    by_rad, by_rad_on, by_rad_off = {}, {}, {}
    for c in cfgs:
        t = (c.get("timing") or {}).get("total_ms")
        if t:
            by_rad.setdefault(c["radiation"], []).append(t)
            (by_rad_on if is_on_axis(c) else by_rad_off).setdefault(
                c["radiation"], []).append(t)
    all_t = [t for v in by_rad.values() for t in v]
    if not all_t:
        note = ('<section id="performance"><div class="sec-head"><h2>Performance'
                '<span class="sec-sub">benchmark timing &amp; resolution convergence'
                '</span></h2></div><p class="note">Benchmark results file contains no '
                'timed configurations — run <code>python tests/validation/'
                'run_validation.py --benchmark</code> to refresh it.</p></section>')
        return note, 0

    def median_rows(rows_by_rad, tag):
        return [(RAD_LABELS.get(r, r), statistics.median(v), RAD_COLORS.get(r),
                 f"{RAD_LABELS.get(r, r)} ({tag}): median {statistics.median(v):.2f} "
                 f"ms, range {min(v):.2f}–{max(v):.2f} ms over {len(v)} configs")
                for r, v in rows_by_rad.items()]

    on_t = [t for v in by_rad_on.values() for t in v]
    off_t = [t for v in by_rad_off.values() for t in v]
    med_on = statistics.median(on_t) if on_t else 0
    med_off = statistics.median(off_t) if off_t else 0
    median_chart = ""
    if by_rad_on:
        median_chart += ('<h4 class="sub-h">on-axis (θ_v = 0)</h4>'
                         + hbar_chart(median_rows(by_rad_on, "on-axis"),
                                      "Median on-axis light-curve time"))
    if by_rad_off:
        median_chart += ('<h4 class="sub-h">off-axis (θ_v/θ_c = 1, 2, 4)</h4>'
                         + hbar_chart(median_rows(by_rad_off, "off-axis"),
                                      "Median off-axis light-curve time"))

    # --- full per-configuration timing, one README-style chart per radiation ---
    rad_order = [r for r in RAD_COLORS if r in by_rad] + \
                [r for r in by_rad if r not in RAD_COLORS]
    timing_details = "".join(
        f'<details class="group" open><summary>'
        f'<code>{esc(RAD_LABELS.get(r, r))}</code>'
        f'<span class="group-meta">{len([c for c in cfgs if c["radiation"] == r])} '
        f'configurations · jet × medium × viewing angle, stacked by stage</span>'
        f'</summary><div class="pad">{grouped_timing_chart(cfgs, r)}</div></details>'
        for r in rad_order)

    # --- resolution convergence: classify every config x dimension at fiducial ---
    DIM_META = {"phi": ("φ", "φ resolution [pts/deg]"),
                "theta": ("θ", "θ resolution [pts/deg]"),
                "t": ("t", "t resolution [pts/decade]")}
    conv_counts = {dim: {"pass": 0, "acceptable": 0, "fail": 0} for dim in FIDUCIAL}
    status_by = {}  # (id(config), dim) -> (status, mean_err, max_err)
    ranked = []
    for c in cfgs:
        for dim in FIDUCIAL:
            conv = c.get(f"{dim}_convergence")
            if not conv:
                continue
            r = classify_convergence(conv, dim)
            if r is None:
                continue
            conv_counts[dim][r[0]] += 1
            status_by[(id(c), dim)] = r
            ranked.append((r[1], r[2], r[0], dim, c))
    conv_rows = [(f"{dim} resolution",
                  [(n["pass"], "st-pass", "✓", "pass"),
                   (n["acceptable"], "st-warn", "△", "acceptable (max err ≥ 15%)"),
                   (n["fail"], "st-fail", "✕", "fail (mean err ≥ 5%)")])
                 for dim, n in conv_counts.items()]
    n_conv = sum(sum(n.values()) for n in conv_counts.values())
    n_conv_fail = sum(n["fail"] for n in conv_counts.values())
    n_conv_acc = sum(n["acceptable"] for n in conv_counts.values())

    # --- full convergence grid: one table per radiation type ---
    STATUS_CELL = {"pass": ("st-pass-text", "✓"), "acceptable": ("st-warn-text", "△"),
                   "fail": ("st-fail-text", "✕")}
    cfg_index = {(c["jet_type"], c["medium"], float(c.get("theta_obs_ratio", 0)),
                  c["radiation"]): c for c in cfgs}
    grid_tables = []
    for r in rad_order:
        head = ("<tr><th rowspan='2'>configuration</th>"
                + "".join(f"<th colspan='3'>θ_v/θ_c = {a:g}</th>" for a in BENCH_ANGLES)
                + "</tr><tr>"
                + "".join(f"<th>{DIM_META[d][0]}</th>" for _ in BENCH_ANGLES
                          for d in FIDUCIAL) + "</tr>")
        rows = []
        for jet, jlabel in BENCH_JETS:
            for med in BENCH_MEDIA:
                cells = []
                for a in BENCH_ANGLES:
                    c = cfg_index.get((jet, med, a, r))
                    for dim in FIDUCIAL:
                        s = c and status_by.get((id(c), dim))
                        if not s:
                            cells.append('<td class="conv-cell">·</td>')
                            continue
                        cls, icon = STATUS_CELL[s[0]]
                        cells.append(
                            f'<td class="conv-cell {cls}" title="{esc(jlabel)} / '
                            f'{esc(med)} / θ_v/θ_c = {a:g} — {DIM_META[dim][0]}: '
                            f'{s[0]}, mean {s[1]:.2%}, max {s[2]:.2%}">{icon}</td>')
                rows.append(f'<tr><td class="tname">{esc(jlabel)} / {esc(med)}</td>'
                            + "".join(cells) + "</tr>")
        n_rad = sum(1 for (*_, rr) in cfg_index if rr == r) * len(FIDUCIAL)
        grid_tables.append(
            f'<details class="group" open><summary><code>{esc(RAD_LABELS.get(r, r))}'
            f'</code><span class="group-meta">{n_rad} dimension checks · hover a '
            f'cell for errors</span></summary><table class="conv-table">'
            f'<thead>{head}</thead><tbody>{"".join(rows)}</tbody></table></details>')
    grid_html = (f'<div class="card"><h3>Convergence status — every configuration '
                 f'(✓ pass · △ acceptable · ✕ fail, at fiducial)</h3>'
                 f'{"".join(grid_tables)}</div>')

    # --- convergence curves: per radiation x dimension, aggregated over configs ---
    curve_panels = []
    for r in rad_order:
        rc = [c for c in cfgs if c["radiation"] == r]
        for dim in FIDUCIAL:
            convs = [c.get(f"{dim}_convergence") for c in rc]
            convs = [cv for cv in convs if cv and cv.get("values")]
            if not convs:
                continue
            xvals = convs[0]["values"]
            bands = [b for b in BAND_COLORS if any(
                b in (cv.get("mean_errors_by_band") or {}) for cv in convs)]
            band_series, y_top = [], 0.0
            for b in bands:
                med = []
                for i in range(len(xvals)):
                    pts = [cv["mean_errors_by_band"][b][i] for cv in convs
                           if b in (cv.get("mean_errors_by_band") or {})
                           and len(cv["mean_errors_by_band"][b]) > i
                           and cv["mean_errors_by_band"][b][i] is not None]
                    med.append(statistics.median(pts) if pts else None)
                band_series.append((b, med))
                y_top = max(y_top, max((v for v in med if v is not None), default=0))
            worst_max = worst_max_errors(convs, len(xvals))
            y_top = max(y_top, max((v for v in worst_max if v is not None),
                                   default=0), 0.16) * 1.15
            curve_panels.append(conv_panel(
                f"{RAD_LABELS.get(r, r)} — {DIM_META[dim][0]}", xvals, band_series,
                worst_max, FIDUCIAL.get(dim), y_top, DIM_META[dim][1]))
    band_lg = "".join(legend_swatch(col, f"{b} (median mean error)")
                      for b, col in BAND_COLORS.items())
    curve_lg = (band_lg + '<span class="lg"><span class="sw worst-sw"></span>'
                'worst max error (any config/band)</span>')
    curves_html = (f'<div class="card"><h3>Error vs resolution '
                   f'(median over configurations)</h3>'
                   f'<div class="legend legend-top">{curve_lg}</div>'
                   f'<div class="fig-grid">{fig_cells(curve_panels)}</div></div>')

    # --- largest fiducial errors, per-config detail curves (non-PASS first) ---
    ranked.sort(key=lambda w: (w[2] == "pass", -w[0]))
    detail_panels = []
    for _me, _mx, _st, dim, c in ranked[:6]:
        conv = c[f"{dim}_convergence"]
        xvals = conv["values"]
        band_series = [(b, (conv.get("mean_errors_by_band") or {}).get(b, []))
                       for b in BAND_COLORS
                       if b in (conv.get("mean_errors_by_band") or {})]
        worst_max = worst_max_errors([conv], len(xvals))
        y_top = max([v for _, ys in band_series for v in ys if v is not None]
                    + [v for v in worst_max if v is not None] + [0.16]) * 1.15
        title = (f'{c["jet_type"]}/{c["medium"]}/'
                 f'{RAD_LABELS.get(c["radiation"], c["radiation"])} '
                 f'θ_v/θ_c={c.get("theta_obs_ratio", 0):g} — {DIM_META[dim][0]}')
        detail_panels.append(conv_panel(title, xvals, band_series, worst_max,
                                        FIDUCIAL.get(dim), y_top, DIM_META[dim][1]))
    n_bad = sum(1 for w in ranked if w[2] != "pass")
    detail_note = (f"{n_bad} checks outside PASS — largest fiducial errors shown"
                   if n_bad else "all checks PASS — largest fiducial errors shown")
    detail_lg = "".join(legend_swatch(col, f"{b} (mean error)")
                        for b, col in BAND_COLORS.items()) + \
        '<span class="lg"><span class="sw worst-sw"></span>max error (worst band)</span>'
    detail_html = (f'<details class="group"><summary><code>largest fiducial errors '
                   f'— per-config detail</code><span class="group-meta">'
                   f'{detail_note}</span></summary><div class="pad">'
                   f'<div class="legend legend-top">{detail_lg}</div>'
                   f'<div class="fig-grid">{fig_cells(detail_panels)}</div>'
                   '</div></details>'
                   if detail_panels else "")

    verdict = (f'<span class="chip st-pass-text">✓ {n_conv}/{n_conv} checks '
               f'converge</span>' if n_conv_fail == 0 and n_conv_acc == 0 else
               f'<span class="chip {"st-fail-text" if n_conv_fail else "st-warn-text"}">'
               f'{"✕" if n_conv_fail else "△"} {n_conv_fail} fail · {n_conv_acc} '
               f'acceptable of {n_conv}</span>')

    tiles = f"""
<div class="tiles">
<div class="tile"><div class="k">median light curve · on-axis</div>
  <div class="v">{fmt_ms(med_on)}</div><div class="s">30-point broadband, single core</div></div>
<div class="tile"><div class="k">median light curve · off-axis</div>
  <div class="v">{fmt_ms(med_off)}</div><div class="s">θ_v/θ_c = 1, 2, 4</div></div>
<div class="tile"><div class="k">configurations</div><div class="v">{len(cfgs)}</div>
  <div class="s">jets × media × radiation × viewing</div></div>
<div class="tile"><div class="k">convergence checks</div><div class="v">{n_conv}</div>
  <div class="s">φ / θ / t resolution × configs</div></div>
<div class="tile"><div class="k">fastest / slowest</div>
  <div class="v">{min(all_t):.3g}–{max(all_t):.3g}<span class="u"> ms</span></div>
  <div class="s">across all configs</div></div>
<div class="tile"><div class="k">hardware</div><div class="v" style="font-size:17px">
  {esc(cpu)}</div><div class="s">run {esc(ts)}</div></div>
</div>"""

    return (f'<section id="performance"><div class="sec-head"><h2>Performance'
            f'<span class="sec-sub">benchmark timing &amp; resolution convergence'
            f'</span></h2>{verdict}</div>{tiles}{how_to_read(PERFORMANCE_HOWTO)}'
            f'<div class="card"><h3>Light-curve time by radiation type (median)</h3>'
            f'{median_chart}</div>'
            f'<div class="card"><h3>Timing by configuration</h3>{timing_details}</div>'
            f'<div class="card"><h3>Resolution convergence at fiducial '
            f'(mean&nbsp;&lt;&nbsp;5%, max&nbsp;&lt;&nbsp;15%)</h3>'
            f'{stacked_status_bar(conv_rows, "Convergence status by grid dimension")}'
            f'</div>{grid_html}{curves_html}{detail_html}</section>', n_conv_fail)


# --------------------------------------------------------------------- page

CSS = """
:root { color-scheme: light dark;
  --surface: #fcfcfb; --plane: #f9f9f7; --ink: #0b0b0b; --ink-2: #52514e;
  --muted: #898781; --grid: #e1e0d9; --ring: rgba(11,11,11,0.10);
  --blue: #5B8ADB; --good: #2e9e4f; --good-text: #1d7a38; --bad: #D46565;
  --warn: #D49B40; --band: #f0efec; }
@media (prefers-color-scheme: dark) { :root {
  --surface: #191a1c; --plane: #0e0f11; --ink: #f5f5f4; --ink-2: #c3c2b7;
  --muted: #8b8a84; --grid: #2b2c2e; --ring: rgba(255,255,255,0.10);
  --blue: #6d99e3; --good: #3fae5f; --good-text: #4cc06e; --bad: #dd7c7c;
  --warn: #ddad5c; --band: #232426; } }
* { box-sizing: border-box; }
body { margin: 0; font: 15px/1.55 -apple-system, "Segoe UI", Roboto,
  "Helvetica Neue", Arial, sans-serif; background: var(--plane); color: var(--ink); }
.wrap { max-width: 1040px; margin: 0 auto; padding: 0 28px 80px; }
.masthead { padding: 34px 0 22px; display: flex; align-items: center; gap: 18px;
  flex-wrap: wrap; }
.masthead svg.logo { width: 52px; height: 56px; flex-shrink: 0; }
.masthead h1 { font-size: 25px; letter-spacing: -0.015em; margin: 0; }
.masthead .meta { color: var(--ink-2); font-size: 13px; margin-top: 3px; }
.verdict { margin-left: auto; font-size: 14px; font-weight: 700; padding: 9px 20px;
  border-radius: 999px; letter-spacing: 0.03em; }
.verdict.ok { background: color-mix(in srgb, var(--good) 13%, transparent);
  color: var(--good-text);
  border: 1px solid color-mix(in srgb, var(--good) 38%, transparent); }
.verdict.bad { background: color-mix(in srgb, var(--bad) 13%, transparent);
  color: var(--bad); border: 1px solid color-mix(in srgb, var(--bad) 42%, transparent); }
nav { position: sticky; top: 0; z-index: 5;
  background: color-mix(in srgb, var(--plane) 86%, transparent);
  backdrop-filter: blur(10px); border-bottom: 1px solid var(--grid);
  margin: 0 -28px; padding: 10px 28px; display: flex; gap: 8px; flex-wrap: wrap; }
nav a { text-decoration: none; color: var(--ink-2); font-size: 13px; font-weight: 600;
  padding: 6px 14px; border-radius: 999px; border: 1px solid transparent; }
nav a:hover { background: var(--surface); border-color: var(--ring); color: var(--ink); }
nav a .n { color: var(--muted); font-weight: 500; }
nav a.archive { color: var(--muted); }
nav .plat-label { margin-left: auto; color: var(--muted); font-size: 13px;
  padding: 7px 4px 7px 14px; }
nav a.plat { color: var(--ink-2); }
nav .plat-cur { color: var(--ink); font-size: 13px; font-weight: 700;
  padding: 7px 10px; }
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
  gap: 12px; margin: 22px 0 8px; }
.tile { background: var(--surface); border: 1px solid var(--ring); border-radius: 12px;
  padding: 16px 18px; box-shadow: 0 1px 2px rgba(0,0,0,0.03); }
.tile .k { font-size: 11px; color: var(--muted); text-transform: uppercase;
  letter-spacing: 0.07em; font-weight: 650; }
.tile .v { font-size: 29px; font-weight: 750; font-variant-numeric: tabular-nums;
  letter-spacing: -0.01em; margin-top: 3px; }
.tile .v .u { font-size: 15px; font-weight: 550; color: var(--ink-2); }
.tile .s { font-size: 12px; color: var(--ink-2); margin-top: 2px; }
.tile.hero-good .v { color: var(--good-text); }
.tile.hero-bad .v { color: var(--bad); }
.card { background: var(--surface); border: 1px solid var(--ring); border-radius: 12px;
  padding: 20px 22px; margin: 14px 0; overflow-x: auto;
  box-shadow: 0 1px 2px rgba(0,0,0,0.03); }
.card h3 { margin: 0 0 14px; font-size: 12.5px; color: var(--ink-2); font-weight: 650;
  text-transform: uppercase; letter-spacing: 0.06em; }
svg { display: block; width: 100%; height: auto; }
.st-pass { fill: var(--good); } .st-fail { fill: var(--bad); }
.st-skip { fill: var(--muted); } .st-warn { fill: var(--warn); }
.seg-label { fill: #fff; font-size: 12px; font-weight: 650; }
.axis-label { fill: var(--ink-2); font-size: 12px; }
.val-label { fill: var(--ink-2); font-size: 11.5px; font-variant-numeric: tabular-nums; }
.val-label .mut { fill: var(--muted); }
.tol-band { fill: var(--band); stroke: var(--grid); }
.exp-tick { stroke: var(--muted); stroke-width: 1.6; }
.fig-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr));
  gap: 14px; margin: 12px 0 6px; }
.fig-cell { border: 1px solid var(--grid); border-radius: 10px;
  padding: 8px 6px 2px; background: var(--surface); }
.fig-frame { fill: none; stroke: var(--grid); }
.fig-title { fill: var(--ink); font-size: 12.5px; font-weight: 650; }
.ax-title { fill: var(--muted); font-size: 11px; }
.tick-label { fill: var(--muted); font-size: 10.5px;
  font-variant-numeric: tabular-nums; }
.grid-line { stroke: var(--grid); stroke-width: 1; }
.phase-band { fill: var(--blue); opacity: 0.08; }
.phase-label { fill: var(--ink-2); font-size: 10px; }
.guide { stroke: var(--warn); stroke-width: 1.6; stroke-dasharray: 5 4; fill: none; }
.guide-sw { background: repeating-linear-gradient(90deg, var(--warn) 0 5px,
  transparent 5px 9px); height: 3px; margin-top: 4px; border-radius: 0; }
.vline-label { font-size: 10.5px; font-weight: 650; }
.axis-line { stroke: var(--muted); stroke-width: 1; }
.bar-total { fill: var(--ink-2); font-size: 8.5px; font-weight: 650;
  font-variant-numeric: tabular-nums; }
.group-label { fill: var(--ink); font-size: 11px; font-weight: 700; }
.thresh-mean { stroke: var(--warn); stroke-width: 1.2; stroke-dasharray: 2 3; }
.thresh-max { stroke: var(--bad); stroke-width: 1.2; stroke-dasharray: 2 3; }
.thresh-label { fill: var(--muted); font-size: 9.5px; }
.fid-tick { stroke: var(--muted); stroke-width: 1.4; stroke-dasharray: 5 4; }
.worst-line { stroke: var(--ink-2); stroke-width: 1.6; stroke-dasharray: 5 4; }
.worst-sw { background: repeating-linear-gradient(90deg, var(--ink-2) 0 5px,
  transparent 5px 9px); height: 3px; margin-top: 4px; border-radius: 0; }
.conv-table { border-collapse: collapse; margin: 8px 0 14px; width: 100%; }
.conv-table th { font-size: 11px; color: var(--muted); font-weight: 600;
  padding: 3px 6px; text-align: center; border-bottom: 1px solid var(--grid);
  white-space: nowrap; }
.conv-table th:first-child { text-align: left; }
.conv-table td.tname { font-size: 12.5px; padding: 3px 10px 3px 0;
  white-space: nowrap; width: 30%; }
.conv-cell { text-align: center; font-size: 12px; font-weight: 700;
  padding: 3px 6px; cursor: default; width: 5.8%; }
.conv-cell:hover { background: var(--band); border-radius: 4px; }
.sub-h { font-size: 12.5px; color: var(--ink-2); font-weight: 650;
  margin: 14px 0 6px; text-transform: uppercase; letter-spacing: 0.04em; }
.tdesc { color: var(--muted); font-size: 12px; line-height: 1.4; margin-top: 1px; }
.howto { font-size: 13px; color: var(--ink-2); margin: -2px 0 16px; }
.howto summary { cursor: pointer; color: var(--muted); font-weight: 600; }
.howto ul { margin: 8px 0 0; padding-left: 22px; }
.howto li { margin: 5px 0; max-width: 76ch; }
.legend { display: flex; gap: 14px; flex-wrap: wrap; margin-top: 12px; }
.legend.legend-top { margin: 2px 0 4px; }
.lg { font-size: 12px; color: var(--ink-2); display: inline-flex; gap: 6px;
  align-items: center; }
.sw { width: 10px; height: 10px; border-radius: 3px; display: inline-block; }
section { margin-top: 44px; }
.sec-head { display: flex; align-items: baseline; gap: 16px; flex-wrap: wrap;
  border-bottom: 1px solid var(--grid); padding-bottom: 9px; margin-bottom: 10px; }
.sec-head h2 { font-size: 19px; margin: 0; letter-spacing: -0.01em; }
.sec-sub { font-size: 13px; color: var(--muted); font-weight: 450; margin-left: 8px; }
.suite-h { font-size: 14px; color: var(--ink-2); margin: 18px 0 6px; }
.note { color: var(--ink-2); font-size: 13px; margin: 6px 0 14px; }
.group { border: 1px solid var(--ring); border-radius: 10px; margin: 10px 0;
  background: var(--surface); box-shadow: 0 1px 2px rgba(0,0,0,0.03); }
.group > summary { cursor: pointer; padding: 11px 16px; display: flex; gap: 12px;
  align-items: baseline; }
.group > summary code { font-size: 13px; font-weight: 600; }
.group-meta { color: var(--muted); font-size: 12.5px; margin-left: auto;
  font-variant-numeric: tabular-nums; display: flex; gap: 10px;
  align-items: baseline; white-space: nowrap; }
.pad { padding: 6px 16px 16px; }
.chip { font-size: 12px; font-weight: 650; white-space: nowrap; }
.st-fail-text { color: var(--bad); } .st-pass-text { color: var(--good-text); }
.st-warn-text { color: var(--warn); }
table { width: 100%; border-collapse: collapse; font-size: 13px; }
th { text-align: left; color: var(--muted); font-weight: 550; font-size: 10.5px;
  text-transform: uppercase; letter-spacing: 0.07em; padding: 4px 16px; }
td { padding: 6px 16px; border-top: 1px solid var(--grid); vertical-align: top; }
tr.case-row:hover td { background: color-mix(in srgb, var(--blue) 5%, transparent); }
.pill { display: inline-flex; gap: 5px; align-items: center; font-weight: 650;
  font-size: 12px; white-space: nowrap; }
.pill.st-pass { color: var(--good-text); } .pill.st-fail { color: var(--bad); }
.pill.st-skip { color: var(--muted); }
.tname { font-family: ui-monospace, "SF Mono", Menlo, monospace; font-size: 12.5px;
  width: 100%; }
.tdur { white-space: nowrap; font-variant-numeric: tabular-nums; color: var(--ink-2); }
.mini-track { display: inline-block; width: 68px; height: 6px; background: var(--grid);
  border-radius: 3px; margin-right: 8px; vertical-align: middle; }
.mini-bar { display: block; height: 6px; min-width: 2px; background: var(--blue);
  border-radius: 3px; }
details pre { white-space: pre-wrap; background: var(--band); border-radius: 8px;
  padding: 10px 12px; font-size: 11.5px; max-height: 320px; overflow: auto; }
details summary { cursor: pointer; }
td details summary { color: var(--blue); font-size: 12px; }
#search { width: 100%; padding: 11px 16px; font-size: 14px; margin: 26px 0 4px;
  border: 1px solid var(--ring); border-radius: 10px; background: var(--surface);
  color: inherit; }
footer { margin-top: 52px; color: var(--muted); font-size: 12px;
  border-top: 1px solid var(--grid); padding-top: 14px; }
"""

JS = """
document.getElementById('search').addEventListener('input', function () {
  const q = this.value.toLowerCase();
  document.querySelectorAll('tr.case-row').forEach(function (r) {
    r.style.display = r.dataset.name.includes(q) ? '' : 'none';
  });
  document.querySelectorAll('details.group').forEach(function (g) {
    const rows = g.querySelectorAll('tr.case-row');
    if (!rows.length) return;
    const any = Array.from(rows).some(function (r) { return r.style.display !== 'none'; });
    g.style.display = any ? '' : 'none';
    if (q) g.open = any;
  });
});
"""


def inline_logo():
    # Inline the SVG as a real element with SMIL animations stripped: GitHub
    # Pages' deployment pipeline rejects pages whose SVG carries <animate>/
    # <animateTransform> (bisected 2026-07-03 — the identical page with a
    # static logo deploys, the animated one fails).
    p = Path(__file__).resolve().parents[1] / "assets" / "logo.svg"
    if p.exists() and p.stat().st_size < 64_000:
        svg = p.read_text(encoding="utf-8")
        svg = re.sub(r"<\?xml[^>]*\?>", "", svg)
        svg = re.sub(r"<animateTransform.*?/>", "", svg, flags=re.S)
        svg = re.sub(r"<animate .*?/>", "", svg, flags=re.S)
        svg = re.sub(r"<svg ", '<svg class="logo" ', svg, count=1)
        svg = re.sub(r'(<svg[^>]*?) width="\d+" height="\d+"', r"\1", svg, count=1)
        return svg
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xml", nargs="*", help="JUnit XML files")
    ap.add_argument("-o", "--output", default="test-report.html")
    ap.add_argument("--label", action="append", default=[], metavar="FILE=LABEL")
    ap.add_argument("--platform-link", action="append", default=[], metavar="NAME=HREF",
                    help="platform switcher entry; HREF '.' marks this report's platform")
    ap.add_argument("--archive-href", default="../",
                    help="relative href of the reports index (default ../)")
    ap.add_argument("--validation", default=None,
                    help="validation regression_results.json")
    ap.add_argument("--benchmark", default=None,
                    help="validation benchmark_history.json")
    args = ap.parse_args()

    labels = dict(kv.split("=", 1) for kv in args.label)
    suites = []
    for path in args.xml:
        if not Path(path).exists():
            print(f"warning: {path} missing, skipped", file=sys.stderr)
            continue
        meta, cases = parse_junit(path)
        suites.append((labels.get(path, Path(path).stem), cases, meta))

    code_rows, physics_cases, all_cases = [], [], []
    for label, cases, _ in suites:
        phys = [c for c in cases if any(m in c["classname"] for m in PHYSICS_MODULES)]
        rest = [c for c in cases if c not in phys]
        all_cases += cases
        physics_cases += phys
        if rest:
            code_rows.append((label, rest))

    val_html, vn_pass, vn_fail = "", 0, 0
    if args.validation and Path(args.validation).exists():
        val_html, vn_pass, vn_fail = validation_section(args.validation)
    bench_html, bn_fail = "", 0
    if args.benchmark and Path(args.benchmark).exists():
        bench_html, bn_fail = benchmark_section(args.benchmark)

    if not suites and not val_html and not bench_html:
        sys.exit("no inputs found")

    n = counts(all_cases)
    total_time = sum(m["time"] for _, _, m in suites)
    n_fail_total = n["failed"] + vn_fail + bn_fail
    ok = n_fail_total == 0
    n_code = sum(len(c) for _, c in code_rows)
    n_phys = len(physics_cases) + vn_pass + vn_fail

    tiles = f"""
<div class="tiles">
<div class="tile {'hero-good' if ok else 'hero-bad'}"><div class="k">verdict</div>
  <div class="v">{"✓ PASS" if ok else f"✕ {n_fail_total}"}</div>
  <div class="s">{len(all_cases)} tests · {vn_pass + vn_fail} validation checks</div></div>
<div class="tile"><div class="k">code correctness</div><div class="v">{n_code}</div>
  <div class="s">C++ unit + Python API tests</div></div>
<div class="tile"><div class="k">physics correctness</div><div class="v">{n_phys}</div>
  <div class="s">closure · golden · invariants{" · validation" if val_html else ""}</div></div>
<div class="tile"><div class="k">failed</div><div class="v">{n_fail_total}</div>
  <div class="s">{n["skipped"]} skipped</div></div>
<div class="tile"><div class="k">suite runtime</div><div class="v">{fmt_t(total_time)}</div>
  <div class="s">unit + physics tests</div></div>
</div>"""

    comp_rows = code_rows + ([("Physics tests (pytest)", physics_cases)]
                             if physics_cases else [])
    nav = ['<a href="#overview">Overview</a>',
           f'<a href="#code">Code correctness <span class="n">{n_code}</span></a>',
           f'<a href="#physics">Physics tests <span class="n">{len(physics_cases)}</span></a>']
    if val_html:
        nav.append(f'<a href="#validation">Physics validation '
                   f'<span class="n">{vn_pass + vn_fail}</span></a>')
    if bench_html:
        nav.append('<a href="#performance">Performance</a>')

    plat_entries = []
    for spec in args.platform_link:
        name, _, href = spec.partition("=")
        plat_entries.append(
            f'<span class="plat-cur">{esc(name)}</span>' if href in (".", "")
            else f'<a class="plat" href="{esc(href)}">{esc(name)}</a>')
    platform_nav = (('<span class="plat-label">platform:</span>' + "".join(plat_entries))
                    if plat_entries else "")

    test_docs = collect_test_docs()
    code_sections = "".join(
        f'<h3 class="suite-h">{esc(label)}</h3>{suite_table(cases, test_docs)}'
        for label, cases in code_rows)

    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>VegasAfterglow tests — {"PASS" if ok else "FAIL"}</title>
<style>{CSS}</style></head><body><div class="wrap">
<div class="masthead">{inline_logo()}<div><h1>VegasAfterglow test &amp; validation report</h1>
<div class="meta">{datetime.now().strftime("%Y-%m-%d %H:%M")} · commit {esc(git_commit())}
 · {esc(platform.system())} {esc(platform.machine())} · Python {platform.python_version()}</div>
</div><span class="verdict {'ok' if ok else 'bad'}">{"✓ ALL PASSING" if ok
    else f"✕ {n_fail_total} FAILING"}</span></div>
<nav>{"".join(nav)}{platform_nav}<a class="archive" href="{esc(args.archive_href)}" title="All published report versions">all versions ↗</a></nav>
<section id="overview">{tiles}
<div class="card"><h3>Outcomes by suite</h3>
{stacked_status_bar(outcome_rows(comp_rows), "Outcomes by suite")}</div>
<div class="card"><h3>Slowest tests</h3>{slowest_chart(all_cases)}</div>
</section>
<input id="search" type="search" placeholder="Filter tests… (file or test name)">
<section id="code"><div class="sec-head"><h2>Code correctness
<span class="sec-sub">unit &amp; API contracts</span></h2></div>{code_sections}</section>
<section id="physics"><div class="sec-head"><h2>Physics tests
<span class="sec-sub">closure relations, exact invariants, golden baselines</span></h2>
</div>{suite_table(physics_cases, test_docs) if physics_cases else "<p class='note'>none</p>"}
</section>
{val_html}
{bench_html}
<footer>Generated by <code>tests/report.py</code> · run <code>make test-report</code>
locally · methodology in TESTING.md</footer>
</div><script>{JS}</script></body></html>"""

    Path(args.output).write_text(doc, encoding="utf-8")
    print(f"report: {args.output}  ({len(all_cases)} tests + {vn_pass + vn_fail} "
          f"validation checks + benchmark{' ✓' if bench_html else ' —'}, "
          f"{n_fail_total} failing)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
