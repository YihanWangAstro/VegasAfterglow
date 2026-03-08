"""Generate benchmark performance SVGs — one per radiation mode.

Called automatically by run_validation.py after benchmark runs.
Outputs to assets/ in the project root.
"""

import json
import math
from pathlib import Path

_ROOT = Path(__file__).parent.parent.parent
_ASSETS = _ROOT / "assets"

STAGES = [
    ("dynamics",     "#D49B40", "Shock Dynamics"),
    ("EAT_grid",     "#4DA5BB", "EAT Grid"),
    ("syn_electrons", "#45AB8A", "Syn. Electrons"),
    ("syn_photons",  "#D46565", "Syn. Photons"),
    ("sync_flux",    "#5B8ADB", "Flux Integration"),
    ("cooling",      "#9B72CF", "IC Cooling"),
    ("ic_photons",   "#D4693A", "IC Photons"),
    ("ssc_flux",     "#52A87C", "SSC Flux"),
]
ALL_STAGE_KEYS = [k for k, *_ in STAGES]

RAD_INFO = {
    "synchrotron":   ("Synchrotron",   "benchmark-sync.svg"),
    "rvs_sync_thin": ("Reverse Shock", "benchmark-rvs.svg"),
    "full_ssc":      ("Full SSC",      "benchmark-ssc.svg"),
}

JET_NAMES = {"tophat": "Tophat", "gaussian": "Gaussian",
             "powerlaw": "Power-law", "two_component": "Two-comp."}
JETS   = ["tophat", "two_component", "gaussian", "powerlaw"]
MEDIAS = ["ISM", "wind"]
ANGLES = [0.0, 1.0, 2.0, 4.0]

FONT = "-apple-system, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif"


def _variant_name(filename, theme):
    if theme == "adaptive":
        return filename
    stem, suffix = filename.rsplit(".", 1)
    return f"{stem}-{theme}.{suffix}"


def _emit_style(emit, theme):
    emit('  <style>')
    if theme == "adaptive":
        emit('    svg { color-scheme: light dark; }')
        emit('    .cp, .cs, .cb, .mi, .mw { fill: #ffffff; }')
        emit('    .gl { stroke: #2D3345; fill: none; }')
        emit('    @media (prefers-color-scheme: light) {')
        emit('      .cp, .cs, .cb, .mi, .mw { fill: #000000; }')
        emit('      .gl { stroke: #d0d7de; }')
        emit('    }')
    elif theme == "light":
        emit('    .cp, .cs, .cb, .mi, .mw { fill: #000000; }')
        emit('    .gl { stroke: #d0d7de; fill: none; }')
    elif theme == "dark":
        emit('    .cp, .cs, .cb, .mi, .mw { fill: #ffffff; }')
        emit('    .gl { stroke: #2D3345; fill: none; }')
    else:
        raise ValueError(f"Unknown theme: {theme}")
    emit('  </style>')


def _make_chart(bench, rad_type, rad_label, out_path, theme="adaptive"):
    cpu_label = bench.get("cpu", "")

    cfg_index = {(c["jet_type"], c["medium"], float(c["theta_obs_ratio"])): c
                 for c in bench["configs"] if c["radiation"] == rad_type}

    groups = []
    for jet in JETS:
        for med in MEDIAS:
            bars = []
            for angle in ANGLES:
                c = cfg_index.get((jet, med, angle))
                if c:
                    sb = c["timing"].get("stage_breakdown", {})
                    bars.append({k: sb.get(k, 0.0) for k in ALL_STAGE_KEYS})
                else:
                    bars.append(None)
            groups.append({"jet": jet, "med": med, "bars": bars})

    SVG_W, SVG_H = 950, 360
    X_LEFT, X_RIGHT = 52, 935
    Y_BOT, Y_TOP = 265, 96
    CHART_W = X_RIGHT - X_LEFT
    CHART_H = Y_BOT - Y_TOP

    N_GROUPS, N_BARS, INTRA, INTER = len(groups), 4, 3, 16
    BAR_W   = int((CHART_W - (N_GROUPS - 1) * INTER - N_GROUPS * (N_BARS - 1) * INTRA) / (N_GROUPS * N_BARS))
    GROUP_W = N_BARS * BAR_W + (N_BARS - 1) * INTRA

    raw_max = max((sum(b.get(k, 0) for k in ALL_STAGE_KEYS) for g in groups for b in g["bars"] if b), default=10.0)
    ideal   = max(raw_max * 1.05 / 4, 1e-6)
    mag     = 10 ** math.floor(math.log10(ideal))
    grid_interval = next(m * mag for m in [1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 7.5, 8, 10] if m >= ideal / mag)
    Y_MAX   = grid_interval * 4

    def gx(gi):      return X_LEFT + gi * (GROUP_W + INTER)
    def bx(gi, bi):  return gx(gi) + bi * (BAR_W + INTRA)
    def bcx(gi, bi): return bx(gi, bi) + BAR_W / 2
    def gcx(gi):     return gx(gi) + GROUP_W / 2
    def fy(ms):      return Y_BOT - ms / Y_MAX * CHART_H

    L = []
    def emit(s=""): L.append(s)

    emit('<?xml version="1.0" encoding="UTF-8"?>')
    emit(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {SVG_W} {SVG_H}"')
    emit(f'     font-family="{FONT}">')
    emit()

    _emit_style(emit, theme)
    emit()

    emit(f'  <rect width="{SVG_W}" height="{SVG_H}" fill="transparent" rx="12"/>')
    emit()

    emit(f'  <text x="475" y="28" text-anchor="middle" font-size="15" font-weight="700"'
         f' class="cp" letter-spacing="0.3">Performance \u00b7 {rad_label}</text>')
    cpu_note = (f"CPU time by stage \u00b7 {cpu_label} \u00b7 single core \u00b7 default resolution"
                if cpu_label else "CPU time by stage \u00b7 single core \u00b7 default resolution")
    emit(f'  <text x="475" y="44" text-anchor="middle" font-size="9" class="cs">{cpu_note}</text>')
    emit(f'  <text x="475" y="57" text-anchor="middle" font-size="9" class="cs">'
         f'each group: \u03b8v/\u03b8c\u202f=\u202f0\u2002\u2502\u20021\u2002\u2502\u20022\u2002\u2502\u20024\u2002(left \u2192 right)</text>')
    emit()

    active_stages = [s for s in STAGES if any(b and b.get(s[0], 0) > 0 for g in groups for b in g["bars"])]
    total_lw = sum(17 + len(label) * 6 for _, _, label in active_stages) + 14 * (len(active_stages) - 1)
    lx = 475 - total_lw // 2
    for key, color, label in active_stages:
        emit(f'  <rect x="{lx}" y="67" width="12" height="12" rx="2" fill="{color}"/>')
        emit(f'  <text x="{lx + 16}" y="78" font-size="10" class="cb">{label}</text>')
        lx += 17 + len(label) * 6 + 14
    emit()

    emit(f'  <line x1="{X_LEFT}" y1="{Y_BOT}" x2="{X_RIGHT}" y2="{Y_BOT}" class="gl" stroke-width="1"/>')
    for n in range(1, 5):
        ms_val = n * grid_interval
        gy = fy(ms_val)
        emit(f'  <line x1="{X_LEFT}" y1="{gy:.1f}" x2="{X_RIGHT}" y2="{gy:.1f}"'
             ' class="gl" stroke-width="1" stroke-dasharray="4,3"/>')
        label_ms = f"{ms_val/1000:.2g} s" if ms_val >= 1000 else f"{ms_val:.4g} ms"
        emit(f'  <text x="{X_LEFT - 4}" y="{gy + 4:.1f}" text-anchor="end" font-size="8" class="cs">{label_ms}</text>')
    emit(f'  <text x="{X_LEFT - 4}" y="{Y_BOT + 4}" text-anchor="end" font-size="8" class="cs">0</text>')
    emit()

    for gi, group in enumerate(groups):
        for bi, bar in enumerate(group["bars"]):
            if bar is None:
                continue
            cx, bot = bcx(gi, bi), 0.0
            for key, color, _ in STAGES:
                ms = bar.get(key, 0.0)
                if ms <= 0:
                    continue
                emit(f'  <rect x="{bx(gi, bi):.1f}" y="{fy(bot + ms):.1f}" width="{BAR_W}"'
                     f' height="{ms / Y_MAX * CHART_H:.1f}" fill="{color}"/>')
                bot += ms
            emit(f'  <text x="{cx:.1f}" y="{fy(bot) - 4:.1f}" text-anchor="middle"'
                 f' font-size="8" font-weight="600" class="cb">{bot:.1f}</text>')
    emit()

    for gi in range(N_GROUPS):
        for bi, ang_lbl in enumerate(["0", "1", "2", "4"]):
            emit(f'  <text x="{bcx(gi, bi):.1f}" y="{Y_BOT + 13}" text-anchor="middle"'
                 f' font-size="7.5" class="cs">{ang_lbl}</text>')
    emit()

    for gi, group in enumerate(groups):
        cx = gcx(gi)
        emit(f'  <text x="{cx:.1f}" y="{Y_BOT + 28}" text-anchor="middle"'
             f' font-size="9.5" font-weight="700" class="cp">{JET_NAMES[group["jet"]]}</text>')
        emit(f'  <text x="{cx:.1f}" y="{Y_BOT + 41}" text-anchor="middle"'
             f' font-size="8.5" class="{"mi" if group["med"] == "ISM" else "mw"}">{group["med"]}</text>')
    emit()

    mid_y = (Y_TOP + Y_BOT) // 2
    emit(f'  <text x="13" y="{mid_y}" text-anchor="middle" font-size="9" class="cs"'
         f' transform="rotate(-90, 13, {mid_y})">Wall time (ms)</text>')
    emit()
    emit('</svg>')

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"  Saved: {out_path}")


def generate(json_path):
    """Generate all benchmark SVGs from the given benchmark_history.json path."""
    with open(json_path, encoding="utf-8") as f:
        bench = json.load(f)
    _ASSETS.mkdir(parents=True, exist_ok=True)
    for rad_type, (rad_label, filename) in RAD_INFO.items():
        for theme in ("adaptive", "light", "dark"):
            out_name = _variant_name(filename, theme)
            _make_chart(bench, rad_type, rad_label, _ASSETS / out_name, theme=theme)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate benchmark SVG charts from benchmark_history.json")
    parser.add_argument("json_path", nargs="?",
                        default=str(_ROOT / "validation" / "benchmark" / "results" / "benchmark_history.json"),
                        help="Path to benchmark_history.json")
    args = parser.parse_args()
    generate(args.json_path)
