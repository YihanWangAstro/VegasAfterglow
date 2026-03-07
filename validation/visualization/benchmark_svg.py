"""Generate benchmark performance SVGs — one per radiation mode.

Called automatically by run_validation.py after benchmark runs.
Outputs to assets/ in the project root.
"""

import json
import os
from pathlib import Path

_ROOT = Path(__file__).parent.parent.parent
_ASSETS = _ROOT / "assets"

STAGES = [
    ("dynamics",      "#D49B40", "Shock Dynamics"),
    ("EAT_grid",      "#4DA5BB", "Equal-Time Grid"),
    ("syn_electrons",  "#45AB8A", "Syn. Electrons"),
    ("syn_photons",    "#D46565", "Syn. Photons"),
    ("sync_flux",      "#5B8ADB", "Flux Integration"),
    ("cooling",        "#9B72CF", "IC Cooling"),
    ("ic_photons",     "#D4693A", "IC Photons"),
    ("ssc_flux",       "#52A87C", "SSC Flux"),
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


def _make_chart(bench, rad_type, rad_label, out_path):
    cpu_label = bench.get("cpu", "")

    cfg_index = {}
    for c in bench["configs"]:
        if c["radiation"] != rad_type:
            continue
        key = (c["jet_type"], c["medium"], float(c["theta_obs_ratio"]))
        cfg_index[key] = c

    groups = []
    for jet in JETS:
        for med in MEDIAS:
            bars = []
            for angle in ANGLES:
                c = cfg_index.get((jet, med, angle))
                if c:
                    sb = c["timing"].get("stage_breakdown", {})
                    bars.append({"total": c["timing"]["total_ms"],
                                 "segs": {k: sb.get(k, 0.0) for k in ALL_STAGE_KEYS}})
                else:
                    bars.append(None)
            groups.append({"jet": jet, "med": med, "bars": bars})

    SVG_W, SVG_H = 950, 360
    X_LEFT  = 52
    X_RIGHT = 935
    CHART_W = X_RIGHT - X_LEFT
    Y_BOT   = 265
    Y_TOP   = 96
    CHART_H = Y_BOT - Y_TOP

    N_GROUPS = len(groups)
    N_BARS   = 4
    INTRA    = 3
    INTER    = 16
    BAR_W    = int((CHART_W - (N_GROUPS - 1) * INTER
                    - N_GROUPS * (N_BARS - 1) * INTRA) / (N_GROUPS * N_BARS))
    GROUP_W  = N_BARS * BAR_W + (N_BARS - 1) * INTRA

    all_totals = [b["total"] for g in groups for b in g["bars"] if b]
    raw_max = max(all_totals) if all_totals else 10.0
    grid_interval = 1.0
    for factor in [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500]:
        if factor * 3 >= raw_max * 1.10:
            grid_interval = factor
            break
    Y_MAX = grid_interval * 3

    def gx(gi):   return X_LEFT + gi * (GROUP_W + INTER)
    def bx(gi, bi): return gx(gi) + bi * (BAR_W + INTRA)
    def bcx(gi, bi): return bx(gi, bi) + BAR_W / 2
    def gcx(gi):  return gx(gi) + GROUP_W / 2
    def fy(ms):   return Y_BOT - ms / Y_MAX * CHART_H

    L = []
    def emit(s=""): L.append(s)

    emit('<?xml version="1.0" encoding="UTF-8"?>')
    emit(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {SVG_W} {SVG_H}"')
    emit(f'     font-family="{FONT}">')
    emit()

    emit('  <!-- Background -->')
    emit(f'  <rect width="{SVG_W}" height="{SVG_H}" fill="#0d1117" rx="12"/>')
    emit()

    emit('  <!-- Title -->')
    emit(f'  <text x="475" y="28" text-anchor="middle" font-size="15" font-weight="700"'
         f' fill="#E6EDF3" letter-spacing="0.3">Performance \u00b7 {rad_label}</text>')
    cpu_note = (f"CPU time by stage \u00b7 {cpu_label} \u00b7 single core \u00b7 default resolution"
                if cpu_label else "CPU time by stage \u00b7 single core \u00b7 default resolution")
    emit(f'  <text x="475" y="44" text-anchor="middle" font-size="9" fill="#A8B8CC">{cpu_note}</text>')
    emit(f'  <text x="475" y="57" text-anchor="middle" font-size="9" fill="#A8B8CC">'
         f'each group: \u03b8v/\u03b8c\u202f=\u202f0\u2002\u2502\u20021\u2002\u2502\u20022\u2002\u2502\u20024\u2002(left \u2192 right)</text>')
    emit()

    active_stages = [s for s in STAGES if any(
        b and b["segs"].get(s[0], 0) > 0 for g in groups for b in g["bars"])]
    emit('  <!-- Legend -->')
    item_widths = [12 + 5 + len(label) * 6 for _, _, label in active_stages]
    total_lw = sum(item_widths) + 14 * (len(active_stages) - 1)
    lx = 475 - total_lw // 2
    for key, color, label in active_stages:
        emit(f'  <rect x="{lx}" y="67" width="12" height="12" rx="2" fill="{color}"/>')
        emit(f'  <text x="{lx + 16}" y="78" font-size="10" fill="#D0DDE8">{label}</text>')
        lx += 12 + 5 + len(label) * 6 + 14
    emit()


    emit('  <!-- Gridlines -->')
    emit(f'  <line x1="{X_LEFT}" y1="{Y_BOT}" x2="{X_RIGHT}" y2="{Y_BOT}"'
         ' stroke="#2D3345" stroke-width="1"/>')
    for n in range(1, 4):
        ms_val = n * grid_interval
        gy = fy(ms_val)
        emit(f'  <line x1="{X_LEFT}" y1="{gy:.1f}" x2="{X_RIGHT}" y2="{gy:.1f}"'
             ' stroke="#2D3345" stroke-width="1" stroke-dasharray="4,3"/>')
        label_ms = f"{ms_val:.0f} ms" if ms_val < 1000 else f"{ms_val/1000:.1f} s"
        emit(f'  <text x="{X_LEFT - 4}" y="{gy + 4:.1f}" text-anchor="end"'
             f' font-size="8" fill="#A8B8CC">{label_ms}</text>')
    emit(f'  <text x="{X_LEFT - 4}" y="{Y_BOT + 4}" text-anchor="end"'
         ' font-size="8" fill="#A8B8CC">0</text>')
    emit()

    emit('  <!-- Stacked bars -->')
    for gi, group in enumerate(groups):
        for bi, bar in enumerate(group["bars"]):
            if bar is None:
                continue
            cx  = bcx(gi, bi)
            bot = 0.0
            for key, color, _ in STAGES:
                ms = bar["segs"].get(key, 0.0)
                if ms <= 0:
                    continue
                ry  = fy(bot + ms)
                rh  = ms / Y_MAX * CHART_H
                emit(f'  <rect x="{bx(gi, bi):.1f}" y="{ry:.1f}" width="{BAR_W}"'
                     f' height="{rh:.1f}" fill="{color}"/>')
                bot += ms
            top_y = fy(bar["total"])
            emit(f'  <text x="{cx:.1f}" y="{top_y - 4:.1f}" text-anchor="middle"'
                 f' font-size="7.5" font-weight="600" fill="#D0DDE8">{bar["total"]:.1f}</text>')
    emit()

    emit('  <!-- Angle labels -->')
    for gi in range(N_GROUPS):
        for bi, ang_lbl in enumerate(["0", "1", "2", "4"]):
            cx = bcx(gi, bi)
            emit(f'  <text x="{cx:.1f}" y="{Y_BOT + 13}" text-anchor="middle"'
                 f' font-size="7.5" fill="#A8B8CC">{ang_lbl}</text>')
    emit()

    emit('  <!-- Group labels -->')
    for gi, group in enumerate(groups):
        cx = gcx(gi)
        emit(f'  <text x="{cx:.1f}" y="{Y_BOT + 28}" text-anchor="middle"'
             f' font-size="9.5" font-weight="700" fill="#E6EDF3">{JET_NAMES[group["jet"]]}</text>')
        med_color = "#5B8ADB" if group["med"] == "ISM" else "#45AB8A"
        emit(f'  <text x="{cx:.1f}" y="{Y_BOT + 41}" text-anchor="middle"'
             f' font-size="8.5" fill="{med_color}">{group["med"]}</text>')
    emit()

    mid_y = (Y_TOP + Y_BOT) // 2
    emit('  <!-- Y-axis label -->')
    emit(f'  <text x="13" y="{mid_y}" text-anchor="middle" font-size="9" fill="#A8B8CC"'
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
        _make_chart(bench, rad_type, rad_label, _ASSETS / filename)
