"""Generate contribution-chart.svg — native SVG, no matplotlib.

Edit DATA below to update fractions. Each tuple: (human, claude_code, codex) — must sum to 1.0.
Run:  python script/developer/generate_contribution_chart.py
"""

from pathlib import Path

# ── DATA ────────────────────────────────────────────────────────────────────
MODULES = [
    "C++ Engine",
    "Python Bindings",
    "Python Wrapper",
    "Tests",
    "Webtool",
    "Docs/Code Comments",
]
ACTIVITIES = ["Architect.", "Coding", "Optim."]

# (human, claude_code, codex)  — each row must sum to 1.0
DATA = {
    "C++ Engine":         {"Architect.": (1.00, 0.00, 0.00), "Coding": (0.90, 0.10, 0.00), "Optim.": (0.80, 0.20, 0.00)},
    "Python Bindings":    {"Architect.": (1.00, 0.00, 0.00), "Coding": (0.70, 0.30, 0.00), "Optim.": (1.00, 0.00, 0.00)},
    "Python Wrapper":     {"Architect.": (1.00, 0.00, 0.00), "Coding": (0.10, 0.90, 0.00), "Optim.": (0.30, 0.70, 0.00)},
    "Tests":              {"Architect.": (0.70, 0.30, 0.00), "Coding": (0.20, 0.80, 0.00), "Optim.": (0.20, 0.80, 0.00)},
    "Webtool":            {"Architect.": (0.50, 0.00, 0.50), "Coding": (0.00, 0.00, 1.00), "Optim.": (0.50, 0.30, 0.20)},
    "Docs/Code Comments": {"Architect.": (0.30, 0.70, 0.00), "Coding": (0.20, 0.80, 0.00), "Optim.": (0.00, 1.00, 0.00)},
}

# ── Colors — architecture.svg palette ───────────────────────────────────────
COLOR_HUMAN  = "#45AB8A"   # green  (Environment column in architecture.svg)
COLOR_CLAUDE = "#D4693A"   # muted orange (desaturated shock front)
COLOR_CODEX  = "#4A90D9"   # blue   (OpenAI brand palette)

FONT = "-apple-system, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif"

PALETTES = {
    "dark": {
        "cp": "#E6EDF3",
        "cs": "#A8B8CC",
        "pct": "#FFFFFF",
        "gl": "#2D3345",
        "arr": "#3D4A5C",
    },
    "light": {
        "cp": "#24292F",
        "cs": "#57606A",
        "pct": "#FFFFFF",
        "gl": "#D0D7DE",
        "arr": "#57606A",
    },
}

# ── Layout ───────────────────────────────────────────────────────────────────
SVG_W, SVG_H = 950, 340
BAR_W    = 32
INTRA    = 8                              # gap between bars in a group
GROUP_W  = 3 * BAR_W + 2 * INTRA         # 112
X_LEFT   = 60                            # left edge of chart area
X_RIGHT  = 930                           # right edge of chart area
CHART_W  = X_RIGHT - X_LEFT              # 870
N_MOD    = len(MODULES)
INTER    = (CHART_W - N_MOD * GROUP_W) / (N_MOD - 1)   # ~39.6
Y_BOT    = 255    # SVG y coordinate for fraction 0.0 (bar bottom)
Y_TOP    = 75     # SVG y coordinate for fraction 1.0 (bar top)
CHART_H  = Y_BOT - Y_TOP                 # 180 px


def gx(m):
    """Left x of module group m."""
    return X_LEFT + m * (GROUP_W + INTER)


def bar_cx(m, a):
    """Center x of bar (module index m, activity index a)."""
    return gx(m) + a * (BAR_W + INTRA) + BAR_W / 2


def fy(frac):
    """SVG y coordinate for a given fraction value (0–1)."""
    return Y_BOT - frac * CHART_H


def _output_path(theme):
    name = "contribution-chart.svg" if theme == "adaptive" else f"contribution-chart-{theme}.svg"
    return Path(__file__).parent.parent.parent / "assets" / name


def _emit_style(emit, theme):
    dark = PALETTES["dark"]
    light = PALETTES["light"]
    emit('  <style>')
    if theme == "adaptive":
        emit('    svg { color-scheme: light dark; }')
        emit(f'    .cp {{ fill: {dark["cp"]}; }}')
        emit(f'    .cs {{ fill: {dark["cs"]}; }}')
        emit(f'    .pct {{ fill: {dark["pct"]}; }}')
        emit(f'    .gl {{ stroke: {dark["gl"]}; fill: none; }}')
        emit(f'    .arr {{ fill: {dark["arr"]}; }}')
        emit('    @media (prefers-color-scheme: light) {')
        emit(f'      .cp {{ fill: {light["cp"]}; }}')
        emit(f'      .cs {{ fill: {light["cs"]}; }}')
        emit(f'      .pct {{ fill: {light["pct"]}; }}')
        emit(f'      .gl {{ stroke: {light["gl"]}; }}')
        emit(f'      .arr {{ fill: {light["arr"]}; }}')
        emit('    }')
    elif theme == "light":
        emit(f'    .cp {{ fill: {light["cp"]}; }}')
        emit(f'    .cs {{ fill: {light["cs"]}; }}')
        emit(f'    .pct {{ fill: {light["pct"]}; }}')
        emit(f'    .gl {{ stroke: {light["gl"]}; fill: none; }}')
        emit(f'    .arr {{ fill: {light["arr"]}; }}')
    elif theme == "dark":
        emit(f'    .cp {{ fill: {dark["cp"]}; }}')
        emit(f'    .cs {{ fill: {dark["cs"]}; }}')
        emit(f'    .pct {{ fill: {dark["pct"]}; }}')
        emit(f'    .gl {{ stroke: {dark["gl"]}; fill: none; }}')
        emit(f'    .arr {{ fill: {dark["arr"]}; }}')
    else:
        raise ValueError(f"Unknown theme: {theme}")
    emit('  </style>')


def _write_chart(theme):
    L = []

    def emit(s=""):
        L.append(s)

    emit('<?xml version="1.0" encoding="UTF-8"?>')
    emit(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {SVG_W} {SVG_H}"')
    emit(f'     font-family="{FONT}">')
    emit()

    _emit_style(emit, theme)
    emit()

    emit('  <defs>')
    emit('    <marker id="arr" markerWidth="7" markerHeight="5" refX="7" refY="2.5" orient="auto">')
    emit('      <polygon points="0 0, 7 2.5, 0 5" class="arr"/>')
    emit('    </marker>')
    emit('  </defs>')
    emit()

    emit(f'  <rect width="{SVG_W}" height="{SVG_H}" fill="transparent" rx="12"/>')
    emit()

    emit('  <text x="475" y="30" text-anchor="middle" font-size="16" font-weight="700"'
         ' class="cp" letter-spacing="0.3">Contribution Overview</text>')
    emit()

    for lx, label, color in [(356, "Human", COLOR_HUMAN), (429, "Claude Code", COLOR_CLAUDE), (542, "CodeX", COLOR_CODEX)]:
        emit(f'  <rect x="{lx}" y="44" width="14" height="14" rx="2" fill="{color}"/>')
        emit(f'  <text x="{lx + 18}" y="56" font-size="11" class="cs">{label}</text>')
    emit()

    emit(f'  <line x1="{X_LEFT}" y1="{Y_BOT}" x2="{X_RIGHT}" y2="{Y_BOT}" class="gl" stroke-width="1"/>')
    for pct in [0.25, 0.50, 0.75]:
        gy = fy(pct)
        emit(f'  <line x1="{X_LEFT}" y1="{gy:.1f}" x2="{X_RIGHT}" y2="{gy:.1f}" class="gl" stroke-width="1" stroke-dasharray="4,3"/>')
        emit(f'  <text x="{X_LEFT - 5}" y="{gy + 4:.1f}" text-anchor="end" font-size="8.5" class="cs">{pct * 100:.0f}%</text>')
    emit(f'  <text x="{X_LEFT - 5}" y="{Y_BOT + 4}" text-anchor="end" font-size="8.5" class="cs">0%</text>')
    emit(f'  <text x="{X_LEFT - 5}" y="{Y_TOP + 4}" text-anchor="end" font-size="8.5" class="cs">100%</text>')
    emit()

    for m, mod in enumerate(MODULES):
        for a, act in enumerate(ACTIVITIES):
            bcx = bar_cx(m, a)
            bot = 0.0
            for frac, color in zip(DATA[mod][act], [COLOR_HUMAN, COLOR_CLAUDE, COLOR_CODEX]):
                if frac > 0:
                    rx = bcx - BAR_W / 2
                    emit(f'  <rect x="{rx:.1f}" y="{fy(bot + frac):.1f}" width="{BAR_W}" height="{frac * CHART_H:.1f}" fill="{color}"/>')
                    if frac >= 0.08:
                        emit(f'  <text x="{bcx:.1f}" y="{fy(bot + frac / 2) + 4:.1f}" text-anchor="middle"'
                             f' font-size="9" font-weight="700" class="pct">{frac * 100:.0f}%</text>')
                    bot += frac
    emit()

    for m, mod in enumerate(MODULES):
        for a, act in enumerate(ACTIVITIES):
            bcx = bar_cx(m, a)
            emit(f'  <text x="{bcx:.1f}" y="{Y_BOT + 15}" text-anchor="middle" font-size="8.5" class="cs">{act}</text>')
    emit()

    for m, mod in enumerate(MODULES):
        mid = (bar_cx(m, 0) + bar_cx(m, len(ACTIVITIES) - 1)) / 2
        emit(f'  <text x="{mid:.1f}" y="{Y_BOT + 33}" text-anchor="middle" font-size="10" font-weight="700" class="cp">{mod}</text>')
    emit()

    arr_y = Y_BOT + 55
    emit(f'  <line x1="{X_LEFT}" y1="{arr_y}" x2="{X_RIGHT}" y2="{arr_y}" class="gl" stroke-width="1" marker-end="url(#arr)"/>')
    emit(f'  <text x="{X_LEFT}" y="{arr_y + 13}" font-size="8.5" class="cs" font-style="italic">Machine Language</text>')
    emit(f'  <text x="{X_RIGHT}" y="{arr_y + 13}" text-anchor="end" font-size="8.5" class="cs" font-style="italic">Natural Language</text>')
    emit()

    emit('</svg>')

    out = _output_path(theme)
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"Saved: {out}")


def main():
    for theme in ("adaptive", "light", "dark"):
        _write_chart(theme)


if __name__ == "__main__":
    main()
