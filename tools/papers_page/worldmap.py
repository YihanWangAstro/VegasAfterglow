"""Static world-map SVGs for the README.

Renders the institute markers from ``data/papers.json`` over a
Natural Earth 110m land outline (public domain, vendored as
``world_land.json``) in equirectangular projection. Light and dark
variants match the README benchmark charts; the weekly workflow
publishes them next to the papers page so the README embed updates
automatically.
"""

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build import aggregate_institutes, load_excludes, is_excluded  # noqa: E402

PALETTES = {
    "light": {"sea": "#f7fafd", "land": "#e3e1da", "coast": "#cfcdc4",
              "primary": "#5B8ADB", "co": "#8b8a84", "text": "#52514e"},
    "dark": {"sea": "#0e0f11", "land": "#26272a", "coast": "#3a3b3e",
             "primary": "#6d99e3", "co": "#8b8a84", "text": "#c3c2b7"},
}
# crop empty polar regions: lat 72N .. 56S
LAT_TOP, LAT_BOT = 72.0, -56.0


def render(institutes, palette, W=1000):
    land = json.loads((HERE / "world_land.json").read_text())
    lw, lh = land["width"], land["height"]
    y0 = (90 - LAT_TOP) / 180 * lh
    y1 = (90 - LAT_BOT) / 180 * lh
    H = round((y1 - y0) / lw * W)
    sx = W / lw
    sy = H / (y1 - y0)
    c = PALETTES[palette]

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
        f'font-family="-apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif">',
        f'<rect width="{W}" height="{H}" rx="10" fill="{c["sea"]}"/>',
        f'<g transform="scale({sx:.4f} {sy:.4f}) translate(0 {-y0:.1f})">',
        f'<g fill="{c["land"]}" stroke="{c["coast"]}" stroke-width="0.7">',
    ]
    parts.extend(f'<path d="{d}"/>' for d in land["paths"])
    parts.append("</g></g>")

    def xy(lat, lng):
        x = (lng + 180) / 360 * W
        y = ((90 - lat) / 180 * lh - y0) * sy
        return x, y

    # co-author tier under primary tier
    for tier, color, base, opac in (("co", c["co"], 3.4, 0.45),
                                    ("primary", c["primary"], 4.6, 0.62)):
        parts.append(f'<g fill="{color}" fill-opacity="{opac}" '
                     f'stroke="{color}" stroke-width="1">')
        for i in institutes:
            n = len(i[tier])
            if not n or i.get("lat") is None or i.get("lng") is None:
                continue
            if tier == "co" and i["primary"]:
                continue
            x, y = xy(i["lat"], i["lng"])
            r = base + 2.2 * math.sqrt(n)
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}"/>')
        parts.append("</g>")

    n_primary = sum(1 for i in institutes if i["primary"])
    n_co = sum(1 for i in institutes if not i["primary"])
    lx, ly = 14, H - 16
    parts.append(
        f'<g font-size="13" fill="{c["text"]}">'
        f'<circle cx="{lx}" cy="{ly - 4}" r="5.5" fill="{c["primary"]}" fill-opacity="0.62"/>'
        f'<text x="{lx + 11}" y="{ly}">first / corresponding author '
        f'({n_primary})</text>'
        f'<circle cx="{lx + 240}" cy="{ly - 4}" r="4.5" fill="{c["co"]}" fill-opacity="0.45"/>'
        f'<text x="{lx + 251}" y="{ly}">co-author ({n_co})</text>'
        "</g>")
    parts.append("</svg>")
    return "".join(parts)


def render_summary(stats, palette, W=900, H=42):
    c = PALETTES[palette]
    accent = c["primary"]
    parts = " &#183;  ".join(
        f'<tspan font-weight="700" fill="{accent}">{v}</tspan> {label}'
        for v, label in stats)
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
        f'font-family="-apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif">'
        f'<text x="{W // 2}" y="{H // 2 + 5}" text-anchor="middle" '
        f'font-size="15.5" fill="{c["text"]}">{parts}</text></svg>')


def main():
    data = json.loads((HERE / "data" / "papers.json").read_text())
    excludes = load_excludes()
    papers = [p for p in data["papers"]
              if p.get("usage") != "cited" and not is_excluded(p, excludes)]
    institutes = aggregate_institutes(papers)
    n_primary = sum(1 for i in institutes if i["primary"])
    countries = {i["country"] for i in institutes if i.get("country")}
    stats = [
        (len(papers), "papers"),
        (n_primary, "first / corresponding author institutes"),
        (len(institutes) - n_primary, "co-author institutes"),
        (len(countries), "countries"),
        (data["generated"], "last updated"),
    ]
    outdir = HERE / "site"
    outdir.mkdir(exist_ok=True)
    for palette in PALETTES:
        (outdir / f"worldmap-{palette}.svg").write_text(render(institutes, palette))
        (outdir / f"papers-summary-{palette}.svg").write_text(
            render_summary(stats, palette))
    print(f"wrote worldmap + papers-summary SVGs "
          f"({len(papers)} papers, {len(institutes)} institutes)")


if __name__ == "__main__":
    main()
