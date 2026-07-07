"""Render the "Papers using VegasAfterglow" page from ``data/papers.json``.

Produces a single self-styled ``site/index.html`` in the same visual
language as the test & validation report (tests/report.py): masthead +
logo, stat tiles, sticky nav, cards. The world map is a Leaflet card
with one marker per institute, sized by paper count; institutes come
from first-author / corresponding-author affiliations.
"""

import html
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
SITE = HERE / "site"
REPO = HERE.parents[1]

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
html { scroll-behavior: smooth; }
section { scroll-margin-top: 64px; }
body { margin: 0; font: 15px/1.55 -apple-system, "Segoe UI", Roboto,
  "Helvetica Neue", Arial, sans-serif; background: var(--plane); color: var(--ink); }
.wrap { max-width: 1040px; margin: 0 auto; padding: 0 28px 80px; }
.masthead { padding: 34px 0 22px; display: flex; align-items: center; gap: 18px;
  flex-wrap: wrap; }
.masthead svg.logo { width: 52px; height: 56px; flex-shrink: 0; }
.masthead h1 { font-size: 25px; letter-spacing: -0.015em; margin: 0; }
.masthead .meta { color: var(--ink-2); font-size: 13px; margin-top: 3px; }
.verdict { margin-left: auto; font-size: 14px; font-weight: 700; padding: 9px 20px;
  border-radius: 999px; letter-spacing: 0.03em;
  background: color-mix(in srgb, var(--blue) 13%, transparent); color: var(--blue);
  border: 1px solid color-mix(in srgb, var(--blue) 38%, transparent); }
nav { position: sticky; top: 0; z-index: 1001;
  background: color-mix(in srgb, var(--plane) 86%, transparent);
  backdrop-filter: blur(10px); border-bottom: 1px solid var(--grid);
  margin: 0 -28px; padding: 10px 28px; display: flex; gap: 8px; flex-wrap: wrap; }
nav a { text-decoration: none; color: var(--ink-2); font-size: 13px; font-weight: 600;
  padding: 6px 14px; border-radius: 999px; border: 1px solid transparent; }
nav a:hover { background: var(--surface); border-color: var(--ring); color: var(--ink); }
nav a .n { color: var(--muted); font-weight: 500; }
nav a.archive { color: var(--muted); margin-left: auto; }
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
  gap: 12px; margin: 22px 0 8px; }
.tile { background: var(--surface); border: 1px solid var(--ring); border-radius: 12px;
  padding: 16px 18px; box-shadow: 0 1px 2px rgba(0,0,0,0.03); }
.tile .k { font-size: 11px; color: var(--muted); text-transform: uppercase;
  letter-spacing: 0.07em; font-weight: 650; }
.tile .v { font-size: 29px; font-weight: 750; font-variant-numeric: tabular-nums;
  letter-spacing: -0.01em; margin-top: 3px; }
.tile .s { font-size: 12px; color: var(--ink-2); margin-top: 2px; }
.card { background: var(--surface); border: 1px solid var(--ring); border-radius: 12px;
  padding: 20px 22px; margin: 14px 0;
  box-shadow: 0 1px 2px rgba(0,0,0,0.03); }
.card h3 { margin: 0 0 14px; font-size: 12.5px; color: var(--ink-2); font-weight: 650;
  text-transform: uppercase; letter-spacing: 0.06em; }
section { margin-top: 44px; }
.sec-head { display: flex; align-items: baseline; gap: 16px; flex-wrap: wrap;
  border-bottom: 1px solid var(--grid); padding-bottom: 9px; margin-bottom: 10px; }
.sec-head h2 { font-size: 19px; margin: 0; letter-spacing: -0.01em; }
.sec-sub { font-size: 13px; color: var(--muted); font-weight: 450; margin-left: 8px; }
footer { margin-top: 52px; color: var(--muted); font-size: 12px;
  border-top: 1px solid var(--grid); padding-top: 14px; }
footer p { margin: 5px 0; }
.leaflet-container { background: var(--band); font: inherit; }
.leaflet-popup-content-wrapper, .leaflet-popup-tip {
  background: var(--surface); color: var(--ink);
  box-shadow: 0 3px 14px rgba(0,0,0,0.3); }
.leaflet-container a.leaflet-popup-close-button { color: var(--muted); }
.leaflet-bar a { background: var(--surface); color: var(--ink);
  border-bottom-color: var(--grid); }
.leaflet-bar a:hover { background: var(--band); }
.leaflet-control-attribution { background: color-mix(in srgb,
  var(--surface) 80%, transparent) !important; color: var(--muted); }
.leaflet-control-attribution a { color: var(--ink-2); }
.inst-pop { display: flex; gap: 10px; align-items: flex-start; max-width: 260px; }
.inst-pop .nm { font-weight: 700; font-size: 13px; }
.inst-pop .ct { color: var(--muted); font-size: 11.5px; }
.inst-pop ul { margin: 6px 0 0; padding-left: 16px; font-size: 11.5px; }
.inst-cols { display: grid; grid-template-columns: 1fr 1fr; gap: 0 26px; }
@media (max-width: 720px) { .inst-cols { grid-template-columns: 1fr; } }
.inst-table { width: 100%; border-collapse: collapse; font-size: 13px;
  align-self: start; }
.inst-table th { text-align: left; color: var(--muted); font-weight: 550;
  font-size: 10.5px; text-transform: uppercase; letter-spacing: 0.07em;
  padding: 4px 12px; }
.inst-table td { padding: 7px 12px; border-top: 1px solid var(--grid);
  vertical-align: middle; }
.inst-table tr:hover td { background: color-mix(in srgb, var(--blue) 6%,
  transparent); }
.inst-table .cnt { font-variant-numeric: tabular-nums; font-weight: 700; }
.inst-table .cnt.co { color: var(--muted); font-weight: 550; }
.legend { display: flex; gap: 14px; flex-wrap: wrap; margin-top: 12px; }
.lg { font-size: 12px; color: var(--ink-2); display: inline-flex; gap: 6px;
  align-items: center; }
.sw { width: 11px; height: 11px; border-radius: 50%; display: inline-block; }
.papers-grid { display: grid; gap: 14px;
  grid-template-columns: repeat(auto-fill, minmax(228px, 1fr)); margin-top: 12px; }
.pcard { background: var(--surface); border: 1px solid var(--ring);
  border-radius: 12px; overflow: hidden; display: flex; flex-direction: column;
  box-shadow: 0 1px 2px rgba(0,0,0,0.03);
  transition: transform 0.16s ease, box-shadow 0.16s ease,
    border-color 0.16s ease; }
.pcard:hover { transform: translateY(-3px);
  box-shadow: 0 8px 22px rgba(0,0,0,0.10);
  border-color: color-mix(in srgb, var(--blue) 45%, var(--ring)); }
.pcard .thumb { display: block; width: 100%; aspect-ratio: 16 / 10;
  object-fit: cover; object-position: left top; background: #fff;
  border-bottom: 1px solid var(--grid);
  transition: transform 0.25s ease; transform-origin: left top; }
.pcard:hover .thumb:not(.placeholder) { transform: scale(1.035); }
.pcard .thumb.placeholder { display: flex; align-items: center;
  justify-content: center; color: var(--muted); font-size: 26px;
  background: var(--band); }
.pcard .body { padding: 12px 14px 13px; display: flex; flex-direction: column;
  gap: 7px; flex: 1; }
.pcard h4 { margin: 0; font-size: 13.5px; line-height: 1.35;
  letter-spacing: -0.005em; display: -webkit-box; -webkit-line-clamp: 3;
  -webkit-box-orient: vertical; overflow: hidden; }
.pcard h4 a { color: inherit; text-decoration: none; }
.pcard h4 a:hover { color: var(--blue); }
.pcard .byline { color: var(--muted); font-size: 11.5px; }
.pcard .foot { margin-top: auto; padding-top: 8px; display: flex; gap: 6px;
  align-items: center; }
.pcard .foot .links { margin-left: auto; display: flex; gap: 5px; }
.pcard .foot .links a { text-decoration: none; font-size: 10.5px;
  font-weight: 700; color: var(--blue); padding: 2px 9px; border-radius: 999px;
  border: 1px solid color-mix(in srgb, var(--blue) 35%, transparent); }
.pcard .foot .links a:hover {
  background: color-mix(in srgb, var(--blue) 10%, transparent); }
"""


def inline_logo():
    p = REPO / "assets" / "logo.svg"
    if p.exists() and p.stat().st_size < 64_000:
        svg = p.read_text(encoding="utf-8")
        svg = re.sub(r"<\?xml[^>]*\?>", "", svg)
        svg = re.sub(r"<animateTransform.*?/>", "", svg, flags=re.S)
        svg = re.sub(r"<animate .*?/>", "", svg, flags=re.S)
        svg = re.sub(r"<svg ", '<svg class="logo" ', svg, count=1)
        return svg
    return ""


def esc(s):
    return html.escape(str(s or ""), quote=True)


def byline(authors, max_names=4):
    if not authors:
        return ""
    if len(authors) > max_names:
        return ", ".join(authors[:max_names]) + f" … (+{len(authors) - max_names})"
    return ", ".join(authors)


def short_summary(abstract, limit=260):
    """First sentences of the abstract, up to ~limit chars."""
    text = re.sub(r"\s+", " ", (abstract or "")).strip()
    text = re.sub(r"^\s*abstract[:.\s]*", "", text, flags=re.IGNORECASE)
    text = re.sub(r"\\[a-zA-Z]+", "", text).replace("$", "")
    text = re.sub(r"[{}~^]", "", text)
    text = re.sub(r"\s+", " ", text).strip()
    if not text:
        return ""
    out = ""
    for m in re.finditer(r"[^.!?]*[.!?]+(?:\s|$)", text):
        if out and len(out) + len(m.group(0)) > limit:
            break
        out += m.group(0)
        if len(out) >= limit:
            break
    return (out or text[:limit]).strip()


def paper_card(p):
    primary = (f"https://arxiv.org/abs/{p['arxiv']}" if p.get("arxiv")
               else f"https://doi.org/{p['doi']}" if p.get("doi") else "#")
    links = []
    if p.get("arxiv"):
        links.append(f'<a href="https://arxiv.org/abs/{esc(p["arxiv"])}">arXiv</a>')
    if p.get("doi"):
        links.append(f'<a href="https://doi.org/{esc(p["doi"])}">DOI</a>')

    figs = p.get("figures") or []
    if figs:
        fig = figs[0]
        title_attr = ("Figure mentioning VegasAfterglow" if fig.get("matched")
                      else "Representative figure")
        thumb = (f'<img class="thumb" src="{esc(fig["file"])}" alt="figure" '
                 f'title="{title_attr}" loading="lazy">')
    else:
        thumb = '<div class="thumb placeholder">✦</div>'

    first = p["authors"][0] if p["authors"] else ""
    etal = " et al." if len(p["authors"]) > 1 else ""
    venue = esc(p.get("journal") or "")
    year = esc(p.get("year") or "")

    return (f'<div class="pcard">'
            f'<a href="{esc(primary)}">{thumb}</a>'
            f'<div class="body">'
            f'<h4 title="{esc(p["title"])}"><a href="{esc(primary)}">{esc(p["title"])}</a></h4>'
            f'<div class="byline">{esc(first)}{etal} · {year}'
            f'{" · " + venue if venue else ""}</div>'
            f'<div class="foot">'
            f'<div class="links">{"".join(links)}</div></div>'
            f'</div></div>')


def aggregate_institutes(papers):
    inst = {}

    def specific(insts):
        keep = [i for i in insts if "system" not in (i["name"] or "").lower()]
        return keep or insts

    def add(i, p, tier):
        # institutes arrive via OpenAlex ids or ROR ids depending on the
        # attribution path; the ROR URL is the shared identity
        key = i.get("ror") or i.get("id") or i.get("name")
        rec = inst.setdefault(key, {**i, "primary": [], "co": []})
        rec[tier].append({"title": p["title"], "year": p.get("year")})

    for p in papers:
        for i in specific(p.get("institutes") or []):
            add(i, p, "primary")
        for i in p.get("co_institutes") or []:
            add(i, p, "co")
    for rec in inst.values():
        # a paper never counts an institute in both tiers
        titles = {x["title"] for x in rec["primary"]}
        rec["co"] = [x for x in rec["co"] if x["title"] not in titles]
    return sorted(inst.values(),
                  key=lambda r: (-len(r["primary"]), -len(r["co"]), r["name"] or ""))


def map_card(institutes):
    markers = []
    for i in institutes:
        if i.get("lat") is None or i.get("lng") is None:
            continue
        n1, n2 = len(i["primary"]), len(i["co"])
        shown = i["primary"][:5] + i["co"][:max(0, 5 - n1)]
        papers_li = "".join(f"<li>{esc(p['title'][:70])}</li>" for p in shown)
        counts = " · ".join(
            ([f"{n1} first/corresponding"] if n1 else [])
            + ([f"{n2} co-author"] if n2 else []))
        popup = (f'<div class="inst-pop"><div>'
                 f'<div class="nm">{esc(i["name"])}</div>'
                 f'<div class="ct">{esc(i.get("city") or "")}'
                 f'{", " if i.get("city") else ""}{esc(i.get("country") or "")}'
                 f' — {counts}</div>'
                 f'<ul>{papers_li}</ul></div></div>')
        markers.append({"lat": i["lat"], "lng": i["lng"], "n1": n1, "n2": n2,
                        "popup": popup, "name": i["name"]})
    return markers


def institute_table(institutes):
    rows = []
    for i in institutes:
        name = (f'<a href="{esc(i["homepage"])}" style="color:inherit">{esc(i["name"])}</a>'
                if i.get("homepage") else esc(i["name"]))
        loc = ", ".join(x for x in (i.get("city"), i.get("country")) if x)
        n1, n2 = len(i["primary"]), len(i["co"])
        rows.append(f"<tr><td>{name}</td><td>{esc(loc)}</td>"
                    f'<td class="cnt">{n1 or "—"}</td>'
                    f'<td class="cnt co">{n2 or "—"}</td></tr>')
    head = ('<tr><th>Institute</th><th>Location</th>'
            '<th>First / Corr.</th><th>Co-author</th></tr>')
    half = (len(rows) + 1) // 2
    cols = "".join(f'<table class="inst-table">{head}{"".join(chunk)}</table>'
                   for chunk in (rows[:half], rows[half:]) if chunk)
    return f'<div class="inst-cols">{cols}</div>'


def load_excludes():
    f = HERE / "exclude_papers.json"
    if not f.exists():
        return set()
    return {str(x).strip().lower() for x in json.loads(f.read_text())}


def is_excluded(p, excludes):
    idents = {p.get("key", ""), p.get("doi", ""), p.get("arxiv", "")}
    return bool({str(i).lower() for i in idents if i} & excludes)


def main():
    SITE.mkdir(parents=True, exist_ok=True)
    data = json.loads((DATA / "papers.json").read_text())
    excludes = load_excludes()
    papers, cited_only = [], []
    for p in data["papers"]:
        if is_excluded(p, excludes):
            cited_only.append((p, "manually excluded"))
        elif p.get("usage") == "cited":
            cited_only.append((p, "cites the methods paper without using the code"))
        else:
            papers.append(p)
    institutes = aggregate_institutes(papers)
    n_primary = sum(1 for i in institutes if i["primary"])
    n_co_only = len(institutes) - n_primary
    countries = {i["country"] for i in institutes if i.get("country")}
    latest = max((p.get("date") or "" for p in papers), default="")
    markers = map_card(institutes)

    tiles = f"""
<div class="tiles">
  <div class="tile"><div class="k">Papers</div><div class="v">{len(papers)}</div>
    <div class="s">using VegasAfterglow</div></div>
  <div class="tile"><div class="k">Institutes</div><div class="v">{n_primary}</div>
    <div class="s">first / corresponding · {n_co_only} co-author</div></div>
  <div class="tile"><div class="k">Countries</div><div class="v">{len(countries)}</div>
    <div class="s">on the map</div></div>
  <div class="tile"><div class="k">Latest paper</div><div class="v"
    style="font-size:20px">{esc(latest) or "—"}</div>
    <div class="s">updated weekly</div></div>
</div>"""

    cards = "\n".join(paper_card(p) for p in papers)


    page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Papers using VegasAfterglow</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"
  integrity="sha384-sHL9NAb7lN7rfvG5lfHpm643Xkcjzp4jFvuavGOndn6pjVqS6ny56CAt3nsEVT4H"
  crossorigin="anonymous">
<style>{CSS}</style></head><body><div class="wrap">
<div class="masthead">{inline_logo()}<div>
  <h1>Papers using VegasAfterglow</h1>
  <div class="meta">updated {esc(data["generated"])} · sources: NASA ADS full-text
  &amp; citations, OpenAlex, curated additions</div></div>
  <div class="verdict">{len(papers)} papers</div></div>
<nav>
  <a href="#map">World map <span class="n">{len(institutes)}</span></a>
  <a href="#papers">Papers <span class="n">{len(papers)}</span></a>
  <a class="archive" href="https://github.com/YihanWangAstro/VegasAfterglow">GitHub ↗</a>
</nav>
{tiles}
<section id="map"><div class="sec-head"><h2>Where VegasAfterglow is used</h2>
<span class="sec-sub">first / corresponding author institutes in blue,
co-author institutes dimmed</span></div>
<div class="card"><div id="mapdiv" style="height:440px;border-radius:10px;
  border:1px solid var(--grid)"></div>
<div class="legend"><span class="lg"><span class="sw"
  style="background:#5B8ADB"></span>first / corresponding author institute</span>
<span class="lg"><span class="sw" style="background:#8b8a84;opacity:.55"></span>
co-author institute</span></div></div>
<div class="card"><h3>Institutes</h3>{institute_table(institutes)}</div>
</section>
<section id="papers"><div class="sec-head"><h2>Papers</h2>
<span class="sec-sub">newest first · hover a thumbnail for its context ·
click through for the full paper</span></div>
<div class="papers-grid">
{cards}
</div>
</section>
<footer>
<p>Institutes are attributed from first-author and corresponding-author
affiliations as recorded by OpenAlex; co-first authorship is not encoded in
bibliographic databases and cannot be attributed automatically.</p>
<p>Discovery: NASA ADS full-text search (<code>full:"VegasAfterglow"</code>) and
citations of the methods paper, OpenAlex citations and full-text search, plus
curated additions (<code>tools/papers_page/extra_papers.json</code>). Figures are
extracted from arXiv source packages, preferring figures whose caption mentions
VegasAfterglow; each figure links back to its paper and remains © its authors.</p>
</footer>
</div>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"
  integrity="sha384-cxOPjt7s7Iz04uaHJceBmS+qpjv2JkIHNVcuOrM+YHwZOmJGBXI00mdUXEq65HTH"
  crossorigin="anonymous"></script>
<script>
const markers = {json.dumps(markers).replace("</", "<\\/")};
const map = L.map('mapdiv', {{worldCopyJump: true, scrollWheelZoom: false}});
L.tileLayer('https://{{s}}.basemaps.cartocdn.com/rastertiles/voyager/{{z}}/{{x}}/{{y}}{{r}}.png', {{
  attribution: '&copy; OpenStreetMap &copy; CARTO', maxZoom: 12 }}).addTo(map);
const group = L.featureGroup();
markers.forEach(m => {{
  const primary = m.n1 > 0;
  L.circleMarker([m.lat, m.lng], {{
    radius: primary ? 7 + 3 * Math.sqrt(m.n1) : 5 + 2 * Math.sqrt(m.n2),
    color: primary ? '#5B8ADB' : '#75746e',
    weight: primary ? 1.5 : 1,
    fillColor: primary ? '#5B8ADB' : '#75746e',
    fillOpacity: primary ? 0.5 : 0.3 }})
   .bindPopup(m.popup)
   .bindTooltip(`${{m.name}} (${{m.n1 + m.n2}})`).addTo(group);
}});
group.addTo(map);
map.setMinZoom(2);
map.fitBounds(group.getBounds().pad(0.03), {{maxZoom: 5}});
</script>
</body></html>"""

    (SITE / "index.html").write_text(page)
    print(f"wrote {SITE / 'index.html'} ({len(papers)} papers, "
          f"{len(institutes)} institutes, {len(markers)} on map)")


if __name__ == "__main__":
    main()
