"""Discover papers that use VegasAfterglow and collect their metadata.

Discovery is a union of independent legs, so papers that use the code
without formally citing the methods paper are still found:

1. OpenAlex works citing the methods paper (no API key required)
2. OpenAlex full-text search for "VegasAfterglow"
3. NASA ADS full-text search ``full:"VegasAfterglow"`` and
   ``citations(bibcode:...)`` when ``ADS_API_TOKEN`` is set (CI)
4. Manual additions from ``extra_papers.json`` (DOIs or arXiv ids) for
   anything spotted by hand, e.g. on ResearchGate

For each paper the first-author and corresponding-author affiliations
(the attribution rule for the map) are resolved through OpenAlex
institution records, which carry ROR ids, homepages, logos, and
coordinates.

Writes ``data/papers.json`` next to this script.
"""

import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

METHODS_PAPER_OPENALEX = "W4414825884"
METHODS_PAPER_BIBCODE = "2026JHEAp..5000490W"
CONTACT = "mailto:yihan.astro@gmail.com"
UA = "VegasAfterglow-papers-page (+https://github.com/YihanWangAstro/VegasAfterglow)"

OPENALEX_FIELDS = (
    "id,title,display_name,publication_year,publication_date,type,doi,ids,"
    "primary_location,locations,authorships,abstract_inverted_index,"
    "cited_by_count"
)


def http_json(url, headers=None, retries=4):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA, **(headers or {})})
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.load(resp)
        except Exception as exc:  # noqa: BLE001 - retry then surface
            if attempt == retries - 1:
                print(f"  ! giving up on {url}: {exc}", file=sys.stderr)
                return None
            time.sleep(4 * (attempt + 1))
    return None


def openalex(path, **params):
    params.setdefault("mailto", CONTACT.removeprefix("mailto:"))
    url = f"https://api.openalex.org/{path}?{urllib.parse.urlencode(params)}"
    out = http_json(url)
    time.sleep(0.3)
    return out


def openalex_paged(filter_expr):
    results, cursor = [], "*"
    while cursor:
        page = openalex("works", filter=filter_expr, select=OPENALEX_FIELDS,
                        **{"per-page": 50, "cursor": cursor})
        if not page:
            break
        results.extend(page["results"])
        cursor = page["meta"].get("next_cursor")
        if not page["results"]:
            break
    return results


def abstract_from_inverted(inv):
    if not inv:
        return ""
    pos = {}
    for word, indices in inv.items():
        for i in indices:
            pos[i] = word
    return " ".join(pos[i] for i in sorted(pos))


def arxiv_id_of(work):
    for loc in [work.get("primary_location")] + (work.get("locations") or []):
        if not loc:
            continue
        for url in (loc.get("landing_page_url"), loc.get("pdf_url")):
            m = re.search(r"arxiv\.org/(?:abs|pdf)/([0-9]{4}\.[0-9]{4,5})", url or "")
            if m:
                return m.group(1)
    m = re.search(r"([0-9]{4}\.[0-9]{4,5})", (work.get("ids") or {}).get("arxiv") or "")
    return m.group(1) if m else None


def norm_title(title):
    return re.sub(r"[^a-z0-9]", "", (title or "").lower())


def arxiv_lookup_by_title(title):
    """Recover an arXiv id for a published article whose record lacks the
    preprint link; requires a near-exact title match to avoid false hits."""
    import difflib
    import xml.etree.ElementTree as ET
    q = urllib.parse.quote(f'ti:"{title}"')
    url = f"http://export.arxiv.org/api/query?search_query={q}&max_results=5"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=30) as resp:
            root = ET.fromstring(resp.read())
    except Exception:  # noqa: BLE001
        return None
    finally:
        time.sleep(3)  # arXiv API politeness
    ns = {"a": "http://www.w3.org/2005/Atom"}
    for entry in root.findall("a:entry", ns):
        cand = (entry.findtext("a:title", "", ns) or "").strip()
        if difflib.SequenceMatcher(None, norm_title(cand), norm_title(title)).ratio() >= 0.92:
            m = re.search(r"abs/([0-9]{4}\.[0-9]{4,5})", entry.findtext("a:id", "", ns) or "")
            if m:
                return m.group(1)
    return None


def attributed_authorships(work):
    """First author + corresponding authors — the affiliation rule for the map.

    Co-first authorship is not encoded in any bibliographic database, so it
    cannot be automated; the page footer documents the rule as implemented.
    """
    auths = work.get("authorships") or []
    first = [a for a in auths if a.get("author_position") == "first"]
    corr = [a for a in auths if a.get("is_corresponding") and a not in first]
    if len(corr) > 3:
        # Degenerate flagging: large collaborations sometimes mark every
        # author as corresponding. Fall back to the first author alone.
        corr = []
    return first + corr


class InstituteResolver:
    """OpenAlex institution records: name, ROR, homepage, logo, coordinates."""

    def __init__(self, cache_path):
        self.cache_path = cache_path
        self.cache = {}
        if cache_path.exists():
            self.cache = json.loads(cache_path.read_text())

    def resolve(self, inst_id):
        short = inst_id.rsplit("/", 1)[-1]
        if short in self.cache:
            return self.cache[short]
        rec = openalex(f"institutions/{short}")
        if not rec:
            return None
        geo = rec.get("geo") or {}
        homepage = rec.get("homepage_url")
        logo = None
        if homepage:
            from urllib.parse import urlparse
            host = re.sub(r"^www\.", "", urlparse(homepage).netloc)
            if host:
                logo = f"https://www.google.com/s2/favicons?domain={host}&sz=64"
        out = {
            "id": short,
            "name": rec.get("display_name"),
            "ror": rec.get("ror"),
            "country": geo.get("country") or (rec.get("country_code") or "").upper(),
            "city": geo.get("city"),
            "lat": geo.get("latitude"),
            "lng": geo.get("longitude"),
            "homepage": homepage,
            "logo": logo,  # homepage favicon: institutions' wikidata images are campus photos, not logos
        }
        self.cache[short] = out
        return out

    def save(self):
        self.cache_path.write_text(json.dumps(self.cache, indent=1))


def ads_query(query, token):
    papers, start = [], 0
    fields = "bibcode,title,abstract,author,aff,year,pubdate,pub,doi,identifier"
    while True:
        url = ("https://api.adsabs.harvard.edu/v1/search/query?"
               + urllib.parse.urlencode({"q": query, "fl": fields, "rows": 200, "start": start}))
        out = http_json(url, headers={"Authorization": f"Bearer {token}"})
        if not out:
            break
        resp = out["response"]
        papers.extend(resp["docs"])
        start += 200
        if start >= resp["numFound"]:
            break
    return papers


def ads_discover(token):
    """ADS legs: full-text + citations. Returns {key: doc} keyed by DOI/arXiv."""
    docs = {}
    for q in (f'full:"VegasAfterglow" -bibcode:{METHODS_PAPER_BIBCODE}',
              f"citations(bibcode:{METHODS_PAPER_BIBCODE})"):
        for d in ads_query(q, token):
            key = None
            for ident in d.get("identifier") or []:
                m = re.match(r"(?:arXiv:)?([0-9]{4}\.[0-9]{4,5})$", ident)
                if m:
                    key = ("arxiv", m.group(1))
                    break
            if not key and d.get("doi"):
                key = ("doi", d["doi"][0].lower())
            if key:
                docs.setdefault(key, d)
    return docs


def main():
    DATA.mkdir(exist_ok=True)
    resolver = InstituteResolver(DATA / "institutes_cache.json")

    # Carry expensive derived state (figures, usage verdicts) across refetches.
    prev = {}
    prev_file = DATA / "papers.json"
    if prev_file.exists():
        for old in json.loads(prev_file.read_text()).get("papers", []):
            prev[old["key"]] = old

    print("OpenAlex: citations of the methods paper ...")
    works = {w["id"]: w for w in openalex_paged(f"cites:{METHODS_PAPER_OPENALEX}")}
    for w in works.values():
        w["_sources"] = {"openalex:cites"}

    print("OpenAlex: full-text search ...")
    for w in openalex_paged("fulltext.search:VegasAfterglow"):
        entry = works.setdefault(w["id"], {**w, "_sources": set()})
        entry["_sources"] = entry.get("_sources", set()) | {"openalex:fulltext"}

    # Manual additions (e.g. found on ResearchGate): list of DOIs / arXiv ids
    extra_file = HERE / "extra_papers.json"
    if extra_file.exists():
        for ident in json.loads(extra_file.read_text()):
            ident = str(ident).strip()
            path = (f"works/doi:{ident}" if "/" in ident else f"works/arxiv:{ident}")
            w = openalex(path, select=OPENALEX_FIELDS)
            if w:
                entry = works.setdefault(w["id"], {**w, "_sources": set()})
                entry["_sources"] = entry.get("_sources", set()) | {"manual"}
            else:
                print(f"  ! extra_papers entry not found on OpenAlex: {ident}")

    # ADS legs (CI only — needs a token)
    ads_extra = {}
    token = os.environ.get("ADS_API_TOKEN")
    if token:
        print("ADS: full-text + citations ...")
        ads_docs = ads_discover(token)
        known_dois = {(w.get("doi") or "").lower().removeprefix("https://doi.org/")
                      for w in works.values()}
        known_arxiv = {arxiv_id_of(w) for w in works.values()}
        for (kind, ident), doc in ads_docs.items():
            if (kind == "doi" and ident in known_dois) or (kind == "arxiv" and ident in known_arxiv):
                continue
            w = openalex(f"works/{kind}:{ident}", select=OPENALEX_FIELDS)
            if w:
                entry = works.setdefault(w["id"], {**w, "_sources": set()})
                entry["_sources"] = entry.get("_sources", set()) | {"ads"}
            else:
                ads_extra[ident] = doc  # rendered from ADS metadata alone
    else:
        print("ADS_API_TOKEN not set — skipping the ADS discovery legs")

    works.pop(f"https://openalex.org/{METHODS_PAPER_OPENALEX}", None)

    # Deduplicate preprint/published pairs by normalized title; prefer published.
    by_title = {}
    for w in works.values():
        key = norm_title(w.get("title"))
        cur = by_title.get(key)
        if cur is None:
            by_title[key] = w
            continue
        w_pre = w.get("type") == "preprint"
        cur_pre = cur.get("type") == "preprint"
        if cur_pre and not w_pre:
            keep, drop = w, cur
        elif w_pre and not cur_pre:
            keep, drop = cur, w
        else:  # same type: keep the more recent record
            keep, drop = ((w, cur) if (w.get("publication_date") or "") > (cur.get("publication_date") or "")
                          else (cur, w))
        keep["_sources"] = keep.get("_sources", set()) | drop.get("_sources", set())
        if not arxiv_id_of(keep) and arxiv_id_of(drop):
            keep.setdefault("locations", []).extend(drop.get("locations") or [])
        by_title[key] = keep

    # Second dedup pass: near-identical titles ("14-yr-old" vs "14-year-old")
    import difflib
    merged = []
    for w in sorted(by_title.values(), key=lambda w: w.get("type") == "preprint"):
        dup = None
        for kept in merged:
            r = difflib.SequenceMatcher(None, norm_title(w.get("title")),
                                        norm_title(kept.get("title"))).ratio()
            if r >= 0.92:
                dup = kept
                break
        if dup is None:
            merged.append(w)
        else:
            dup["_sources"] = dup.get("_sources", set()) | w.get("_sources", set())
            if not arxiv_id_of(dup) and arxiv_id_of(w):
                dup.setdefault("locations", []).extend(w.get("locations") or [])
            if not dup.get("abstract_inverted_index") and w.get("abstract_inverted_index"):
                dup["abstract_inverted_index"] = w["abstract_inverted_index"]
            if not attributed_authorships(dup) or not any(
                    a.get("institutions") for a in attributed_authorships(dup)):
                if any(a.get("institutions") for a in attributed_authorships(w)):
                    dup["authorships"] = w["authorships"]

    papers = []
    for w in sorted(merged, key=lambda w: w.get("publication_date") or "", reverse=True):
        attributed = attributed_authorships(w)
        insts, seen = [], set()
        for a in attributed:
            for inst in a.get("institutions") or []:
                if inst.get("id") and inst["id"] not in seen:
                    seen.add(inst["id"])
                    rec = resolver.resolve(inst["id"])
                    if rec:
                        insts.append(rec)
        # co-author institutes: every other author's affiliations (dim tier)
        co_insts = []
        for a in w.get("authorships") or []:
            if a in attributed:
                continue
            for inst in a.get("institutions") or []:
                if inst.get("id") and inst["id"] not in seen and len(co_insts) < 40:
                    seen.add(inst["id"])
                    rec = resolver.resolve(inst["id"])
                    if rec:
                        co_insts.append(rec)
        arxiv = arxiv_id_of(w)
        if not arxiv and w.get("title"):
            arxiv = arxiv_lookup_by_title(w["title"])
        authors = [a["author"]["display_name"] for a in (w.get("authorships") or [])]
        journal = ((w.get("primary_location") or {}).get("source") or {}).get("display_name") or ""
        doi = (w.get("doi") or "").removeprefix("https://doi.org/")
        carried = prev.get(w["id"].rsplit("/", 1)[-1], {})
        papers.append({
            **{k: carried[k] for k in ("figures", "usage") if k in carried},
            "key": w["id"].rsplit("/", 1)[-1],
            "title": re.sub(r"<[^>]+>", "", w.get("title") or ""),
            "year": w.get("publication_year"),
            "date": w.get("publication_date"),
            "type": w.get("type"),
            "authors": authors,
            "journal": journal,
            "doi": doi,
            "arxiv": arxiv,
            "abstract": abstract_from_inverted(w.get("abstract_inverted_index")),
            # source-attributed institutes (figures.py) are authoritative and
            # survive refetches; fresh OpenAlex fills in only for new papers
            "institutes": carried.get("institutes") or insts,
            "co_institutes": carried.get("co_institutes") or co_insts,
            "sources": sorted(w.get("_sources", set())),
        })

    for ident, d in ads_extra.items():
        papers.append({
            "key": f"ads_{re.sub(r'[^A-Za-z0-9]', '_', d['bibcode'])}",
            "title": (d.get("title") or [""])[0],
            "year": int(d["year"]) if d.get("year") else None,
            "date": d.get("pubdate", "").replace("-00", "-01"),
            "type": "article",
            "authors": d.get("author") or [],
            "journal": d.get("pub") or "",
            "doi": (d.get("doi") or [""])[0],
            "arxiv": ident if "." in ident and "/" not in ident else None,
            "abstract": d.get("abstract") or "",
            "institutes": [],
            "sources": ["ads"],
        })

    resolver.save()
    out = {"generated": time.strftime("%Y-%m-%d"), "papers": papers}
    (DATA / "papers.json").write_text(json.dumps(out, indent=1))
    n_inst = len({i["id"] for p in papers for i in p["institutes"]})
    print(f"wrote {len(papers)} papers, {n_inst} institutes -> {DATA / 'papers.json'}")


if __name__ == "__main__":
    main()
