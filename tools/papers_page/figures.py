"""Figures and usage classification from arXiv source packages.

For every paper in ``data/papers.json`` with an arXiv id, download the
e-print source once and derive two things:

**Usage verdict** — papers that merely cite the methods paper (e.g. in a
software round-up) must not appear on the page. With the bibliography
stripped, a paper counts as ``used`` when a figure caption mentions
VegasAfterglow, the body mentions it at least twice, or a sentence
mentioning it (by name or through the methods-paper citation key)
carries a usage verb (fit, model, use, ...). One passing name-drop with
no usage context is ``cited`` and gets excluded at build time. Papers
without retrievable sources stay ``unverified`` (kept, flagged).

**Figures** — figure environments whose caption mentions VegasAfterglow,
exported as web-ready PNGs; papers with no caption match fall back to
their first figure, labeled ``representative``.

PDF/EPS conversion prefers ``pdftoppm`` (poppler); on macOS it falls
back to ``sips``. Writes images to ``site/figures/`` and records
``figures`` + ``usage`` back into ``data/papers.json``.
"""

import json
import os
import re
import shutil
import urllib.parse
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.request
from pathlib import Path

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parent))
from authors import attributed_affiliations, affiliations_from_pdf_text  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
SITE_FIGS = HERE / "site" / "figures"

MAX_FIGS_PER_PAPER = 2
TARGET_WIDTH = 880
CAPTION_RE = re.compile(r"vegas\s*afterglow|vegasafterglow", re.IGNORECASE)
USAGE_VERB_RE = re.compile(
    r"\b(us(?:e|ed|ing)|adopt\w*|employ\w*|fit\w*|model\w*|perform\w*|"
    r"r(?:un|an)|carr\w*|generat\w*|simulat\w*|infer\w*|constrain\w*|"
    r"calculat\w*|comput\w*|estimat\w*|mcmc|posterior|package|code|framework)\b",
    re.IGNORECASE)


GROBID_URL = os.environ.get("GROBID_URL", "").rstrip("/")


def fetch_pdf(arxiv_id):
    url = f"https://arxiv.org/pdf/{arxiv_id}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "VegasAfterglow-papers-page"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            return resp.read()
    except Exception as exc:  # noqa: BLE001
        print(f"  ! {arxiv_id}: pdf download failed: {exc}")
        return None
    finally:
        time.sleep(3)


def pdf_text_affiliations(pdf_blob):
    """Affiliations from rendered first pages (fallback when the LaTeX
    source cannot be parsed and GROBID is unavailable)."""
    try:
        from pypdf import PdfReader
        import io as _io
        reader = PdfReader(_io.BytesIO(pdf_blob))
        text = "\n".join((page.extract_text() or "") for page in reader.pages[:2])
        return affiliations_from_pdf_text(text)
    except Exception as exc:  # noqa: BLE001
        print(f"  ! pdf text fallback failed: {exc}")
        return []


def grobid_header(pdf_blob):
    """Structured authors from GROBID's header model: name, corresponding
    flag, affiliation strings. Returns None when GROBID is unreachable."""
    if not GROBID_URL:
        return None
    import xml.etree.ElementTree as ET
    boundary = "----vegasafterglowboundary"
    body = ((f"--{boundary}\r\n"
             'Content-Disposition: form-data; name="input"; filename="p.pdf"\r\n'
             "Content-Type: application/pdf\r\n\r\n").encode()
            + pdf_blob + f"\r\n--{boundary}--\r\n".encode())
    req = urllib.request.Request(
        f"{GROBID_URL}/api/processHeaderDocument", data=body, method="POST",
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}",
                 "User-Agent": "VegasAfterglow-papers-page"})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            tei = resp.read()
        root = ET.fromstring(tei)
    except Exception as exc:  # noqa: BLE001
        print(f"  ! grobid failed: {exc}")
        return None
    ns = {"t": "http://www.tei-c.org/ns/1.0"}
    authors = []
    for a in root.findall(".//t:sourceDesc//t:biblStruct//t:author", ns):
        pers = a.find("t:persName", ns)
        name = " ".join(x.text for x in pers.iter() if x.text) if pers is not None else ""
        affs = []
        for aff in a.findall("t:affiliation", ns):
            orgs = [o.text.strip() for o in aff.findall(".//t:orgName", ns)
                    if o.text and o.text.strip()]
            place = [aff.findtext(".//t:settlement", "", ns),
                     aff.findtext(".//t:country", "", ns)]
            parts = orgs + [x for x in place if x]
            if orgs:
                affs.append(", ".join(dict.fromkeys(parts)))
        authors.append({"name": name, "corresp": a.get("role") == "corresp",
                        "affs": affs})
    return authors or None


def fetch_eprint(arxiv_id, dest):
    url = f"https://arxiv.org/e-print/{arxiv_id}"
    req = urllib.request.Request(url, headers={"User-Agent": "VegasAfterglow-papers-page"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        dest.write_bytes(resp.read())
    time.sleep(3)  # arXiv politeness


def tex_files(root):
    return sorted(root.rglob("*.tex"), key=lambda p: p.stat().st_size, reverse=True)


def strip_comments(tex):
    return re.sub(r"(?<!\\)%.*", "", tex)


def strip_bibliography(tex):
    tex = re.sub(r"\\begin\{thebibliography\}.*?\\end\{thebibliography\}", "", tex, flags=re.S)
    tex = re.sub(r"\\bibliography\{[^}]*\}", "", tex)
    return tex


def va_cite_keys(root):
    """Bib keys of the methods paper, from .bbl/.bib/thebibliography entries."""
    keys = set()
    entry_res = [
        (re.compile(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}(.{0,600})", re.S), 1, 2),
        (re.compile(r"@\w+\{([^,]+),(.{0,600})", re.S), 1, 2),
    ]
    for path in list(root.rglob("*.bbl")) + list(root.rglob("*.bib")) + list(root.rglob("*.tex")):
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        for rx, kg, bg in entry_res:
            for m in rx.finditer(text):
                if CAPTION_RE.search(m.group(bg)):
                    keys.add(m.group(kg).strip())
    return keys


def classify_usage(body_tex, caption_hit, cite_keys):
    """'used' vs 'cited' from the bibliography-stripped body text."""
    if caption_hit:
        return "used"
    mentions = list(CAPTION_RE.finditer(body_tex))
    sentences = []
    for m in mentions:
        lo = body_tex.rfind(".", 0, m.start()) + 1
        hi = body_tex.find(".", m.end())
        sentences.append(body_tex[lo:hi if hi != -1 else m.end() + 300])
    for key in cite_keys:
        for m in re.finditer(re.escape(key), body_tex):
            lo = body_tex.rfind(".", 0, m.start()) + 1
            hi = body_tex.find(".", m.end())
            sentences.append(body_tex[lo:hi if hi != -1 else m.end() + 300])
    if len(mentions) >= 2:
        return "used"
    if any(USAGE_VERB_RE.search(s) for s in sentences):
        return "used"
    return "cited"


def find_figures(tex):
    """Yield (caption, [graphics files]) per figure environment, in order."""
    for m in re.finditer(r"\\begin\{figure\*?\}(.*?)\\end\{figure\*?\}", tex, re.S):
        body = m.group(1)
        cap = re.search(r"\\caption(?:\[[^\]]*\])?\{", body)
        caption = ""
        if cap:
            depth, i = 1, cap.end()
            while i < len(body) and depth:
                depth += {"{": 1, "}": -1}.get(body[i], 0)
                i += 1
            caption = body[cap.end():i - 1]
        graphics = re.findall(r"\\includegraphics\*?(?:\[[^\]]*\])?\{([^}]+)\}", body)
        # plotone/plottwo (AASTeX), \fig inside \gridline, overpic environments
        graphics += re.findall(r"\\plot(?:one|two)\{([^}]+)\}", body)
        graphics += re.findall(r"\\fig\{([^}]+)\}", body)
        graphics += re.findall(r"\\begin\{overpic\}(?:\[[^\]]*\])?\{([^}]+)\}", body)
        if graphics:
            yield caption, graphics


def clean_caption(caption, limit=420):
    txt = caption
    txt = re.sub(r"\\(?:cite[pt]?|ref|label|footnote)\*?(?:\[[^\]]*\])*\{[^}]*\}", "", txt)
    txt = re.sub(r"\\(?:text|math)(?:bf|it|rm|sc|sf|tt)\{([^}]*)\}", r"\1", txt)
    txt = re.sub(r"\\(?:emph|textsc|textbf|textit|texttt|mbox)\{([^}]*)\}", r"\1", txt)
    txt = re.sub(r"\\[a-zA-Z]+\s*", " ", txt)
    txt = txt.replace("{", "").replace("}", "").replace("~", " ")
    txt = re.sub(r"\s+", " ", txt).strip()
    if len(txt) > limit:
        txt = txt[:limit].rsplit(" ", 1)[0] + " …"
    return txt


def resolve_graphic(root, name):
    cands = [name] + [f"{name}{ext}" for ext in (".pdf", ".png", ".jpg", ".jpeg", ".eps", ".ps")]
    for c in cands:
        hits = [p for p in root.rglob(Path(c).name) if p.is_file()]
        exact = [p for p in hits if str(p.relative_to(root)).endswith(c)]
        if exact or hits:
            return (exact or hits)[0]
    return None


def flatten_white(img):
    """Composite any transparency onto white — PIL's plain RGB convert
    composites onto black, which blanks matplotlib transparent-background
    figures."""
    from PIL import Image
    if img.mode in ("RGBA", "LA", "P"):
        rgba = img.convert("RGBA")
        flat = Image.new("RGB", rgba.size, (255, 255, 255))
        flat.paste(rgba, mask=rgba.split()[3])
        return flat
    return img.convert("RGB")


def to_png(src, out_png):
    suffix = src.suffix.lower()
    if suffix in (".png", ".jpg", ".jpeg"):
        try:
            from PIL import Image
            img = Image.open(src)
            img = flatten_white(img)
            if img.width > TARGET_WIDTH:
                img = img.resize((TARGET_WIDTH, int(img.height * TARGET_WIDTH / img.width)))
            img.save(out_png)
            return True
        except Exception:
            shutil.copy(src, out_png.with_suffix(suffix))
            return True
    if suffix in (".pdf", ".eps", ".ps"):
        pdf = src
        if suffix in (".eps", ".ps"):
            if not shutil.which("gs"):
                return False
            pdf = out_png.with_suffix(".tmp.pdf")
            if subprocess.run(["gs", "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=pdfwrite",
                               f"-sOutputFile={pdf}", str(src)], capture_output=True).returncode:
                return False
        if shutil.which("pdftoppm"):
            ok = subprocess.run(["pdftoppm", "-png", "-singlefile", "-r", "110",
                                 str(pdf), str(out_png.with_suffix(""))],
                                capture_output=True).returncode == 0
        elif shutil.which("gs"):
            ok = subprocess.run(["gs", "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m",
                                 "-r110", "-dFirstPage=1", "-dLastPage=1",
                                 f"-sOutputFile={out_png}", str(pdf)],
                                capture_output=True).returncode == 0
        elif sys.platform == "darwin":
            ok = subprocess.run(["sips", "-s", "format", "png", str(pdf),
                                 "--out", str(out_png)], capture_output=True).returncode == 0
        else:
            ok = False
        pdf_tmp = out_png.with_suffix(".tmp.pdf")
        if pdf_tmp.exists():
            pdf_tmp.unlink()
        if ok:
            try:
                from PIL import Image
                img = flatten_white(Image.open(out_png))
                if img.width > TARGET_WIDTH:
                    img = img.resize((TARGET_WIDTH, int(img.height * TARGET_WIDTH / img.width)))
                img.save(out_png)
            except Exception:
                pass
        return ok and out_png.exists()
    return False



def ror_resolve(aff, cache):
    """Institute record from the ROR affiliation matcher (name, geo, homepage)."""
    key = "aff:" + aff.lower()[:120]
    if key in cache:
        return cache[key]
    # Street addresses and postcodes derail the ROR matcher ("University of
    # Nevada, 4505 S. Maryland Pkwy., Las Vegas" resolves to Reno): query on
    # the digit-free segments only.
    segs = [t.strip() for t in aff.split(",")]
    query = ", ".join(t for t in segs if t and not re.search(r"\d", t)) or aff
    url = ("https://api.ror.org/v2/organizations?affiliation="
           + urllib.parse.quote(query))
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "VegasAfterglow-papers-page"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            items = json.load(resp).get("items") or []
    except Exception:  # noqa: BLE001
        return None
    finally:
        time.sleep(1)
    chosen = next((i for i in items if i.get("chosen")), None)
    # ROR's matcher can pick the wrong campus of a multi-campus system
    # ("University of Nevada, ..., Las Vegas" -> Reno); prefer a strong
    # candidate whose city actually appears in the affiliation string.
    aff_l = aff.lower()
    city_hits = []
    for it in items[:5]:
        loc = (it["organization"].get("locations") or [{}])[0]
        city = ((loc.get("geonames_details") or {}).get("name") or "").lower()
        if city and city in aff_l and it.get("score", 0) >= 0.9:
            city_hits.append(it)
    if city_hits:
        chosen = max(city_hits, key=lambda it: it.get("score", 0))
    if not chosen:
        cache[key] = None
        return None
    org = chosen["organization"]
    name = next((n["value"] for n in org.get("names", [])
                 if "ror_display" in (n.get("types") or [])), None)
    loc = (org.get("locations") or [{}])[0].get("geonames_details") or {}
    homepage = next((l["value"] for l in org.get("links", [])
                     if l.get("type") == "website"), None)
    logo = None
    if homepage:
        host = re.sub(r"^www\.", "", urllib.parse.urlparse(homepage).netloc)
        if host:
            logo = f"https://www.google.com/s2/favicons?domain={host}&sz=64"
    rec = {
        "id": org["id"].rsplit("/", 1)[-1],
        "name": name or aff[:60],
        "ror": org["id"],
        "country": loc.get("country_name"),
        "city": loc.get("name"),
        "lat": loc.get("lat"),
        "lng": loc.get("lng"),
        "homepage": homepage,
        "logo": logo,
    }
    cache[key] = rec
    return rec


def extract_for_paper(paper):
    """Returns (figures, usage_verdict, first_author_aff_strings)."""
    arxiv_id = paper.get("arxiv")
    if not arxiv_id:
        return [], "unverified", []
    workdir = Path(tempfile.mkdtemp(prefix="vae_fig_"))
    try:
        blob = workdir / "eprint"
        try:
            fetch_eprint(arxiv_id, blob)
        except Exception as exc:  # noqa: BLE001
            print(f"  ! {arxiv_id}: e-print download failed: {exc}")
            return [], "unverified", []
        src = workdir / "src"
        src.mkdir()
        try:
            with tarfile.open(blob) as tar:
                tar.extractall(src, filter="data")
        except tarfile.ReadError:
            # single-file submissions: gzipped tex, or PDF-only (no source)
            import gzip
            try:
                text = gzip.decompress(blob.read_bytes())
            except OSError:
                return [], "unverified", []
            if text.lstrip()[:5] == b"%PDF-":
                return [], "unverified", []
            (src / "main.tex").write_bytes(text)

        figures, body_parts, affs = [], [], []
        for tex_path in tex_files(src):
            tex = strip_comments(tex_path.read_text(errors="replace"))
            figures.extend(find_figures(tex))
            body_parts.append(strip_bibliography(tex))
            if not affs:
                affs = attributed_affiliations(tex)
        body = "\n".join(body_parts)
        matched = [(cap, g) for cap, g in figures if CAPTION_RE.search(cap)]
        usage = classify_usage(body, bool(matched), va_cite_keys(src))
        if not figures or usage == "cited":
            return [], usage, affs

        chosen = [(c, g, True) for c, g in matched[:MAX_FIGS_PER_PAPER]]
        if not chosen:
            chosen = [(figures[0][0], figures[0][1], False)]

        out = []
        for idx, (caption, graphics, was_match) in enumerate(chosen):
            gfile = resolve_graphic(src, graphics[0])
            if not gfile:
                continue
            out_png = SITE_FIGS / f"{paper['key']}_{idx}.png"
            if to_png(gfile, out_png):
                out.append({"file": f"figures/{out_png.name}", "caption": clean_caption(caption),
                            "matched": was_match})
        return out, usage, affs
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def main():
    SITE_FIGS.mkdir(parents=True, exist_ok=True)
    data = json.loads((DATA / "papers.json").read_text())
    cache_file = DATA / "institutes_cache.json"
    cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}
    for paper in data["papers"]:
        existing = paper.get("figures")
        have = (paper.get("usage") and existing is not None
                and all((HERE / "site" / f["file"]).exists() for f in existing))
        if have:
            continue
        print(f"analyzing: {paper.get('arxiv') or paper['key']} — {paper['title'][:60]}")
        paper["figures"], paper["usage"], affs = extract_for_paper(paper)
        co_affs = []
        pdf_blob = None
        if paper.get("arxiv") and (GROBID_URL or not affs):
            pdf_blob = fetch_pdf(paper["arxiv"])
        if pdf_blob and GROBID_URL:
            gr = grobid_header(pdf_blob)
            if gr:
                primary_idx = {0} | {i for i, a in enumerate(gr) if a["corresp"]}
                for i in sorted(primary_idx):
                    for aff in gr[i]["affs"]:
                        if aff not in affs:
                            affs.append(aff)
                co_affs = [aff for i, a in enumerate(gr) if i not in primary_idx
                           for aff in a["affs"]]
        if not affs and pdf_blob:
            affs = pdf_text_affiliations(pdf_blob)
            if affs:
                print(f"  affiliations from PDF text: {affs}")
        if "manual" in (paper.get("sources") or []):
            paper["usage"] = "used"  # curated entries are trusted
        elif paper["usage"] == "cited":
            print("  -> cites VegasAfterglow but does not use it (excluded)")
        if affs:
            # The source text is authoritative for the attribution rule
            # (first author, "equal contribution" co-firsts, corresponding
            # authors); bibliographic records don't encode co-firsts and
            # over-flag corresponding authors.
            insts, seen = [], set()
            for aff in affs:
                rec = ror_resolve(aff, cache)
                if rec and rec["id"] not in seen:
                    seen.add(rec["id"])
                    insts.append(rec)
            if insts:
                paper["institutes"] = insts
                print(f"  attributed institutes: {[i['name'] for i in insts]}")
        if co_affs:
            primary_rors = {i.get("ror") for i in paper.get("institutes") or []}
            co, seen = [], set()
            for aff in co_affs:
                rec = ror_resolve(aff, cache)
                if (rec and rec["id"] not in seen
                        and rec.get("ror") not in primary_rors and len(co) < 40):
                    seen.add(rec["id"])
                    co.append(rec)
            if co:
                paper["co_institutes"] = co
    cache_file.write_text(json.dumps(cache, indent=1))
    (DATA / "papers.json").write_text(json.dumps(data, indent=1))
    n = sum(len(p.get("figures") or []) for p in data["papers"])
    verdicts = {}
    for p in data["papers"]:
        verdicts[p.get("usage")] = verdicts.get(p.get("usage"), 0) + 1
    print(f"total figures: {n}; verdicts: {verdicts}")


if __name__ == "__main__":
    main()
