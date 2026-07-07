"""Attributed-affiliation extraction from LaTeX sources.

Implements the attribution rule for the map: the affiliations of the
first author, every co-first author (marked "equal contribution"), and
every corresponding author (``\\correspondingauthor``, ``\\email``, or
``\\thanks{...e-mail...}``). Each selected author may carry several
affiliations.

Handles the two dominant layouts in astro papers:

* AASTeX — one ``\\author{}`` per author followed by its
  ``\\affiliation{}``/``\\affil{}`` lines.
* MNRAS/A&A — a single ``\\author{}`` block with ``$^{1,2}$``
  superscripts and a numbered affiliation list.
"""

import re

EQUAL_RE = re.compile(r"equal(?:ly)?\s+contrib|contributed\s+equally|co[- ]first",
                      re.IGNORECASE)
MAIL_RE = re.compile(r"e-?mail|@", re.IGNORECASE)


def _macro_defs(tex):
    r"""User \newcommand/\def bodies, for \affiliation{\myinst} patterns."""
    defs = {}
    for m in re.finditer(r"\\(?:newcommand|renewcommand|def)\s*\{?\\(\w+)\}?\s*\{", tex):
        depth, i = 1, m.end()
        while i < len(tex) and depth:
            depth += {"{": 1, "}": -1}.get(tex[i], 0)
            i += 1
        body = tex[m.end():i - 1]
        if 0 < len(body) < 400:
            defs[m.group(1)] = body
    return defs


def _expand(s, defs):
    for _ in range(2):
        new = re.sub(r"\\(\w+)",
                     lambda m: defs.get(m.group(1), m.group(0)), s)
        if new == s:
            break
        s = new
    return s


def _clean(aff):
    aff = re.sub(r"\\[a-zA-Z]+\s*", " ", aff)
    aff = aff.replace("{", "").replace("}", "").replace("~", " ")
    aff = re.sub(r"\s+", " ", aff).strip(" ;,")
    return aff


def _norm_name(name):
    name = re.sub(r"\\[a-zA-Z]+", "", name)
    return re.sub(r"[^a-z]", "", name.lower())


def _aastex(tex, defs):
    """One \\author per person, each followed by its affiliations."""
    marks = list(re.finditer(r"\\author\*?(?:\[[^\]]*\])?\{", tex))
    if not marks:
        return None
    authors = []
    for k, m in enumerate(marks):
        # author name: balanced-brace read
        depth, i = 1, m.end()
        while i < len(tex) and depth:
            depth += {"{": 1, "}": -1}.get(tex[i], 0)
            i += 1
        name = tex[m.end():i - 1]
        end = marks[k + 1].start() if k + 1 < len(marks) else min(i + 4000, len(tex))
        seg = tex[i:end]
        affs = [_clean(_expand(a, defs)) for a in
                re.findall(r"\\affil(?:iation)?\{([^}]+)\}", seg)]
        notes = " ".join(re.findall(
            r"\\(?:altaffiliation|thanks|footnote)\{([^}]*)\}", seg))
        has_email = bool(re.search(r"\\email\{", seg))
        authors.append({"name": name, "affs": [a for a in affs if len(a) > 5],
                        "equal": bool(EQUAL_RE.search(notes)),
                        "email": has_email or bool(MAIL_RE.search(notes))})
    if not any(a["affs"] for a in authors):
        return None

    corresponding_names = " ".join(
        re.findall(r"\\correspondingauthor\{([^}]*)\}", tex))
    selected = set()
    for idx, a in enumerate(authors):
        if idx == 0 or a["equal"]:
            selected.add(idx)
        elif corresponding_names and _norm_name(a["name"]) and \
                _norm_name(a["name"]) in _norm_name(corresponding_names):
            selected.add(idx)
    if not corresponding_names:
        # \email marks the corresponding author unless everyone carries one
        emailed = [i for i, a in enumerate(authors) if a["email"]]
        if 0 < len(emailed) <= 3 and len(emailed) < len(authors):
            selected.update(emailed)
    if sum(1 for i in selected if authors[i]["equal"]) > 4:
        # runaway equal-contribution match: fall back to first + corresponding
        selected = {i for i in selected if i == 0 or not authors[i]["equal"]}

    affs = []
    for i in sorted(selected):
        for a in authors[i]["affs"][:4]:
            if a not in affs:
                affs.append(a)
    return affs or None


def _superscript_block(tex, defs):
    """MNRAS/A&A: one author block, $^{1,2}$ superscripts, numbered list."""
    m = re.search(r"\\author\*?(?:\[[^\]]*\])?\{", tex)
    if not m:
        return None
    depth, i = 1, m.end()
    while i < len(tex) and depth:
        depth += {"{": 1, "}": -1}.get(tex[i], 0)
        i += 1
    block = tex[m.end():i - 1]

    # numbered affiliations: $^{1}$Dept..., ^1 Dept..., \textsuperscript{1}Dept
    aff_map = {}
    for mm in re.finditer(r"(?:\$\^\{?(\d+)\}?\$|\\textsuperscript\{(\d+)\})"
                          r"\s*([^\n\\$]{10,240})", block):
        n = mm.group(1) or mm.group(2)
        text = _clean(_expand(mm.group(3), defs))
        # affiliation definitions contain location-ish text, not names
        if len(text) > 12 and n not in aff_map:
            aff_map[n] = text
    if not aff_map:
        return None

    # author entries: Name$^{1,2}$ possibly with \thanks{E-mail...}
    entries = list(re.finditer(
        r"([A-Z][^,$\\]{2,60}?)\s*\$\^\{?([\d,\s]+)\}?\$(\s*\\thanks\{[^}]*\})?",
        block))
    if not entries:
        return None
    wanted_nums = []

    def add_nums(s):
        for n in re.split(r"[,\s]+", s.strip()):
            if n and n not in wanted_nums:
                wanted_nums.append(n)

    add_nums(entries[0].group(2))  # first author
    for e in entries[1:]:
        thanks = e.group(3) or ""
        if MAIL_RE.search(thanks) or EQUAL_RE.search(thanks):
            add_nums(e.group(2))
    # corresponding marked on the first author is already included
    affs = [aff_map[n] for n in wanted_nums if n in aff_map]
    return affs or None


def attributed_affiliations(tex):
    """Affiliation strings for first / co-first / corresponding authors."""
    defs = _macro_defs(tex)
    return _aastex(tex, defs) or _superscript_block(tex, defs) or []


def affiliations_from_pdf_text(text):
    """Fallback: numbered affiliations from a rendered first page.

    Takes the affiliation(s) numbered 1 — the first author's — from the
    header text before the abstract. Used when the LaTeX source could not
    be parsed; layout-agnostic but blind to co-first/corresponding marks.
    """
    head = re.split(r"\bA\s*B\s*S\s*T\s*R\s*A\s*C\s*T\b|\bAbstract\b", text)[0]
    inst_re = re.compile(
        r"Universit|Institut|Observator|Department|Center|Centre|"
        r"Laborator|College|School|Academy|Facility|Space")
    affs = []
    for m in re.finditer(r"(?m)^\s*(\d{1,2})\s*([A-Z][^\n]{15,220})", head):
        n, t = int(m.group(1)), m.group(2).strip().rstrip(",;")
        if inst_re.search(t):
            affs.append((n, t))
    affs.sort(key=lambda x: x[0])
    ones = [t for n, t in affs if n == 1]
    if ones:
        return ones[:2]
    return [t for _, t in affs[:1]]
