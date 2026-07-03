#!/usr/bin/env python3
"""Generate the reports/ index page listing the latest and per-release reports.

Usage:
    python tests/report_index.py versions.txt out/index.html

versions.txt holds one version directory name per line (e.g. "v2.0.1"),
typically harvested from the gh-pages branch. Stdlib only.
"""
import re
import sys
from pathlib import Path


def version_key(v):
    """Sort key for 'v2.10.1'-style names: numeric fields compare numerically."""
    return [int(part) if part.isdigit() else part
            for part in re.split(r"(\d+)", v)]

TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>VegasAfterglow validation reports</title>
<style>
:root {{ color-scheme: light dark;
  --ink: #0b0b0b; --ink-2: #52514e; --muted: #898781; --blue: #2a63c2; }}
@media (prefers-color-scheme: dark) {{ :root {{
  --ink: #f5f5f4; --ink-2: #c3c2b7; --muted: #8b8a84; --blue: #6d99e3; }} }}
body {{ font: 16px/1.6 -apple-system, "Segoe UI", Roboto, "Helvetica Neue",
  Arial, sans-serif; max-width: 640px; margin: 10vh auto; padding: 0 24px;
  color: var(--ink); }}
h1 {{ font-size: 22px; letter-spacing: -0.01em; }}
h2 {{ font-size: 15px; color: var(--ink-2); margin-top: 32px;
  text-transform: uppercase; letter-spacing: 0.05em; }}
a {{ color: var(--blue); }}
li {{ margin: 6px 0; }}
.note {{ color: var(--muted); font-size: 14px; margin-top: 32px; }}
.tag {{ font-size: 12px; color: var(--blue); border: 1px solid currentColor;
  border-radius: 999px; padding: 1px 9px; margin-left: 8px; }}
</style></head><body>
<h1>VegasAfterglow test &amp; validation reports</h1>
<p><a href="latest/"><strong>Latest report</strong></a> — refreshed on every
release and by the weekly scheduled validation run.</p>
<h2>Release snapshots</h2>
{body}
<p class="note">Each release report is an immutable snapshot of the full test,
physics validation, and performance results for that version. Source:
<a href="https://github.com/YihanWangAstro/VegasAfterglow">VegasAfterglow</a>.</p>
</body></html>
"""


def main():
    versions_file, out_path = sys.argv[1], sys.argv[2]
    versions = []
    p = Path(versions_file)
    if p.exists():
        versions = sorted({v.strip() for v in p.read_text(encoding="utf-8").split() if v.strip()},
                          key=version_key, reverse=True)
    if versions:
        stable = next((v for v in versions if "-" not in v), None)
        body = "<ul>" + "".join(
            f'<li><a href="{v}/">VegasAfterglow {v}</a>'
            + (' <span class="tag">latest stable</span>' if v == stable else "")
            + "</li>" for v in versions
        ) + "</ul>"
    else:
        body = '<p class="note">No release snapshots published yet.</p>'
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(TEMPLATE.format(body=body), encoding="utf-8")
    print(f"index: {out} ({len(versions)} release snapshots)")


if __name__ == "__main__":
    main()
