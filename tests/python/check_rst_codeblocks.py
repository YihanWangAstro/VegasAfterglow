#!/usr/bin/env python3
"""Check that Python code blocks in RST files have valid syntax.

Extracts ``.. code-block:: python`` blocks, normalises common documentation
conventions (trailing ``...`` in function calls, bodyless function stubs),
then runs ``ast.parse()`` on each block.

Blocks whose first non-blank line is ``# noqa: syntax`` are skipped.

Exit code 0 if all blocks parse; non-zero otherwise.
"""

import ast
import re
import sys
from pathlib import Path

_DIRECTIVE_RE = re.compile(r"^(\s*)\.\. code-block:: python\s*$")

# Trailing ``...`` as a positional arg after keyword args: ``, ...)``
_TRAILING_ELLIPSIS_RE = re.compile(r",\s*\.\.\.\s*\)")

# Mid-call ``...`` as placeholder between keyword args: ``, ...,``
_MID_ELLIPSIS_RE = re.compile(r",\s*\.\.\.\s*,")

# Standalone ``...`` as sole argument: ``func(...)``
_SOLE_ELLIPSIS_RE = re.compile(r"\(\s*\.\.\.\s*\)")


def extract_blocks(text: str):
    """Yield (line_number, dedented_code) for each Python code block.

    Handles code blocks nested inside RST list items by tracking the
    directive's indent level and stopping when indentation drops back.
    """
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        m = _DIRECTIVE_RE.match(lines[i])
        if not m:
            i += 1
            continue

        directive_indent = len(m.group(1))
        i += 1

        # Skip blank lines between directive and body
        while i < len(lines) and not lines[i].strip():
            i += 1

        if i >= len(lines):
            break

        # Determine content indent from first non-blank line
        first_content = lines[i]
        content_indent = len(first_content) - len(first_content.lstrip())
        if content_indent <= directive_indent:
            # No actual content after directive
            continue

        # Collect body lines: keep going while lines are blank or indented
        # deeper than the directive
        body_lines = []
        code_start_line = i + 1  # 1-based line number
        while i < len(lines):
            line = lines[i]
            if not line.strip():
                # Blank line — could be internal to the block or the end.
                # Peek ahead: if next non-blank line is still indented enough,
                # this blank line is part of the block.
                j = i + 1
                while j < len(lines) and not lines[j].strip():
                    j += 1
                if j < len(lines):
                    next_indent = len(lines[j]) - len(lines[j].lstrip())
                    if next_indent >= content_indent:
                        body_lines.append("")
                        i += 1
                        continue
                # End of block
                break
            else:
                line_indent = len(line) - len(line.lstrip())
                if line_indent < content_indent:
                    # Indentation dropped — block ended
                    break
                body_lines.append(line[content_indent:])
                i += 1

        code = "\n".join(body_lines).rstrip()
        if code:
            yield code_start_line, code


def normalise(code: str) -> str:
    """Apply documentation-convention normalisation so ``ast.parse`` succeeds
    on idiomatic doc snippets."""
    # Replace trailing ``, ...)`` → ``)``
    code = _TRAILING_ELLIPSIS_RE.sub(")", code)
    # Replace mid-call ``, ...,`` → ``,`` (drop placeholder)
    code = _MID_ELLIPSIS_RE.sub(",", code)
    # Replace ``(...)`` → ``()``
    code = _SOLE_ELLIPSIS_RE.sub("()", code)

    # If the block ends with a function/class header that has no body, append
    # ``...`` (valid Python stub).
    lines = code.split("\n")
    for line in reversed(lines):
        if line.strip():
            last_real = line.strip()
            header_indent = len(line) - len(line.lstrip())
            if last_real.endswith(":") and (
                last_real.startswith(("def ", "class ", "async def "))
                or last_real == ":"
            ):
                code += "\n" + " " * (header_indent + 4) + "..."
            break

    return code


def check_file(path: Path) -> list[tuple[int, str]]:
    """Return list of (line, error_message) for blocks that fail parsing."""
    text = path.read_text(encoding="utf-8")
    errors = []
    for line_no, code in extract_blocks(text):
        # Skip blocks with noqa directive
        for line in code.split("\n"):
            if line.strip():
                if line.strip() == "# noqa: syntax":
                    break
                # First non-blank line is not noqa — proceed
                norm = normalise(code)
                try:
                    ast.parse(norm)
                except SyntaxError as exc:
                    errors.append((line_no, str(exc)))
                break
    return errors


def main() -> int:
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} FILE [FILE ...]", file=sys.stderr)
        return 2

    failed = False
    for arg in sys.argv[1:]:
        path = Path(arg)
        if not path.is_file():
            continue
        for line_no, msg in check_file(path):
            print(f"{path}:{line_no}: SyntaxError: {msg}")
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
