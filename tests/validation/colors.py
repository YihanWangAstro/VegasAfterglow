"""ANSI color helpers for terminal output."""

import sys

_USE_COLOR = hasattr(sys.stdout, "isatty") and sys.stdout.isatty()

_W = 70  # standard header width


def _c(code, text):
    return f"\033[{code}m{text}\033[0m" if _USE_COLOR else str(text)


def _bold(t):       return _c("1", t)
def _green(t):      return _c("32", t)
def _red(t):        return _c("31", t)
def _yellow(t):     return _c("33", t)
def _cyan(t):       return _c("36", t)
def _dim(t):        return _c("2", t)
def _bold_cyan(t):  return _c("1;36", t)
def _bold_green(t): return _c("1;32", t)
def _bold_red(t):   return _c("1;31", t)


def _header(title):
    """Major section header with full-width bars."""
    bar = "=" * _W
    return f"\n{_bold_cyan(bar)}\n{_bold_cyan(title)}\n{_bold_cyan(bar)}"


def _subheader(title):
    """Sub-section header with === markers."""
    return f"\n{_bold_cyan(f'=== {title} ===')}\n"


def _bar():
    """Closing bar."""
    return _bold_cyan("=" * _W)
