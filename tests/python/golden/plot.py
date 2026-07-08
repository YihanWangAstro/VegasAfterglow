"""Plot each golden baseline against a fresh recomputation.

Run from the repo root:

    python tests/python/golden/plot.py            # all cases -> golden/plots/
    python tests/python/golden/plot.py rs_thick   # selected cases only
    python tests/python/golden/plot.py --out /tmp/golden_plots

One figure per golden case: for every flux component with signal, the stored
baseline (solid) is overlaid with the current run (dashed) per frequency, and
a deviation panel below shows |current - baseline| normalized by the
test_golden.py acceptance term rtol*|baseline| + atol*peak — bins above the
1.0 line are exactly the bins the golden test would fail.
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

GOLDEN_DIR = os.path.dirname(os.path.abspath(__file__))
if GOLDEN_DIR not in sys.path:
    sys.path.insert(0, GOLDEN_DIR)
import regenerate as regen  # noqa: E402

RTOL = regen.RTOL
ATOL_PEAK = regen.ATOL_PEAK
COMPONENTS = regen.COMPONENTS
NU_LABELS = ("radio 1e9 Hz", "optical 1e14 Hz", "X-ray 1e17 Hz", "1e22 Hz")
NU_COLORS = ("#8b8a84", "#45AB8A", "#5B8ADB", "#D4693A")


def render_case(name, current=None):
    """Build the baseline-vs-current figure for one golden case.

    ``current`` takes precomputed components (e.g. from the test fixture cache);
    when omitted the case is recomputed fresh. Returns (figure, worst normalized
    deviation).
    """
    baseline = np.load(os.path.join(GOLDEN_DIR, f"{name}.npz"))
    if current is None:
        model = regen.build_model(regen.CONFIGS[name])
        current = regen.compute_components(model)

    components = [c for c in COMPONENTS if np.any(baseline[c]) or np.any(np.asarray(current[c]))]
    n = len(components)
    fig, axes = plt.subplots(
        2, n, figsize=(4.2 * n, 6.5), sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1.0]}, squeeze=False)
    t = baseline["t"]
    worst = 0.0

    for col, comp in enumerate(components):
        ref = np.asarray(baseline[comp])
        cur = np.asarray(current[comp])
        ax, axd = axes[0][col], axes[1][col]

        # normalized deviation: >1 means the golden test would fail that bin
        tol = RTOL * np.abs(ref) + ATOL_PEAK * np.max(np.abs(ref))
        dev = np.abs(cur - ref) / (tol if np.any(ref) else 1.0)
        if np.any(ref):
            worst = max(worst, float(dev.max()))

        for i, (label, color) in enumerate(zip(NU_LABELS, NU_COLORS)):
            if not np.any(ref[i]) and not np.any(cur[i]):
                continue
            # wide translucent baseline under a thin opaque current line: when
            # they agree the thin line rides inside the halo; any divergence
            # shows as the thin line leaving it (dashes would just overlap)
            ax.plot(t, ref[i], color=color, lw=4.0, alpha=0.30, label=label)
            ax.plot(t, cur[i], color=color, lw=1.2)
            axd.plot(t, dev[i], color=color, lw=1.2)

        ax.set_xscale("log")
        ax.set_yscale("log")
        # reverse-shock post-crossing decay spans hundreds of decades; anchor
        # the floor to the faintest frequency's peak so every series keeps its
        # structure while the exponential tail is trimmed
        both = np.maximum(ref, cur)
        row_peaks = both.max(axis=1)
        if np.any(row_peaks > 0):
            ax.set_ylim(row_peaks[row_peaks > 0].min() * 1e-6, row_peaks.max() * 3)
        ax.set_title(comp)
        ax.grid(alpha=0.25, which="both")
        if col == 0:
            ax.set_ylabel("flux density [baseline halo, current line]")
            axd.set_ylabel("|diff| / (rtol|ref| + atol peak)")
            ax.legend(fontsize=8, loc="best")

        axd.axhline(1.0, color="#c04040", lw=1.0, ls=":")
        axd.set_xscale("log")
        axd.set_yscale("log")
        axd.set_ylim(1e-8, 10)
        axd.grid(alpha=0.25, which="both")
        axd.set_xlabel("t [s]")

    fig.suptitle(f"{name} — worst normalized deviation {worst:.2g} "
                 f"({'PASS' if worst <= 1 else 'FAIL'} at rtol={RTOL:g}, atol={ATOL_PEAK:g}*peak)")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig, worst


def plot_case(name, out_dir):
    fig, worst = render_case(name)
    out = os.path.join(out_dir, f"{name}.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out, worst


def main():
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("names", nargs="*", help="golden case names (default: all)")
    parser.add_argument("--out", default=os.path.join(GOLDEN_DIR, "plots"),
                        help="output directory (default: tests/python/golden/plots)")
    args = parser.parse_args()

    names = args.names or sorted(regen.CONFIGS)
    unknown = [n for n in names if n not in regen.CONFIGS]
    if unknown:
        parser.error(f"unknown golden case(s): {', '.join(unknown)}")

    os.makedirs(args.out, exist_ok=True)
    for name in names:
        out, worst = plot_case(name, args.out)
        print(f"{name:28s} worst {worst:8.2g}  -> {out}")


if __name__ == "__main__":
    main()
