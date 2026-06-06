"""``Fitter.draw_fit`` implementation and its private helpers.

Kept out of the ``Fitter`` class itself so the fitter module stays focused on
data setup, sampler orchestration, and post-fit prediction. The public entry
point is ``draw_fit(fitter, ...)``, which the ``Fitter.draw_fit``
method delegates to.
"""

import logging
from typing import Literal, Optional, Tuple

import numpy as np

from .colors import (
    _auto_filter_name,
    _broad_band,
    _disambiguate_filter_colors,
    _filter_color,
)
from .style import _ERRBAR_STYLE, _MARKER_STYLE, _journal_style_ticks

logger = logging.getLogger(__name__)


def _classify_point_data(point_t, point_nu, point_flux, point_err, point_labels=None):
    """Split the Fitter's point-data lists into light curves and a spectrum count.

    Each entry in ``point_t`` etc. is one ``add_*`` call. A light curve has a
    single unique nu (multiple times); a spectrum has a single unique t
    (multiple frequencies). Mixed entries (rare) are split by unique nu.

    Returns ``(lc_list, spectrum_count)`` where ``lc_list`` is a list of
    ``(nu_hz, t_arr, f_arr, e_arr, label_or_None)`` tuples.
    """
    if point_labels is None:
        point_labels = [None] * len(point_t)
    lc_list = []
    spectrum_count = 0
    for t, nu, f, e, lbl in zip(point_t, point_nu, point_flux, point_err, point_labels):
        t = np.asarray(t)
        nu = np.asarray(nu)
        f = np.asarray(f)
        e = np.asarray(e)
        nu_unique = np.unique(nu)
        t_unique = np.unique(t)
        if nu_unique.size == 1 and t_unique.size >= 1:
            lc_list.append((float(nu_unique[0]), t, f, e, lbl))
        elif t_unique.size == 1 and nu_unique.size > 1:
            spectrum_count += 1
        else:
            for nu0 in nu_unique:
                mask = nu == nu0
                lc_list.append((float(nu0), t[mask], f[mask], e[mask], lbl))
    return lc_list, spectrum_count


def _flux_log10_range(flux, err):
    """``(log10_lo, log10_hi)`` covering the data ± error bars.

    Falls back to point values where ``f - err`` goes non-positive, and
    returns ``(0.0, 0.0)`` if nothing is finite/positive.
    """
    f = np.asarray(flux, dtype=float)
    e = np.asarray(err, dtype=float)
    lower = np.where(f - e > 0, f - e, f)
    upper = f + e
    mask = (
        (f > 0) & np.isfinite(f) & np.isfinite(lower) & np.isfinite(upper) & (lower > 0)
    )
    if not mask.any():
        return 0.0, 0.0
    return float(np.log10(lower[mask].min())), float(np.log10(upper[mask].max()))


def _stack_entries(entries, gap_decades):
    """Rank-based shifts with bottom-up overlap relaxation, re-centred on zero.

    ``entries`` is a list of ``(sort_key, kind_key, lo, hi)`` tuples where
    ``lo``/``hi`` bound the log10-flux range of one observed band on a single
    axis. Returns ``{kind_key: shift_decades}``.
    """
    entries = sorted(entries, key=lambda x: x[0])
    n = len(entries)
    if n == 0:
        return {}

    offset = (n - 1) / 2.0
    shifts = [(i - offset) * gap_decades for i in range(n)]

    # Bottom-up: push later entries up just enough so adjacent shifted ranges
    # touch instead of overlap.
    for i in range(1, n):
        prev_top = shifts[i - 1] + entries[i - 1][3]
        curr_bot = shifts[i] + entries[i][2]
        if curr_bot < prev_top:
            delta = prev_top - curr_bot
            for j in range(i, n):
                shifts[j] += delta

    median_shift = float(np.median(shifts))
    return {entries[i][1]: shifts[i] - median_shift for i in range(n)}


def _auto_shifts_by_frequency(lc_list, band_list, gap_decades):
    """Per-band log10 shifts: rank-based spacing, relaxed to avoid overlap.

    Light curves and band-integrated entries are stacked **independently**
    when both are present, because ``draw_fit`` puts them on different
    y-axes (``F_nu`` vs ``F``) with different unit systems — interleaving them
    in a single rank order would make the overlap relaxation compare
    log10-flux values that live on different scales.

    Returns ``{kind_key: shift_decades}`` where ``kind_key`` is either
    ``("lc", nu_hz)`` or ``("band", (nu_min, nu_max))``.
    """
    lc_entries = [
        (float(nu), ("lc", float(nu)), *_flux_log10_range(f, e))
        for nu, _t, f, e, _lbl in lc_list
    ]
    band_entries = [
        (
            float(np.sqrt(b.nu_min * b.nu_max)),
            ("band", (b.nu_min, b.nu_max)),
            *_flux_log10_range(b.flux, b.err),
        )
        for b in band_list
    ]
    return {
        **_stack_entries(lc_entries, gap_decades),
        **_stack_entries(band_entries, gap_decades),
    }


def _match_decade_span(ax_L, ax_R):
    """Force both log-y axes to span the same number of decades.

    Power-law slopes alpha appear with the same visual slope on both axes when
    log10(y1_max/y1_min) == log10(y2_max/y2_min).
    """
    yL0, yL1 = ax_L.get_ylim()
    yR0, yR1 = ax_R.get_ylim()
    if yL0 <= 0 or yR0 <= 0 or yL1 <= 0 or yR1 <= 0:
        return
    spanL = np.log10(yL1 / yL0)
    spanR = np.log10(yR1 / yR0)
    target = max(spanL, spanR)

    def _expand(lo, hi, target):
        mid = 0.5 * (np.log10(lo) + np.log10(hi))
        return 10 ** (mid - target / 2), 10 ** (mid + target / 2)

    ax_L.set_ylim(*_expand(yL0, yL1, target))
    ax_R.set_ylim(*_expand(yR0, yR1, target))


def _find_log_crossings(t, y, y_target):
    """Times where the log-log curve ``y(t)`` crosses ``y_target``.

    Linear interpolation in (log10 t, log10 y); ignores non-positive / non-finite
    samples. Returns an empty list when ``y_target <= 0`` or no crossing exists.
    """
    if y_target <= 0:
        return []
    t = np.asarray(t)
    y = np.asarray(y)
    good = np.isfinite(t) & np.isfinite(y) & (t > 0) & (y > 0)
    if good.sum() < 2:
        return []
    log_t = np.log10(t[good])
    log_y = np.log10(y[good])
    target = np.log10(y_target)
    diffs = log_y - target
    crossings = []
    for i in range(len(diffs) - 1):
        d0, d1 = diffs[i], diffs[i + 1]
        if d0 * d1 > 0 or d0 == d1:
            continue
        f = d0 / (d0 - d1)
        crossings.append(10 ** (log_t[i] + f * (log_t[i + 1] - log_t[i])))
    return crossings


def _draw_freq_panel(
    ax_bot, model, t_lo, t_hi, params, redshift, lc_color_map, band_color_map, band_obs
):
    """Draw nu_a / nu_m / nu_c in the observer frame on the bottom panel.

    Uses the theta-column closest to the line of sight (``theta_v``) so the
    Doppler boost reflects what the off-axis observer sees: on-axis this is the
    jet axis; off-axis it's the ring at ``theta ~ theta_v`` -- the region the
    beaming cone passes through and which dominates the observed emission.
    """
    det = model.details(t_min=t_lo, t_max=t_hi)
    theta_v = float(getattr(params, "theta_v", 0.0))
    j_obs = int(np.argmin(np.abs(np.asarray(det.theta) - theta_v)))

    def _slice(arr):
        return np.asarray(arr)[0, j_obs, :]

    nu_a, nu_m, nu_c = _slice(det.fwd.nu_a), _slice(det.fwd.nu_m), _slice(det.fwd.nu_c)
    t_obs = _slice(det.fwd.t_obs)
    boost = _slice(det.fwd.Doppler) / (1.0 + redshift)

    # (curve, label, color) triples shared by the line plot below and the
    # crossing-marker loop so the marker edge can encode which break frequency
    # is doing the crossing.
    break_curves = [
        (nu_a * boost, r"$\nu_a$", "firebrick"),
        (nu_m * boost, r"$\nu_m$", "yellowgreen"),
        (nu_c * boost, r"$\nu_c$", "royalblue"),
    ]

    for curve, label, color in break_curves:
        good = np.isfinite(curve) & (curve > 0) & np.isfinite(t_obs) & (t_obs > 0)
        if good.any():
            ax_bot.plot(t_obs[good], curve[good], label=label, color=color, lw=1)
    ax_bot.set_xscale("log")
    ax_bot.set_yscale("log")

    # Crossing markers: face = observed-band color, edge = break-curve color,
    # so a quick glance reads off "ν_m crossed the X-ray band here". Bands
    # (axhspan) get no markers -- their edges are usually wider than the curve
    # features, so square markers there just clutter the panel.
    for nu, band_color in lc_color_map.items():
        ax_bot.axhline(nu, color=band_color, linestyle="--", alpha=0.7, lw=1)
        for curve, _label, curve_color in break_curves:
            times = _find_log_crossings(t_obs, curve, nu)
            if times:
                ax_bot.plot(
                    times,
                    [nu] * len(times),
                    marker="o",
                    color=band_color,
                    **{
                        **_MARKER_STYLE,
                        "markeredgecolor": curve_color,
                        "markeredgewidth": 0.5,
                    },
                )
    for b in band_obs:
        ax_bot.axhspan(
            b.nu_min,
            b.nu_max,
            color=band_color_map[(b.nu_min, b.nu_max)],
            alpha=0.15,
        )

    ax_bot.set_ylabel(r"$\nu$ [Hz]")
    ax_bot.yaxis.set_label_coords(-0.13, 0.5)  # match upper-panel y-label x
    if ax_bot.get_legend_handles_labels()[0]:
        ax_bot.legend(loc="best", fontsize=9, ncol=3)


def draw_fit(
    fitter,
    best_params: Optional[np.ndarray] = None,
    *,
    ci: float = 0.68,
    n_samples: int = 100,
    obs_noise: Literal["frac", "abs", "none"] = "none",
    t_range: Optional[Tuple[float, float]] = None,
    n_t: int = 500,
    shifts: Optional[dict] = None,
    auto_shift_gap: float = 1.0,
    show_nu_panel: bool = True,
    resolution: Optional[Tuple[float, float, float]] = None,
    fig=None,
    axes=None,
):
    """Diagnostic plot of observation data overlaid with the fitted model.

    See :py:meth:`VegasAfterglow.fitting.fitter.Fitter.draw_fit` for the
    full parameter documentation -- this free function holds the
    implementation so the ``Fitter`` class can stay focused on fitting.
    """
    import matplotlib.pyplot as plt

    if best_params is None:
        fitter._require_fitted()
        best_params = fitter.result.top_k_params[0]
    params = fitter._to_params(best_params)

    # Classify point-data entries; spectrum entries get a warning and skip.
    lc_list, spectrum_count = _classify_point_data(
        fitter._point_t,
        fitter._point_nu,
        fitter._point_flux,
        fitter._point_err,
        fitter._point_labels,
    )
    if spectrum_count > 0:
        logger.warning(
            "draw_fit: %d spectrum entries skipped "
            "(v1 supports light-curve and band-integrated data only)",
            spectrum_count,
        )

    if not lc_list and not fitter._band_obs:
        raise ValueError(
            "draw_fit: no observation data — "
            "call add_flux_density(...) or add_flux(...) first."
        )

    has_lc = bool(lc_list)
    has_band = bool(fitter._band_obs)
    dual = has_lc and has_band

    # Figure setup: top panel grows with the number of stacked bands so
    # auto-shift always has room; hspace=0 glues the panels together.
    if axes is not None:
        ax_top, ax_bot = axes
        if fig is None:
            fig = ax_top.figure
    else:
        n_bands = len(lc_list) + len(fitter._band_obs)
        top_inches = max(3.0, 0.25 * n_bands + 2.2)
        bot_inches = 2.16
        fig_width = 4.8
        if show_nu_panel:
            fig, (ax_top, ax_bot) = plt.subplots(
                2,
                1,
                sharex=True,
                height_ratios=(top_inches, bot_inches),
                figsize=(fig_width, top_inches + bot_inches + 0.4),
                dpi=300,
                gridspec_kw={"hspace": 0.0},
            )
            plt.setp(ax_top.get_xticklabels(), visible=False)
        else:
            fig, ax_top = plt.subplots(figsize=(fig_width, top_inches + 0.5), dpi=300)
            ax_bot = None
    ax_top_R = ax_top.twinx() if dual else None

    # Observed data extent: drives both the visible x-limits and the model
    # grid (the grid extends past on each side so the curves reach the edges).
    all_t = [t_arr for _, t_arr, _, _, _ in lc_list] + [
        np.asarray(b.t) for b in fitter._band_obs
    ]
    all_t_concat = np.concatenate(all_t) if all_t else np.array([1e2, 1e7])
    t_data_min = float(all_t_concat.min())
    t_data_max = float(all_t_concat.max())
    if t_range is None:
        t_lo = t_data_min * 1e-1
        t_hi = t_data_max * 1e2
    else:
        t_lo, t_hi = float(t_range[0]), float(t_range[1])
    t_grid = np.geomspace(t_lo, t_hi, n_t)

    # Auto-shift: rank-based uniform spacing along frequency. Keeps shifts
    # bounded by (n-1)/2 * auto_shift_gap so legend labels stay readable.
    shift_map = _auto_shifts_by_frequency(lc_list, fitter._band_obs, auto_shift_gap)
    if shifts:
        shift_map.update(shifts)

    # Resolve a filter/instrument name per entry: user-supplied label wins;
    # otherwise auto-match the frequency against the canonical filter table
    # (covers the common ``nu=filter('r')`` pattern). May still be None, in
    # which case downstream code falls back to a broad-band name.
    lc_names = {nu: (lbl or _auto_filter_name(nu)) for nu, *_, lbl in lc_list}
    band_names = {
        (b.nu_min, b.nu_max): (
            b.name or _auto_filter_name(float(np.sqrt(b.nu_min * b.nu_max)))
        )
        for b in fitter._band_obs
    }

    # Colors via the canonical-filter registry; same-color filters at
    # different frequencies (e.g. 'r' and 'VT_R') get a lightness spread so
    # they remain visually distinct.
    color_items = [
        (("lc", nu), _filter_color(lc_names[nu], nu), nu) for nu, *_ in lc_list
    ] + [
        (
            ("band", (b.nu_min, b.nu_max)),
            _filter_color(
                band_names[(b.nu_min, b.nu_max)], float(np.sqrt(b.nu_min * b.nu_max))
            ),
            float(np.sqrt(b.nu_min * b.nu_max)),
        )
        for b in fitter._band_obs
    ]
    colors = _disambiguate_filter_colors(color_items)
    lc_color_map = {nu: colors[("lc", nu)] for nu, *_ in lc_list}
    band_color_map = {
        (b.nu_min, b.nu_max): colors[("band", (b.nu_min, b.nu_max))]
        for b in fitter._band_obs
    }

    def _shifted_label(base, shift):
        return f"{base} $\\times 10^{{{shift:.0f}}}$" if abs(shift) > 1e-9 else base

    if obs_noise not in ("frac", "abs", "none"):
        raise ValueError(
            f"draw_fit: obs_noise must be one of 'frac', 'abs', 'none'; got {obs_noise!r}"
        )

    def _broaden(central, lower, upper, sigma):
        """Posterior-predictive widening: combine model spread and observation
        noise in quadrature so the band represents *"where would a new
        observation fall?"* rather than *"where is the underlying model?"*.
        ``sigma`` may be a scalar or an array broadcastable with ``central``."""
        half_lo = central - lower
        half_hi = upper - central
        return (
            central - np.sqrt(half_lo * half_lo + sigma * sigma),
            central + np.sqrt(half_hi * half_hi + sigma * sigma),
        )

    def _obs_sigma(central, flux, err):
        """Per-mode observation noise for the broadening step. Returns ``None``
        when ``obs_noise='none'``; ``'abs'`` returns the per-band median
        ``err`` as a scalar; ``'frac'`` returns ``central * median(err/flux)``
        evaluated only over positive flux samples (median is robust to outliers
        and skips marginal detections that would inflate the ratio)."""
        if obs_noise == "none":
            return None
        err = np.asarray(err)
        if obs_noise == "abs":
            return float(np.median(err))
        flux = np.asarray(flux)
        positive = flux > 0
        if not positive.any():
            # No usable points for the fractional ratio; fall back to abs σ.
            return float(np.median(err))
        return central * float(np.median(err[positive] / flux[positive]))

    # When posterior samples are available and the user did not opt out
    # (n_samples=0), show the posterior median trajectory with a shaded
    # ``ci`` credible band. Otherwise fall back to the MAP curve.
    use_band = n_samples > 0 and fitter._has_posterior_samples()

    with fitter._override_resolution(resolution):
        # ---- Light curves on the left axis ----
        if lc_list:
            lc_nus = np.array([nu for nu, *_ in lc_list], dtype=float)
            if use_band:
                cb = fitter.flux_density_credible(
                    t_grid, lc_nus, ci=ci, n_samples=n_samples
                )
                lc_central, lc_lower, lc_upper = cb.median, cb.lower, cb.upper
            else:
                lc_central = np.asarray(
                    fitter.flux_density_grid(best_params, t_grid, lc_nus).total
                )
                lc_lower = lc_upper = None

            for i, (nu, t_data, f_data, e_data, _user_label) in enumerate(lc_list):
                shift = float(shift_map.get(("lc", nu), 0.0))
                scale = 10.0**shift
                color = lc_color_map[nu]
                label = _shifted_label(lc_names[nu] or _broad_band(nu), shift)
                ax_top.errorbar(
                    t_data,
                    f_data * scale,
                    yerr=e_data * scale,
                    color=color,
                    ecolor=color,
                    **_ERRBAR_STYLE,
                )
                ax_top.plot(
                    t_grid, lc_central[i] * scale, "-", color=color, label=label, lw=1
                )
                if lc_lower is not None:
                    band_lo, band_hi = lc_lower[i], lc_upper[i]
                    sigma = _obs_sigma(lc_central[i], f_data, e_data)
                    if sigma is not None:
                        band_lo, band_hi = _broaden(
                            lc_central[i], band_lo, band_hi, sigma
                        )
                    ax_top.fill_between(
                        t_grid,
                        band_lo * scale,
                        band_hi * scale,
                        color=color,
                        alpha=0.2,
                        lw=0,
                    )

        # ---- Band-integrated on the right axis (or the left, if no LC data) ----
        band_ax = ax_top_R if dual else ax_top
        for b in fitter._band_obs:
            key = (b.nu_min, b.nu_max)
            shift = float(shift_map.get(("band", key), 0.0))
            scale = 10.0**shift
            color = band_color_map[key]
            if use_band:
                cb = fitter.flux_credible(
                    t_grid, key, ci=ci, n_samples=n_samples, num_points=b.num_points
                )
                central, lower, upper = cb.median, cb.lower, cb.upper
            else:
                central = np.asarray(
                    fitter.flux(best_params, t_grid, key, num_points=b.num_points).total
                )
                lower = upper = None
            nu_center = float(np.sqrt(b.nu_min * b.nu_max))
            label = _shifted_label(band_names[key] or _broad_band(nu_center), shift)
            band_ax.errorbar(
                np.asarray(b.t),
                np.asarray(b.flux) * scale,
                yerr=np.asarray(b.err) * scale,
                color=color,
                ecolor=color,
                **_ERRBAR_STYLE,
            )
            band_ax.plot(t_grid, central * scale, "--", color=color, label=label, lw=1)
            if lower is not None:
                sigma = _obs_sigma(central, b.flux, b.err)
                if sigma is not None:
                    lower, upper = _broaden(central, lower, upper, sigma)
                band_ax.fill_between(
                    t_grid,
                    lower * scale,
                    upper * scale,
                    color=color,
                    alpha=0.2,
                    lw=0,
                )

    ax_top.set_xscale("log")
    ax_top.set_yscale("log")
    # Pin y-label x-coords so long rotated labels can't extend across the
    # hspace=0 panel boundary into the bottom panel's label region.
    if has_lc:
        ax_top.set_ylabel(r"$F_\nu$ [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$]")
        ax_top.yaxis.set_label_coords(-0.13, 0.5)
    if dual:
        ax_top_R.set_yscale("log")
        ax_top_R.set_ylabel(r"$F$ [erg cm$^{-2}$ s$^{-1}$]")
        ax_top_R.yaxis.set_label_coords(1.13, 0.5)
        _match_decade_span(ax_top, ax_top_R)
    elif has_band and not has_lc:
        ax_top.set_ylabel(r"$F$ [erg cm$^{-2}$ s$^{-1}$]")
        ax_top.yaxis.set_label_coords(-0.13, 0.5)

    if t_range is None:
        ax_top.set_xlim(t_data_min * 1e-1, t_data_max * 1e2)
    else:
        ax_top.set_xlim(t_lo, t_hi)

    # Combined legend (gather handles from twin axes when present).
    h1, l1 = ax_top.get_legend_handles_labels()
    h2, l2 = ax_top_R.get_legend_handles_labels() if ax_top_R else ([], [])
    ax_top.legend(h1 + h2, l1 + l2, loc="best", fontsize=9)

    # --- bottom panel ----------------------------------------------------
    if ax_bot is not None:
        _draw_freq_panel(
            ax_bot,
            fitter._build_model(params),
            t_lo,
            t_hi,
            params,
            fitter.z,
            lc_color_map,
            band_color_map,
            fitter._band_obs,
        )
        ax_bot.set_xscale("log")
        ax_bot.set_xlabel(r"$t_{\rm obs}$ [s]")
    else:
        ax_top.set_xlabel(r"$t_{\rm obs}$ [s]")

    # Journal-paper polish: inward major + minor ticks on all four sides.
    _journal_style_ticks(ax_top, ax_top_R, ax_bot)

    return fig, (ax_top, ax_bot)
