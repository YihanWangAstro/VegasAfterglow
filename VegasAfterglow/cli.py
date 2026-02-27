"""Command-line light curve generator for VegasAfterglow.

Usage::

    vegasgen                          # defaults: tophat, ISM, on-axis
    vegasgen --jet gaussian --z 0.5   # override params
    vegasgen --nu R J --plot          # filter names + quick plot
    vegasgen -o lc.csv                # save to file
"""

import argparse
import sys

import numpy as np

from .units import _NAMED_BANDS, _c_A
from .units import keV as _keV


def _lumi_dist_from_z(z):
    """Approximate luminosity distance from redshift (flat LCDM, H0=67.4)."""
    c_km_s = 299792.458
    H0 = 67.4  # km/s/Mpc
    Mpc_cm = 3.0856775814913673e24
    # Hubble-law approximation: d_L = c*z/H0 * (1+z), adequate for quick use
    return c_km_s * z / H0 * (1 + z) * Mpc_cm


def _parse_frequency(value):
    """Parse a frequency value: try as float (Hz), then as filter name."""
    try:
        return float(value)
    except ValueError:
        pass
    from . import units

    try:
        return units.filter(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"'{value}' is not a valid frequency (Hz) or filter name"
        )


def _parse_nu_entry(value):
    """Parse a --nu entry as point frequency, named band, or [nu_min,nu_max].

    Returns either a float (point frequency in Hz) or a tuple
    (nu_min, nu_max, label) for a frequency band.
    """
    # Named band
    if value in _NAMED_BANDS:
        nu_min, nu_max = _NAMED_BANDS[value]
        return (nu_min, nu_max, value)
    # Bracket band: [nu_min,nu_max]
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1]
        parts = inner.split(",")
        if len(parts) != 2:
            raise argparse.ArgumentTypeError(
                f"'{value}' — band format must be [nu_min,nu_max]"
            )
        try:
            nu_min, nu_max = float(parts[0]), float(parts[1])
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"'{value}' — band edges must be numeric Hz values"
            )
        return (nu_min, nu_max, value)
    # Point frequency or filter name
    return _parse_frequency(value)


def _format_nu_label(nu):
    """Format a frequency for CSV/JSON column headers."""
    return _format_energy(nu)


def _broad_band(nu):
    """Return the broad-band name for a frequency."""
    from .units import eV as _eV_Hz

    E_eV = nu / _eV_Hz

    if nu < 1e12:  # < 1 THz
        return "Radio"
    if nu < 4e14:  # 1 THz - 400 THz  (λ > 750 nm)
        return "IR"
    if nu < 7.5e14:  # 400 - 750 THz  (400-750 nm)
        return "Optical"
    if E_eV < 100:  # 750 THz - 0.1 keV  (λ > 12 nm)
        return "UV"
    if E_eV < 1e5:  # 0.1 - 100 keV
        return "X-ray"
    if E_eV < 1e9:  # 100 keV - 1 GeV
        return r"$\gamma$-ray"
    if E_eV < 1e12:  # 1 GeV - 1 TeV
        return "GeV"
    return "TeV"


def _format_nu_latex(nu, label=None):
    """Format a frequency as a LaTeX label with broad-band prefix for plots."""
    from . import units

    band = _broad_band(nu)

    # Show filter name only if the user actually typed a filter name
    if label is not None:
        try:
            float(label)
        except ValueError:
            if label in units._VEGA_FILTERS:
                return rf"{band} (${label}$-band)"
            if label in units._ST_FILTERS:
                return rf"{band} ({label})"
            if label in units._SURVEY_FILTERS:
                return rf"{band} ({label})"

    return rf"{band} ({_format_energy(nu, latex=True)})"


_FLUX_SCALES = {"mJy": 1e-26, "Jy": 1e-23, "uJy": 1e-29, "cgs": 1.0}
_TIME_SCALES = {"s": 1.0, "day": 86400.0, "hr": 3600.0, "min": 60.0}


def parse_args(argv=None):
    from . import __version__

    p = argparse.ArgumentParser(
        prog="vegasgen",
        description="Quick GRB afterglow light curve generator (VegasAfterglow).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  vegasgen                                  # default tophat/ISM/on-axis\n"
            "  vegasgen --jet gaussian --z 0.5           # Gaussian jet at z=0.5\n"
            "  vegasgen --nu R J F606W --plot            # filter-name frequencies + plot\n"
            "  vegasgen --medium wind --A_star 0.1 --k_m 1.5 -o lc.csv\n"
            "  vegasgen --duration 100 --rvs --plot      # thick shell + reverse shock\n"
            "  vegasgen --rvs --plot                     # forward + reverse shock\n"
            "  vegasgen --rvs --rvs_eps_B 0.1 --plot     # RS with different microphysics\n"
            "  vegasgen --ssc --rvs --rvs_ssc --plot     # SSC on both shocks\n"
        ),
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s (VegasAfterglow {__version__})",
    )

    # -- Jet ----------------------------------------------------------------
    jet = p.add_argument_group("jet")
    jet.add_argument(
        "--jet", choices=["tophat", "gaussian", "powerlaw"], default="tophat"
    )
    jet.add_argument(
        "--theta_c", type=float, default=0.1, help="half-opening angle [rad]"
    )
    jet.add_argument(
        "--E_iso", type=float, default=1e52, help="isotropic-equivalent energy [erg]"
    )
    jet.add_argument("--Gamma0", type=float, default=300, help="initial Lorentz factor")
    jet.add_argument(
        "--k_e", type=float, default=2, help="energy power-law index (powerlaw jet)"
    )
    jet.add_argument(
        "--k_g", type=float, default=2, help="Gamma power-law index (powerlaw jet)"
    )
    jet.add_argument(
        "--spreading", action="store_true", help="enable lateral spreading"
    )
    jet.add_argument(
        "--duration", type=float, default=1, help="engine activity duration [s]"
    )

    # -- Medium -------------------------------------------------------------
    med = p.add_argument_group("medium")
    med.add_argument("--medium", choices=["ism", "wind"], default="ism")
    med.add_argument(
        "--n_ism", type=float, default=1.0, help="ISM number density [cm^-3]"
    )
    med.add_argument("--A_star", type=float, default=0.1, help="wind parameter A*")
    med.add_argument(
        "--k_m", type=float, default=2, help="wind density power-law index"
    )

    # -- Observer -----------------------------------------------------------
    obs = p.add_argument_group("observer")
    obs.add_argument("--z", type=float, default=0.01, help="redshift")
    obs.add_argument(
        "--lumi_dist",
        type=float,
        default=None,
        help="luminosity distance [cm] (auto from z if omitted)",
    )
    obs.add_argument("--theta_obs", type=float, default=0, help="viewing angle [rad]")

    # -- Radiation ----------------------------------------------------------
    rad = p.add_argument_group("radiation")
    rad.add_argument(
        "--eps_e", type=float, default=0.1, help="electron energy fraction"
    )
    rad.add_argument(
        "--eps_B", type=float, default=1e-3, help="magnetic energy fraction"
    )
    rad.add_argument("--p", type=float, default=2.3, help="electron spectral index")
    rad.add_argument(
        "--xi_e", type=float, default=1, help="electron acceleration fraction"
    )
    rad.add_argument(
        "--ssc", action="store_true", help="enable synchrotron self-Compton"
    )
    rad.add_argument(
        "--kn", action="store_true", help="enable Klein-Nishina corrections"
    )
    rad.add_argument("--cmb_cooling", action="store_true", help="enable CMB IC cooling")

    # -- Reverse shock ------------------------------------------------------
    rvs = p.add_argument_group("reverse shock")
    rvs.add_argument("--rvs", action="store_true", help="enable reverse shock")
    rvs.add_argument(
        "--rvs_eps_e",
        type=float,
        default=None,
        help="RS electron energy fraction (default: same as --eps_e)",
    )
    rvs.add_argument(
        "--rvs_eps_B",
        type=float,
        default=None,
        help="RS magnetic energy fraction (default: same as --eps_B)",
    )
    rvs.add_argument(
        "--rvs_p",
        type=float,
        default=None,
        help="RS electron spectral index (default: same as --p)",
    )
    rvs.add_argument(
        "--rvs_xi_e",
        type=float,
        default=None,
        help="RS electron acceleration fraction (default: same as --xi_e)",
    )
    rvs.add_argument(
        "--rvs_ssc", action="store_true", help="RS synchrotron self-Compton"
    )
    rvs.add_argument(
        "--rvs_kn", action="store_true", help="RS Klein-Nishina corrections"
    )

    # -- Frequencies --------------------------------------------------------
    freq = p.add_argument_group("frequencies")
    freq.add_argument(
        "--nu",
        nargs="+",
        default=["1e9", "5e14", "1e18"],
        help=("frequencies, filters, or bands " "(e.g. 1e9 R F606W XRT [7e16,2e18])"),
    )
    freq.add_argument(
        "--num_nu_band",
        type=int,
        default=15,
        help="frequency sampling points per band integration (default: 15)",
    )

    # -- Time grid ----------------------------------------------------------
    tg = p.add_argument_group("time grid")
    tg.add_argument("--t_min", type=float, default=1, help="start time [s]")
    tg.add_argument("--t_max", type=float, default=1e8, help="end time [s]")
    tg.add_argument("--num_t", type=int, default=100, help="number of time points")

    # -- Resolution ---------------------------------------------------------
    res = p.add_argument_group("resolution")
    res.add_argument(
        "--res",
        nargs=3,
        type=float,
        default=[0.15, 1.0, 10],
        metavar=("PHI", "THETA", "T"),
        help="resolution (phi_ppd, theta_ppd, t_ppd)",
    )

    # -- Output -------------------------------------------------------------
    out = p.add_argument_group("output")
    out.add_argument(
        "-o", "--output", default=None, help="output file (default: stdout)"
    )
    out.add_argument("--format", choices=["csv", "json"], default="csv")
    out.add_argument("--flux_unit", choices=list(_FLUX_SCALES), default="mJy")
    out.add_argument("--time_unit", choices=list(_TIME_SCALES), default="s")
    out.add_argument("--plot", action="store_true", help="show/save a quick plot")
    out.add_argument(
        "--font", default=None, help="plot font family (e.g. 'Helvetica', 'Palatino')"
    )

    return p.parse_args(argv)


def build_jet(args):
    from . import GaussianJet, PowerLawJet, TophatJet

    common = dict(
        theta_c=args.theta_c,
        E_iso=args.E_iso,
        Gamma0=args.Gamma0,
        spreading=args.spreading,
        duration=args.duration,
    )
    if args.jet == "tophat":
        return TophatJet(**common)
    if args.jet == "gaussian":
        return GaussianJet(**common)
    if args.jet == "powerlaw":
        return PowerLawJet(**common, k_e=args.k_e, k_g=args.k_g)


def build_medium(args):
    from . import ISM, Wind

    if args.medium == "ism":
        return ISM(n_ism=args.n_ism)
    return Wind(A_star=args.A_star, k_m=args.k_m)


def build_observer(args):
    from . import Observer

    d_L = args.lumi_dist if args.lumi_dist is not None else _lumi_dist_from_z(args.z)
    return Observer(lumi_dist=d_L, z=args.z, theta_obs=args.theta_obs)


def build_radiation(args):
    from . import Radiation

    fwd_rad = Radiation(
        eps_e=args.eps_e,
        eps_B=args.eps_B,
        p=args.p,
        xi_e=args.xi_e,
        ssc=args.ssc,
        kn=args.kn,
        cmb_cooling=args.cmb_cooling,
    )
    rvs_rad = None
    if args.rvs:
        rvs_rad = Radiation(
            eps_e=args.rvs_eps_e if args.rvs_eps_e is not None else args.eps_e,
            eps_B=args.rvs_eps_B if args.rvs_eps_B is not None else args.eps_B,
            p=args.rvs_p if args.rvs_p is not None else args.p,
            xi_e=args.rvs_xi_e if args.rvs_xi_e is not None else args.xi_e,
            ssc=args.rvs_ssc,
            kn=args.rvs_kn,
        )
    return fwd_rad, rvs_rad


def _format_energy(nu, latex=False):
    """Format a frequency as a human-readable energy/wavelength/frequency string.

    Uses wavelength (nm) for optical/IR/UV, photon energy for X-ray and above,
    and frequency (GHz/MHz/Hz) for radio.
    """
    E_keV = nu / _keV
    lam_nm = _c_A / nu / 10

    def _v(val, unit):
        s = f"{val:.3g}"
        return rf"${s}$ {unit}" if latex else f"{s} {unit}"

    if E_keV >= 1e9:
        return _v(E_keV / 1e9, "TeV")
    if E_keV >= 1e6:
        return _v(E_keV / 1e6, "GeV")
    if E_keV >= 1e3:
        return _v(E_keV / 1e3, "MeV")
    if E_keV >= 0.1:
        return _v(E_keV, "keV")
    if 100 < lam_nm < 10000:
        return _v(lam_nm, "nm")
    if nu >= 1e12:
        return _v(nu / 1e12, "THz")
    if nu >= 1e9:
        s = f"{nu / 1e9:g}"
        return rf"${s}$ GHz" if latex else f"{s} GHz"
    if nu >= 1e6:
        s = f"{nu / 1e6:.0f}"
        return rf"${s}$ MHz" if latex else f"{s} MHz"
    s = f"{nu:.0f}"
    return rf"${s}$ Hz" if latex else f"{s} Hz"


def _format_band(nu_min, nu_max, name=None, latex=False):
    """Format a frequency band label (plain text or LaTeX)."""
    lo = _format_energy(nu_min, latex=latex)
    hi = _format_energy(nu_max, latex=latex)
    sep = "\u2013" if latex else "-"
    range_str = f"{lo}{sep}{hi}"
    if name and name in _NAMED_BANDS:
        return f"{name} ({range_str})" if latex else f"{name}({range_str})"
    return range_str


def parse_frequencies(nu_args):
    """Parse --nu entries into point frequencies and bands.

    Returns (point_nus, point_labels, bands) where:
      - point_nus: ndarray of point frequencies [Hz] (may be empty)
      - point_labels: list of original strings for each point frequency
      - bands: list of (nu_min, nu_max, label) tuples (may be empty)
    """
    point_nus = []
    point_labels = []
    bands = []
    for v in nu_args:
        entry = _parse_nu_entry(v)
        if isinstance(entry, tuple):
            bands.append(entry)
        else:
            point_nus.append(entry)
            point_labels.append(v)
    nus = np.array(point_nus) if point_nus else np.array([])
    return nus, point_labels, bands


def _has_data(arr):
    """Check if a flux array has actual data (not 0-d or empty)."""
    return arr.ndim >= 1 and arr.size > 0


def _get_components(result):
    """Return list of (name, flux_array) for all active flux components.

    Sub-components are only included when there are two or more (i.e. when
    the decomposition is non-trivial).  A single sub-component is identical
    to total and would just add clutter.
    """
    components = [("total", result.total)]
    subs = []
    if _has_data(result.fwd.sync):
        subs.append(("fwd_sync", result.fwd.sync))
    if _has_data(result.fwd.ssc):
        subs.append(("fwd_ssc", result.fwd.ssc))
    if _has_data(result.rvs.sync):
        subs.append(("rvs_sync", result.rvs.sync))
    if _has_data(result.rvs.ssc):
        subs.append(("rvs_ssc", result.rvs.ssc))
    if len(subs) > 1:
        components.extend(subs)
    return components


def write_csv(times, nus, point_result, bands, band_results, args, file):
    t_scale = _TIME_SCALES[args.time_unit]
    f_scale = _FLUX_SCALES[args.flux_unit]

    has_points = point_result is not None and len(nus) > 0
    has_bands = len(bands) > 0

    if has_points:
        point_comps = _get_components(point_result)
    if has_bands:
        band_comps = _get_components(band_results[0])

    # Header
    file.write("# VegasAfterglow light curve\n")
    file.write(
        f"# jet={args.jet} theta_c={args.theta_c} E_iso={args.E_iso:.1e}"
        f" Gamma0={args.Gamma0} medium={args.medium}"
        f" z={args.z} theta_obs={args.theta_obs}\n"
    )
    file.write(f"# eps_e={args.eps_e} eps_B={args.eps_B} p={args.p}\n")
    if has_points and has_bands:
        file.write(
            f"# t_unit={args.time_unit}"
            f" flux_density_unit={args.flux_unit}"
            f" band_flux_unit=erg/cm2/s\n"
        )
    elif has_bands:
        file.write(f"# t_unit={args.time_unit} band_flux_unit=erg/cm2/s\n")
    else:
        file.write(f"# t_unit={args.time_unit} flux_unit={args.flux_unit}\n")

    # Column headers
    cols = [f"t({args.time_unit})"]
    if has_points:
        for nu in nus:
            nu_label = _format_nu_label(nu)
            for comp_name, _ in point_comps:
                cols.append(f"F_{comp_name}({nu_label})")
    if has_bands:
        for nu_min, nu_max, name in bands:
            band_label = _format_band(nu_min, nu_max, name)
            for comp_name, _ in band_comps:
                cols.append(f"F_{comp_name}({band_label})")
    file.write(",".join(cols) + "\n")

    # Data
    for j in range(len(times)):
        row = [f"{times[j] / t_scale:.6e}"]
        if has_points:
            for i in range(len(nus)):
                for _, comp_flux in point_comps:
                    row.append(f"{comp_flux[i, j] / f_scale:.6e}")
        if has_bands:
            for b_idx in range(len(bands)):
                for _, comp_flux in _get_components(band_results[b_idx]):
                    row.append(f"{comp_flux[j]:.6e}")
        file.write(",".join(row) + "\n")


def write_json(times, nus, point_result, bands, band_results, args, file):
    import json

    t_scale = _TIME_SCALES[args.time_unit]
    f_scale = _FLUX_SCALES[args.flux_unit]

    has_points = point_result is not None and len(nus) > 0
    has_bands = len(bands) > 0

    data = {
        "parameters": {
            "jet": args.jet,
            "theta_c": args.theta_c,
            "E_iso": args.E_iso,
            "Gamma0": args.Gamma0,
            "medium": args.medium,
            "z": args.z,
            "theta_obs": args.theta_obs,
            "eps_e": args.eps_e,
            "eps_B": args.eps_B,
            "p": args.p,
        },
        "units": {"time": args.time_unit},
        "times": (times / t_scale).tolist(),
    }

    if has_points:
        point_comps = _get_components(point_result)
        data["units"]["flux_density"] = args.flux_unit
        data["frequencies_Hz"] = nus.tolist()
        data["flux_density"] = {
            _format_nu_label(nus[i]): {
                comp_name: (comp_flux[i] / f_scale).tolist()
                for comp_name, comp_flux in point_comps
            }
            for i in range(len(nus))
        }

    if has_bands:
        data["units"]["band_flux"] = "erg/cm2/s"
        data["bands"] = {}
        for b_idx, (nu_min, nu_max, name) in enumerate(bands):
            band_label = _format_band(nu_min, nu_max, name)
            band_comps = _get_components(band_results[b_idx])
            data["bands"][band_label] = {
                "nu_min_Hz": nu_min,
                "nu_max_Hz": nu_max,
                "flux": {
                    comp_name: comp_flux.tolist() for comp_name, comp_flux in band_comps
                },
            }

    json.dump(data, file, indent=2)
    file.write("\n")


def write_output(times, nus, point_result, bands, band_results, args):
    writer = write_json if args.format == "json" else write_csv
    if args.output:
        with open(args.output, "w") as f:
            writer(times, nus, point_result, bands, band_results, args, f)
        print(f"Saved to {args.output}", file=sys.stderr)
    else:
        writer(times, nus, point_result, bands, band_results, args, sys.stdout)


_FLUX_LABELS = {
    "mJy": r"mJy",
    "Jy": r"Jy",
    "uJy": r"$\mu$Jy",
    "cgs": r"erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$",
}

_TIME_LABELS = {"s": "s", "day": "days", "hr": "hr", "min": "min"}


def _freq_colors(nus):
    """Map frequencies to colors, adapting to the frequency range.

    When the range spans > 3 decades, colors reflect absolute position on the
    radio-to-X-ray spectrum.  For narrower ranges, colors are evenly spaced
    across the palette so nearby frequencies remain visually distinct.
    """
    from matplotlib.colors import LinearSegmentedColormap

    palette = ["#E03530", "#E8872E", "#D4C43A", "#2AB07E", "#2878B5", "#7B3FA0"]
    cmap = LinearSegmentedColormap.from_list("freq", palette)

    log_nus = np.log10(nus)
    log_range = log_nus.max() - log_nus.min()

    if len(nus) == 1:
        t_vals = np.array([0.5])
    elif log_range > 3:
        # Wide range: map to absolute position on 10^7 – 10^22 Hz scale
        t_vals = np.clip((log_nus - 7) / (22 - 7), 0, 1)
    else:
        # Narrow range: spread evenly across the full palette
        t_vals = np.linspace(0, 1, len(nus))

    hex_colors = []
    for t in t_vals:
        rgba = cmap(t)
        hex_colors.append(
            f"#{int(rgba[0]*255):02x}{int(rgba[1]*255):02x}{int(rgba[2]*255):02x}"
        )
    return hex_colors


def _setup_plot_style(font=None):
    """Configure matplotlib for publication-quality output.

    Uses mathtext with STIX fonts (no LaTeX subprocess) for fast rendering.
    """
    import matplotlib as mpl

    _SANS_SERIF = {
        "helvetica",
        "arial",
        "liberation sans",
        "dejavu sans",
        "gill sans",
        "futura",
        "optima",
        "verdana",
        "tahoma",
    }

    if font is None:
        family = "sans-serif"
        font_list = ["Helvetica"]
        math_font = "Helvetica"
    elif font.lower() in _SANS_SERIF:
        family = "sans-serif"
        font_list = [font]
        math_font = font
    else:
        family = "serif"
        font_list = [font]
        math_font = font

    style = {
        "text.usetex": False,
        "mathtext.fontset": "custom",
        "mathtext.rm": math_font,
        "mathtext.it": f"{math_font}:italic",
        "mathtext.bf": f"{math_font}:bold",
        "font.family": family,
        f"font.{family}": font_list,
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "legend.fontsize": 5.5,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "axes.linewidth": 0.5,
        "lines.linewidth": 1.0,
        "axes.grid": True,
        "axes.grid.which": "both",
        "grid.color": "0.7",
        "grid.linestyle": ":",
        "grid.linewidth": 0.3,
        "grid.alpha": 0.4,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
    }

    mpl.rcParams.update(style)


def _sci_tex(val):
    """Format a number as LaTeX scientific notation: 1e52 -> 10^{52}."""
    if val == 0:
        return "0"
    exp = int(np.floor(np.log10(abs(val))))
    coeff = val / 10**exp
    if abs(coeff - 1) < 0.01:
        return rf"10^{{{exp}}}"
    return rf"{coeff:.1f}\times10^{{{exp}}}"


def _build_param_text(args):
    """Build a compact two-line parameter summary for figure title area."""
    jet_names = {"tophat": "Top-hat", "gaussian": "Gaussian", "powerlaw": "Power-law"}
    line1 = (
        rf"{jet_names.get(args.jet, args.jet)}, {args.medium.upper()}, "
        rf"$E_{{\mathrm{{iso}}}}={_sci_tex(args.E_iso)}$ erg, "
        rf"$\Gamma_0={args.Gamma0:g}$, "
        rf"$z={args.z:g}$"
    )
    line2 = (
        rf"$\theta_c={args.theta_c:g}$, "
        rf"$\theta_{{\mathrm{{obs}}}}={args.theta_obs:g}$ rad, "
        rf"$\epsilon_e={args.eps_e:g}$, "
        rf"$\epsilon_B={_sci_tex(args.eps_B)}$, "
        rf"$p={args.p:g}$"
    )
    if args.rvs:
        rvs_eps_e = args.rvs_eps_e if args.rvs_eps_e is not None else args.eps_e
        rvs_eps_B = args.rvs_eps_B if args.rvs_eps_B is not None else args.eps_B
        rvs_p = args.rvs_p if args.rvs_p is not None else args.p
        line2 += (
            rf", $\epsilon_{{e,r}}={rvs_eps_e:g}$, "
            rf"$\epsilon_{{B,r}}={_sci_tex(rvs_eps_B)}$, "
            rf"$p_r={rvs_p:g}$"
        )
    return line1 + "\n" + line2


_COMP_STYLES = {
    "total": "-",
    "fwd_sync": "--",
    "fwd_ssc": ":",
    "rvs_sync": "-.",
    "rvs_ssc": (0, (1, 2)),  # loosely dotted
}

_COMP_LABELS = {
    "fwd_sync": "FS sync",
    "fwd_ssc": "FS SSC",
    "rvs_sync": "RS sync",
    "rvs_ssc": "RS SSC",
}


def _add_legend(ax, plotted_comps, handles=None, labels=None):
    """Add legend with band entries + component linestyle key."""
    from matplotlib.lines import Line2D

    if handles is None or labels is None:
        handles, labels = ax.get_legend_handles_labels()
    if plotted_comps:
        for comp_name in ["fwd_sync", "fwd_ssc", "rvs_sync", "rvs_ssc"]:
            if comp_name in plotted_comps:
                handles.append(
                    Line2D(
                        [0],
                        [0],
                        color="0.4",
                        linestyle=_COMP_STYLES[comp_name],
                        linewidth=0.9,
                    )
                )
                labels.append(_COMP_LABELS[comp_name])

    ax.legend(
        handles,
        labels,
        loc="best",
        frameon=False,
        borderpad=0.4,
        handlelength=1.5,
    )


def _smart_ylim(ax, flux_arrays, scale=1.0):
    """Set smart log y-axis limits from flux arrays, capped at 12 decades."""
    all_pos = []
    for arr in flux_arrays:
        scaled = arr / scale
        pos = scaled[scaled > 0]
        if pos.size > 0:
            all_pos.append(pos)
    if all_pos:
        combined = np.concatenate(all_pos)
        f_max, f_min = np.max(combined), np.min(combined)
        y_top = f_max * 10
        y_bot = f_min / 10
        if y_top / y_bot > 1e12:
            y_bot = y_top * 1e-12
        ax.set_ylim(bottom=y_bot, top=y_top)


def _plot_point_ax(
    ax, times, nus, nu_labels, result, t_scale, f_scale, args, colors=None
):
    """Plot point-frequency flux density on the given axes.

    Returns the set of plotted sub-component names (for legend building).
    """
    components = _get_components(result)
    if colors is None:
        colors = _freq_colors(nus)
    plotted_comps = set()

    for i, nu in enumerate(nus):
        color = colors[i]
        freq_label = _format_nu_latex(nu, nu_labels[i])
        for comp_name, comp_flux in components:
            f = comp_flux[i] / f_scale
            mask = f > 0
            if not np.any(mask):
                continue
            is_total = comp_name == "total"
            ax.plot(
                (times / t_scale)[mask],
                f[mask],
                color=color,
                linestyle=_COMP_STYLES[comp_name],
                linewidth=1.2 if is_total else 0.9,
                alpha=1.0 if is_total else 0.7,
                label=freq_label if is_total else None,
                zorder=2 if is_total else 1,
            )
            if not is_total:
                plotted_comps.add(comp_name)

    ax.set_xscale("log")
    ax.set_yscale("log")
    _smart_ylim(ax, [cf for _, cf in components], f_scale)
    ax.set_ylabel(rf"$F_\nu$ ({_FLUX_LABELS[args.flux_unit]})")
    return plotted_comps


def _plot_band_ax(ax, times, bands, band_results, t_scale, colors=None):
    """Plot band-integrated flux on the given axes.

    Returns the set of plotted sub-component names (for legend building).
    """
    band_freqs = np.array([np.sqrt(nu_min * nu_max) for nu_min, nu_max, _ in bands])
    if colors is None:
        colors = _freq_colors(band_freqs)
    plotted_comps = set()

    for b_idx, (nu_min, nu_max, name) in enumerate(bands):
        color = colors[b_idx]
        label = _format_band(nu_min, nu_max, name, latex=True)
        components = _get_components(band_results[b_idx])
        for comp_name, comp_flux in components:
            f = comp_flux
            mask = f > 0
            if not np.any(mask):
                continue
            is_total = comp_name == "total"
            ax.plot(
                (times / t_scale)[mask],
                f[mask],
                color=color,
                linestyle=_COMP_STYLES[comp_name],
                linewidth=1.2 if is_total else 0.9,
                alpha=1.0 if is_total else 0.7,
                label=label if is_total else None,
                zorder=2 if is_total else 1,
            )
            if not is_total:
                plotted_comps.add(comp_name)

    ax.set_xscale("log")
    ax.set_yscale("log")
    all_flux = []
    for br in band_results:
        for _, cf in _get_components(br):
            all_flux.append(cf)
    _smart_ylim(ax, all_flux)
    ax.set_ylabel(r"$F$ (erg cm$^{-2}$ s$^{-1}$)")
    return plotted_comps


def plot_lightcurve(times, nus, nu_labels, point_result, bands, band_results, args):
    try:
        saving = args.output and args.output.lower().endswith(
            (".png", ".pdf", ".jpg", ".svg")
        )
        if saving:
            import matplotlib

            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "Error: --plot requires matplotlib. Install with: pip install matplotlib",
            file=sys.stderr,
        )
        sys.exit(1)

    _setup_plot_style(font=args.font)

    t_scale = _TIME_SCALES[args.time_unit]
    f_scale = _FLUX_SCALES[args.flux_unit]

    has_points = point_result is not None and len(nus) > 0
    has_bands = len(bands) > 0

    if has_points and has_bands:
        # Assign colors across all entries together for visual distinction
        all_freqs = list(nus) + [
            np.sqrt(nu_min * nu_max) for nu_min, nu_max, _ in bands
        ]
        all_colors = _freq_colors(np.array(all_freqs))
        pt_colors = all_colors[: len(nus)]
        bd_colors = all_colors[len(nus) :]

        fig, ax_left = plt.subplots(figsize=(3.5, 3.0))
        ax_right = ax_left.twinx()

        pt_comps = _plot_point_ax(
            ax_left,
            times,
            nus,
            nu_labels,
            point_result,
            t_scale,
            f_scale,
            args,
            colors=pt_colors,
        )
        bd_comps = _plot_band_ax(
            ax_right,
            times,
            bands,
            band_results,
            t_scale,
            colors=bd_colors,
        )

        ax_left.set_xlabel(rf"$t_\mathrm{{obs}}$ ({_TIME_LABELS[args.time_unit]})")

        # Combined legend from both axes
        h1, l1 = ax_left.get_legend_handles_labels()
        h2, l2 = ax_right.get_legend_handles_labels()
        _add_legend(ax_left, pt_comps | bd_comps, handles=h1 + h2, labels=l1 + l2)

        ax_left.set_title(_build_param_text(args), fontsize=7, color="0.3", pad=6)
        fig.subplots_adjust(left=0.16, right=0.84, bottom=0.18, top=0.87)
    elif has_bands:
        fig, ax = plt.subplots(figsize=(3.5, 3.0))
        bd_comps = _plot_band_ax(ax, times, bands, band_results, t_scale)
        ax.set_xlabel(rf"$t_\mathrm{{obs}}$ ({_TIME_LABELS[args.time_unit]})")
        _add_legend(ax, bd_comps)
        ax.set_title(_build_param_text(args), fontsize=7, color="0.3", pad=6)
        fig.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.87)
    else:
        fig, ax = plt.subplots(figsize=(3.5, 3.0))
        pt_comps = _plot_point_ax(
            ax,
            times,
            nus,
            nu_labels,
            point_result,
            t_scale,
            f_scale,
            args,
        )
        ax.set_xlabel(rf"$t_\mathrm{{obs}}$ ({_TIME_LABELS[args.time_unit]})")
        _add_legend(ax, pt_comps)
        ax.set_title(_build_param_text(args), fontsize=7, color="0.3", pad=6)
        fig.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.87)

    if args.output and args.output.lower().endswith((".png", ".pdf", ".jpg", ".svg")):
        dpi = 300 if args.output.lower().endswith(".png") else 150
        fig.savefig(args.output, dpi=dpi)
        print(f"Saved plot to {args.output}", file=sys.stderr)
    else:
        plt.show()


def main():
    args = parse_args()

    jet = build_jet(args)
    medium = build_medium(args)
    observer = build_observer(args)
    fwd_rad, rvs_rad = build_radiation(args)

    from . import Model

    model = Model(
        jet, medium, observer, fwd_rad, rvs_rad=rvs_rad, resolutions=tuple(args.res)
    )

    times = np.logspace(np.log10(args.t_min), np.log10(args.t_max), args.num_t)
    nus, nu_labels, bands = parse_frequencies(args.nu)

    n_point = len(nus)
    n_band = len(bands)
    shock_desc = "FS+RS" if args.rvs else "FS"
    parts = []
    if n_point:
        parts.append(f"{n_point} freq")
    if n_band:
        parts.append(f"{n_band} band{'s' if n_band > 1 else ''}")
    desc = ", ".join(parts) if parts else "0 freq"
    print(
        f"Computing {args.jet}/{args.medium} ({shock_desc}) light curve "
        f"({desc}, {args.num_t} time points)...",
        file=sys.stderr,
    )

    # Point-frequency flux density
    point_result = None
    if n_point:
        point_result = model.flux_density_grid(times, nus)

    # Band-integrated flux
    band_results = []
    for nu_min, nu_max, _ in bands:
        band_results.append(model.flux(times, nu_min, nu_max, args.num_nu_band))

    if args.plot:
        plot_lightcurve(times, nus, nu_labels, point_result, bands, band_results, args)
    else:
        write_output(times, nus, point_result, bands, band_results, args)


if __name__ == "__main__":
    main()
