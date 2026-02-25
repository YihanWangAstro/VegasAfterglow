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


def _format_nu_label(nu):
    """Format a frequency for CSV/JSON column headers."""
    named = {1e9: "radio", 5e14: "optical", 1e18: "xray"}
    if nu in named:
        return named[nu]
    if nu >= 1e9:
        return f"{nu:.2e}Hz"
    return f"{nu}Hz"


def _broad_band(nu):
    """Return the broad-band name for a frequency."""
    _h_cgs = 6.62607015e-27
    _eV_cgs = 1.602176634e-12
    E_eV = nu * _h_cgs / _eV_cgs

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
    """Format a frequency as a LaTeX label with broad-band prefix for plots.

    Parameters
    ----------
    nu : float
        Frequency in Hz.
    label : str or None
        Original user input string. If it's a filter name, show as filter;
        if numeric, show wavelength/energy.
    """
    from . import units

    _c_A = 2.99792458e18  # speed of light [A/s]
    _h_cgs = 6.62607015e-27
    _eV_cgs = 1.602176634e-12
    band = _broad_band(nu)

    # Show filter name only if the user actually typed a filter name
    if label is not None:
        try:
            float(label)
        except ValueError:
            # Not a number — user typed a filter name
            if label in units._VEGA_FILTERS:
                return rf"{band} (${label}$-band)"
            if label in units._ST_FILTERS:
                return rf"{band} ({label})"

    E_keV = nu * _h_cgs / _eV_cgs / 1e3
    lam_nm = _c_A / nu / 10  # wavelength in nm

    if nu < 1e6:
        return rf"{band} (${nu:.0f}$ Hz)"
    if nu < 1e9:
        return rf"{band} (${nu / 1e6:.0f}$ MHz)"
    if nu < 1e12:
        val = nu / 1e9
        s = f"{val:g}"
        return rf"{band} (${s}$ GHz)"
    if 100 < lam_nm < 10000:
        s = f"{lam_nm:.3g}"
        return rf"{band} (${s}$ nm)"
    E_GeV = E_keV / 1e6
    E_TeV = E_keV / 1e9
    if E_TeV >= 1:
        s = f"{E_TeV:.3g}"
        return rf"{band} (${s}$ TeV)"
    if E_GeV >= 1:
        s = f"{E_GeV:.3g}"
        return rf"{band} (${s}$ GeV)"
    if E_keV >= 1e3:
        s = f"{E_keV / 1e3:.3g}"
        return rf"{band} (${s}$ MeV)"
    if E_keV >= 0.1:
        s = f"{E_keV:.3g}"
        return rf"{band} (${s}$ keV)"
    # Fallback: scientific notation
    exp = int(np.floor(np.log10(nu)))
    coeff = nu / 10**exp
    if abs(coeff - 1) < 0.01:
        return rf"{band} ($10^{{{exp}}}$ Hz)"
    return rf"{band} (${coeff:.1f}\times10^{{{exp}}}$ Hz)"


_FLUX_SCALES = {"mJy": 1e-26, "Jy": 1e-23, "uJy": 1e-29, "cgs": 1.0}
_TIME_SCALES = {"s": 1.0, "day": 86400.0, "hr": 3600.0, "min": 60.0}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="vegasgen",
        description="Quick GRB afterglow light curve generator (VegasAfterglow).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  vegasgen                                  # default tophat/ISM/on-axis\n"
            "  vegasgen --jet gaussian --z 0.5           # Gaussian jet at z=0.5\n"
            "  vegasgen --nu R J F606W --plot            # filter-name frequencies + plot\n"
            "  vegasgen --medium wind --A_star 0.1 -o lc.csv\n"
        ),
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

    # -- Medium -------------------------------------------------------------
    med = p.add_argument_group("medium")
    med.add_argument("--medium", choices=["ism", "wind"], default="ism")
    med.add_argument(
        "--n_ism", type=float, default=1.0, help="ISM number density [cm^-3]"
    )
    med.add_argument("--A_star", type=float, default=0.1, help="wind parameter A*")

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

    # -- Frequencies --------------------------------------------------------
    freq = p.add_argument_group("frequencies")
    freq.add_argument(
        "--nu",
        nargs="+",
        default=["1e9", "5e14", "1e18"],
        help="frequencies in Hz or filter names (e.g. 1e9 R F606W)",
    )

    # -- Time grid ----------------------------------------------------------
    tg = p.add_argument_group("time grid")
    tg.add_argument("--t_min", type=float, default=100, help="start time [s]")
    tg.add_argument("--t_max", type=float, default=1e8, help="end time [s]")
    tg.add_argument("--num_t", type=int, default=200, help="number of time points")

    # -- Resolution ---------------------------------------------------------
    res = p.add_argument_group("resolution")
    res.add_argument(
        "--res",
        nargs=3,
        type=float,
        default=[0.15, 0.5, 10],
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
    return Wind(A_star=args.A_star)


def build_observer(args):
    from . import Observer

    d_L = args.lumi_dist if args.lumi_dist is not None else _lumi_dist_from_z(args.z)
    return Observer(lumi_dist=d_L, z=args.z, theta_obs=args.theta_obs)


def build_radiation(args):
    from . import Radiation

    return Radiation(
        eps_e=args.eps_e,
        eps_B=args.eps_B,
        p=args.p,
        xi_e=args.xi_e,
        ssc=args.ssc,
        kn=args.kn,
        cmb_cooling=args.cmb_cooling,
    )


def parse_frequencies(nu_args):
    """Return (nus_array, nu_labels) where nu_labels are the original strings."""
    nus = np.array([_parse_frequency(v) for v in nu_args])
    return nus, list(nu_args)


def write_csv(times, nus, result, args, file):
    t_scale = _TIME_SCALES[args.time_unit]
    f_scale = _FLUX_SCALES[args.flux_unit]

    # Header
    file.write("# VegasAfterglow light curve\n")
    file.write(
        f"# jet={args.jet} theta_c={args.theta_c} E_iso={args.E_iso:.1e}"
        f" Gamma0={args.Gamma0} medium={args.medium}"
        f" z={args.z} theta_obs={args.theta_obs}\n"
    )
    file.write(f"# eps_e={args.eps_e} eps_B={args.eps_B} p={args.p}\n")
    file.write(f"# t_unit={args.time_unit} flux_unit={args.flux_unit}\n")

    # Column headers
    cols = [f"t({args.time_unit})"]
    for nu in nus:
        cols.append(f"F_total({_format_nu_label(nu)})")
    file.write(",".join(cols) + "\n")

    # Data — result.total shape is (n_nu, n_t) from flux_density_grid
    flux = result.total
    for j in range(len(times)):
        row = [f"{times[j] / t_scale:.6e}"]
        for i in range(len(nus)):
            row.append(f"{flux[i, j] / f_scale:.6e}")
        file.write(",".join(row) + "\n")


def write_json(times, nus, result, args, file):
    import json

    t_scale = _TIME_SCALES[args.time_unit]
    f_scale = _FLUX_SCALES[args.flux_unit]
    flux = result.total

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
        "units": {"time": args.time_unit, "flux": args.flux_unit},
        "frequencies_Hz": nus.tolist(),
        "times": (times / t_scale).tolist(),
        "flux": {
            _format_nu_label(nus[i]): (flux[i] / f_scale).tolist()
            for i in range(len(nus))
        },
    }
    json.dump(data, file, indent=2)
    file.write("\n")


def write_output(times, nus, result, args):
    writer = write_json if args.format == "json" else write_csv
    if args.output:
        with open(args.output, "w") as f:
            writer(times, nus, result, args, f)
        print(f"Saved to {args.output}", file=sys.stderr)
    else:
        writer(times, nus, result, args, sys.stdout)


_FLUX_LABELS = {
    "mJy": r"mJy",
    "Jy": r"Jy",
    "uJy": r"$\mu$Jy",
    "cgs": r"erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$",
}

_TIME_LABELS = {"s": "s", "day": "days", "hr": "hr", "min": "min"}

# VegasAfterglow signature palette — vibrant, high-contrast, colorblind-friendly
_COLORS = [
    "#E63946",  # vivid red
    "#1D3557",  # deep navy
    "#2A9D8F",  # teal green
    "#E9C46A",  # warm gold
    "#F77F00",  # bright orange
    "#7209B7",  # deep violet
    "#4CC9F0",  # sky blue
]


def _setup_plot_style(font=None):
    """Configure matplotlib for publication-quality output."""
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
        family = "serif"
        font_list = ["Times New Roman"]
    elif font.lower() in _SANS_SERIF:
        family = "sans-serif"
        font_list = [font]
    else:
        family = "serif"
        font_list = [font]

    mpl.rcParams.update(
        {
            "text.usetex": True,
            "font.family": family,
            f"font.{family}": font_list,
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.minor.visible": True,
            "ytick.minor.visible": True,
            "xtick.major.size": 5,
            "ytick.major.size": 5,
            "xtick.minor.size": 3,
            "ytick.minor.size": 3,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.5,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.03,
        }
    )


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
        rf"$E_{{\rm iso}}={_sci_tex(args.E_iso)}$ erg, "
        rf"$\Gamma_0={args.Gamma0:g}$, "
        rf"$z={args.z:g}$"
    )
    line2 = (
        rf"$\theta_c={args.theta_c:g}$, "
        rf"$\theta_{{\rm obs}}={args.theta_obs:g}$ rad, "
        rf"$\epsilon_e={args.eps_e:g}$, "
        rf"$\epsilon_B={_sci_tex(args.eps_B)}$, "
        rf"$p={args.p:g}$"
    )
    return line1 + "\n" + line2


def plot_lightcurve(times, nus, nu_labels, result, args):
    try:
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
    flux = result.total

    # Single-column journal figure: 3.5 in wide, golden ratio height
    fig, ax = plt.subplots(figsize=(3.5, 2.6))

    for i, nu in enumerate(nus):
        color = _COLORS[i % len(_COLORS)]
        f = flux[i] / f_scale
        mask = f > 0
        ax.plot(
            (times / t_scale)[mask],
            f[mask],
            color=color,
            label=_format_nu_latex(nu, nu_labels[i]),
            zorder=2,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(rf"$t_\mathrm{{obs}}$ ({_TIME_LABELS[args.time_unit]})")
    ax.set_ylabel(rf"$F_\nu$ ({_FLUX_LABELS[args.flux_unit]})")

    ax.legend(
        loc="best",
        frameon=True,
        fancybox=False,
        edgecolor="0.7",
        framealpha=0.9,
        borderpad=0.4,
        handlelength=1.5,
    )

    # Parameter text above the plot
    ax.set_title(
        _build_param_text(args),
        fontsize=7,
        color="0.3",
        pad=6,
    )

    fig.tight_layout(pad=0.3)

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
    radiation = build_radiation(args)

    from . import Model

    model = Model(jet, medium, observer, radiation, resolutions=tuple(args.res))

    times = np.logspace(np.log10(args.t_min), np.log10(args.t_max), args.num_t)
    nus, nu_labels = parse_frequencies(args.nu)

    print(
        f"Computing {args.jet}/{args.medium} light curve "
        f"({len(nus)} bands, {args.num_t} time points)...",
        file=sys.stderr,
    )

    result = model.flux_density_grid(times, nus)

    if args.plot:
        plot_lightcurve(times, nus, nu_labels, result, args)
    else:
        write_output(times, nus, result, args)


if __name__ == "__main__":
    main()
