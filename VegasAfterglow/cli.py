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
    """Format a frequency for column headers."""
    named = {1e9: "radio", 5e14: "optical", 1e18: "xray"}
    if nu in named:
        return named[nu]
    if nu >= 1e9:
        return f"{nu:.2e}Hz"
    return f"{nu}Hz"


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
    return np.array([_parse_frequency(v) for v in nu_args])


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


def plot_lightcurve(times, nus, result, args):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "Error: --plot requires matplotlib. Install with: pip install matplotlib",
            file=sys.stderr,
        )
        sys.exit(1)

    t_scale = _TIME_SCALES[args.time_unit]
    f_scale = _FLUX_SCALES[args.flux_unit]
    flux = result.total

    fig, ax = plt.subplots(figsize=(8, 5))
    for i, nu in enumerate(nus):
        ax.plot(times / t_scale, flux[i] / f_scale, label=_format_nu_label(nu))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(f"Time ({args.time_unit})")
    ax.set_ylabel(f"Flux density ({args.flux_unit})")
    ax.set_title(
        f"{args.jet} jet, {args.medium.upper()}, "
        f"z={args.z}, θ_obs={args.theta_obs:.2f} rad"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if args.output and args.output.lower().endswith((".png", ".pdf", ".jpg", ".svg")):
        fig.savefig(args.output, dpi=150)
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
    nus = parse_frequencies(args.nu)

    print(
        f"Computing {args.jet}/{args.medium} light curve "
        f"({len(nus)} bands, {args.num_t} time points)...",
        file=sys.stderr,
    )

    result = model.flux_density_grid(times, nus)

    if args.plot:
        plot_lightcurve(times, nus, result, args)
    else:
        write_output(times, nus, result, args)


if __name__ == "__main__":
    main()
