#!/usr/bin/env python
"""Generate documentation figures for VegasAfterglow docs.

Usage:
    python docs/source/gen_figures.py

Saves all PNGs to assets/ (the git-tracked image directory).
build_docs.sh copies them to docs/source/_static/images/ at build time.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Resolve project root (two levels up from this script)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
ASSETS_DIR = os.path.join(PROJECT_ROOT, "assets")

from VegasAfterglow import (
    TophatJet, PowerLawJet, TwoComponentJet,
    ISM, Wind, Observer, Radiation, Model,
)
from VegasAfterglow.units import uas


def sky_image_single():
    """Single-frame on-axis sky image (matches sky_image.rst 'Single Frame')."""
    model = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=200),
        medium=ISM(n_ism=1),
        observer=Observer(lumi_dist=1e26, z=0.1, theta_obs=0),
        fwd_rad=Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3),
    )

    img = model.sky_image([1e6], nu_obs=1e9, fov=500 * uas, npixel=128)

    fig, ax = plt.subplots(dpi=200)
    extent = img.extent / uas

    im = ax.imshow(
        img.image[0].T,
        origin="lower",
        extent=extent,
        cmap="inferno",
        norm=LogNorm(),
    )
    ax.set_xlabel(r"$\Delta x$ ($\mu$as)")
    ax.set_ylabel(r"$\Delta y$ ($\mu$as)")
    ax.set_title(r"$t_{\rm obs} = 10^6$ s, $\nu = 1$ GHz")
    fig.colorbar(im, label=r"Surface brightness (erg/cm$^2$/s/Hz/sr)")
    plt.tight_layout()
    fig.savefig(os.path.join(ASSETS_DIR, "sky_image_single.png"), bbox_inches="tight")
    plt.close(fig)
    print("  -> sky_image_single.png")


def sky_image_offaxis():
    """Off-axis 3-panel sky image evolution (matches sky_image.rst 'Off-Axis Observer')."""
    model = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=200),
        medium=ISM(n_ism=1),
        observer=Observer(lumi_dist=1e26, z=0.1, theta_obs=0.4),
        fwd_rad=Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3),
    )

    times = np.array([1e5, 1e6, 1e7])
    imgs = model.sky_image(times, nu_obs=1e9, fov=5000 * uas, npixel=128)
    extent = imgs.extent / uas

    vmin = imgs.image[imgs.image > 0].min()
    vmax = imgs.image.max()

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), dpi=200,
                             constrained_layout=True)
    for i, (ax, t) in enumerate(zip(axes, times)):
        im = ax.imshow(
            imgs.image[i].T,
            origin="lower",
            extent=extent,
            cmap="inferno",
            norm=LogNorm(vmin=vmin, vmax=vmax),
        )
        exp = int(np.floor(np.log10(t)))
        ax.set_title(rf"$t_{{\rm obs}} = 10^{exp}$ s")
        ax.set_xlabel(r"$\Delta x$ ($\mu$as)")
        if i == 0:
            ax.set_ylabel(r"$\Delta y$ ($\mu$as)")
    fig.colorbar(im, ax=axes, label=r"erg/cm$^2$/s/Hz/sr", shrink=0.85)
    fig.savefig(os.path.join(ASSETS_DIR, "sky_image_offaxis.png"), bbox_inches="tight")
    plt.close(fig)
    print("  -> sky_image_offaxis.png")


def sky_image_flux_comparison():
    """Image-integrated flux vs direct calculation (matches sky_image.rst 'Flux from Image')."""
    model = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=200),
        medium=ISM(n_ism=1),
        observer=Observer(lumi_dist=1e26, z=0.1, theta_obs=0),
        fwd_rad=Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3),
    )

    t_obs = np.logspace(3, 8, 30)
    nu_obs = 1e9

    img = model.sky_image(t_obs, nu_obs=nu_obs, fov=2000 * uas, npixel=128)
    flux_from_image = img.image.sum(axis=(1, 2)) * img.pixel_solid_angle

    flux_direct = model.flux_density_grid(t_obs, np.array([nu_obs])).total[0, :]

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(5, 5), dpi=200, sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )

    ax1.loglog(t_obs, flux_direct, "k-", label="flux_density_grid")
    ax1.loglog(t_obs, flux_from_image, "o", ms=4, color="C1", label="sky_image (integrated)")
    ax1.set_ylabel(r"Flux density (erg/cm$^2$/s/Hz)")
    ax1.legend()

    ratio = flux_from_image / flux_direct
    ax2.semilogx(t_obs, ratio, "o-", ms=4, color="C1")
    ax2.axhline(1, color="k", ls="--", lw=0.8)
    ax2.set_ylabel("image / direct")
    ax2.set_xlabel("Observer time (s)")
    ax2.set_ylim(0.95, 1.05)
    plt.tight_layout()
    fig.savefig(os.path.join(ASSETS_DIR, "sky_image_flux_comparison.png"), bbox_inches="tight")
    plt.close(fig)
    print("  -> sky_image_flux_comparison.png")


def reverse_shock_lc():
    """Forward + reverse shock light curves (matches models.rst 'Reverse Shock Emission')."""
    jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300, duration=100)
    fwd_rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3)
    rvs_rad = Radiation(eps_e=1e-2, eps_B=1e-4, p=2.4)

    model = Model(
        jet=jet,
        medium=ISM(n_ism=1),
        observer=Observer(lumi_dist=1e26, z=0.1, theta_obs=0),
        fwd_rad=fwd_rad,
        rvs_rad=rvs_rad,
        resolutions=(0.1, 0.25, 10),
    )

    times = np.logspace(2, 8, 200)
    bands = np.array([1e9, 1e14, 1e17])
    results = model.flux_density_grid(times, bands)

    fig, ax = plt.subplots(figsize=(5, 3.6), dpi=200)
    colors = ["C0", "C1", "C2"]
    for i, nu in enumerate(bands):
        exp = int(np.floor(np.log10(nu)))
        base = nu / 10**exp
        label_fwd = rf"${base:.0f} \times 10^{{{exp}}}$ Hz (fwd)"
        label_rvs = rf"${base:.0f} \times 10^{{{exp}}}$ Hz (rvs)"
        fwd = results.fwd.sync[i, :]
        rvs = results.rvs.sync[i, :]
        ax.loglog(times, fwd, color=colors[i], label=label_fwd)
        ax.loglog(times, rvs, color=colors[i], ls="--", label=label_rvs)

    # Clip y-axis to reasonable range around the data
    fwd_max = max(results.fwd.sync[i, :].max() for i in range(len(bands)))
    ax.set_ylim(fwd_max * 1e-6, fwd_max * 10)
    ax.set_xlabel("Observer time (s)")
    ax.set_ylabel(r"Flux density (erg/cm$^2$/s/Hz)")
    ax.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    fig.savefig(os.path.join(ASSETS_DIR, "reverse_shock_lc.png"), bbox_inches="tight")
    plt.close(fig)
    print("  -> reverse_shock_lc.png")


def ssc_lc():
    """Synchrotron + SSC light curves (matches models.rst 'Self-Synchrotron Compton')."""
    model = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300),
        medium=ISM(n_ism=1),
        observer=Observer(lumi_dist=1e26, z=0.1, theta_obs=0),
        fwd_rad=Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3, ssc=True, kn=True),
    )

    times = np.logspace(2, 8, 200)
    bands = np.array([1e9, 1e14, 1e17])
    results = model.flux_density_grid(times, bands)

    fig, ax = plt.subplots(figsize=(5, 3.6), dpi=200)
    colors = ["C0", "C1", "C2"]
    for i, nu in enumerate(bands):
        exp = int(np.floor(np.log10(nu)))
        base = nu / 10**exp
        label_sync = rf"${base:.0f} \times 10^{{{exp}}}$ Hz (sync)"
        label_ssc = rf"${base:.0f} \times 10^{{{exp}}}$ Hz (SSC)"
        ax.loglog(times, results.fwd.sync[i, :], color=colors[i], label=label_sync)
        ax.loglog(times, results.fwd.ssc[i, :], color=colors[i], ls="--", label=label_ssc)
    ax.set_xlabel("Observer time (s)")
    ax.set_ylabel(r"Flux density (erg/cm$^2$/s/Hz)")
    ax.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    fig.savefig(os.path.join(ASSETS_DIR, "ssc_lc.png"), bbox_inches="tight")
    plt.close(fig)
    print("  -> ssc_lc.png")


def basic_lightcurves():
    """Light curves + spectra (matches basic_usage.rst 'Setting up a simple afterglow model')."""
    medium = ISM(n_ism=1)
    jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
    obs = Observer(lumi_dist=1e26, z=0.1, theta_obs=0)
    rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3, xi_e=1)
    model = Model(jet=jet, medium=medium, observer=obs, fwd_rad=rad, resolutions=(0.1, 0.25, 10))

    # --- Light curves ---
    times = np.logspace(2, 8, 200)
    bands = np.array([1e9, 1e14, 1e17])
    results = model.flux_density_grid(times, bands)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.6), dpi=200)

    for i, nu in enumerate(bands):
        exp = int(np.floor(np.log10(nu)))
        base = nu / 10**exp
        ax1.loglog(times, results.total[i, :],
                   label=rf"${base:.1f} \times 10^{{{exp}}}$ Hz")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel(r"Flux Density (erg/cm$^2$/s/Hz)")
    ax1.legend(fontsize=8)
    ax1.set_title("Multi-wavelength Light Curves")

    # --- Spectra ---
    frequencies = np.logspace(5, 22, 200)
    epochs = np.array([1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8])
    results_spec = model.flux_density_grid(epochs, frequencies)
    colors = plt.cm.viridis(np.linspace(0, 1, len(epochs)))

    for i, t in enumerate(epochs):
        exp = int(np.floor(np.log10(t)))
        base = t / 10**exp
        ax2.loglog(frequencies, results_spec.total[:, i], color=colors[i],
                   label=rf"${base:.1f} \times 10^{{{exp}}}$ s")
    for i, band in enumerate(bands):
        ax2.axvline(band, ls="--", color=f"C{i}", alpha=0.5)
    ax2.set_xlabel("Frequency (Hz)")
    ax2.set_ylabel(r"Flux Density (erg/cm$^2$/s/Hz)")
    ax2.legend(fontsize=6, ncol=2)
    ax2.set_title("Synchrotron Spectra")

    plt.tight_layout()
    fig.savefig(os.path.join(ASSETS_DIR, "basic_lc_spec.png"), bbox_inches="tight")
    plt.close(fig)
    print("  -> basic_lc_spec.png")


def basic_bolometric():
    """Bolometric flux (matches basic_usage.rst 'Calculate bolometric flux')."""
    medium = ISM(n_ism=1)
    jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
    obs = Observer(lumi_dist=1e26, z=0.1, theta_obs=0)
    rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3, xi_e=1)
    model = Model(jet=jet, medium=medium, observer=obs, fwd_rad=rad, resolutions=(0.1, 0.25, 10))

    times = np.logspace(2, 8, 100)
    flux_bat = model.flux(times, 3.6e18, 3.6e19, 20)
    flux_v = model.flux(times, 4.6e14, 5.6e14, 20)

    fig, ax = plt.subplots(figsize=(5, 3.6), dpi=200)
    ax.loglog(times, flux_bat.total, label="Swift/BAT (15-150 keV)", linewidth=2)
    ax.loglog(times, flux_v.total, label="V-band optical", linewidth=2)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(r"Integrated Flux (erg/cm$^2$/s)")
    ax.legend()
    ax.set_title("Broadband Light Curves")
    plt.tight_layout()
    fig.savefig(os.path.join(ASSETS_DIR, "basic_bolometric.png"), bbox_inches="tight")
    plt.close(fig)
    print("  -> basic_bolometric.png")


def introspection_jet():
    """Power-law jet profiles (matches introspection.rst 'Jet Property Introspection')."""
    jet = PowerLawJet(theta_c=0.1, E_iso=1e52, Gamma0=300, k_e=2.0, k_g=1.5)
    medium = ISM(n_ism=1)
    obs = Observer(lumi_dist=1e26, z=0.1, theta_obs=0)
    rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3)
    model = Model(jet=jet, medium=medium, observer=obs, fwd_rad=rad)

    phi = 0.0
    theta = np.linspace(0, 0.5, 100)
    E_iso_profile = model.jet_E_iso(phi, theta)
    Gamma0_profile = model.jet_Gamma0(phi, theta)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), dpi=200)
    ax1.semilogy(np.degrees(theta), E_iso_profile)
    ax1.set_xlabel("Polar Angle (degrees)")
    ax1.set_ylabel(r"$E_{\rm iso}$ (erg)")
    ax1.set_title("Jet Energy Profile")
    ax1.grid(True, alpha=0.3)

    ax2.semilogy(np.degrees(theta), Gamma0_profile)
    ax2.set_xlabel("Polar Angle (degrees)")
    ax2.set_ylabel(r"$\Gamma_0$")
    ax2.set_title("Jet Lorentz Factor Profile")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(os.path.join(ASSETS_DIR, "introspection_jet.png"), bbox_inches="tight")
    plt.close(fig)
    print("  -> introspection_jet.png")


def introspection_medium():
    """Wind medium density profile (matches introspection.rst 'Medium Density Introspection')."""
    jet = PowerLawJet(theta_c=0.1, E_iso=1e52, Gamma0=300, k_e=2.0, k_g=1.5)
    wind = Wind(A_star=1.0, n_ism=0.1, n0=1e3, k_m=2)
    obs = Observer(lumi_dist=1e26, z=0.1, theta_obs=0)
    rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3)
    model = Model(jet=jet, medium=wind, observer=obs, fwd_rad=rad)

    phi = 0.0
    theta = 0.1
    r = np.logspace(13, 20, 200)
    rho_profile = model.medium(phi, theta, r)
    n_profile = rho_profile / 1.67e-24  # approximate number density (assuming pure hydrogen)

    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
    ax.loglog(r, n_profile)
    ax.set_xlabel("Radius (cm)")
    ax.set_ylabel(r"Number Density (cm$^{-3}$)")
    ax.set_title("Medium Density Profile")
    ax.grid(True, alpha=0.3)
    ax.axhline(1e3, color="red", ls="--", alpha=0.7, label=r"Inner density ($n_0$)")
    ax.axhline(0.1, color="blue", ls="--", alpha=0.7, label=r"Outer ISM density ($n_{\rm ism}$)")
    ax.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(ASSETS_DIR, "introspection_medium.png"), bbox_inches="tight")
    plt.close(fig)
    print("  -> introspection_medium.png")


def introspection_twocomp():
    """Two-component jet profiles (matches introspection.rst 'Two-Component Jet Analysis')."""
    jet = TwoComponentJet(
        theta_c=0.05, E_iso=1e53, Gamma0=300,
        theta_w=0.15, E_iso_w=1e52, Gamma0_w=100,
    )
    medium = ISM(n_ism=1)
    obs = Observer(lumi_dist=1e26, z=0.1, theta_obs=0)
    rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3)
    model = Model(jet=jet, medium=medium, observer=obs, fwd_rad=rad)

    theta = np.linspace(0, 0.3, 200)
    E_iso_profile = model.jet_E_iso(0, theta)
    Gamma0_profile = model.jet_Gamma0(0, theta)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 7), dpi=200)
    ax1.semilogy(np.degrees(theta), E_iso_profile)
    ax1.axvline(np.degrees(0.05), color="red", ls="--", alpha=0.7, label="Core boundary")
    ax1.axvline(np.degrees(0.15), color="blue", ls="--", alpha=0.7, label="Wide component boundary")
    ax1.set_ylabel(r"$E_{\rm iso}$ (erg)")
    ax1.set_title("Two-Component Jet: Energy Profile")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.semilogy(np.degrees(theta), Gamma0_profile)
    ax2.axvline(np.degrees(0.05), color="red", ls="--", alpha=0.7, label="Core boundary")
    ax2.axvline(np.degrees(0.15), color="blue", ls="--", alpha=0.7, label="Wide component boundary")
    ax2.set_xlabel("Polar Angle (degrees)")
    ax2.set_ylabel(r"$\Gamma_0$")
    ax2.set_title("Two-Component Jet: Lorentz Factor Profile")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(os.path.join(ASSETS_DIR, "introspection_twocomp.png"), bbox_inches="tight")
    plt.close(fig)
    print("  -> introspection_twocomp.png")


if __name__ == "__main__":
    os.makedirs(ASSETS_DIR, exist_ok=True)
    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 11,
    })

    print("Generating documentation figures...")
    sky_image_single()
    sky_image_offaxis()
    sky_image_flux_comparison()
    reverse_shock_lc()
    ssc_lc()
    basic_lightcurves()
    basic_bolometric()
    introspection_jet()
    introspection_medium()
    introspection_twocomp()
    print("Done. All figures saved to assets/")
