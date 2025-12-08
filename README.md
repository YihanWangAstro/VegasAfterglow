# VegasAfterglow

<img align="left" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/logo.svg" alt="VegasAfterglow Logo" width="350"/>

[![C++ Version](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![PyPI version](https://img.shields.io/pypi/v/VegasAfterglow.svg)](https://pypi.org/project/VegasAfterglow/)
[![Build Status](https://github.com/YihanWangAstro/VegasAfterglow/actions/workflows/PyPI-build.yml/badge.svg)](https://github.com/YihanWangAstro/VegasAfterglow/actions/workflows/PyPI-build.yml)
[![License](https://img.shields.io/badge/License-BSD--3--Clause-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-Linux%20|%20macOS%20|%20Windows-lightgrey.svg)]()
[![Python Version](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![Documentation](https://img.shields.io/badge/Documentation-Online-brightgreen.svg)](https://yihanwangastro.github.io/VegasAfterglow/docs/index.html)

<div align="left">

**[Latest Release Notes](CHANGELOG.md#v030---2025-09-06)** | **[Full Changelog](CHANGELOG.md)** | **[Install Now](#installation)**
</div>

**VegasAfterglow** is a high-performance C++ framework with a user-friendly Python interface designed for the comprehensive modeling of afterglows. It achieves exceptional computational efficiency, enabling the generation of multi-wavelength light curves in milliseconds and facilitating robust Markov Chain Monte Carlo (MCMC) parameter inference in seconds to minutes. The framework incorporates advanced models for shock dynamics (both forward and reverse shocks), diverse radiation mechanisms (synchrotron with self-absorption, and inverse Compton scattering with Klein-Nishina corrections), and complex structured jet configurations. For lightweight afterglow modeling, one can also consider the [PyFRS](https://github.com/leiwh/PyFRS) package.
<br clear="left"/>

---

## Table of Contents

- [VegasAfterglow](#vegasafterglow)
  - [Table of Contents](#table-of-contents)
  - [Features](#features)
  - [Performance Highlights](#performance-highlights)
  - [Installation](#installation)
    - [Python Installation](#python-installation)
    - [C++ Installation](#c-installation)
  - [Usage](#usage)
    - [Quick Start](#quick-start)
    - [Light Curve \& Spectrum Calculation](#light-curve--spectrum-calculation)
    - [Internal Quantities Evolution](#internal-quantities-evolution)
    - [MCMC Parameter Fitting](#mcmc-parameter-fitting)
  - [Documentation](#documentation)
  - [Contributing](#contributing)
  - [License](#license)
  - [Acknowledgments \& Citation](#acknowledgments--citation)

---

## Features

<h3 align="center">Shock Dynamics</h3>

<img align="right" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/shock_dynamics.svg" width="450"/>

- **Forward and Reverse Shock Modeling:** Simulates both shocks via shock crossing dynamics with arbitrary magnetization levels and shell thicknesses.
- **Relativistic and Non-Relativistic Regimes:** Accurately models shock evolution across all velocity regimes.
- **Adiabatic and Radiative Blast Waves:** Supports smooth transition between adiabatic and radiative blast waves.
- **Ambient Medium:** Supports uniform Interstellar Medium (ISM), stellar wind environments, and user-defined density profiles.
- **Energy and Mass Injection:** Supports user-defined profiles for continuous energy and/or mass injection into the blast wave.

<br clear="right"/>

<h3 align="center">Jet Structure & Geometry</h3>

<img align="right" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/jet_geometry.svg" width="450"/>

- **Structured Jet Profiles:** Allows user-defined angular profiles for energy distribution, initial Lorentz factor, and magnetization.
- **Arbitrary Viewing Angles:** Supports off-axis observers at any viewing angle relative to the jet axis.
- **Jet Spreading:** Includes lateral expansion dynamics for realistic jet evolution (experimental).
- **Non-Axisymmetric Jets:** Capable of modeling complex, non-axisymmetric jet structures.

<br clear="right"/>

<h3 align="center">Radiation Mechanisms</h3>

<img align="right" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/radiation_mechanisms.svg" width="450"/>

- **Synchrotron Radiation:** Calculates synchrotron emission from shocked electrons.
- **Synchrotron Self-Absorption (SSA):** Includes SSA effects, crucial at low frequencies.
- **Inverse Compton (IC) Scattering:** Models IC processes, including:
  - Synchrotron Self-Compton (SSC) from both forward and reverse shocks.
  - Pairwise IC between forward and reverse shock electron and photon populations (experimental).
  - Includes Klein-Nishina corrections for accurate synchrotron and IC emission.

<br clear="right"/>

---

## Performance Highlights

<img align="right" src="https://github.com/YihanWangAstro/VegasAfterglow/raw/main/assets/convergence_plot.png" width="400"/>

VegasAfterglow delivers exceptional computational performance through deep optimization of its core algorithms:

- **Ultra-fast Light Curve Computation:** Generates a 100-point single-frequency light curve (forward shock & synchrotron only) from a structured jet viewed off-axis in approximately 1 millisecond on an Apple M2 chip with a single core.

- **Rapid MCMC Exploration:** Enables parameter estimation with 10,000 MCMC steps for 8 parameters on 20 data points across multi-wavelength light curves and spectra on an 8-core Apple M2 chip in:
  - ~50 seconds for on-axis structured jet scenarios

This level of performance is achieved through optimized algorithm implementation and efficient memory access patterns, facilitating comprehensive Bayesian inference on standard laptop hardware in seconds to minutes rather than hours or days. The accelerated convergence speed enables rapid iteration through different physical models and makes VegasAfterglow suitable for both detailed analysis of individual GRB events and large-scale population studies.

<br clear="right"/>

---

## Installation

VegasAfterglow is available as a Python package with C++ source code also provided for direct use.

### Python Installation

To install VegasAfterglow using pip:

```bash
pip install VegasAfterglow
```

This is the recommended method for most users. VegasAfterglow requires Python 3.8 or higher.

<details>
<summary><b>Alternative: Install from Source</b> <i>(click to expand/collapse)</i></summary>
<br>

For cases where pip installation is not viable or when the development version is required:

1. Clone this repository:

```bash
git clone https://github.com/YihanWangAstro/VegasAfterglow.git
```

2. Navigate to the directory and install the Python package:

```bash
cd VegasAfterglow
pip install .
```

Standard development environments typically include the necessary prerequisites (C++20 compatible compiler). For build-related issues, refer to the prerequisites section in [C++ Installation](#c-installation).
</details>

### C++ Installation

For advanced users who need to compile and use the C++ library directly:

<details>
<summary><b>Instructions for C++ Installation</b> <i>(click to expand/collapse)</i></summary>
<br>

1. Clone the repository (if not previously done):

```bash
git clone https://github.com/YihanWangAstro/VegasAfterglow.git
cd VegasAfterglow
```

2. Compile and run tests:

```bash
make tests
```

Upon successful compilation, you can create custom C++ problem generators using the VegasAfterglow interfaces. For implementation details, refer to the [Creating Custom Problem Generators with C++](#creating-custom-problem-generators-with-c) section or examine the example problem generators in `tests/demo/`.

<details>
<summary><b>Build Prerequisites</b> <i>(click to expand for dependency information)</i></summary>
<br>

The following development tools are required:

- **C++20 compatible compiler**:
  - **Linux**: GCC 10+ or Clang 13+
  - **macOS**: Apple Clang 13+ (with Xcode 13+) or GCC 10+ (via Homebrew)
  - **Windows**: MSVC 19.29+ (Visual Studio 2019 16.10+) or MinGW-w64 with GCC 10+

- **Build tools**:
  - Make (GNU Make 4.0+ recommended)

</details>
</details>

---

## Usage

### Quick Start

We provide basic example scripts (`script/quick.ipynb`, `script/details.ipynb` and `script/mcmc.ipynb`) that demonstrate how to set up and run afterglow simulations. This section shows how to calculate light curves and spectra for a simple GRB afterglow model without the need for observational data and perform MCMC parameter fitting with observational data. The notebook can be run using either Jupyter Notebook or VSCode with the Jupyter extension.

To avoid conflicts when updating the repository in the future, make a copy of the example notebook in the same directory and work with the copy instead of the original.

### Light Curve & Spectrum Calculation

The example below walks through the main components needed to model a GRB afterglow, from setting up the physical parameters to producing light curves and spectra via `script/quick.ipynb`.

<details>
<summary><b>Model Setup</b> <i>(click to expand/collapse)</i></summary>
<br>

First, let's set up the physical components of our afterglow model, including the environment, jet, observer, and radiation parameters:

```python
import numpy as np
import matplotlib.pyplot as plt
from VegasAfterglow import ISM, TophatJet, Observer, Radiation, Model

# 1. Define the circumburst environment (constant density ISM)
medium = ISM(n_ism=1) #in cgs unit

# 2. Configure the jet structure (top-hat with opening angle, energy, and Lorentz factor)
jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300) #in cgs unit

# 3. Set observer parameters (distance, redshift, viewing angle)
obs = Observer(lumi_dist=1e26, z=0.1, theta_obs=0) #in cgs unit

# 4. Define radiation microphysics parameters
rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3)

# 5. Combine all components into a complete afterglow model
model = Model(jet=jet, medium=medium, observer=obs, fwd_rad=rad)
```

</details>

<details>
<summary><b>Light Curve Calculation</b> <i>(click to expand/collapse)</i></summary>
<br>

Now, let's compute and plot multi-wavelength light curves to see how the afterglow evolves over time:

```python
# 1. Create logarithmic time array from 10² to 10⁸ seconds (100s to ~3yrs)
times = np.logspace(2, 8, 100)

# 2. Define observing frequencies (radio, optical, X-ray bands in Hz)
bands = np.array([1e9, 1e14, 1e17])

# 3. Calculate the afterglow emission at each time and frequency
# NOTE: times array must be in ascending order, frequencies can be in random order
results = model.flux_density_grid(times, bands)

# 4. Visualize the multi-wavelength light curves
plt.figure(figsize=(4.8, 3.6),dpi=200)

# 5. Plot each frequency band
for i, nu in enumerate(bands):
    exp = int(np.floor(np.log10(nu)))
    base = nu / 10**exp
    plt.loglog(times, results.total[i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz')

def add_note(plt):
    plt.annotate('jet break',xy=(3e4, 1e-26), xytext=(3e3, 5e-28), arrowprops=dict(arrowstyle='->'))
    plt.annotate(r'$\nu_m=\nu_a$',xy=(8e5, 2e-25), xytext=(7.5e4, 5e-24), arrowprops=dict(arrowstyle='->'))
    plt.annotate(r'$\nu=\nu_a$',xy=(4e6, 4e-25), xytext=(7.5e5, 5e-24), arrowprops=dict(arrowstyle='->'))

add_note(plt)
plt.xlabel('Time (s)')
plt.ylabel('Flux Density (erg/cm²/s/Hz)')
plt.legend()
plt.title('Light Curves')
plt.savefig('quick-lc.png',dpi=300,bbox_inches='tight')
```

<div align="center">
<img src="assets/quick-lc.png" alt="Afterglow Light Curves" width="600"/>

Running the light curve script will produce this figure showing the afterglow evolution across different frequencies.
</div>
</details>

<details>
<summary><b>Spectrum Analysis</b> <i>(click to expand/collapse)</i></summary>
<br>

We can also examine how the broadband spectrum evolves at different times after the burst:

```python
# 1. Define broad frequency range (10⁵ to 10²² Hz)
frequencies = np.logspace(5, 22, 100)

# 2. Select specific time epochs for spectral snapshots
epochs = np.array([1e2, 1e3, 1e4, 1e5 ,1e6, 1e7, 1e8])

# 3. Calculate spectra at each epoch
# NOTE: epochs array must be in ascending order, frequencies can be in random order
results = model.flux_density_grid(epochs, frequencies)

# 4. Plot broadband spectra at each epoch
plt.figure(figsize=(4.8, 3.6),dpi=200)
colors = plt.cm.viridis(np.linspace(0,1,len(epochs)))

for i, t in enumerate(epochs):
    exp = int(np.floor(np.log10(t)))
    base = t / 10**exp
    plt.loglog(frequencies, results.total[:,i], color=colors[i], label=fr'${base:.1f} \times 10^{{{exp}}}$ s')

# 5. Add vertical lines marking the bands from the light curve plot
for i, band in enumerate(bands):
    exp = int(np.floor(np.log10(band)))
    base = band / 10**exp
    plt.axvline(band,ls='--',color='C'+str(i))

plt.xlabel('frequency (Hz)')
plt.ylabel('flux density (erg/cm²/s/Hz)')
plt.legend(ncol=2)
plt.title('Synchrotron Spectra')
plt.savefig('quick-spec.png',dpi=300,bbox_inches='tight')
```

<div align="center">
<img src="assets/quick-spec.png" alt="Broadband Spectra" width="600"/>

The spectral analysis code will generate this visualization showing spectra at different times, with vertical lines indicating the frequencies calculated in the light curve example.
</div>
</details>

<details>
<summary><b>Time-Frequency Pairs Calculation</b> <i>(click to expand/collapse)</i></summary>
<br>

If you want to calculate flux at specific time-frequency pairs (t_i, nu_i) instead of a grid (t_i, nu_j), you can use the alternative series interfaces:

```python
# Define time and frequency arrays (must be the same length)
times = np.logspace(2, 8, 200)
frequencies = np.logspace(9, 17, 200)

# For time-frequency pairs (times array must be in ascending order)
results = model.flux_density(times, frequencies)

# The returned results is a FluxDict object with different flux components
print("Result attributes:", dir(results))  # Shows .fwd, .rvs, .total attributes
print("Total flux shape:", results.total.shape)  # Same shape as input arrays
print("Forward shock shape:", results.fwd.sync.shape)  # Forward shock synchrotron component
```

**Key differences:**
- `flux_density_grid()`: Calculates flux on a time-frequency grid (NxM output from N times and M frequencies)
- `flux_density()`: Calculates flux at paired time-frequency points (N output from N time-frequency pairs), requires ascending order time arrays
- `flux_density_exposures()`: Same as above but with exposure time averaging for realistic observational scenarios

**Return value structure:**
All flux calculation methods return a `FluxDict` object with:
- `.total`: Combined flux from all components
- `.fwd`: Forward shock flux (has `.sync` and `.ssc` attributes)
- `.rvs`: Reverse shock flux (has `.sync` and `.ssc` attributes)

</details>



### Internal Quantities Evolution

The example below walks through how you can check the evolution of internal quantities under various reference frames via `script/details.ipynb`.

<details>
<summary><b>Model Setup</b> <i>(click to expand/collapse)</i></summary>
<br>

Same as for light curve generation, let's set up the physical components of our afterglow model, including the environment, jet, observer, and radiation parameters:

```python
import numpy as np
import matplotlib.pyplot as plt
from VegasAfterglow import ISM, TophatJet, Observer, Radiation, Model

medium = ISM(n_ism=1)

jet = TophatJet(theta_c=0.3, E_iso=1e52, Gamma0=100)

obs = Observer(lumi_dist=1e26, z=0.1, theta_obs=0.)

rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3)

model = Model(jet=jet, medium=medium, observer=obs, fwd_rad=rad, resolutions=(0.1,5,10))
```

</details>

<details>
<summary><b>Get the simulation quantities</b> <i>(click to expand/collapse)</i></summary>
<br>

Now, let's get the internal simulation quantities:

```python

# Get the simulation details over a time range
details = model.details(t_min=1e0, t_max=1e8)

# Print the available attributes
print("Simulation details attributes:", dir(details))
print("Forward shock attributes:", dir(details.fwd))
```
You will get a `SimulationDetails` object with the following structure:

**Main grid coordinates:**
- `details.phi`: 1D numpy array of azimuthal angles in `radians`
- `details.theta`: 1D numpy array of polar angles in `radians`
- `details.t_src`: 3D numpy array of source frame times on coordinate (phi_i, theta_j, t_k) grid in `seconds`

**Forward shock details (accessed via `details.fwd`):**
- `details.fwd.t_comv`: 3D numpy array of comoving times for the forward shock in `seconds`
- `details.fwd.t_obs`: 3D numpy array of observer times for the forward shock in `seconds`
- `details.fwd.Gamma`: 3D numpy array of downstream Lorentz factors for the forward shock
- `details.fwd.Gamma_th`: 3D numpy array of thermal Lorentz factors for the forward shock
- `details.fwd.r`: 3D numpy array of lab frame radii in `cm`
- `details.fwd.B_comv`: 3D numpy array of downstream comoving magnetic field strengths for the forward shock in `Gauss`
- `details.fwd.theta`: 3D numpy array of polar angles for the forward shock in `radians`
- `details.fwd.N_p`: 3D numpy array of downstream shocked proton number per solid angle for the forward shock
- `details.fwd.N_e`: 3D numpy array of downstream synchrotron electron number per solid angle for the forward shock
- `details.fwd.gamma_a`: 3D numpy array of comoving frame self-absorption Lorentz factors for the forward shock
- `details.fwd.gamma_m`: 3D numpy array of comoving frame injection Lorentz factors for the forward shock
- `details.fwd.gamma_c`: 3D numpy array of comoving frame cooling Lorentz factors for the forward shock
- `details.fwd.gamma_M`: 3D numpy array of comoving frame maximum Lorentz factors for the forward shock
- `details.fwd.nu_a`: 3D numpy array of comoving frame self-absorption frequencies for the forward shock in `Hz`
- `details.fwd.nu_m`: 3D numpy array of comoving frame injection frequencies for the forward shock in `Hz`
- `details.fwd.nu_c`: 3D numpy array of comoving frame cooling frequencies for the forward shock in `Hz`
- `details.fwd.nu_M`: 3D numpy array of comoving frame maximum frequencies for the forward shock in `Hz`
- `details.fwd.I_nu_max`: 3D numpy array of comoving frame synchrotron maximum specific intensities for the forward shock in `erg/cm²/s/Hz`
- `details.fwd.Doppler`: 3D numpy array of Doppler factors for the forward shock

**Reverse shock details (accessed via `details.rvs`, if reverse shock is enabled):**
- Similar attributes as forward shock but for the reverse shock component

</details>

<details>
<summary><b>Checking the evolution of various parameters</b> <i>(click to expand/collapse)</i></summary>
<br>

To analyze the temporal evolution of physical parameters across different reference frames, we can visualize how key quantities evolve in the source, comoving, and observer frames. The following analysis demonstrates the comprehensive tracking of shock dynamics and microphysical parameters throughout the afterglow evolution:

**Multi-parameter evolution visualization:**
This code creates a comprehensive multi-panel figure displaying the temporal evolution of fundamental shock parameters (Lorentz factor, magnetic field, particle numbers, radius, and peak synchrotron power) across all three reference frames:

```python
attrs =['Gamma', 'B_comv', 'N_p','r','N_e','I_nu_max']
ylabels = [r'$\Gamma$', r'$B^\prime$ [G]', r'$N_p$', r'$r$ [cm]', r'$N_e$', r'$I_{\nu, \rm max}^\prime$ [erg/s/Hz]']

frames = ['t_src', 't_comv', 't_obs']
titles = ['source frame', 'comoving frame', 'observer frame']
colors = ['C0', 'C1', 'C2']
xlabels = [r'$t_{\rm src}$ [s]', r'$t^\prime$ [s]', r'$t_{\rm obs}$ [s]']
plt.figure(figsize= (4.2*len(frames), 3*len(attrs)))

#plot the evolution of various parameters for phi = 0 and theta = 0 (so the first two indexes are 0)
for i, frame in enumerate(frames):
    for j, attr in enumerate(attrs):
        plt.subplot(len(attrs), len(frames) , j * len(frames) + i + 1)
        if j == 0:
            plt.title(titles[i])
        value = getattr(details.fwd, attr)
        if frame == 't_src':
            t = getattr(details, frame)
        else:
            t = getattr(details.fwd, frame)
        plt.loglog(t[0, 0, :], value[0, 0, :], color='k',lw=2.5)
        plt.loglog(t[0, 0, :], value[0, 0, :], color=colors[i])

        plt.xlabel(xlabels[i])
        plt.ylabel(ylabels[j])

plt.tight_layout()
plt.savefig('shock_quantities.png', dpi=300,bbox_inches='tight')
```

<div align="center">
<img src="assets/shock_quantities.png" alt="Shock evolution" width="1000"/>
</div>

**Electron energy distribution analysis:**
This visualization focuses specifically on the characteristic electron energies (self-absorption, injection, and cooling) across all three reference frames:

```python
frames = ['t_src', 't_comv', 't_obs']
xlabels = [r'$t_{\rm src}$ [s]', r'$t^\prime$ [s]', r'$t_{\rm obs}$ [s]']
plt.figure(figsize= (4.2*len(frames), 3.6))

for i, frame in enumerate(frames):
    plt.subplot(1, len(frames), i + 1)
    if frame == 't_src':
        t = getattr(details, frame)
    else:
        t = getattr(details.fwd, frame)
    plt.loglog(t[0, 0, :], details.fwd.gamma_a[0, 0, :],label=r'$\gamma_a^\prime$',c='firebrick')
    plt.loglog(t[0, 0, :], details.fwd.gamma_m[0, 0, :],label=r'$\gamma_m^\prime$',c='yellowgreen')
    plt.loglog(t[0, 0, :], details.fwd.gamma_c[0, 0, :],label=r'$\gamma_c^\prime$',c='royalblue')
    plt.loglog(t[0, 0, :], details.fwd.gamma_a[0, 0, :]*details.fwd.Doppler[0,0,:],label=r'$\gamma_a$',ls='--',c='firebrick')
    plt.loglog(t[0, 0, :], details.fwd.gamma_m[0, 0, :]*details.fwd.Doppler[0,0,:],label=r'$\gamma_m$',ls='--',c='yellowgreen')
    plt.loglog(t[0, 0, :], details.fwd.gamma_c[0, 0, :]*details.fwd.Doppler[0,0,:],label=r'$\gamma_c$',ls='--',c='royalblue')
    plt.xlabel(xlabels[i])
    plt.ylabel(r'$\gamma_e^\prime$')
    plt.legend(ncol=2)
plt.tight_layout()
plt.savefig('electron_quantities.png', dpi=300,bbox_inches='tight')
```

<div align="center">
<img src="assets/electron_quantities.png" alt="Shock evolution" width="1000"/>
</div>

**Synchrotron frequency evolution:**
This analysis tracks the evolution of characteristic synchrotron frequencies:

```python
frames = ['t_src', 't_comv', 't_obs']
xlabels = [r'$t_{\rm src}$ [s]', r'$t^\prime$ [s]', r'$t_{\rm obs}$ [s]']
plt.figure(figsize= (4.2*len(frames), 3.6))

for i, frame in enumerate(frames):
    plt.subplot(1, len(frames), i + 1)
    if frame == 't_src':
        t = getattr(details, frame)
    else:
        t = getattr(details.fwd, frame)
    plt.loglog(t[0, 0, :], details.fwd.nu_a[0, 0, :],label=r'$\nu_a^\prime$',c='firebrick')
    plt.loglog(t[0, 0, :], details.fwd.nu_m[0, 0, :],label=r'$\nu_m^\prime$',c='yellowgreen')
    plt.loglog(t[0, 0, :], details.fwd.nu_c[0, 0, :],label=r'$\nu_c^\prime$',c='royalblue')
    plt.loglog(t[0, 0, :], details.fwd.nu_a[0, 0, :]*details.fwd.Doppler[0,0,:],label=r'$\nu_a$',ls='--',c='firebrick')
    plt.loglog(t[0, 0, :], details.fwd.nu_m[0, 0, :]*details.fwd.Doppler[0,0,:],label=r'$\nu_m$',ls='--',c='yellowgreen')
    plt.loglog(t[0, 0, :], details.fwd.nu_c[0, 0, :]*details.fwd.Doppler[0,0,:],label=r'$\nu_c$',ls='--',c='royalblue')
    plt.xlabel(xlabels[i])
    plt.ylabel(r'$\nu$ [Hz]')
    plt.legend(ncol=2)
plt.tight_layout()
plt.savefig('photon_quantities.png', dpi=300,bbox_inches='tight')
```
<div align="center">
<img src="assets/photon_quantities.png" alt="Shock evolution" width="1000"/>
</div>

**Doppler factor spatial distribution:**
This polar plot visualizes the spatial distribution of the Doppler factor across the jet structure, showing how relativistic beaming varies with angular position and radial distance:

```python
plt.figure(figsize=(6,6))
ax = plt.subplot(111, polar=True)

theta = details.fwd.theta[0,:,:]
r     = details.fwd.r[0,:,:]
D     = details.fwd.Doppler[0,:,:]

# Polar contour plot
scale = 3.0
c = ax.contourf(theta*scale, r, np.log10(D), levels=30, cmap='viridis')

ax.set_rscale('log')
true_ticks = np.linspace(0, 0.3, 6)
ax.set_xticks(true_ticks * scale)
ax.set_xticklabels([f"{t:.2f}" for t in true_ticks])
ax.set_xlim(0,0.3*scale)
ax.set_ylabel(r'$\theta$ [rad]')
ax.set_xlabel(r'$r$ [cm]')

plt.colorbar(c, ax=ax, label=r'$\log_{10} D$')

plt.tight_layout()
plt.savefig('doppler.png', dpi=300,bbox_inches='tight')
```
<div align="center">
<img src="assets/doppler.png" alt="Shock evolution" width="600"/>
</div>

**Equal arrival time (EAT) surface visualization:**
This final visualization maps the equal arrival time surfaces in polar coordinates, illustrating how light from different parts of the jet reaches the observer at the same time:

```python
plt.figure(figsize=(6,6))
ax = plt.subplot(111, polar=True)

theta = details.fwd.theta[0,:,:]
r     = details.fwd.r[0,:,:]
t_obs = details.fwd.t_obs[0,:,:]

scale = 3.0
c = ax.contourf(theta*scale, r, np.log10(t_obs), levels=30, cmap='viridis')

ax.set_rscale('log')
true_ticks = np.linspace(0, 0.3, 6)
ax.set_xticks(true_ticks * scale)
ax.set_xticklabels([f"{t:.2f}" for t in true_ticks])
ax.set_xlim(0,0.3*scale)
ax.set_ylabel(r'$\theta$ [rad]')
ax.set_xlabel(r'$r$ [cm]')

plt.colorbar(c, ax=ax, label=r'$\log_{10} (t_{\rm obs}/s)$')

plt.tight_layout()
plt.savefig('EAT.png', dpi=300,bbox_inches='tight')
```
<div align="center">
<img src="assets/EAT.png" alt="Shock evolution" width="600"/>

</div>
</details>


### MCMC Parameter Fitting

We provide some example data files in the `data` folder. Remember to keep your copy in the same directory as the original to ensure all data paths work correctly.

<details>
<summary><b>1. Preparing Data and Configuring the Model</b> <i>(click to expand/collapse)</i></summary>
<br>

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import corner
from VegasAfterglow import ObsData, Setups, Fitter, ParamDef, Scale
```

VegasAfterglow provides flexible options for loading observational data through the `ObsData` class. You can add light curves (specific flux vs. time) and spectra (specific flux vs. frequency) in multiple ways.

```python
# Create an instance to store observational data
data = ObsData()

# Method 1: Add data directly from lists or numpy arrays

# For light curves
t_data = [1e3, 2e3, 5e3, 1e4, 2e4]  # Time in seconds
flux_data = [1e-26, 8e-27, 5e-27, 3e-27, 2e-27]  # Specific flux in erg/cm²/s/Hz
flux_err = [1e-28, 8e-28, 5e-28, 3e-28, 2e-28]  # Specific flux error in erg/cm²/s/Hz
data.add_flux_density(nu=4.84e14, t=t_data, f_nu=flux_data, err=flux_err)  # All quantities in CGS units
# You can also assign weights to each data point to account for systematic uncertainties or correlations. You don't need to worry about the weights' normalization, the code will normalize them automatically.
#data.add_flux_density(nu=4.84e14, t=t_data, f_nu=flux_data, err=flux_err, weights=np.ones(len(t_data)))

# For spectra
nu_data = [...]  # Frequencies in Hz
spectrum_data = [...] # Specific flux values in erg/cm²/s/Hz
spectrum_err = [...]   # Specific flux errors in erg/cm²/s/Hz
data.add_spectrum(t=3000, nu=nu_data, f_nu=spectrum_data, err=spectrum_err, weights=np.ones(len(nu_data)))  # All quantities in CGS units
```

```python
# Method 2: Load from CSV files

data = ObsData()
# Define your bands and files
bands = [2.4e17, 4.84e14, 1.4e14]  # Example: X-ray, optical R-band
lc_files = ["data/ep.csv", "data/r.csv", "data/vt-r.csv"]

# Load light curves from files
for nu, fname in zip(bands, lc_files):
    df = pd.read_csv(fname)
    data.add_flux_density(nu=nu, t=df["t"], f_nu=df["Fv_obs"], err=df["Fv_err"])  # All quantities in CGS units

times = [3000] # Example: time in seconds
spec_files = ["data/ep-spec.csv"]

# Load spectra from files
for t, fname in zip(times, spec_files):
    df = pd.read_csv(fname)
    data.add_spectrum(t=t, nu=df["nu"], f_nu=df["Fv_obs"], err=df["Fv_err"])  # All quantities in CGS units
```

> **Note:** The `ObsData` interface is designed to be flexible. You can mix and match different data sources, and add multiple light curves at different frequencies as well as multiple spectra at different times.

The `Setups` class defines the global properties and environment for your model. These settings remain fixed during the MCMC process. Check the [documentation](https://yihanwangastro.github.io/VegasAfterglow/docs/index.html) for all available options.

```python
cfg = Setups()

# Source properties
cfg.lumi_dist = 3.364e28    # Luminosity distance [cm]
cfg.z = 1.58               # Redshift

# Physical model configuration
cfg.medium = "wind"        # Ambient medium: "wind", "ism", etc. (see documentation)
cfg.jet = "powerlaw"       # Jet structure: "powerlaw", "gaussian", "tophat", etc. (see documentation)
```

These settings affect how the model is calculated but are not varied during the MCMC process.
</details>

<details>
<summary><b>2. Defining Parameters and Running MCMC</b> <i>(click to expand/collapse)</i></summary>
<br>

The `ParamDef` class is used to define the parameters for MCMC exploration. Each parameter requires a name, prior range, and sampling scale:

```python
mc_params = [
    ParamDef("E_iso",      1e50,  1e54,  Scale.LOG),       # Isotropic energy [erg]
    ParamDef("Gamma0",        5,  1000,  Scale.LOG),       # Lorentz factor at the core
    ParamDef("theta_c",     0.0,   0.5,  Scale.LINEAR),    # Core half-opening angle [rad]
    ParamDef("k_e",           2,     2,  Scale.FIXED),     # Energy power law index
    ParamDef("k_g",           2,     2,  Scale.FIXED),     # Lorentz factor power law index
    ParamDef("theta_v",     0.0,   0.0,  Scale.FIXED),     # Viewing angle [rad]
    ParamDef("p",             2,     3,  Scale.LINEAR),    # Shocked electron power law index
    ParamDef("eps_e",      1e-2,   0.5,  Scale.LOG),       # Electron energy fraction
    ParamDef("eps_B",      1e-4,   0.5,  Scale.LOG),       # Magnetic field energy fraction
    ParamDef("A_star",     1e-3,     1,  Scale.LOG),       # Wind parameter
    ParamDef("xi_e",       1e-3,     1,  Scale.LOG),       # Electron acceleration fraction
]
```

**Scale Types:**

- `Scale.LOG`: Sample in logarithmic space (log10) - ideal for parameters spanning multiple orders of magnitude
- `Scale.LINEAR`: Sample in linear space - appropriate for parameters with narrower ranges
- `Scale.FIXED`: Keep parameter fixed at the initial value - use for parameters you don't want to vary

**Parameter Choices:**
The parameters you include depend on your model configuration (See documentation for all options):

- For "wind" medium: use `A_star` parameter
- For "ISM" medium: use `n_ism` parameter instead
- Different jet structures may require different parameters

Initialize the `Fitter` class with your data and configuration, then run the MCMC process:

```python
# Create the fitter object
fitter = Fitter(data, cfg)

# Option 1: MCMC with emcee (faster, recommended for quick fitting)
result = fitter.fit(
    mc_params,
    resolution=(0.3, 1, 10),       # Grid resolution (phi, theta, t)
    sampler="emcee",               # MCMC sampler
    nsteps=10000,                  # Number of steps per walker
    nburn=1000,                    # Burn-in steps to discard
    thin=1,                        # Save every nth sample
    npool=8,                       # Number of parallel processes
    top_k=10,                      # Number of best-fit parameters to return
    outdir="bilby_output",         # Output directory (default)
    label="afterglow_fit",         # Run label (default: "afterglow")
)

# Option 2: Nested sampling with dynesty (slower but computes Bayesian evidence)
result = fitter.fit(
    mc_params,
    resolution=(0.3, 1, 10),       # Grid resolution (phi, theta, t)
    sampler="dynesty",             # Nested sampling algorithm
    nlive=500,                     # Number of live points
    dlogz=0.1,                     # Stopping criterion (evidence tolerance)
    sample="rwalk",                # Sampling method
    npool=8,                       # Number of parallel processes
    top_k=10,                      # Number of best-fit parameters to return
    outdir="bilby_output",         # Output directory (default)
    label="afterglow_fit",         # Run label (default: "afterglow")
)
```

**Important Notes:**
- Parameters with `Scale.LOG` are sampled as `log10_<name>` (e.g., `log10_E_iso`)
- The sampler works in log10 space for LOG-scale parameters, then transforms back
- Use `npool` to parallelize likelihood evaluations across multiple CPU cores

The `result` object contains:

- `samples`: The posterior samples (shape: [n_samples, 1, n_params])
- `labels`: Parameter names (with `log10_` prefix for LOG-scale params)
- `latex_labels`: LaTeX-formatted labels for plotting (e.g., `$\log_{10}(E_{\rm iso})$`)
- `top_k_params`: Top-k maximum likelihood parameter values
- `top_k_log_probs`: Log probabilities for top-k samples
- `bilby_result`: Full bilby Result object (for advanced diagnostics)

</details>

<details>
<summary><b>3. Analyzing Results and Generating Predictions</b> <i>(click to expand/collapse)</i></summary>
<br>

Check the top-k parameters and their uncertainties:

```python
# Print top-k parameters (maximum likelihood)
top_k_data = []
for i in range(result.top_k_params.shape[0]):
    row = {'Rank': i+1, 'chi^2': f"{-2*result.top_k_log_probs[i]:.2f}"}
    for name, val in zip(result.labels, result.top_k_params[i]):
        row[name] = f"{val:.4f}"
    top_k_data.append(row)

top_k_df = pd.DataFrame(top_k_data)
print("Top-k parameters:")
print(top_k_df.to_string(index=False))
```

Use the best-fit parameters to generate model predictions

```python
# Define time and frequency ranges for predictions
t_out = np.logspace(2, 9, 150)

nu_out = np.logspace(16,20,150)

best_params = result.top_k_params[0]

# Generate model light curves at the specified bands using the best-fit parameters
lc = fitter.flux_density_grid(best_params, t_out, band)

# Generate model spectra at the specified times using the best-fit parameters
spec = fitter.flux_density_grid(best_params, times, nu_out)
```

Now you can plot the best-fit model:

```python
# Function to plot model light curves along with observed data
def draw_bestfit(t,lc_fit, nu, spec_fit):
    fig =plt.figure(figsize=(4.5, 7.5))

    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    shift = [1,1,200]
    colors = ['blue', 'orange', 'green']
    for i, file, sft, c in zip(range(len(lc_files)), lc_files, shift, colors ):
        df = pd.read_csv(file)
        ax1.errorbar(df["t"], df["Fv_obs"]*sft, df["Fv_err"]*sft, fmt='o',markersize=4,label=file, color=c,markeredgecolor='k', markeredgewidth=.4)
        ax1.plot(t, np.array(lc_fit[i,:])*sft, color=c,lw=1)

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('t [s]')
    ax1.set_ylabel(r'$F_\nu$ [erg/cm$^2$/s/Hz]')
    ax1.legend()

    for i, file, sft, c in zip(range(len(spec_files)), spec_files, shift, colors ):
        df = pd.read_csv(file)
        ax2.errorbar(df["nu"], df["Fv_obs"]*sft, df["Fv_err"]*sft, fmt='o',markersize=4,label=file, color=c,markeredgecolor='k', markeredgewidth=.4)
        ax2.plot(nu, np.array(spec_fit[:,i])*sft, color=c,lw=1)

    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel(r'$\nu$ [Hz]')
    ax2.set_ylabel(r'$F_\nu$ [erg/cm$^2$/s/Hz]')
    ax2.legend()
    plt.tight_layout()

draw_bestfit(t_out, lc, nu_out, spec)
```

Corner plots are essential for visualizing parameter correlations and posterior distributions:

```python
def plot_corner(flat_chain, labels, filename="corner_plot.png"):
    fig = corner.corner(
        flat_chain,
        labels=labels,
        quantiles=[0.16, 0.5, 0.84],  # For median and ±1σ
        show_titles=True,
        title_kwargs={"fontsize": 14},
        label_kwargs={"fontsize": 14},
        truths=np.median(flat_chain, axis=0),  # Show median values
        truth_color='red',
        bins=30,
        smooth=1,
        fill_contours=True,
        levels=[0.16, 0.5, 0.68],  # 1σ and 2σ contours
        color='k'
    )
    fig.savefig(filename, dpi=300, bbox_inches='tight')

# Create the corner plot
flat_chain = result.samples.reshape(-1, result.samples.shape[-1])
plot_corner(flat_chain, result.latex_labels)
```

</details>

---

## Documentation

Comprehensive documentation is available at **[Documentation](https://yihanwangastro.github.io/VegasAfterglow/docs/index.html)** including:

- **Installation Guide**: Detailed instructions for setting up VegasAfterglow
- **Examples**: Practical examples showing common use cases
- **Python API Reference**: Complete documentation of the Python interface
- **C++ API Reference**: Detailed documentation of C++ classes and functions
- **Contributing Guide**: Information for developers who wish to contribute

The documentation is regularly updated with the latest features and improvements (not yet officially released).

For a complete history of changes and new features, see our [**Changelog**](CHANGELOG.md).

---

## Contributing

If you encounter any issues, have questions about the code, or want to request new features:

1. **GitHub Issues** - The most straightforward and fastest way to get help:
   - Open an issue at [Issues](https://github.com/YihanWangAstro/VegasAfterglow/issues)
   - You can report bugs, suggest features, or ask questions
   - This allows other users to see the problem/solution as well
   - Can be done anonymously if preferred

2. **Pull Requests** - If you've implemented a fix or feature:
   - Fork the repository
   - Create a branch for your changes
   - Submit a pull request with your changes

3. **Email** - For private questions or discussions:
   - Contact the maintainers directly via email
   - Yihan Wang: yihan.astro@gmail.com
   - Bing Zhang: bing.zhang@unlv.edu
   - Connery Chen: connery.chen@unlv.edu

We value all contributions and aim to respond to issues promptly. All communications are extremely welcome!

---

## License

VegasAfterglow is released under the **BSD-3-Clause License**.

The BSD 3-Clause License is a permissive open source license that allows you to:

- Freely use, modify, and distribute the software in source and binary forms
- Use the software for commercial purposes
- Integrate the software into proprietary applications

**Requirements:**

- You must include the original copyright notice and the license text
- You cannot use the names of the authors or contributors to endorse derived products
- The license provides no warranty or liability protection

For the full license text, see the [LICENSE](LICENSE) file in the repository.

---

## Acknowledgments & Citation

We would like to thank the contributors who helped improve VegasAfterglow. **Special thanks to Weihua Lei, Shaoyu Fu, Liang-Jun Chen, Iris Yin, Cuiyuan Dai and Binbin Zhang** for their invaluable work as beta testers, providing feedback and helping with bug fixes during development. We also thank the broader community for their suggestions and support.

If you use VegasAfterglow in your research, please cite the relevant paper(s):

[https://ui.adsabs.harvard.edu/abs/2025arXiv250710829W/abstract](https://ui.adsabs.harvard.edu/abs/2025arXiv250710829W/abstract)
