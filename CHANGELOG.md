# Changelog

All notable changes to VegasAfterglow will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.1.0] - 2025-12-08

### Added

#### **Bilby MCMC Integration**
- **Multi-Sampler Support**: Full integration with bilby framework for Bayesian inference
  - `emcee`: Affine-invariant MCMC ensemble sampler (fast but does not compute evidence)
  - `dynesty`: Dynamic nested sampling (slow but computes Bayesian evidence)
  - Support for all bilby samplers (`nestle`, `cpnest`, `pymultinest`, `ultranest`, etc.)
- **Enhanced Result Object**: `FitResult` now includes LaTeX-formatted labels
  - `latex_labels`: Properly formatted labels for corner plots (e.g., `$\log_{10}(E_{\rm iso})$`)
  - Automatic formatting for LOG-scale parameters with `log10_` prefix
- **Parallelization**: Multi-core support via `npool` parameter
  - Efficient parallel likelihood evaluations across CPU cores
  - Works with both emcee and dynesty samplers

#### **Flexible Parameter Interface**
- **LOG-Scale Parameter Naming**: Parameters with `Scale.LOG` automatically prefixed with `log10_`
  - Example: `E_iso` with `Scale.LOG` becomes `log10_E_iso` in sampler space
  - Automatic transformation between log10 and physical values
- **Generic Sampler Configuration**: Pass sampler-specific parameters via `**sampler_kwargs`
  - Emcee defaults: `nsteps=5000`, `nburn=1000`, `thin=1`, `nwalkers=2*ndim`
  - Dynesty defaults: `nlive=500`, `dlogz=0.1`, `sample="rwalk"`
  - Extensible to any bilby-supported sampler

### Documentation

#### **Comprehensive MCMC Guides**
- **README Updates**: Complete bilby integration examples
  - Emcee positioned as primary option (Option 1)
  - Dynesty as secondary option (Option 2) for evidence calculation
  - All corner plot examples updated to use `result.latex_labels`
- **Enhanced RST Documentation**:
  - `mcmc_fitting.rst`: Complete `Fitter.fit()` interface reference with all parameters
  - `quickstart.rst`: Updated with current sampler options and examples
  - `python_api.rst`: Refreshed with bilby integration examples
  - References to bilby documentation for additional samplers

---

## [v1.0.3] - 2025-09-29

### Changed

#### **Default Physics Settings**
- **Self-Absorption Heating**: Disabled self-absorption heating by default for improved numerical stability
- **Strong Absorption Subsegmentation**: Enhanced subsegmentation (5/2 ratio) for improved handling of strong self-absorption regions

#### **MCMC Framework Improvements**
- **Enhanced Move Strategy**: Implemented new MCMC move strategy for better parameter space exploration
- **Wider Initial Spread**: Increased initial parameter spread for more robust sampling
- **General Medium Support**: Extended `-k` parameter support for general medium configurations

### Fixed

#### **Model Interface**
- **Flux Unit Consistency**: Fixed model flux unit handling for consistent calculations across interfaces
- **Reverse Shock Dynamics**: Implemented hardcut Gamma treatment for more stable reverse shock evolution

### Performance

#### **Inverse Compton Optimization**
- **Early Break Implementation**: Added early break condition for inverse Compton calculations to improve computational efficiency

### Documentation

#### **Enhanced Guidelines**
- **MCMC Fitting Guidelines**: Added comprehensive guidelines for MCMC parameter fitting
- **Interface Documentation**: Updated documentation for interface consistency and parameter usage
- **Code Examples**: Enhanced examples and troubleshooting documentation

### Development

#### **Code Quality**
- **Interface Consistency**: Major cleanup for interface name consistency across Python bindings
- **Code Cleanup**: General code cleanup and optimization throughout the codebase

---

## [v1.0.2] - 2025-09-15

### Changed

#### **⚠️ BREAKING: Python Interface Method Name Updates**
- **Standardized Method Names**: Updated method names across Python interfaces for consistency and clarity
  - **PyModel Methods**:
    - `specific_flux()` → `flux_density_grid()` (multi-dimensional flux density calculations)
    - `specific_flux_series()` → `flux_density()` (time-series flux density calculations)
    - `specific_flux_series_with_expo()` → `flux_density_exposures()` (exposure-averaged flux density)
  - **VegasMC Methods**:
    - `specific_flux()` → `flux_density_grid()` (matches PyModel interface)

#### **⚠️ BREAKING: Data Input Parameter Name Cleanup**
- **ObsData Interface**: Simplified parameter names by removing explicit "_cgs" suffixes for cleaner API
  - **Method**: `add_flux_density(nu, t, f_nu, err, weights=None)`
    - Previous: `add_flux_density(nu_cgs, t_cgs, Fnu_cgs, Fnu_err, weights=None)`
    - Updated: `add_flux_density(nu, t, f_nu, err, weights=None)`
  - **Method**: `add_flux(nu_min, nu_max, num_points, t, flux, err, weights=None)`
    - Previous: `add_flux(nu_min, nu_max, num_points, t_cgs, F, F_err, weights=None)`
    - Updated: `add_flux(nu_min, nu_max, num_points, t, flux, err, weights=None)`
  - **Method**: `add_spectrum(t, nu, f_nu, err, weights=None)`
    - Previous: `add_spectrum(t_cgs, nu_cgs, Fnu_cgs, Fnu_err, weights=None)`
    - Updated: `add_spectrum(t, nu, f_nu, err, weights=None)`

### Added

#### **Enhanced Documentation**
- **CGS Unit Documentation**: Added comprehensive CGS unit documentation across all interfaces
  - Added docstrings to pybind11 methods specifying CGS units for all physical quantities
  - Updated all examples in README, documentation, and notebooks with clear CGS unit comments
  - Consistent unit documentation: `nu [Hz]`, `t [s]`, `f_nu [erg/cm²/s/Hz]`, `flux [erg/cm²/s]`

#### **Interface Consistency Verification**
- **Parameter Binding Validation**: Verified all pybind11 interface parameters match documentation
  - All ObsData methods confirmed to use simplified parameter names
  - Method names verified to be consistent across all interfaces
  - Complete synchronization between C++ bindings and Python documentation

### Migration Notes

**Method Name Changes**: Update your method calls as follows:

**For PyModel/Model objects:**
```python
# Old method names → New method names
results = model.specific_flux(t, nu)              → results = model.flux_density_grid(t, nu)
flux = model.specific_flux_series(t, nu)          → flux = model.flux_density(t, nu)
flux = model.specific_flux_series_with_expo(...)  → flux = model.flux_density_exposures(...)
```

**For VegasMC objects:**
```python
# Old method names → New method names
results = vegasmc.specific_flux(params, t, nu)    → results = vegasmc.flux_density_grid(params, t, nu)
```

**For ObsData methods:**
```python
# Old parameter names → New parameter names
data.add_flux_density(nu_cgs=1e14, t_cgs=times, Fnu_cgs=flux, Fnu_err=errors)
# becomes:
data.add_flux_density(nu=1e14, t=times, f_nu=flux, err=errors)  # All quantities in CGS units

data.add_spectrum(t_cgs=1000, nu_cgs=freqs, Fnu_cgs=spectrum, Fnu_err=errors)
# becomes:
data.add_spectrum(t=1000, nu=freqs, f_nu=spectrum, err=errors)  # All quantities in CGS units

data.add_flux(nu_min=1e14, nu_max=1e15, num_points=5, t_cgs=times, F=flux, F_err=errors)
# becomes:
data.add_flux(nu_min=1e14, nu_max=1e15, num_points=5, t=times, flux=flux, err=errors)  # All quantities in CGS units
```

**Important**: All physical quantities are still expected in CGS units as before - only the parameter names have been simplified for a cleaner interface.

---

## [v1.0.1] - 2025-09-15

### Added

#### **Enhanced MCMC Capabilities**
- **Frequency Integrated Flux Support**: Added support for frequency-integrated observations in MCMC fitting framework
  - New `MultiBandData.add_light_curve()` overload for broadband flux measurements over frequency ranges
  - Enables modeling of bolometric observations and filter-integrated measurements
  - Improved handling of observational data where effective frequency depends on spectral shape

#### **Documentation & Examples**
- **Comprehensive API Documentation**: Added extensive documentation with detailed physics explanations
  - Complete parameter descriptions with units and typical ranges
  - Detailed method documentation with usage examples
  - Enhanced code examples demonstrating advanced features

### Changed

#### **⚠️ BREAKING: Return Object Interface Redesign**
- **Named Member Access**: Replaced dictionary-style access with structured object interfaces for better usability and IDE support
  - **specific_flux()** methods now return `FluxDict` objects with named members:
    - `results.total` (total flux array)
    - `results.fwd.sync` (forward shock synchrotron)
    - `results.fwd.ssc` (forward shock SSC)
    - `results.rvs.sync` (reverse shock synchrotron)
    - `results.rvs.ssc` (reverse shock SSC)
  - **details()** method now returns `SimulationDetails` objects with named members:
    - `details.phi`, `details.theta`, `details.t_src` (coordinate arrays)
    - `details.fwd.*` (forward shock quantities: `t_obs`, `Gamma`, `r`, `B_comv`, etc.)
    - `details.rvs.*` (reverse shock quantities, if enabled)
  - **Enhanced Type Safety**: All return objects have well-defined attributes for better development experience

#### **API Consistency Improvements**
- **Parameter Naming Standardization**: Major cleanup of parameter names across all interfaces for consistency
  - **ConfigParams**: `fwd_SSC` → `fwd_ssc`, `rvs_SSC` → `rvs_ssc`, `IC_cooling` → `ssc_cooling`, `KN` → `kn`
  - **PyRadiation**: `IC_cooling` → `ssc_cooling`, `SSC` → `ssc`, `KN` → `kn`
  - **Model Constructor**: `forward_rad` → `fwd_rad`, `reverse_rad` → `rvs_rad`
  - All parameter names now follow the consistent snake_case convention

#### **Enhanced Code Documentation**
- **Comprehensive Comments**: Added detailed documentation to all major classes and structures
  - Complete parameter descriptions with physics context
  - Method documentation with detailed explanations
  - Consistent documentation style following project conventions

### Fixed

#### **Interface Consistency**
- **Documentation Updates**: Updated all documentation and examples to reflect new parameter names
  - README.md examples updated with new parameter conventions
  - Tutorial documentation fully synchronized with API changes
  - All code examples verified for consistency

### Migration Notes

**Parameter Name Changes**: This release standardizes parameter naming conventions. Update your code as follows:

**For MCMC Configuration (ConfigParams/Setups):**
```python
# Old names → New names
cfg.fwd_SSC = True    → cfg.fwd_ssc = True
cfg.rvs_SSC = True    → cfg.rvs_ssc = True
cfg.IC_cooling = True → cfg.ssc_cooling = True
cfg.KN = True         → cfg.kn = True
```

**For Radiation Physics (PyRadiation/Radiation):**
```python
# Old names → New names
Radiation(eps_e=0.1, eps_B=0.01, p=2.3, IC_cooling=True, SSC=True, KN=True)
# becomes:
Radiation(eps_e=0.1, eps_B=0.01, p=2.3, ssc_cooling=True, ssc=True, kn=True)
```

**For Model Constructor:**
```python
# Old names → New names
Model(jet=jet, medium=medium, observer=obs, forward_rad=rad, reverse_rad=rad_rvs)
# becomes:
Model(jet=jet, medium=medium, observer=obs, fwd_rad=rad, rvs_rad=rad_rvs)
```

**For Return Object Access (BREAKING CHANGE):**
```python
# OLD: Dictionary-style access (no longer supported)
results = model.specific_flux(times, frequencies)
sync_flux = results['sync']  #  No longer works

# NEW: Named member access
results = model.specific_flux(times, frequencies)
total_flux = results.total              # New way
sync_flux = results.fwd.sync           # New way
ssc_flux = results.fwd.ssc             # New way
rvs_sync = results.rvs.sync            # New way (if reverse shock enabled)

# OLD: Dictionary-style details access
details = model.details(t_min, t_max)
gamma = details['Gamma_fwd']     # No longer works

# NEW: Named member access
details = model.details(t_min, t_max)
gamma = details.fwd.Gamma              # New way
time_obs = details.fwd.t_obs           # New way
coordinates = details.theta            # New way
```

---

## [v1.0.0] - 2025-09-06

### Major Features & Enhancements

#### **Advanced MCMC Framework**
- **More model availability**: All built-in jet and medium models now supported in MCMC fitting, including reverse shock, inverse Compton, magnetar injection and etc. See documentations for detials.

#### **Adaptive Mesh Generation**
- **Dynamic Grid Optimization**: New adaptive angular grid generation based on jet properties and viewing angles. Grid points distributed according to Doppler boosting factors for optimal efficiency
- **Performance Gains**: ~5x faster convergence for reverse shocks.

#### **New Jet Models**
- **StepPowerLawJet**: Uniform core with sharp transition to power-law wings for realistic jet structures
- **Enhanced TwoComponentJet**: Separate narrow and wide components with independent energy and Lorentz factor profiles
- **Improved PowerLawJet**: Split power-law indices for energy (`k_e`) and Lorentz factor (`k_g`) angular dependence

### **Performance & Computational Improvements**

#### **Shock Physics Enhancements**
- **Variable Naming Standardization**: Major refactoring for clarity (`EAT_fwd` → `t_obs_fwd`, `Gamma_rel` → `Gamma_th`)
- **Reverse Shock Optimization**: Major code refactoring for reverse shock dynamics. A unified model for shock crossing and post-crossing evolution.
- **Better Crossing Dynamics**: Track the internal energy evolution during shock crossing for improved accuracy with more accurate shock heating and adiabatic cooling. **Note: The enhanced adiabatic cooling and detailed shock heating treatment leads to even weaker reverse shock emission compared to previous versions.**

#### **Numerical & Memory Optimizations**
- **Memory Efficiency**: Reduced memory footprint through optimized array operations and memory access patterns
- **Grid Pre-computation**: Enhanced caching strategies for frequently used calculations

###  **Enhanced Python Interface**

#### **Data Management Improvements**
- **MultiBandData Redesign**: Unified handling of light curves and spectra with flexible weighting
- **Series Calculations**: New methods for calculating flux at specific time-frequency pairs
- **Memory-Efficient Storage**: Optimized data structures for large multi-wavelength datasets

###  **Development & Build System**

#### **Code Quality Enhancements**
- **Enhanced Template System**: Improved compile-time type checking and memory management
- **Better Documentation**: Comprehensive API documentation with detailed parameter descriptions
- **Updated Examples**: New MCMC tutorials with real data fitting demonstrations

### **API Changes & Migration**

#### **Parameter Interface Updates**
- **PowerLawJet**: Split single `k` parameter into separate `k_e` (energy) and `k_g` (Lorentz factor) indices for more flexible modeling
- **TwoComponentJet**: Parameter names standardized to `(theta_c, E_iso, Gamma0, theta_w, E_iso_w, Gamma0_w)` for core and wide components
- **StepPowerLawJet**: New jet model with parameters `(theta_c, E_iso, Gamma0, E_iso_w, Gamma0_w, k_e, k_g)` for core and power-law wing components
- **Wind Medium**: Extended with optional `n_ism` and `n_0` parameters for stratified medium modeling: `Wind(A_star, n_ism=0, n_0=inf)`
- **Medium Class**: Simplified from `Medium(rho, mass)` to `Medium(rho)` - removed separate mass parameter
- **Model Methods**:
  - Removed `specific_flux_sorted_series()` method
  - Added `specific_flux_series_with_expo(t, nu, expo_time, num_points=10)` for exposure time averaging
  - Changed `details(t_obs)` to `details(t_min, t_max)` interface
- **Resolution Parameters**: Default resolution changed from `(0.3, 3.0, 5.0)` to `(0.3, 1, 10)` for optimal performance/accuracy balance
- **MCMC Parameters**: Major restructuring with new parameters for all jet types, medium configurations, and reverse shock physics
  - Added: `theta_v, n_ism, n0, A_star, k_e, k_g, duration, E_iso_w, Gamma0_w, theta_w, L0, t0, q`
  - Added reverse shock parameters: `p_r, eps_e_r, eps_B_r, xi_e_r`
  - Removed: `k_jet` parameter (replaced by `k_e, k_g`)
- **MCMC Data Handling**: Enhanced `MultiBandData` with optional `weights` parameter for both light curves and spectra
- **MCMC Model Interface**: Replaced separate `light_curves()` and `spectra()` methods with unified `specific_flux()` method

---

**Migration Notes**: This release includes some API changes. Most existing code will work with minimal modifications. See the documentation for detailed migration guidance.

**⚠️ Important Physics Changes**: With the improved reverse shock dynamics, we now find that the reverse shock emission is even weaker than what we reported in our previous code paper, which itself was already weaker than the analytical scalings. This further reduction is due to enhanced adiabatic cooling and a more detailed treatment of shock heating. Users fitting reverse shock data may therefore need to re-evaluate their models and parameters.


## [v0.2.8] - 2025-08-15

### Improved

- **Inverse Compton Performance**: Major performance optimization (~10x speedup) for inverse Compton scattering calculations
- **Photon Interface**: Changed parameter name from `P_nu_max` to `I_nu_max` for better consistency in photon interface
- **Spectrum Smoothing**: Enhanced spectrum smoothing at `nu_max` for better numerical stability

### Documentation

- Updated documentation for detailed simulation quantities evolution
- Enhanced examples and API documentation

## [v0.2.7] - 2025-07-31

### Added

- **Internal Quantities Evolution Interface**: New Python interface to check the evolution of internal simulation quantities under various reference frames
- Comprehensive documentation for detailed simulation quantities evolution

### Fixed

- **Performance Issue**: Fixed jet edge detection function that was mistakenly set to 90 degrees, making the computational domain unnecessarily large regardless of the jet profile

### Documentation

- Updated documentation for detailed simulation quantities evolution
- Enhanced examples showing how to track shock dynamics and microphysical parameters

## [v0.2.6] - 2025-07-19

### Added

- **Built-in Two-Component Jet**: Added support for built-in two-component jet configurations

## [v0.2.5] - 2025-07-19

### Improved

- **Jet Edge Detection**: Improved jet edge detection algorithm for user-defined jet profiles

### Documentation

- Updated README with latest features
- Enhanced documentation with additional examples
- Updated API documentation

## [v0.2.4] - 2025-06-28

### Added

- **Magnetar Spin-down Documentation**: Added comprehensive documentation for magnetar spin-down energy injection

### Documentation

- Updated documentation for magnetar energy injection features
- Enhanced README with additional usage examples
- Improved API reference documentation

## [v0.2.3] - 2025-06-27

### Changed

- **API Breaking Change**: Modified Python-level user-defined medium/jet unit system for consistency
- **Parameter Naming**: Changed jet duration parameter name from `T0` (confusing to observer frame) to `duration`

### Improved

- **Data Handling**: Sorted series flux density for better data organization

## [v0.2.2] - 2025-06-23

### Changed

- **Default Resolution**: Changed default resolution settings when unspecified for better performance

### Improved

- **Reverse Shock**: Enhanced reverse shock smoothing algorithms

## [v0.2.1] - 2025-06-22

### Improved

- **Reverse Shock Modeling**:
  - Enhanced reverse shock smoothing
  - Improved post-crossing reverse shock calibration
  - Better single electron peak power calibration

## [v0.1.9] - 2025-06-19

### Changed

- **Python Support**: Removed support for Python 3.7 (minimum requirement now Python 3.8+)

### Fixed

- **Reverse Shock**: Significant corrections to reverse shock calculations

### Documentation

- Updated documentation and examples
- Enhanced API reference

## [v0.1.8] - 2025-06-09

### Improved

- **Reverse Shock**: Refined reverse shock calculations and algorithms
- **Code Quality**: General code cleanup and optimization

## [v0.1.7] - 2025-05-24

### Fixed

- **macOS Compatibility**: Fixed macOS-specific bugs
- **Deep Newtonian Regime**: Refinement for deep Newtonian regime calculations
  - Shocked electron calculations moved from energy space to momentum space
  - Enhanced self-absorption refinement
  - Interface improvements

### Improved

- **Code Quality**: General updates and optimizations

## [v0.1.6] - 2025-05-15

### Added

- **Energy Injection**: Added missing Python bindings header for energy injection functionality

### Improved

- **Magnetar Integration**: Enhanced C-level magnetar binding for better performance

## [v0.1.5] - 2025-05-14

### Added

- **Inverse Compton Enhancements**:
  - Updated IC spectrum calculations
  - Added magnetar injection Python interface
  - Enhanced IC code with cleanup and optimization

### Changed

- **Python Support**:
  - Removed Python 3.6 support
  - Added Python 3.7 support

### Fixed

- **Unit System**: Corrected unit handling in various calculations
- **Build System**: Various build fixes and improvements

### Documentation

- Updated examples in documentation
- Enhanced README with new features
- Improved API documentation

## [v0.1.4] - 2025-05-11

### Added

- **Enhanced Build System**: Improved CMake configuration and cross-platform support

### Fixed

- **Windows Compatibility**: Resolved various Windows build issues
- **External Dependencies**: Fixed external library integration

## [v0.1.3] - 2025-05-07

### Improved

- **Build System**: Enhanced CMake build configuration
- **Documentation**: Updated README and logo integration

### Fixed

- **Cross-platform Issues**: Resolved various build issues across platforms

## [v0.1.2] - 2025-05-06

### Added

- **Enhanced Documentation**: Added comprehensive documentation system
- **Logo Integration**: Added project logo and branding

### Improved

- **Code Structure**: Enhanced code organization and commenting
- **Build System**: Improved build workflow and packaging

## [v0.1.1] - 2025-05-06

### Fixed

- **Build System**: Fixed wheel building and packaging issues
- **Documentation**: Updated README with installation instructions

### Improved

- **Layout**: Enhanced README layout and organization

## [v0.1.0] - 2025-05-05

### Added

- **Initial Release**: First public release of VegasAfterglow
- **Core Features**:
  - High-performance C++ framework with Python interface
  - Comprehensive GRB afterglow modeling
  - Forward and reverse shock dynamics
  - Synchrotron and inverse Compton radiation
  - Structured jet configurations
  - MCMC parameter fitting capabilities
- **Radiation Mechanisms**:
  - Synchrotron emission with self-absorption
  - Inverse Compton scattering with Klein-Nishina corrections
  - Synchrotron Self-Compton (SSC)
- **Physical Models**:
  - Multiple ambient medium types (ISM, wind, user-defined)
  - Various jet structures (top-hat, power-law, Gaussian, user-defined)
  - Relativistic and non-relativistic shock evolution
  - Energy and mass injection
- **Performance**: Ultra-fast light curve computation (millisecond timescales)
- **Cross-platform Support**: Linux, macOS, and Windows compatibility
- **Python Interface**: User-friendly Python bindings for easy integration

### Documentation

- Comprehensive API documentation
- Installation guides for Python and C++
- Usage examples and tutorials
- Quick start guide with Jupyter notebooks

---

## Version History Summary

| Version | Release Date | Key Features |
|---------|--------------|--------------|
| v1.0.3  | 2025-09-29   | Self-absorption heating disabled by default, MCMC improvements, flux unit fixes |
| v1.0.2  | 2025-09-15   | Python interface method name updates, enhanced documentation, API consistency |
| v1.0.1  | 2025-09-15   | Frequency integrated flux support, return object interface redesign |
| v1.0.0  | 2025-09-06   | **Major Release**: Advanced MCMC framework, adaptive mesh generation, new jet models |
| v0.2.8  | 2025-08-15   | Inverse Compton performance optimization (~10x), interface improvements |
| v0.2.7  | 2025-07-31   | Internal quantities evolution interface, performance fixes |
| v0.2.6  | 2025-07-19   | Two-component jet support |
| v0.2.5  | 2025-07-19   | Enhanced jet edge detection |
| v0.2.4  | 2025-06-28   | Magnetar spin-down documentation |
| v0.2.3  | 2025-06-27   | Unit system updates, parameter naming |
| v0.2.2  | 2025-06-23   | Default resolution improvements |
| v0.2.1  | 2025-06-22   | Reverse shock enhancements |
| v0.1.9  | 2025-06-19   | Python 3.7 support removed, reverse shock fixes |
| v0.1.8  | 2025-06-09   | Reverse shock refinements |
| v0.1.7  | 2025-05-24   | macOS fixes, deep Newtonian regime improvements |
| v0.1.6  | 2025-05-15   | Magnetar integration, energy injection bindings |
| v0.1.5  | 2025-05-14   | Inverse Compton enhancements, Python support updates |
| v0.1.4  | 2025-05-11   | Build system improvements, Windows compatibility |
| v0.1.3  | 2025-05-07   | Enhanced build system, cross-platform fixes |
| v0.1.2  | 2025-05-06   | Documentation system, logo integration |
| v0.1.1  | 2025-05-06   | Build fixes, documentation updates |
| v0.1.0  | 2025-05-05   | Initial public release |

---

For detailed information about each release, see the individual version sections above.
For the latest updates and development progress, visit our [GitHub repository](https://github.com/YihanWangAstro/VegasAfterglow).
