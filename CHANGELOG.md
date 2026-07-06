# Changelog

All notable changes to VegasAfterglow will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Roadmap

### v2.1.0 — Planned

- [ ] Proton processes
- [ ] Time-evolving microphysical parameters (ε_e, ε_B, p, ξ_e)
- [ ] Polarization
- [ ] External inverse Compton

Have a feature request? [Open an issue](https://github.com/YihanWangAstro/VegasAfterglow/issues) — community input helps shape the roadmap.

---

## [v2.0.6] - 2026-07-05

### Changed

- **MCMC fitting throughput up 3-4x**: the grid, observer, and default-resolution improvements below combine to make typical likelihood evaluations 3-4x faster than v2.0.5 (on-axis tophat fit: 2.2 ms → 0.6 ms per evaluation; gaussian jet with reverse shock: 29 ms → 7 ms), which carries one-to-one into MCMC throughput on any machine. The walker retuning below adds ~25% on top
- **Angular and temporal grids sharpened**: off-axis axisymmetric jets now integrate azimuth over the half range (the flux is mirror-symmetric about the jet-observer plane), roughly doubling off-axis evaluation speed at identical accuracy. Azimuthal and polar refinement respond to the actual Doppler structure, resolving viewing angles at and just outside a sharp jet edge that were previously under-resolved (errors of several percent reduced to below one percent), and the forward-shock time lattice concentrates points around the deceleration turnover at no extra cost
- **Lower default resolution, per mode**: with the sharper grids, the default drops from (0.1, 0.25, 10) to (0.06, 0.15, 6) for forward-shock-only runs and (0.06, 0.2, 10) with a reverse shock (its components need denser polar and time lattices) — 15-35% faster at unchanged worst-case accuracy, calibrated per emission component against the validation gates (mean < 5%, max < 15%) across jet families and viewing angles. An explicit `resolutions` is honored as given. Golden baselines regenerated
- **Magnetized reverse-shock accuracy**: shells with `sigma0 > 0` automatically tighten the dynamics ODE tolerance one decade — their jump conditions are far more tolerance-sensitive than the unmagnetized path, and at the previous uniform tolerance a sigma ~ 1 shell silently carried ~1% reverse-shock flux error. New golden baselines pin the reverse shock at sigma = 1 and 10 alongside the existing 0.1
- **Reverse-shock evaluations reuse the forward-shock observer geometry**: both shocked regions ride the same contact discontinuity, so the equal-arrival-time grids are computed once per call instead of twice (~5% faster with a reverse shock, bit-identical results). The series flux path also processes each time interval as a block with its boundary spectra carried between adjacent intervals, and axisymmetric geometry factors assemble their logarithms from factorized small passes (~10% faster on multi-band series data)
- **MCMC fitting: twice the walkers, half the steps**: the automatic emcee walker count roughly doubles — sized to fill the likelihood thread pool with several evaluation waves per step on 8-16 core machines (~25% higher sampling throughput; `nwalkers=...` still overrides). Since each step now draws twice the samples, **half the previous `nsteps` yields the same posterior sample volume** — halve `nsteps` (and `nburn`) in existing scripts to keep runtimes unchanged with better ensemble mixing; example scripts and docs updated accordingly. The default move mixture also shifts to `DEMove 0.9 / DESnookerMove 0.1` (measured ~20% shorter autocorrelation times than the previous 0.7/0.3 weighting)
- **Observer-layer performance**: off-axis flux evaluations are 10-20% faster (exact results preserved) — the equal-arrival-time grids take their logarithms in vectorized whole-tensor passes, flux contributions batch their final exponential per grid column, and the spectrum evaluation skips provably-trivial work (the inverse-Compton correction when IC cooling is off, and the optically-thick branch's far-field exponential)
- **Opt-in fast math improved** (`AFTERGLOW_FAST_MATH`, still off by default): the polynomial `log2`/`exp2` kernels are upgraded to degree-5 minimax fits (relative error at the 1e-7 level, from 4e-4), gain full domain handling — non-positive, subnormal, infinite, and NaN inputs behave like the standard library — and now cover every `log2`/`exp2` call site in the library. Enabling the option speeds up off-axis structured-jet evaluations by ~1.2x; golden baselines are always generated with the default exact-libm build, and the regeneration script enforces this
- Removed the hard radiative-efficiency truncation (`eps_rad = 0` deep in slow cooling); the efficiency now follows the analytic form everywhere, consistently across forward, reverse, and simple shock solvers
- Generic `Ejecta` and `Medium` constructors now validate their inputs (callability, finite on-axis values) like the typed factories
- Validation suite moved under the unified test tree at `tests/validation/`; the release validation PDF was replaced by the HTML report published at [reports/latest](https://yihanwangastro.github.io/VegasAfterglow/reports/latest/), with per-release snapshots archived in the [reports index](https://yihanwangastro.github.io/VegasAfterglow/reports/)
- `run_validation.py` slimmed to running and checking the suites (report flags removed); the `[test]` extra now needs only pytest, pytest-cov, and matplotlib

### Removed

- **`cmb_cooling` option removed** from `Radiation`, `Fitter`, and the CLI: inverse Compton cooling off the CMB is negligible for GRB afterglow shocks (the comoving magnetic energy density dwarfs the CMB energy density in all relevant regimes) and the flag added API and code complexity without practical use

### Fixed

- **Radiative losses in the blast-wave dynamics were underestimated**: the cooling Lorentz factor in the dynamics-side radiative efficiency omitted the 8π from `B'² = 8π ε_B e` (the synchrotron module had it right), placing the cooling break ~25× too high — the fireball evolved too adiabatically. With the fix, radiative deceleration strengthens and late-time light curves dim by up to ~10-30% at `eps_e = 0.1` (up to ~70% at `eps_e = 0.3`); the effect grows with `eps_e`. The efficiency also gains the same trans-relativistic cooling correction the synchrotron module uses. Present in all releases since v1.0.0; golden baselines regenerated, all closure-relation and crossing-scaling checks pass unchanged
- **Reverse-shock rates are evaluated on the physical domain** (`Gamma3 <= Gamma4`, `0 <= m3 <= m4`, `x3, U3 >= 0`): adaptive-step overshoot at the seed scale of the shocked-shell variables could previously escalate into a runaway that stalled the solver ("ODE exceeded steps" warnings) and corrupted that grid row's light curve; healthy solves are unaffected
- **Reverse-shock crossing end** is now located by bisection on the ODE dense output instead of at stored grid times: the frozen crossing state was previously up to one time-grid step late, biasing the post-crossing reverse-shock light curve near the crossing peak (golden baselines regenerated)
- **Reverse-shock crossing rate** now computed via cancellation-free Lorentz-factor identities: the direct `(beta4 - beta3)/(1 - beta3)` form lost up to all precision at crossing onset where the shocked and unshocked ejecta move at nearly equal speeds
- **Magnetized shock jump conditions**: the downstream four-velocity cubic is now solved by root deflation; the previous direct evaluation lost precision for relativistic shocks with high relative Lorentz factors
- Relativistic kinematics helpers now use `(Gamma - 1)(Gamma + 1)` in place of `Gamma^2 - 1` and avoid `1 - beta` forms, making results independent of platform FMA behavior
- **Reverse-shock relative Lorentz factor** now computed via a cancellation-free identity: near equal shell/shock velocities the old form lost precision to catastrophic cancellation, seeding numerical noise that could destabilize resolution-convergence of reverse-shock light curves (platform/SIMD dependent)
- **Reverse shock with magnetized ejecta (`sigma0 > 0`)**: fixed a forward-shock flux runaway when reverse-shock emission was enabled — the shell-crossing rate is now gated on shell penetration and capped at the fast magnetosonic speed, so light curves evolve correctly across sigma0 from 0 to >10 (unmagnetized results unchanged)

### Added

- **Unified test framework** (see `TESTING.md`): `make test` runs the C++ unit tests, the Python suite, and the full validation suite, and writes `test-report.html` — a single self-contained page with every test outcome and description, physics-validation figures (shock evolution, characteristic frequencies, spectral regimes, measured-vs-expected checks), full resolution-convergence results, and per-configuration performance timing; `make test-quick` runs just the fast tiers
- Comprehensive test coverage: closure relations against standard afterglow theory, exact invariants, golden-baseline regression, parameter-space corner sweeps (including a magnetized-ejecta + reverse-shock regression), and expanded C++ physics assertions — every test carries a one-line description of what it asserts
- CI workflow running the full Python test suite on Linux, macOS, and Windows, with C++ tests on Linux/macOS and a unified report artifact per OS
- Weekly scheduled validation run that refreshes the published report

---

## [v2.0.5] - 2026-06-01

### Added

#### ► **`Fitter.draw_fit()` — one-call diagnostic plot with credible bands**

- Two-panel data + model figure in one call. Top: data + posterior-median light curves with a shaded 68% credible band. Bottom: observer-frame ν_a / ν_m / ν_c, marking where each break crosses an observed band.
- Band style via ``obs_noise``: ``'none'`` (default) for the model-curve band, ``'frac'`` for a posterior-predictive band with σ ∝ flux, ``'abs'`` for constant σ. ``ci=`` / ``n_samples=`` tune the band; ``n_samples=0`` falls back to MAP.
- Filter / instrument names from ``add_flux_density(label=...)`` (or auto-detected from ``nu=filter("r")``) drive a curated color palette: SDSS, Johnson, 2MASS, HST WFC3, SVOM VT, Swift, Einstein Probe, Fermi.

#### ► **`fitter.save(path)` / `Fitter.load(path)` — one-line persistence**

- ``fitter.save(path)`` writes a full snapshot (constructor args, data, parameter defs, samples) to bilby-native HDF5 / JSON; ``Fitter.load(path)`` reconstructs the configured ``Fitter`` + ``FitResult`` in one call, ready for predictions.
- Custom ``jet`` / ``medium`` / ``extinction`` callables can't round-trip — pass them back via ``Fitter.load(path, jet=...)``. ``Fitter.fit(...)`` also caches its result on ``self.result``.

### Changed

#### ► **Synchrotron break shapes use Granot & Sari (2002) Table 2**

- Regime-dependent smoothing parameters at the ν_m / ν_c / ν_a breaks for ISM (k=0), replacing the uniform `s = 1` approximation.
- Light-curve amplitudes above ν_c shift by a p-dependent factor: ~−1% at p=2.05, ~+6% at p=2.3, ~+19% at p=2.6. MCMC fits calibrated against earlier releases may need a re-run.

#### ► **Continuous synchrotron formula across regime crossings**

- The optically-thin spectrum is now a single double-smoothed expression that asymptotes to both slow- and fast-cooling segments — removes the ~3% light-curve jump at the ν_m–ν_c crossing.
- The ν_a join is also smoothed continuously as ν_a moves through ν_m and ν_c, removing similar regime-switch discontinuities.

### Removed

#### ► **Dropped Python 3.8 support**

- Minimum Python version is now 3.9. Python 3.8 reached end-of-life in October 2024, and several runtime dependencies (`bilby>=2.7`, `numpy>=2.0`) already required 3.9+.

### Fixed

#### ► **Inverse Compton correction self-consistency**

- The IC steepening factor `(1+Y_c)/(1+Y(ν))` is applied only to the optically-thin branch above ν_c, not to the full spectrum, so the SSA-dominated branch is no longer over-suppressed in fast cooling or when ν_a > ν_c.
- ν_a now incorporates IC cooling in segments where ν_a > ν_c. Default behaviour with `ssc=False` is unchanged.

#### ► **Spectrum below ν_c independent of cooling treatment**

- The optically-thin synchrotron spectrum below ν_c no longer depends on the ν_c position; the ~4% upward leak in models with stronger IC cooling is gone. KN-cooling spectra are now correctly bounded above by the no-IC spectrum everywhere.

---

## [v2.0.4] - 2026-05-17

### Added

#### ► **`FitResult.summary()` / `save()` / `load()`**

- `result.summary()` returns a top-K best-fit table (rank, chi², parameter values) matching the boilerplate every fit notebook re-implements today — drops ~10 lines per notebook. Renders cleanly in both `print(result.summary())` and Jupyter last-line auto-display.
- `result.save(path)` / `FitResult.load(path)` persist a fit to disk in bilby's native HDF5 (`.h5`/`.hdf5`) or JSON (`.json`) format. Files are interoperable with bilby tooling — `bilby.read_in_result(path).plot_corner()` works on any saved file, and any pre-existing bilby Result file loads as a `FitResult`. Skips re-running the MCMC when reopening a notebook.

---

## [v2.0.3] - 2026-05-09

### Fixed

#### ► **Extinction now applied in post-fit `flux_density_grid`**

- `Fitter.flux_density_grid(best_params, t, nu)` now applies the same per-nu host-galaxy extinction used during the fit's chi-squared evaluation, so plotted best-fit light curves match the data instead of being systematically brighter at optical wavelengths
- `Fitter.flux(...)` (band-integrated) and `Fitter.model(...)` (raw `Model` escape hatch) remain consistent with the fit-time band-integrated path and explicitly do not apply extinction; docstrings updated to reflect this

---

## [v2.0.2] - 2026-05-01

### Added

#### ► **Host-Galaxy Dust Extinction**

- New `extinction` keyword on `Fitter` applies host-galaxy attenuation to model fluxes before chi-squared. Pass `"smc"`, `"lmc"`, or `"mw"` for built-in Pei (1992) laws, or a callable `f(lam_cm, params) -> k(lam)` for a custom law (e.g., a fitted `R_V`)
- Fit host extinction via `ParamDef("A_V", 0, 3, Scale.linear)`; attenuation is `f_obs = f_model · 10^(-0.4 · A_V · k(λ_rest))`
- Built-in laws zero out below the Lyman limit (912 Å), so X-ray data is unaffected — Galactic line-of-sight dust must still be removed from the data upstream
- See [Host-Galaxy Extinction](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/model_configurations.html#host-galaxy-extinction) and [Custom Extinction Profiles](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/advanced.html#custom-extinction-profiles)

---

## [v2.0.0] - 2026-03-01

**[Validation Report](https://yihanwangastro.github.io/VegasAfterglow/reports/latest/)**

### Added

#### ► **Redesigned MCMC Fitting Module**

The `Fitter` class has been rebuilt on top of the `Model` API, replacing the internal `VegasMC` batch evaluator. This brings full Python-level flexibility while retaining C++ performance.

- **Simplified data workflow**: Add observations directly to the fitter with `fitter.add_flux_density()`, `fitter.add_spectrum()`, and `fitter.add_flux()`
- **Custom jet and medium profiles**: Pass `jet` or `medium` factory functions to `Fitter` to use arbitrary angular energy/Lorentz factor profiles or density structures in MCMC fitting (see [Advanced Fitting](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/advanced.html#advanced-fitting))
- **`@gil_free` decorator**: Compile custom jet/medium profile functions to native code for full multi-threaded performance, eliminating the GIL bottleneck that slows Python callbacks during parallel MCMC evaluation (requires `numba`)
- **Custom priors**: Supply bilby `Prior` objects via the `priors` argument to go beyond uniform priors — works with all samplers (see [Custom Priors](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/advanced.html#custom-priors))
- **Custom likelihood functions**: Override the default Gaussian log-likelihood with `log_likelihood_fn` for non-standard noise models or upper limits (see [Custom Likelihood](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/advanced.html#custom-likelihood))
- **`logscale_screen` standalone utility**: Logarithmic data subsampling now available as a module-level function
- **Thread-based parallelism**: Likelihood evaluations parallelized via `ThreadPoolExecutor` with the GIL released during C++ computation, providing near-native throughput

For the full MCMC fitting guide, including advanced customization examples, see the [MCMC Parameter Fitting](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/index.html) documentation.

#### ► **Unit Conversion Constants**
- New `VegasAfterglow.units` module with multiplicative constants for common astronomical units (days, GHz, keV, mJy, Mpc, degrees, etc.)
- Multiply to convert inputs: `5*GHz`, `t_data*day`, `f_data*mJy`; divide to convert outputs: `flux.total / mJy`
- Magnitude converters for AB, Vega, and ST (HST STMAG) systems with built-in filter data for Johnson-Cousins (UBVRI), 2MASS (JHKs), Swift UVOT, SDSS (griz), and HST WFC3/ACS filters
- Survey/telescope filters: SVOM VT (VT_B, VT_R) and WFST wide-band (w)
- Named instrument bands via `units.band()`: look up frequency ranges for Swift XRT/BAT, Einstein Probe FXT/WXT, SVOM MXT/ECLAIRs, and Fermi LAT/GBM — use with `Fitter.add_flux()` or `Model.flux()`

#### ► **Smooth Synchrotron Spectra**
- New default synchrotron spectral model with smooth transitions at break frequencies (ν_m, ν_c, ν_a)
- Produces more physical light curves without artificial kinks at spectral breaks
- The original sharp broken power-law model remains available as `PowerLawSyn`

#### ► **CMB Inverse Compton Cooling**
- Electron cooling now includes inverse Compton scattering off the cosmic microwave background (CMB), which contributes to the Compton-Y parameter alongside synchrotron self-Compton
- Can be enabled independently of SSC via the ``cmb_cooling`` flag, or combined with SSC for full IC cooling
- More relevant for blazar jets than GRB afterglows, but included for completeness as the CMB energy density scales as (1+z)⁴

#### ► **Command-Line Interface (`vegasgen`)**

- New `vegasgen` command generates multi-band light curves directly from the terminal — no Python scripting required
- All physical parameters (jet type, medium, observer, radiation) configurable via command-line arguments with sensible defaults
- Frequencies can be specified as Hz values or standard filter names (e.g., `--nu R J F606W 1e18`), with support for Johnson-Cousins, 2MASS, Swift UVOT, and HST filters
- Output to CSV or JSON, or generate publication-quality plots with `--plot` (Times New Roman, LaTeX labels, colorblind-friendly palette)
- See the [CLI documentation](https://vegasafterglow.readthedocs.io/en/latest/using_cli.html) for the full argument reference

#### ► **Smarter Adaptive Grid Generation**
- Grid generation now accounts for the ambient medium density profile, producing better time sampling near the deceleration time
- Sharp features in jet angular profiles (e.g., core-wing boundaries) are automatically detected and resolved as grid anchor points
- Improved phi-direction sampling for off-axis observers

#### ► **Smoother Reverse Shock Evolution**
- Engine shutdown transitions smoothly rather than cutting off abruptly, producing more stable ODE integration and smoother light curves
- More accurate initial conditions with proper swept-mass and thermal energy integration

#### ► **Lightweight Core Installation**
- `pip install VegasAfterglow` now requires only `numpy` — no MCMC dependencies
- MCMC fitting available via `pip install VegasAfterglow[mcmc]` (bilby, emcee, dynesty)
- Core physics imports always work; MCMC types give a clear error message when extras are missing

#### ► **Redback Integration**
- Full documentation and tutorial for using VegasAfterglow models within the Redback transient inference framework
- New `redback-inference.ipynb` example notebook

#### ► **Validation and Testing**

- **Benchmark tests**: Verify numerical convergence across resolution parameters and measure computation speed for all jet/medium/radiation configurations. Default resolution `(0.1, 0.25, 10)` validated to converge with mean error < 5%
- **Regression tests**: Check that simulation outputs match theoretical predictions from GRB afterglow theory — shock dynamics power-law scaling, characteristic frequency evolution, spectral shape indices, and coasting/Blandford-McKee/Sedov-Taylor phases
- **Python API tests**: Verify model creation, flux calculations, output correctness, and input validation across all jet types and configurations
- **One-command validation**: Run `python validation/run_validation.py --all` to execute the full suite and generate a PDF report with convergence plots and diagnostics
- **Interactive PDF reports**: Clickable navigation links in validation reports — summary grids link to detail pages, detail pages link back to summaries. Separate clickable regions for ISM and Wind medium columns in regression test summaries
- **GitHub Actions CI**: Automated test workflow for Ubuntu and macOS across Python 3.9–3.11 (see `.github/workflows/test.yml`)
- See the [Validation & Testing](README.md#validation--testing) section in the README for full details
- Powered by [Claude Code](https://claude.ai/claude-code)

#### ► **Built-in Profiling**

- See where computation time is spent (dynamics, electrons, photons, flux) via `Model.profile_data()` in Python
- Useful for identifying bottlenecks

#### ► **Per-Cell Spectrum Access in Simulation Details**

- `details.fwd.sync_spectrum[i,j,k](nu_comv)` returns the comoving synchrotron specific intensity at any frequency for a given grid cell
- `details.fwd.ssc_spectrum[i,j,k](nu_comv)` returns the comoving SSC specific intensity (when SSC is enabled)
- `details.fwd.Y_spectrum[i,j,k](gamma)` returns the Compton Y parameter as a function of electron Lorentz factor
- All evaluators accept NumPy arrays for vectorized queries

#### ► **Python Object Introspection**

- All key Python objects now have informative `__repr__` output: `Observer`, `Radiation`, `Magnetar`, `Model`, `FluxDict`, `Flux`, `ShockDetails`, `SimulationDetails`
- Read-only properties on `Model` for inspecting state after construction: `model.observer`, `model.fwd_rad`, `model.rvs_rad`, `model.resolutions`, `model.rtol`, `model.axisymmetric`
- Read-only properties on `Observer` (`lumi_dist`, `z`, `theta_obs`, `phi_obs`) and `Radiation` (`eps_e`, `eps_B`, `p`, `xi_e`, `ssc`, `kn`)

#### ► **Improved Inverse Compton Spectrum Computation**

- Adaptive frequency grid concentrates resolution near synchrotron break frequencies (ν_a, ν_m, ν_c, ν_M), improving accuracy of the seed photon spectrum
- Redesigned electron energy integration produces smoother SSC light curves and spectra

#### ► **Sky Image Generation**

- New `model.sky_image()` method generates spatially resolved afterglow images at any observer time and frequency, supporting all jet types and viewing angles
- Batch evaluation over multiple observer times for efficient multi-frame image sequences and movies
- See the [sky image example notebook](script/sky-image.ipynb) and [documentation](https://vegasafterglow.readthedocs.io/en/latest/examples/sky_image.html) for usage

#### ► **Interactive Web Tool**

- Browser-based real-time, parameter-tunable afterglow modeling tool for light curves, spectra, and sky images, with mobile-phone compatible UI and observational data input/upload support (CSV/TXT/XLS/XLSX) — no installation required: [vegasafterglow.vercel.app](https://www.vegasafterglow.com)

#### ► **Smooth Electron Distribution for IC Integration**

- Electron energy distribution now uses smooth broken power-law transitions at cooling and injection breaks, consistent with the synchrotron photon spectrum

### Performance

#### ► **~5x Faster Computation for Top-Hat and Two-Component Jets**
- Automatic symmetry detection eliminates redundant computation for symmetric jet structures
- Top-hat jets now compute a single representative grid point and broadcast results across the full angular grid
- Two-component jets similarly benefit by computing only unique angular groups
- Applies to all stages: shock dynamics, electron distribution, synchrotron and IC photon generation

#### ► **Reduced Memory Usage in Flux Computation**
- Eliminated large intermediate arrays during observer flux integration
- Direct accumulation replaces the previous two-pass (store-then-sum) approach

### Changed

#### ► **⚠️ BREAKING: MCMC Fitting Interface Redesigned**

The MCMC fitting module has been completely redesigned. Existing fitting scripts will need to be updated.

- **`VegasMC` removed**: The C++ batch evaluator has been replaced by the new `Fitter` class, which uses the Python `Model` API internally with equivalent performance
- **`ObsData` removed**: Use `fitter.add_flux_density()`, `fitter.add_spectrum()`, and `fitter.add_flux()` to add observational data directly to the `Fitter`
- **`Setups` removed**: All model configuration is now passed directly to the `Fitter` constructor as keyword arguments: `Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="ism", rvs_shock=True, ...)`
- **`jet` and `medium` accept callables**: In addition to built-in type strings, you can now pass factory functions that return custom `Ejecta` or `Medium` objects

See the [MCMC Parameter Fitting](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/index.html) documentation for migration examples.

#### ► **⚠️ BREAKING: `ssc_cooling` Parameter Removed**
- The `ssc_cooling` parameter has been removed from `Radiation` and `Fitter`
- IC cooling is now automatically enabled when `ssc=True`
- Update calls: `Radiation(..., ssc_cooling=True, ssc=True)` → `Radiation(..., ssc=True)`

#### ► **Default Resolution**
- Default grid resolution changed from `(0.3, 1, 10)` to `(0.1, 0.25, 10)` for better accuracy out of the box
- Thanks to the improved adaptive grid algorithm (medium-aware time sampling, jet edge anchoring), the code converges with fewer grid points than before — so the new default is both more accurate and comparably fast

#### ► **Default MCMC Sampler**
- Dynesty (nested sampling) is now the recommended default sampler, replacing emcee
- Emcee move strategy updated to `DEMove + DESnookerMove` for more efficient sampling

#### ► **Inverse Compton Cooling**
- Expanded Klein-Nishina regime classification (5 regimes, up from 3) for more accurate IC cooling in the strong-KN limit

#### ► **Wind Medium Defaults**
- `Wind(A_star, n_ism=None, n0=None, k_m=2)`: `n_ism` and `n0` now default to `None` instead of `0` and `inf`, handled internally — cleaner intent for the common pure-wind case
- Wind density power-law index renamed from `k` to `k_m` to match the MCMC fitting parameter name

#### ► **Python >= 3.8 Required**
- Minimum Python version raised from 3.7 to 3.8

### Fixed
- **Reverse shock ν_c continuity**: Eliminated small jump in cooling frequency across the shock crossing boundary
- **Reverse shock initial conditions**: Smoother transition for thick-shell reverse shock crossing
- **Grid sizing bug**: Fixed incorrect grid dimensions when using merged grid sizes in auto_grid

---

## [v1.1.0] - 2025-12-08

### Added

#### ► **Bilby MCMC Integration**
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

#### ► **Flexible Parameter Interface**
- **LOG-Scale Parameter Naming**: Parameters with `Scale.log` automatically prefixed with `log10_`
  - Example: `E_iso` with `Scale.log` becomes `log10_E_iso` in sampler space
  - Automatic transformation between log10 and physical values
- **Generic Sampler Configuration**: Pass sampler-specific parameters via `**sampler_kwargs`
  - Emcee defaults: `nsteps=5000`, `nburn=1000`, `thin=1`, `nwalkers=2*ndim`
  - Dynesty defaults: `nlive=500`, `dlogz=0.1`, `sample="rwalk"`
  - Extensible to any bilby-supported sampler

### Documentation

#### ► **Comprehensive MCMC Guides**
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

#### ► **Default Physics Settings**
- **Self-Absorption Heating**: Disabled self-absorption heating by default for improved numerical stability
- **Strong Absorption Subsegmentation**: Enhanced subsegmentation (5/2 ratio) for improved handling of strong self-absorption regions

#### ► **MCMC Framework Improvements**
- **Enhanced Move Strategy**: Implemented new MCMC move strategy for better parameter space exploration
- **Wider Initial Spread**: Increased initial parameter spread for more robust sampling
- **General Medium Support**: Extended `-k` parameter support for general medium configurations

### Fixed

#### ► **Model Interface**
- **Flux Unit Consistency**: Fixed model flux unit handling for consistent calculations across interfaces
- **Reverse Shock Dynamics**: Implemented hardcut Gamma treatment for more stable reverse shock evolution

### Performance

#### ► **Inverse Compton Optimization**
- **Early Break Implementation**: Added early break condition for inverse Compton calculations to improve computational efficiency

### Documentation

#### ► **Enhanced Guidelines**
- **MCMC Fitting Guidelines**: Added comprehensive guidelines for MCMC parameter fitting
- **Interface Documentation**: Updated documentation for interface consistency and parameter usage
- **Code Examples**: Enhanced examples and troubleshooting documentation

### Development

#### ► **Code Quality**
- **Interface Consistency**: Major cleanup for interface name consistency across Python bindings
- **Code Cleanup**: General code cleanup and optimization throughout the codebase

---

## [v1.0.2] - 2025-09-15

### Changed

#### ► **⚠️ BREAKING: Python Interface Method Name Updates**
- **Standardized Method Names**: Updated method names across Python interfaces for consistency and clarity
  - **PyModel Methods**:
    - `specific_flux()` → `flux_density_grid()` (multi-dimensional flux density calculations)
    - `specific_flux_series()` → `flux_density()` (time-series flux density calculations)
    - `specific_flux_series_with_expo()` → `flux_density_exposures()` (exposure-averaged flux density)
  - **VegasMC Methods**:
    - `specific_flux()` → `flux_density_grid()` (matches PyModel interface)

#### ► **⚠️ BREAKING: Data Input Parameter Name Cleanup**
- **ObsData Interface**: Simplified parameter names by removing explicit "_cgs" suffixes for cleaner API
  - **Method**: `add_flux_density(nu, t, f_nu, err, weights=None)`
    - Previous: `add_flux_density(nu_cgs, t_cgs, Fnu_cgs, Fnu_err, weights=None)`
    - Updated: `add_flux_density(nu, t, f_nu, err, weights=None)`
  - **Method**: `add_flux(band, t, flux, err, num_points=15, weights=None)`
    - Previous: `add_flux(nu_min, nu_max, num_points, t_cgs, F, F_err, weights=None)`
    - Updated: `add_flux(band, t, flux, err, num_points=15, weights=None)` — `band` is a `(nu_min, nu_max)` tuple; use `units.band("XRT")` for named bands
  - **Method**: `add_spectrum(t, nu, f_nu, err, weights=None)`
    - Previous: `add_spectrum(t_cgs, nu_cgs, Fnu_cgs, Fnu_err, weights=None)`
    - Updated: `add_spectrum(t, nu, f_nu, err, weights=None)`

### Added

#### ► **Enhanced Documentation**
- **CGS Unit Documentation**: Added comprehensive CGS unit documentation across all interfaces
  - Added docstrings to pybind11 methods specifying CGS units for all physical quantities
  - Updated all examples in README, documentation, and notebooks with clear CGS unit comments
  - Consistent unit documentation: `nu [Hz]`, `t [s]`, `f_nu [erg/cm²/s/Hz]`, `flux [erg/cm²/s]`

#### ► **Interface Consistency Verification**
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
data.add_flux(band=(1e14, 1e15), t=times, flux=flux, err=errors)  # or band=band("XRT")
```

**Important**: All physical quantities are still expected in CGS units as before - only the parameter names have been simplified for a cleaner interface.

---

## [v1.0.1] - 2025-09-15

### Added

#### ► **Enhanced MCMC Capabilities**
- **Frequency Integrated Flux Support**: Added support for frequency-integrated observations in MCMC fitting framework
  - New `MultiBandData.add_light_curve()` overload for broadband flux measurements over frequency ranges
  - Enables modeling of bolometric observations and filter-integrated measurements
  - Improved handling of observational data where effective frequency depends on spectral shape

#### ► **Documentation & Examples**
- **Comprehensive API Documentation**: Added extensive documentation with detailed physics explanations
  - Complete parameter descriptions with units and typical ranges
  - Detailed method documentation with usage examples
  - Enhanced code examples demonstrating advanced features

### Changed

#### ► **⚠️ BREAKING: Return Object Interface Redesign**
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

#### ► **API Consistency Improvements**
- **Parameter Naming Standardization**: Major cleanup of parameter names across all interfaces for consistency
  - **ConfigParams**: `fwd_SSC` → `fwd_ssc`, `rvs_SSC` → `rvs_ssc`, `KN` → `kn`
  - **PyRadiation**: `SSC` → `ssc`, `KN` → `kn` (IC cooling is now automatically enabled when `ssc=True`)
  - **Model Constructor**: `forward_rad` → `fwd_rad`, `reverse_rad` → `rvs_rad`
  - All parameter names now follow the consistent snake_case convention

#### ► **Enhanced Code Documentation**
- **Comprehensive Comments**: Added detailed documentation to all major classes and structures
  - Complete parameter descriptions with physics context
  - Method documentation with detailed explanations
  - Consistent documentation style following project conventions

### Fixed

#### ► **Interface Consistency**
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
cfg.KN = True         → cfg.kn = True
```

**For Radiation Physics (PyRadiation/Radiation):**
```python
# Old names → New names
Radiation(eps_e=0.1, eps_B=0.01, p=2.3, SSC=True, KN=True)
# becomes:
Radiation(eps_e=0.1, eps_B=0.01, p=2.3, ssc=True, kn=True)
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

#### ► **Advanced MCMC Framework**
- **More model availability**: All built-in jet and medium models now supported in MCMC fitting, including reverse shock, inverse Compton, magnetar injection and etc. See documentations for detials.

#### ► **Adaptive Mesh Generation**
- **Dynamic Grid Optimization**: New adaptive angular grid generation based on jet properties and viewing angles. Grid points distributed according to Doppler boosting factors for optimal efficiency
- **Performance Gains**: ~5x faster convergence for reverse shocks.

#### ► **New Jet Models**
- **StepPowerLawJet**: Uniform core with sharp transition to power-law wings for realistic jet structures
- **Enhanced TwoComponentJet**: Separate narrow and wide components with independent energy and Lorentz factor profiles
- **Improved PowerLawJet**: Split power-law indices for energy (`k_e`) and Lorentz factor (`k_g`) angular dependence

### **Performance & Computational Improvements**

#### ► **Shock Physics Enhancements**
- **Variable Naming Standardization**: Major refactoring for clarity (`EAT_fwd` → `t_obs_fwd`, `Gamma_rel` → `Gamma_th`)
- **Reverse Shock Optimization**: Major code refactoring for reverse shock dynamics. A unified model for shock crossing and post-crossing evolution.
- **Better Crossing Dynamics**: Track the internal energy evolution during shock crossing for improved accuracy with more accurate shock heating and adiabatic cooling. **Note: The enhanced adiabatic cooling and detailed shock heating treatment leads to even weaker reverse shock emission compared to previous versions.**

#### ► **Numerical & Memory Optimizations**
- **Memory Efficiency**: Reduced memory footprint through optimized array operations and memory access patterns
- **Grid Pre-computation**: Enhanced caching strategies for frequently used calculations

###  **Enhanced Python Interface**

#### ► **Data Management Improvements**
- **MultiBandData Redesign**: Unified handling of light curves and spectra with flexible weighting
- **Series Calculations**: New methods for calculating flux at specific time-frequency pairs
- **Memory-Efficient Storage**: Optimized data structures for large multi-wavelength datasets

###  **Development & Build System**

#### ► **Code Quality Enhancements**
- **Enhanced Template System**: Improved compile-time type checking and memory management
- **Better Documentation**: Comprehensive API documentation with detailed parameter descriptions
- **Updated Examples**: New MCMC tutorials with real data fitting demonstrations

### **API Changes & Migration**

#### ► **Parameter Interface Updates**
- **PowerLawJet**: Split single `k` parameter into separate `k_e` (energy) and `k_g` (Lorentz factor) indices for more flexible modeling
- **TwoComponentJet**: Parameter names standardized to `(theta_c, E_iso, Gamma0, theta_w, E_iso_w, Gamma0_w)` for core and wide components
- **StepPowerLawJet**: New jet model with parameters `(theta_c, E_iso, Gamma0, E_iso_w, Gamma0_w, k_e, k_g)` for core and power-law wing components
- **Wind Medium**: Extended with optional `n_ism` and `n_0` parameters for stratified medium modeling: `Wind(A_star, n_ism=0, n_0=inf)`
- **Medium Class**: Simplified from `Medium(rho, mass)` to `Medium(rho)` - removed separate mass parameter
- **Model Methods**:
  - Removed `specific_flux_sorted_series()` method
  - Added `specific_flux_series_with_expo(t, nu, expo_time, num_points=10)` for exposure time averaging
  - Changed `details(t_obs)` to `details(t_min, t_max)` interface
- **Resolution Parameters**: Default resolution changed from `(0.3, 3.0, 5.0)` to `(0.1, 0.25, 10)` for optimal performance/accuracy balance
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
| v2.0.0  | 2026-03-01   | Redesigned MCMC module, ~5x speedup for symmetric jets, smooth synchrotron spectra, custom jet/medium/prior/likelihood, lightweight install |
| v1.1.0  | 2025-12-08   | Bilby MCMC integration, multi-sampler support, flexible parameter interface |
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
