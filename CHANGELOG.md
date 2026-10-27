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

## [Unreleased]

### Fixed

- **Wind-medium start state.** The forward-shock solver used to start the blast wave as if nothing had been swept up yet, but in a stellar-wind medium the starting point already holds a few percent of the deceleration mass. Fluxes came out 3-10% too bright, more so in denser winds, and the answer changed with the requested time window. The start state now accounts for that mass. ISM results move by less than 1e-4; reverse-shock runs are unchanged.
  - The solver started at `Gamma = Gamma0` at `t0 = min(first grid time, 0.1 s, 0.1 t_dec)`, where the swept mass `m2 ∝ r` is already ~10% of the deceleration mass for a wind. The start state now solves `Gamma (m_jet + m2) + Gamma_eff U + E_rad = Gamma0 m_jet + m2`, with the thermal and radiated shares of the pre-start history evaluated at the mid-history Lorentz factor.
- **Silent zero flux near the deceleration time.** If you requested flux in a very narrow window near the deceleration time (e.g. one epoch), the result could be exactly zero without warning. These requests now return the correct flux.
  - v2.0.6 regression: the window ended inside the deceleration refinement band of the time lattice, so the grid's last node evaluated to NaN and every shock ODE row was silently skipped, returning identically zero flux. Now guarded, with regression tests.
- **Reported chi-squared from emcee fits.** `summary()` and the top-k table previously reported chi-squared that was too high by a constant offset; they now show the true data chi-squared. If you compared models using these numbers, re-check them. Dynesty results were unaffected.
  - The sampler stored the posterior log-probability, so reported chi-squared included the prior normalization (−2·log-prior; ~10 for typical uniform priors — e.g. chi²/DOF reading 2.1 when the data chi²/DOF was 1.2). The prior is now subtracted so reported values reflect the data likelihood alone.
- **Fast dynamics for custom media.** A custom `Medium` used to always take the slow per-angle solve; you now set `isotropic=True` on an angle-independent profile to get the fast dynamics, or `False` to keep the full per-angle solve. Update your `Medium` construction to pass an explicit `isotropic` value.
  - A custom `Medium` was always treated as anisotropic, forcing the dynamics to solve every angular cell separately (~10-30x slower than the built-in media for identical physics). `Medium` now requires an explicit `isotropic` declaration: `True` keeps the fast broadcast dynamics (the claim is cross-checked against an angular probe grid and rejected if the profile visibly varies with angle), `False` solves every angular cell. With a numba-compiled density (`NativeFunc`), isotropic custom media run at parity with the built-ins.
- **Exceptions from custom profile callbacks.** A bug (even a typo) in a custom Python profile used to kill the whole process; it now surfaces as a normal Python exception you can read and debug.
  - An error raised inside a custom Python profile previously terminated the interpreter from inside the solver. Callback errors now surface as a regular Python exception carrying the original traceback.
- **`NativeFunc` signature check.** A numba cfunc whose signature doesn't match a `(phi, theta)` or `(phi, theta, r/t)` profile is rejected at construction with a clear message, instead of silently producing wrong physics.
  - A cfunc whose signature did not match its use was previously called with misaligned arguments, producing wrong physics; it is now rejected at construction with a clear message.
- **Credible-band central line.** `flux_density_credible` / `flux_credible` now return a `.median` that is the pointwise median of the sampled curves and always sits inside the band. For a best-fit overlay, evaluate `flux_density_grid` at `result.top_k_params[0]` instead.
  - `.median` previously returned the model at the componentwise-median parameter vector, which for correlated posteriors is not representative and can fall outside the band (the documentation already described the pointwise median). `.median` is now the 50th percentile of the sampled curves, always inside `[lower, upper]`.
- **`draw_fit` central line.** The plotted line now corresponds to the best-fit trajectory and the chi-squared quoted from the fit, rather than the credible-band median.
  - With posterior samples available, the plotted line was the credible-band median rather than the best-fit trajectory, so the curve did not correspond to the quoted chi-squared. The line is now always `best_params` (defaulting to `result.top_k_params[0]`).
- **Fits no longer hang on a misspelt or mis-typed jet/medium.** `Fitter(jet="tophatt")`, or passing a built jet object such as `TophatJet(...)` instead of a name or factory function, used to pass construction, fail silently inside every likelihood call, and leave Dynesty stuck at "Generating initial points from the prior". These now raise an error that lists the valid names, and an error inside a likelihood call is logged the first time it happens instead of being swallowed.
  - The jet/medium argument is checked in the constructor (registered name or callable factory). Both sampler paths still return `-inf` for a failed evaluation, but log the first few distinct errors.

### Added

- **Calibration uncertainty per dataset.** `add_flux_density`, `add_spectrum` and `add_flux` take a `sys_err`, the fractional flux-scale uncertainty of that dataset (for example `0.1` for 10%). The whole dataset can then shift up or down together by about that amount while the relative precision between its points is kept, so there is no need to inflate error bars. `summary()` and `result.calibration_offsets` report how far each dataset moved. See [Calibration Systematics](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/basic_setup.html#calibration-systematics) in the docs.
  - The points of one call share a single multiplicative calibration factor with a Gaussian prior of width `sys_err` in ln-flux, which is integrated out analytically: `chi2 = sum(p r^2) - (sum p r)^2 / (1/sys_err^2 + sum p)` per dataset. No extra sampler dimension; `log_likelihood_fn` still receives one chi-squared.
- **Upper limits in fits.** Use `add_flux_density_limits`, `add_spectrum_limits` or `add_flux_limits` for non-detections, either as n-sigma limits (3-sigma by default) or with the noise level behind them. A limit costs nothing while the model stays well below it, a little when the model touches it, and more and more the further the model rises above it, so a non-detection rules out models only as strongly as its noise allows. Limits count as data, `draw_fit` draws them as arrows, and they are saved and reloaded with the fit. The old trick of entering a zero flux no longer works with the log-flux likelihood and is replaced by these methods. See [Upper Limits](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/basic_setup.html#upper-limits) in the docs.
  - Each limit `U` with noise `sigma` contributes `-2 ln Phi((U - F_model)/sigma)` in linear flux: about 0 far below the limit, `2 ln 2` at the limit, and `((F_model - U)/sigma)^2` beyond.
- **Chi-squared confidence envelopes.** `flux_density_confidence` / `flux_confidence` shade the likelihood-consistent uncertainty band around a best-fit curve; `draw_fit` now uses this band, while `flux_density_credible` / `flux_credible` remain available for pointwise posterior credible bands.
  - They shade the envelope of curves from all posterior samples inside the joint confidence region `chi2 <= chi2_min + dchi2(cl, n_free)` — the band containing the best-fit curve by construction.
- **Jet spreading in fits (experimental).** `Fitter(spreading=True)` turns on lateral jet spreading for the built-in jet types, so fits can use the same spreading dynamics as the light-curve `Model`. Fits previously always ran without spreading; a custom jet factory could enable it, and still sets its own flag.
  - Intended for the top-hat, Gaussian and power-law jets with a core. The flag is saved and restored with the fit. Combining it with a custom jet factory raises an error.


### Changed

- **Log-flux likelihood.** Chi-squared is evaluated on `ln F` with `sigma_ln = err / F_obs` propagated from your linear errors, so bright points no longer dominate data spanning decades. Fluxes and errors must be strictly positive; chi-squared values differ from earlier releases, so re-baseline any absolute comparisons. For a different likelihood, pass `log_likelihood_fn` to `fit()`.
  - `chi2 = sum(((ln F_obs - ln F_model) / sigma_ln)^2)`. The Gaussian normalization `ln(2 pi sigma_ln^2)` is parameter-independent and omitted, so reported chi-squared stays a pure residual sum; it cancels in posterior sampling and same-data Bayes factors.

- **Synchrotron spectral shape.** The Granot & Sari (2002) regime-blended break sharpnesses introduced in v2.0.5 tied the spectral shape to the evolving break frequencies, imprinting visible jumps on light curves, and are calibrated for ISM environments only. Spectral breaks are back to the uniform s = 1 smoothing, the slow/fast-cooling forms switch by break order, and the self-absorption turnover keeps its original p-dependent sharpness. Fluxes near spectral breaks shift relative to v2.0.5/v2.0.6, typically at the few-percent level.
- **Cubic light-curve interpolation.** The straight log-log chord is exact on power-law stretches but cut the corners where spectral breaks sweep through the band; a cubic with slopes from the neighbouring nodes now tracks those transitions, making light curves smooth across grid nodes and 2-5x more accurate near break crossings at the same resolution. Node slopes are limited so the curve stays monotone between nodes and never overshoots them; pure power-law segments are unchanged.
- **Forward-shock time resolution.** Forward-shock-only runs are 30-40% faster at equal light-curve quality: with the cubic time interpolation carrying the smoothness and crossing accuracy, the default time resolution drops from 10 to 5 grid points per decade. Reverse-shock runs keep 10 points per decade (the shock-crossing evolution needs the denser sampling); pass `resolutions` to override either.
- **Flux evaluation speed for fitting.** Narrow-window and band-integrated flux requests, densely sampled light curves, and few-epoch spectra now run noticeably faster with no change to results.
  - Evaluations are clamped to the observed window instead of paying for the full time range (up to ~2x for fit-style calls), repeated-frequency light-curve points share their boundary evaluations (~2x for densely sampled light curves), and grid requests evaluate only the time nodes their samples actually touch (4-5x for spectra at a few epochs).
- **Lane-batched shock dynamics.** The blast-wave ODE now integrates four angular rows in lockstep: structured-jet forward-shock dynamics run ~2x faster and forward+reverse shock pairs 1.3-1.8x, on macOS, Linux, and Windows alike.
  - One shock per SIMD lane through plain auto-vectorizable arithmetic (no intrinsics or platform branches); rows sharing a batch adopt the most demanding row's step size, which shifts fluxes by at most ~1e-3 near the reverse-shock crossing and far less elsewhere.
- **Grid construction speed.** The adaptive azimuthal grid no longer re-derives jet-profile terms at every sample, and the per-group time lattices fill without expression-template overhead — grid setup is roughly 2x faster off-axis and ~30% faster for structured jets.
- **Band-integration default.** The fitter's default `num_points` drops from 15 to 5, making band-flux calls ~2.7x faster with negligible error for optical and X-ray bands. For radio bands or bands wider than two decades pass 9-17, and for bolometric requests pass 33+.
  - The band integrator uses Boole's rule on a log frequency grid, most accurate when `num_points = 4k + 1`. At the new default, integration error stays below a few 1e-4 for optical and X-ray instrument bands. Radio bands (spectral breaks often in-band) and wide bands need more points.

The SSC (inverse Compton) pipeline was rebuilt this release. It was unchanged between v2.0.5 and v2.0.6, so the comparisons below hold against both.

- **SSC spectral accuracy.** Previous releases systematically overestimated SSC flux, worsening toward high frequencies; SSC light curves above the seed peak now come out ~15% fainter for typical parameters — the earlier values were the biased ones.
  - The seed-photon integral now treats each spectral bin as the local power law instead of a linear trapezoid, removing the overestimate (against a converged reference: optical ~20% → ~3%, X-ray ~25% → ~3%, TeV ~45% → ~4%).
- **Radio-band SSC.** Radio-band SSC was under-predicted by ~45% in all previous releases. Low-frequency SSC now brightens by up to ~1.9x and its resolution-dependent noise (up to 15% scatter) disappears.
  - The electron integral was truncated at min(gamma_m, gamma_c)/3, cutting the soft low-energy tail exp(-gamma_m/gamma) whose suppression is cancelled by the 1/gamma^2 Compton weight for a full decade below the injection minimum. The integration now extends to where that tail genuinely dies.
- **SSC speed.** Inverse-Compton spectra compute ~2.8x faster than v2.0.5, Klein-Nishina SSC ~2.2x, with no loss of accuracy.
  - IC spectra are computed only over the frequency band the observer actually samples, determined per time slice from that slice's Doppler range (up to a further ~1.7x for on-axis Thomson SSC, parameter-dependent); the electron, seed, and output grids share a commensurate log lattice, so the KN cross-section is evaluated once per lattice node (instead of once per electron-seed pair) and the scattering convolution reduces to integer-indexed accumulation.
- **SSC spectral-regime baselines.** Two new golden baselines pin the SSC seed spectrum across every ordering of its spectral breaks, including strongly self-absorbed regimes.
  - The dense-medium baselines cover every ordering of the synchrotron self-absorption, injection, and cooling breaks that shapes the SSC seed spectrum.

---

## [v2.0.6] - 2026-07-05

### Changed

- **MCMC fitting is 3-4x faster than v2.0.5.** Your fits run several times faster with no change to results, thanks to the grid and observer improvements below.
  - On-axis tophat fit: 2.2 ms → 0.6 ms per likelihood; gaussian jet with reverse shock: 29 ms → 7 ms. Walker retuning adds ~25% on top.
- **Sharper, faster angular and time grids.** Off-axis fits are about 2x faster and viewing angles near sharp jet edges are now resolved correctly.
  - Off-axis azimuthal integration runs over the mirror-symmetric half range (~2x off-axis speed); refinement now tracks the Doppler structure, fixing under-resolved viewing angles near sharp jet edges (several percent → below one percent); time grid concentrated around the deceleration turnover.
- **Lower default resolution, tuned per mode.** Default runs are 15-35% faster at the same worst-case accuracy; if you pass an explicit `resolutions` it is still honored as given.
  - Defaults are (0.06, 0.15, 6) for forward-shock-only runs and (0.06, 0.2, 10) with a reverse shock — calibrated per emission component across jet families and viewing angles.
- **Faster observer layer, same results.** Off-axis and multi-band evaluations run faster with no change to output.
  - Off-axis flux evaluations are 10-20% faster, reverse-shock runs reuse the forward-shock equal-arrival-time grid (~5%), and multi-band series evaluation is ~10% faster — none of these change results.
- **MCMC: twice the walkers, half the steps.** The automatic walker count now doubles for ~25% higher throughput. **Halve `nsteps` and `nburn` in existing scripts** to get the same posterior volume at unchanged runtime.
  - Walker count doubles to fill the likelihood thread pool (`nwalkers=...` still overrides). Default moves are now `DEMove 0.9 / DESnookerMove 0.1` (~20% shorter autocorrelation times).
- **More accurate magnetized reverse shocks.** Magnetized shells (`sigma0 > 0`) are now solved more tightly, removing a ~1% flux error.
  - `sigma0 > 0` shells automatically tighten the ODE tolerance (the previous uniform tolerance carried ~1% flux error at sigma ~ 1).
- **Opt-in fast math is more accurate and broader.** Still off by default; enabling `AFTERGLOW_FAST_MATH` gives ~1.2x on off-axis structured jets.
  - The fast math kernels are ~1000x more accurate, handle edge-case inputs correctly, and now cover all `log2`/`exp2` call sites.
- **Radiative efficiency now analytic everywhere.** The previous hard truncation is gone, and `Ejecta`/`Medium` constructors now reject invalid inputs.
  - Radiative efficiency follows the analytic form everywhere (hard truncation removed); generic `Ejecta`/`Medium` constructors validate their inputs.
- **Validation report moved online.** The validation suite is unified under `tests/validation/`, and the release PDF is replaced by the published HTML report at [reports/latest](https://yihanwangastro.github.io/VegasAfterglow/reports/latest/).

### Removed

- **`cmb_cooling` removed.** Drop `cmb_cooling` from your calls; IC cooling off the CMB is negligible for GRB afterglow shocks.
  - Removed from `Radiation`, `Fitter`, and the CLI.

### Fixed

- **Fixed: radiative losses were underestimated in all releases since v1.0.0.** Late-time light curves now dim more; fits calibrated against earlier releases should be re-run.
  - The dynamics-side cooling Lorentz factor missed the 8π in `B'² = 8π ε_B e`, placing the cooling break ~25x too high, so the fireball evolved too adiabatically. Late-time light curves now dim by ~10-30% at `eps_e = 0.1` (up to ~70% at `eps_e = 0.3`).
- **More robust reverse-shock solves.** Reverse-shock runs no longer produce solver runaways or crossing-time bias, giving cleaner grids and flux near the crossing peak.
  - Rates evaluated on the physical domain (no more solver runaways corrupting grid rows); crossing end located by bisection instead of at grid times (removes up to one-grid-step bias near the crossing peak); magnetized (`sigma0 > 0`) shell crossing gated on shell penetration, fixing a forward-shock flux runaway.
- **Reproducible results across platforms.** Numerical results are now consistent across platforms and compilers.
  - Cancellation-free forms for the crossing rate, relative Lorentz factor, and relativistic kinematics helpers, plus a more robust magnetized jump-condition solve.

### Added

- **New `radiative_fireball` switch to select adiabatic dynamics.** Set `radiative_fireball=False` (CLI `--adiabatic`) for the adiabatic approximation used by most afterglow codes — useful for cross-code comparisons or low-`eps_e` sources. The default keeps the radiative dynamics.
  - Available on `Model`, `Fitter`, and the CLI.
- **Unified test framework.** Run `make test` (see `TESTING.md`) to exercise everything at once and get a self-contained HTML report.
  - C++ units, Python suite, and full physics validation in one run: closure relations, golden baselines, parameter-corner sweeps, and per-configuration timing.
- **Continuous integration across platforms.** CI now runs on Linux/macOS/Windows, plus a weekly validation run refreshing the published report.

---

## [v2.0.5] - 2026-06-01

### Added

#### ► **`Fitter.draw_fit()` — one-call diagnostic plot with credible bands**

- **One call gives you a full data-plus-model diagnostic figure.** Use `Fitter.draw_fit()` to get a two-panel plot of your data against the fit, with an uncertainty band, in a single step.
  - Top panel: data + posterior-median light curves with a shaded 68% credible band. Bottom panel: observer-frame ν_a / ν_m / ν_c, marking where each break crosses an observed band.
- **Control how the uncertainty band is drawn.** Pick the band style that matches how you want to show scatter, and tune its width, or fall back to the MAP curve.
  - Band style via ``obs_noise``: ``'none'`` (default) for the model-curve band, ``'frac'`` for a posterior-predictive band with σ ∝ flux, ``'abs'`` for constant σ. ``ci=`` / ``n_samples=`` tune the band; ``n_samples=0`` falls back to MAP.
- **Your filter/instrument labels drive plot colors automatically.** Name your data with ``add_flux_density(label=...)`` (or via ``nu=filter("r")``) and each instrument gets a consistent, curated color.
  - Curated color palette covers SDSS, Johnson, 2MASS, HST WFC3, SVOM VT, Swift, Einstein Probe, Fermi.

#### ► **`fitter.save(path)` / `Fitter.load(path)` — one-line persistence**

- **Save and reload an entire fit in one line.** ``fitter.save(path)`` stores everything about a fit to disk and ``Fitter.load(path)`` brings it back ready for predictions, so you can stop and resume work without re-running.
  - ``fitter.save(path)`` writes a full snapshot (constructor args, data, parameter defs, samples) to bilby-native HDF5 / JSON; ``Fitter.load(path)`` reconstructs the configured ``Fitter`` + ``FitResult`` in one call.
- **Custom callables must be passed back in when loading.** If you used custom ``jet`` / ``medium`` / ``extinction`` functions, supply them again on load; also, ``fit(...)`` now keeps its result handy for you.
  - Custom ``jet`` / ``medium`` / ``extinction`` callables can't round-trip — pass them back via ``Fitter.load(path, jet=...)``. ``Fitter.fit(...)`` also caches its result on ``self.result``.

### Changed

#### ► **Synchrotron break shapes use Granot & Sari (2002) Table 2**

- **Spectral breaks now use more accurate, published smoothing shapes.** Break shapes for ISM environments are more realistic than before; no action needed unless you compare against earlier releases.
  - Regime-dependent smoothing parameters at the ν_m / ν_c / ν_a breaks for ISM (k=0), replacing the uniform `s = 1` approximation.
- **Light curves above ν_c shift slightly — re-run older fits.** Flux amplitudes above ν_c change by a small, p-dependent amount, so MCMC fits calibrated on earlier releases may need to be re-run.
  - Amplitude shift above ν_c: ~−1% at p=2.05, ~+6% at p=2.3, ~+19% at p=2.6.

#### ► **Continuous synchrotron formula across regime crossings**

- **Removed a small jump in the light curve at the ν_m–ν_c crossing.** The optically-thin spectrum is now smooth across the slow-to-fast-cooling transition, so the ~3% discontinuity is gone.
  - The optically-thin spectrum is now a single double-smoothed expression that asymptotes to both slow- and fast-cooling segments.
- **The ν_a break is now smooth as well.** Similar regime-switch jumps at the self-absorption break are removed as ν_a moves through ν_m and ν_c.
  - The ν_a join is smoothed continuously as ν_a moves through ν_m and ν_c, removing similar regime-switch discontinuities.

### Removed

#### ► **Dropped Python 3.8 support**

- **Python 3.9 is now the minimum.** If you are on Python 3.8, upgrade to 3.9 or newer to keep using the package.
  - Python 3.8 reached end-of-life in October 2024, and several runtime dependencies (`bilby>=2.7`, `numpy>=2.0`) already required 3.9+.

### Fixed

#### ► **Inverse Compton correction self-consistency**

- **More accurate IC cooling in self-absorbed and fast-cooling models.** The IC steepening correction is no longer over-applied, so the self-absorption-dominated part of the spectrum is no longer over-suppressed.
  - The IC steepening factor `(1+Y_c)/(1+Y(ν))` is applied only to the optically-thin branch above ν_c, not to the full spectrum, so the SSA-dominated branch is no longer over-suppressed in fast cooling or when ν_a > ν_c.
- **ν_a now accounts for IC cooling where relevant.** When ν_a > ν_c, the self-absorption frequency includes IC cooling; default behaviour with `ssc=False` is unchanged.
  - ν_a incorporates IC cooling in segments where ν_a > ν_c.

#### ► **Spectrum below ν_c independent of cooling treatment**

- **Fixed a small spurious flux leak below ν_c in strong-IC-cooling models.** The spectrum below ν_c no longer depends on where ν_c sits, removing a ~4% upward bias in models with stronger IC cooling.
  - The optically-thin synchrotron spectrum below ν_c no longer depends on the ν_c position; the ~4% upward leak is gone. KN-cooling spectra are now correctly bounded above by the no-IC spectrum everywhere.

---

## [v2.0.4] - 2026-05-17

### Added

#### ► **`FitResult.summary()` / `save()` / `load()`**

- **Get a ready-made best-fit table without hand-writing it.** `result.summary()` returns a top-K best-fit table so you no longer need to re-implement that boilerplate in every fit notebook. It displays cleanly whether you `print(result.summary())` or let Jupyter auto-display it.
  - The table gives rank, chi², and parameter values, matching the boilerplate every fit notebook re-implements today — drops ~10 lines per notebook. Renders cleanly in both `print(result.summary())` and Jupyter last-line auto-display.
- **Save a fit to disk and reload it instead of re-running the MCMC.** `result.save(path)` and `FitResult.load(path)` persist a fit and read it back, so reopening a notebook skips the sampler. Files also work with existing bilby tooling.
  - Saves in bilby's native HDF5 (`.h5`/`.hdf5`) or JSON (`.json`) format. Files are interoperable with bilby tooling — `bilby.read_in_result(path).plot_corner()` works on any saved file, and any pre-existing bilby Result file loads as a `FitResult`. Skips re-running the MCMC when reopening a notebook.

---

## [v2.0.3] - 2026-05-09

### Fixed

#### ► **Extinction now applied in post-fit `flux_density_grid`**

- **Plotted best-fit light curves now match your optical data.** `Fitter.flux_density_grid(best_params, t, nu)` applies the same host-galaxy extinction used during the fit, so curves are no longer systematically too bright at optical wavelengths. Re-generate any plots made with this function.
  - Applies the same per-nu host-galaxy extinction used in the fit's chi-squared evaluation, so plotted best-fit light curves match the data instead of being systematically brighter at optical wavelengths.
- **`flux(...)` and `model(...)` deliberately skip extinction — use `flux_density_grid` if you need it applied.** These two paths stay consistent with the fit-time band-integrated calculation and do not attenuate; their docstrings now say so.
  - `Fitter.flux(...)` (band-integrated) and `Fitter.model(...)` (raw `Model` escape hatch) remain consistent with the fit-time band-integrated path and explicitly do not apply extinction; docstrings updated to reflect this.

---

## [v2.0.2] - 2026-05-01

### Added

#### ► **Host-Galaxy Dust Extinction**

- **Attenuate model fluxes for host-galaxy dust before chi-squared.** A new `extinction` keyword on `Fitter` lets you apply host-galaxy dimming using a built-in law or your own.
  - Pass `"smc"`, `"lmc"`, or `"mw"` for built-in Pei (1992) laws, or a callable `f(lam_cm, params) -> k(lam)` for a custom law (e.g., a fitted `R_V`).
- **Fit the extinction amount as a free parameter.** Add `A_V` to your fit to let the sampler constrain host dust.
  - Fit host extinction via `ParamDef("A_V", 0, 3, Scale.linear)`; attenuation is `f_obs = f_model · 10^(-0.4 · A_V · k(λ_rest))`.
- **Built-in laws leave X-ray data untouched, but you must still remove Galactic dust yourself.** The built-in laws stop attenuating below the Lyman limit, so X-ray points are unaffected; Milky Way line-of-sight extinction must be corrected in the data beforehand.
  - Built-in laws zero out below the Lyman limit (912 Å), so X-ray data is unaffected — Galactic line-of-sight dust must still be removed from the data upstream.
- **See the docs for details.** See [Host-Galaxy Extinction](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/model_configurations.html#host-galaxy-extinction) and [Custom Extinction Profiles](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/advanced.html#custom-extinction-profiles).

---

## [v2.0.0] - 2026-03-01

**[Validation Report](https://yihanwangastro.github.io/VegasAfterglow/reports/latest/)**

### Added

#### ► **Redesigned MCMC Fitting Module**

The `Fitter` class has been rebuilt on top of the `Model` API, replacing the internal `VegasMC` batch evaluator. This brings full Python-level flexibility while retaining C++ performance.

- **Add data straight to the fitter.** You can register observations on the `Fitter` directly, with no separate data container to build.
  - Add observations with `fitter.add_flux_density()`, `fitter.add_spectrum()`, and `fitter.add_flux()`.
- **Fit with your own jet and medium shapes.** Supply your own angular energy/Lorentz-factor profiles or density structures to MCMC fitting instead of only the built-in types.
  - Pass `jet` or `medium` factory functions to `Fitter` to use arbitrary angular energy/Lorentz factor profiles or density structures (see [Advanced Fitting](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/advanced.html#advanced-fitting)).
- **Keep custom profiles fast during parallel fitting.** If you write a custom jet/medium function in Python, decorate it so it runs at native speed and doesn't stall multi-threaded MCMC.
  - The `@gil_free` decorator compiles custom jet/medium profile functions to native code for full multi-threaded performance, eliminating the GIL bottleneck that slows Python callbacks during parallel MCMC evaluation (requires `numba`).
- **Use non-uniform priors.** Pass bilby `Prior` objects to move beyond uniform priors; this works with every sampler.
  - Supply bilby `Prior` objects via the `priors` argument (see [Custom Priors](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/advanced.html#custom-priors)).
- **Plug in your own likelihood.** Replace the default Gaussian log-likelihood to handle non-standard noise models or upper limits.
  - Override the default Gaussian log-likelihood with `log_likelihood_fn` (see [Custom Likelihood](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/advanced.html#custom-likelihood)).
- **Logarithmic subsampling is now a plain function.** You can call `logscale_screen` directly without going through the fitter.
  - Logarithmic data subsampling is now available as a module-level `logscale_screen` utility.
- **Parallel fitting runs at near-native speed.** Likelihood evaluations are spread across threads with the GIL released during the C++ work.
  - Likelihood evaluations are parallelized via `ThreadPoolExecutor` with the GIL released during C++ computation, providing near-native throughput.

For the full MCMC fitting guide, including advanced customization examples, see the [MCMC Parameter Fitting](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/index.html) documentation.

#### ► **Unit Conversion Constants**
- **Convert units by simple multiplication.** A new module gives you constants for common astronomical units so you don't hand-code conversion factors.
  - New `VegasAfterglow.units` module with multiplicative constants for common astronomical units (days, GHz, keV, mJy, Mpc, degrees, etc.).
- **Multiply to convert inputs, divide to convert outputs.** Scale your data into and out of code units inline.
  - Multiply to convert inputs: `5*GHz`, `t_data*day`, `f_data*mJy`; divide to convert outputs: `flux.total / mJy`.
- **Work directly in magnitudes.** Converters for AB, Vega, and ST systems ship with filter data for many standard filter sets.
  - Magnitude converters for AB, Vega, and ST (HST STMAG) systems with built-in filter data for Johnson-Cousins (UBVRI), 2MASS (JHKs), Swift UVOT, SDSS (griz), and HST WFC3/ACS filters.
- **Extra survey/telescope filters included.** SVOM VT and WFST wide-band filters are available.
  - Survey/telescope filters: SVOM VT (VT_B, VT_R) and WFST wide-band (w).
- **Look up instrument bands by name.** Fetch a named band's frequency range instead of typing frequencies, then pass it to your flux calls.
  - Named instrument bands via `units.band()`: look up frequency ranges for Swift XRT/BAT, Einstein Probe FXT/WXT, SVOM MXT/ECLAIRs, and Fermi LAT/GBM — use with `Fitter.add_flux()` or `Model.flux()`.

#### ► **Smooth Synchrotron Spectra**
- **Smoother synchrotron spectra by default.** The new default model transitions smoothly across break frequencies instead of showing sharp corners.
  - New default synchrotron spectral model with smooth transitions at break frequencies (ν_m, ν_c, ν_a).
- **More physical light curves.** Light curves no longer show artificial kinks at spectral breaks.
  - Produces more physical light curves without artificial kinks at spectral breaks.
- **Old sharp model still available.** If you need the previous behavior, request it explicitly.
  - The original sharp broken power-law model remains available as `PowerLawSyn`.

#### ► **CMB Inverse Compton Cooling**
- **Electron cooling can now include the CMB.** Cooling can account for inverse Compton scattering off the cosmic microwave background, adding to the Compton-Y parameter alongside SSC.
  - Electron cooling now includes inverse Compton scattering off the CMB, which contributes to the Compton-Y parameter alongside synchrotron self-Compton.
- **Turn CMB cooling on independently or with SSC.** Enable it alone via a flag, or combine it with SSC for full IC cooling.
  - Can be enabled independently of SSC via the ``cmb_cooling`` flag, or combined with SSC for full IC cooling.
- **Mainly matters for blazar jets.** It is included for completeness and is usually minor for GRB afterglows, but grows strongly at high redshift.
  - More relevant for blazar jets than GRB afterglows, but included for completeness as the CMB energy density scales as (1+z)⁴.

#### ► **Command-Line Interface (`vegasgen`)**

- **Generate light curves from the terminal.** A new `vegasgen` command produces multi-band light curves without writing any Python.
  - New `vegasgen` command generates multi-band light curves directly from the terminal — no Python scripting required.
- **Configure the physics from the command line.** Jet, medium, observer, and radiation settings are all command-line arguments with sensible defaults.
  - All physical parameters (jet type, medium, observer, radiation) configurable via command-line arguments with sensible defaults.
- **Ask for frequencies or filter names.** Pass raw Hz values or standard filter names interchangeably.
  - Frequencies can be specified as Hz values or standard filter names (e.g., `--nu R J F606W 1e18`), with support for Johnson-Cousins, 2MASS, Swift UVOT, and HST filters.
- **Save data or make publication plots.** Output to CSV or JSON, or add `--plot` for a ready-to-use figure.
  - Output to CSV or JSON, or generate publication-quality plots with `--plot` (Times New Roman, LaTeX labels, colorblind-friendly palette).
- See the [CLI documentation](https://vegasafterglow.readthedocs.io/en/latest/using_cli.html) for the full argument reference

#### ► **Smarter Adaptive Grid Generation**
- **Better time sampling near the deceleration time.** Grids now adapt to the ambient medium density profile automatically.
  - Grid generation now accounts for the ambient medium density profile, producing better time sampling near the deceleration time.
- **Sharp jet features are resolved automatically.** Features like core-wing boundaries are detected and pinned as grid points without manual tuning.
  - Sharp features in jet angular profiles (e.g., core-wing boundaries) are automatically detected and resolved as grid anchor points.
- **Better sampling for off-axis observers.** Phi-direction sampling is improved.
  - Improved phi-direction sampling for off-axis observers.

#### ► **Smoother Reverse Shock Evolution**
- **Smoother light curves at engine shutdown.** The engine now switches off gradually, giving more stable integration and cleaner light curves.
  - Engine shutdown transitions smoothly rather than cutting off abruptly, producing more stable ODE integration and smoother light curves.
- **More accurate reverse-shock start.** Initial conditions now properly account for swept mass and thermal energy.
  - More accurate initial conditions with proper swept-mass and thermal energy integration.

#### ► **Lightweight Core Installation**
- **A plain install now pulls in only numpy.** The base package no longer drags in MCMC dependencies.
  - `pip install VegasAfterglow` now requires only `numpy` — no MCMC dependencies.
- **Install MCMC support as an extra when you need it.** Use the `[mcmc]` extra to add the fitting stack.
  - MCMC fitting available via `pip install VegasAfterglow[mcmc]` (bilby, emcee, dynesty).
- **Core physics always imports; missing extras fail clearly.** If MCMC extras are absent, MCMC types raise a clear error instead of a cryptic one.
  - Core physics imports always work; MCMC types give a clear error message when extras are missing.

#### ► **Redback Integration**
- **Use VegasAfterglow models inside Redback.** Full documentation and a tutorial show how to run the models in the Redback transient inference framework.
  - Full documentation and tutorial for using VegasAfterglow models within the Redback transient inference framework.
- **New Redback example notebook.** A ready-to-run notebook is included.
  - New `redback-inference.ipynb` example notebook.

#### ► **Validation and Testing**

- **Benchmark tests confirm convergence and speed.** They check that results converge across resolution settings and measure runtime for every configuration.
  - Verify numerical convergence across resolution parameters and measure computation speed for all jet/medium/radiation configurations. Default resolution `(0.1, 0.25, 10)` validated to converge with mean error < 5%.
- **Regression tests check results against theory.** Outputs are compared to known GRB afterglow predictions.
  - Check that simulation outputs match theoretical predictions from GRB afterglow theory — shock dynamics power-law scaling, characteristic frequency evolution, spectral shape indices, and coasting/Blandford-McKee/Sedov-Taylor phases.
- **Python API tests cover the user-facing surface.** They verify model creation, flux calculations, outputs, and input validation across all configurations.
  - Verify model creation, flux calculations, output correctness, and input validation across all jet types and configurations.
- **Run the whole suite with one command.** A single command runs everything and produces a PDF report with convergence plots and diagnostics.
  - Run `python validation/run_validation.py --all` to execute the full suite and generate a PDF report with convergence plots and diagnostics.
- **Validation PDFs are navigable.** Reports include clickable links between summary grids and detail pages.
  - Clickable navigation links in validation reports — summary grids link to detail pages, detail pages link back to summaries. Separate clickable regions for ISM and Wind medium columns in regression test summaries.
- **Tests run automatically in CI.** A GitHub Actions workflow tests Ubuntu and macOS across Python 3.9–3.11.
  - Automated test workflow for Ubuntu and macOS across Python 3.9–3.11 (see `.github/workflows/test.yml`).
- See the [Validation & Testing](README.md#validation--testing) section in the README for full details
- Powered by [Claude Code](https://claude.ai/claude-code)

#### ► **Built-in Profiling**

- **See where your run spends its time.** A profiling method reports time split across dynamics, electrons, photons, and flux.
  - See where computation time is spent (dynamics, electrons, photons, flux) via `Model.profile_data()` in Python.
- **Useful for finding bottlenecks.** Helps you spot the slow stage in a run.

#### ► **Per-Cell Spectrum Access in Simulation Details**

- **Read the synchrotron spectrum of any grid cell.** Query the comoving synchrotron intensity at any frequency for a chosen cell.
  - `details.fwd.sync_spectrum[i,j,k](nu_comv)` returns the comoving synchrotron specific intensity at any frequency for a given grid cell.
- **Read the SSC spectrum of any grid cell.** When SSC is on, query its comoving intensity per cell.
  - `details.fwd.ssc_spectrum[i,j,k](nu_comv)` returns the comoving SSC specific intensity (when SSC is enabled).
- **Read the Compton Y parameter per cell.** Query Y as a function of electron Lorentz factor for a chosen cell.
  - `details.fwd.Y_spectrum[i,j,k](gamma)` returns the Compton Y parameter as a function of electron Lorentz factor.
- **Query many frequencies at once.** All these evaluators accept NumPy arrays for vectorized queries.

#### ► **Python Object Introspection**

- **Objects now print readable summaries.** Printing a key object shows an informative description instead of a bare memory address.
  - All key Python objects now have informative `__repr__` output: `Observer`, `Radiation`, `Magnetar`, `Model`, `FluxDict`, `Flux`, `ShockDetails`, `SimulationDetails`.
- **Inspect a Model after building it.** Read-only properties let you check a model's configuration once constructed.
  - Read-only properties on `Model` for inspecting state after construction: `model.observer`, `model.fwd_rad`, `model.rvs_rad`, `model.resolutions`, `model.rtol`, `model.axisymmetric`.
- **Inspect Observer and Radiation settings too.** Their key parameters are exposed as read-only properties.
  - Read-only properties on `Observer` (`lumi_dist`, `z`, `theta_obs`, `phi_obs`) and `Radiation` (`eps_e`, `eps_B`, `p`, `xi_e`, `ssc`, `kn`).

#### ► **Improved Inverse Compton Spectrum Computation**

- **More accurate seed photon spectrum.** The frequency grid now concentrates resolution near synchrotron breaks.
  - Adaptive frequency grid concentrates resolution near synchrotron break frequencies (ν_a, ν_m, ν_c, ν_M), improving accuracy of the seed photon spectrum.
- **Smoother SSC light curves and spectra.** The electron energy integration was redesigned.
  - Redesigned electron energy integration produces smoother SSC light curves and spectra.

#### ► **Sky Image Generation**

- **Generate resolved afterglow images.** A new method produces spatial images at any observer time and frequency, for all jet types and viewing angles.
  - New `model.sky_image()` method generates spatially resolved afterglow images at any observer time and frequency, supporting all jet types and viewing angles.
- **Make image sequences and movies efficiently.** Evaluate over many observer times in one batch.
  - Batch evaluation over multiple observer times for efficient multi-frame image sequences and movies.
- See the [sky image example notebook](script/sky-image.ipynb) and [documentation](https://vegasafterglow.readthedocs.io/en/latest/examples/sky_image.html) for usage

#### ► **Interactive Web Tool**

- **Model afterglows in your browser.** A no-install web tool lets you tune parameters in real time for light curves, spectra, and sky images, works on mobile, and accepts your own data.
  - Browser-based real-time, parameter-tunable afterglow modeling tool for light curves, spectra, and sky images, with mobile-phone compatible UI and observational data input/upload support (CSV/TXT/XLS/XLSX) — no installation required: [vegasafterglow.vercel.app](https://www.vegasafterglow.com).

#### ► **Smooth Electron Distribution for IC Integration**

- **Smoother electron distribution for IC.** The electron distribution now uses smooth broken power-law transitions at its breaks, matching the synchrotron photon spectrum.
  - Electron energy distribution now uses smooth broken power-law transitions at cooling and injection breaks, consistent with the synchrotron photon spectrum.

### Performance

#### ► **~5x Faster Computation for Top-Hat and Two-Component Jets**
- **Symmetric jets skip redundant work.** The code detects symmetry and avoids recomputing identical angular slices.
  - Automatic symmetry detection eliminates redundant computation for symmetric jet structures.
- **Top-hat jets compute one point and broadcast.** A single representative grid point is evaluated and reused across the angular grid.
  - Top-hat jets now compute a single representative grid point and broadcast results across the full angular grid.
- **Two-component jets benefit similarly.** Only unique angular groups are computed.
  - Two-component jets similarly benefit by computing only unique angular groups.
- **The speedup covers every stage.** It applies from shock dynamics through photon generation.
  - Applies to all stages: shock dynamics, electron distribution, synchrotron and IC photon generation.

#### ► **Reduced Memory Usage in Flux Computation**
- **Flux integration uses less memory.** Large temporary arrays during observer flux integration are gone.
  - Eliminated large intermediate arrays during observer flux integration.
- **Single-pass accumulation.** Flux is accumulated directly instead of stored and summed in two passes.
  - Direct accumulation replaces the previous two-pass (store-then-sum) approach.

### Changed

#### ► **⚠️ BREAKING: MCMC Fitting Interface Redesigned**

The MCMC fitting module has been completely redesigned. Existing fitting scripts will need to be updated.

- **`VegasMC` is gone — use `Fitter`.** Rewrite scripts that used the old batch evaluator to use the new `Fitter`; performance is equivalent.
  - The C++ batch evaluator `VegasMC` has been replaced by the `Fitter` class, which uses the Python `Model` API internally with equivalent performance.
- **`ObsData` is gone — add data on the fitter.** Replace `ObsData` usage with the fitter's data-adding methods.
  - Use `fitter.add_flux_density()`, `fitter.add_spectrum()`, and `fitter.add_flux()` to add observational data directly to the `Fitter`.
- **`Setups` is gone — pass config to the constructor.** Move all model configuration into `Fitter(...)` keyword arguments.
  - All model configuration is now passed directly to the `Fitter` constructor as keyword arguments: `Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="ism", rvs_shock=True, ...)`.
- **`jet` and `medium` now accept callables.** Alongside the built-in type strings, you can pass factory functions for custom structures.
  - In addition to built-in type strings, you can now pass factory functions that return custom `Ejecta` or `Medium` objects.

See the [MCMC Parameter Fitting](https://vegasafterglow.readthedocs.io/en/latest/mcmc_fitting/index.html) documentation for migration examples.

#### ► **⚠️ BREAKING: `ssc_cooling` Parameter Removed**
- **Turning on SSC now enables its cooling automatically — drop `ssc_cooling`.** If your code passes `ssc_cooling`, remove it; `ssc=True` alone now gives full IC cooling: `Radiation(..., ssc_cooling=True, ssc=True)` → `Radiation(..., ssc=True)`.
  - `ssc_cooling` removed from `Radiation` and `Fitter`; IC cooling is enabled automatically whenever `ssc=True`.

#### ► **Default Resolution**
- **More accurate default resolution.** The default grid resolution changed to give better accuracy out of the box; re-run scripts that relied on the old default if you need to reproduce prior results.
  - Default grid resolution changed from `(0.3, 1, 10)` to `(0.1, 0.25, 10)`.
- **New default is more accurate and comparably fast.** The improved adaptive grid converges with fewer points, so the finer default is not slower.
  - Thanks to the improved adaptive grid algorithm (medium-aware time sampling, jet edge anchoring), the code converges with fewer grid points than before — so the new default is both more accurate and comparably fast.

#### ► **Default MCMC Sampler**
- **Dynesty is now the recommended default.** Nested sampling replaces emcee as the suggested sampler.
  - Dynesty (nested sampling) is now the recommended default sampler, replacing emcee.
- **Emcee sampling is more efficient.** Its move strategy was updated.
  - Emcee move strategy updated to `DEMove + DESnookerMove` for more efficient sampling.

#### ► **Inverse Compton Cooling**
- **More accurate IC cooling in the strong-KN limit.** The Klein-Nishina regime classification was expanded.
  - Expanded Klein-Nishina regime classification (5 regimes, up from 3) for more accurate IC cooling in the strong-KN limit.

#### ► **Wind Medium Defaults**
- **Cleaner defaults for pure-wind media.** `n_ism` and `n0` now default to `None` and are handled internally, matching the common pure-wind case.
  - `Wind(A_star, n_ism=None, n0=None, k_m=2)`: `n_ism` and `n0` now default to `None` instead of `0` and `inf`, handled internally.
- **Wind index renamed to `k_m`.** Update code that used `k`; the name now matches the MCMC fitting parameter.
  - Wind density power-law index renamed from `k` to `k_m` to match the MCMC fitting parameter name.

#### ► **Python >= 3.8 Required**
- **Python 3.8 is now the minimum.** Upgrade if you are still on 3.7.
  - Minimum Python version raised from 3.7 to 3.8.

### Fixed
- **Fixed a small cooling-frequency jump in the reverse shock.** Reverse-shock light curves no longer show a discontinuity at shock crossing.
  - Eliminated small jump in cooling frequency (ν_c) across the shock crossing boundary.
- **Smoother thick-shell reverse shock crossing.** Reverse-shock initial conditions transition more smoothly.
  - Smoother transition for thick-shell reverse shock crossing.
- **Fixed a grid sizing bug.** Grid dimensions are now correct when merged grid sizes are used.
  - Fixed incorrect grid dimensions when using merged grid sizes in `auto_grid`.

---

## [v1.1.0] - 2025-12-08

### Added

#### ► **Bilby MCMC Integration**
- **Run your fits through bilby with multiple samplers.** You can now drive VegasAfterglow fits through the bilby framework and pick the sampler that fits your needs — a fast one for quick posteriors or a slower one when you need Bayesian evidence.
  - Full bilby integration: `emcee` (affine-invariant MCMC ensemble sampler, fast but does not compute evidence), `dynesty` (dynamic nested sampling, slow but computes Bayesian evidence), and all other bilby samplers (`nestle`, `cpnest`, `pymultinest`, `ultranest`, etc.).
- **Corner plots now get properly formatted labels for free.** `FitResult` carries ready-made LaTeX labels, so your corner plots show clean math-formatted parameter names without manual setup.
  - `latex_labels` provides properly formatted labels for corner plots (e.g., `$\log_{10}(E_{\rm iso})$`), with automatic formatting for LOG-scale parameters via the `log10_` prefix.
- **Use multiple CPU cores to speed up fits.** Set `npool` to spread likelihood evaluations across cores and cut wall-clock time.
  - Multi-core support via the `npool` parameter for efficient parallel likelihood evaluations; works with both emcee and dynesty samplers.

#### ► **Flexible Parameter Interface**
- **Log-scale parameters are handled automatically.** If you mark a parameter as `Scale.log`, the sampler works in log10 space and converts back to physical values for you — no manual conversions.
  - Parameters with `Scale.log` are automatically prefixed with `log10_` (e.g., `E_iso` with `Scale.log` becomes `log10_E_iso` in sampler space), with automatic transformation between log10 and physical values.
- **Tune any sampler with its own settings.** Pass sampler-specific options through `**sampler_kwargs`; sensible defaults apply if you pass nothing.
  - Generic sampler configuration via `**sampler_kwargs`. Emcee defaults: `nsteps=5000`, `nburn=1000`, `thin=1`, `nwalkers=2*ndim`. Dynesty defaults: `nlive=500`, `dlogz=0.1`, `sample="rwalk"`. Extensible to any bilby-supported sampler.

### Documentation

#### ► **Comprehensive MCMC Guides**
- **README now shows full bilby examples.** The README walks through bilby fitting end to end, with emcee as the default choice and dynesty when you need evidence.
  - Complete bilby integration examples: emcee positioned as primary option (Option 1), dynesty as secondary option (Option 2) for evidence calculation, and all corner plot examples updated to use `result.latex_labels`.
- **Expanded RST reference docs.** The online docs now cover the full fitting interface and current sampler options with worked examples.
  - `mcmc_fitting.rst`: complete `Fitter.fit()` interface reference with all parameters; `quickstart.rst`: updated with current sampler options and examples; `python_api.rst`: refreshed with bilby integration examples; plus references to bilby documentation for additional samplers.

---

## [v1.0.3] - 2025-09-29

### Changed

#### ► **Default Physics Settings**
- **More stable runs by default near strong absorption.** Self-absorption heating is now off by default, which improves numerical stability; you can re-enable it if you need it.
  - Disabled self-absorption heating by default for improved numerical stability.
- **Better handling of strongly self-absorbed regions.** Strong-absorption regions are resolved more finely, giving more reliable results there.
  - Enhanced subsegmentation (5/2 ratio) for improved handling of strong self-absorption regions.

#### ► **MCMC Framework Improvements**
- **Better parameter-space exploration.** The sampler explores the parameter space more effectively, so your chains mix better.
  - Implemented a new MCMC move strategy for better parameter space exploration.
- **More robust sampling from the start.** Walkers now start with a wider initial spread, reducing the chance of getting stuck early.
  - Increased initial parameter spread for more robust sampling.
- **General medium configurations now supported.** You can use the `-k` parameter for general medium setups, not just the previously supported cases.
  - Extended `-k` parameter support for general medium configurations.

### Fixed

#### ► **Model Interface**
- **Consistent flux units across interfaces.** Model flux is now reported in consistent units no matter which interface you use, so results line up.
  - Fixed model flux unit handling for consistent calculations across interfaces.
- **More stable reverse-shock evolution.** Reverse shock calculations are steadier and less prone to numerical trouble.
  - Implemented hardcut Gamma treatment for more stable reverse shock evolution.

### Performance

#### ► **Inverse Compton Optimization**
- **Faster inverse Compton calculations.** IC computations skip unnecessary work and finish sooner.
  - Added an early break condition for inverse Compton calculations to improve computational efficiency.

### Documentation

#### ► **Enhanced Guidelines**
- **New MCMC fitting guidelines.** Added guidance to help you set up and run MCMC parameter fits.
  - Added comprehensive guidelines for MCMC parameter fitting.
- **Updated interface docs.** Documentation now reflects the current interfaces and parameter usage.
  - Updated documentation for interface consistency and parameter usage.
- **More examples and troubleshooting.** Expanded code examples and troubleshooting material.
  - Enhanced examples and troubleshooting documentation.

### Development

#### ► **Code Quality**
- **Consistent interface names across Python bindings.** Interface names were cleaned up for consistency; if you relied on older names, check your calls.
  - Major cleanup for interface name consistency across Python bindings.
- **General code cleanup.** Internal cleanup and optimization across the codebase; no action needed.
  - General code cleanup and optimization throughout the codebase.

---

## [v1.0.2] - 2025-09-15

### Changed

#### ► **⚠️ BREAKING: Python Interface Method Name Updates**
- **Flux methods were renamed — update your calls.** The Python flux methods now have clearer, consistent names; rename any calls in your scripts (see Migration Notes below for the exact mapping).
  - PyModel: `specific_flux()` → `flux_density_grid()` (multi-dimensional flux density calculations), `specific_flux_series()` → `flux_density()` (time-series flux density calculations), `specific_flux_series_with_expo()` → `flux_density_exposures()` (exposure-averaged flux density). VegasMC: `specific_flux()` → `flux_density_grid()` (matches PyModel interface).

#### ► **⚠️ BREAKING: Data Input Parameter Name Cleanup**
- **Data-input parameter names are simpler — update your `ObsData` calls.** The `_cgs` suffixes are gone from the `ObsData` methods (values are still CGS); rename the keyword arguments where you add data.
  - `add_flux_density(nu_cgs, t_cgs, Fnu_cgs, Fnu_err, weights=None)` → `add_flux_density(nu, t, f_nu, err, weights=None)`. `add_flux(nu_min, nu_max, num_points, t_cgs, F, F_err, weights=None)` → `add_flux(band, t, flux, err, num_points=15, weights=None)` — `band` is a `(nu_min, nu_max)` tuple; use `units.band("XRT")` for named bands. `add_spectrum(t_cgs, nu_cgs, Fnu_cgs, Fnu_err, weights=None)` → `add_spectrum(t, nu, f_nu, err, weights=None)`.

### Added

#### ► **Enhanced Documentation**
- **Every interface now documents its CGS units.** You can see the expected units for each physical quantity directly in docstrings and examples, so you no longer have to guess.
  - Added CGS-unit docstrings to pybind11 methods for all physical quantities; updated README, documentation, and notebook examples with unit comments; consistent notation `nu [Hz]`, `t [s]`, `f_nu [erg/cm²/s/Hz]`, `flux [erg/cm²/s]`.

#### ► **Interface Consistency Verification**
- **Parameter names now match the docs across the board.** The Python bindings and documentation were checked against each other, so the names you read in the docs are the names the code actually uses.
  - Verified all pybind11 interface parameters against documentation; all ObsData methods confirmed to use the simplified parameter names; method names verified consistent across interfaces; complete synchronization between C++ bindings and Python documentation.

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
- **You can now fit frequency-integrated observations.** MCMC fitting accepts broadband flux measured over a frequency range, so you can model bolometric and filter-integrated data.
  - New `MultiBandData.add_light_curve()` overload for broadband flux measurements over frequency ranges; enables modeling of bolometric observations and filter-integrated measurements, with improved handling of observational data where the effective frequency depends on spectral shape.

#### ► **Documentation & Examples**
- **Expanded API documentation.** The docs now cover parameters, units, typical ranges, and worked examples, so you can look up how to use each feature.
  - Added extensive documentation with detailed physics explanations: complete parameter descriptions with units and typical ranges, detailed method documentation with usage examples, and enhanced code examples demonstrating advanced features.

### Changed

#### ► **⚠️ BREAKING: Return Object Interface Redesign**
- **Results are now accessed by name instead of dictionary keys — update how you read outputs.** `specific_flux()` and `details()` return structured objects with named attributes (with IDE autocomplete), so replace `results['sync']`-style access with `results.fwd.sync`-style access.
  - `specific_flux()` now returns `FluxDict` objects with named members: `results.total` (total flux array), `results.fwd.sync` (forward shock synchrotron), `results.fwd.ssc` (forward shock SSC), `results.rvs.sync` (reverse shock synchrotron), `results.rvs.ssc` (reverse shock SSC). `details()` now returns `SimulationDetails` objects with named members: `details.phi`, `details.theta`, `details.t_src` (coordinate arrays), `details.fwd.*` (forward shock quantities: `t_obs`, `Gamma`, `r`, `B_comv`, etc.), and `details.rvs.*` (reverse shock quantities, if enabled). All return objects have well-defined attributes for better development experience.

#### ► **API Consistency Improvements**
- **Several parameter names changed to snake_case — rename them in your code.** SSC/KN and forward/reverse radiation arguments were standardized; if you set any of the old names, switch to the new ones.
  - `ConfigParams`: `fwd_SSC` → `fwd_ssc`, `rvs_SSC` → `rvs_ssc`, `KN` → `kn`. `PyRadiation`: `SSC` → `ssc`, `KN` → `kn` (IC cooling is now automatically enabled when `ssc=True`). `Model` constructor: `forward_rad` → `fwd_rad`, `reverse_rad` → `rvs_rad`. All parameter names now follow the consistent snake_case convention.

#### ► **Enhanced Code Documentation**
- **More in-code documentation.** Major classes and structures now carry detailed comments explaining parameters and methods.
  - Added comprehensive comments to all major classes and structures: complete parameter descriptions with physics context, method documentation with detailed explanations, and a consistent documentation style following project conventions.

### Fixed

#### ► **Interface Consistency**
- **Docs and examples match the new parameter names.** README and tutorials were updated so copy-pasted examples work with this release.
  - Updated all documentation and examples to reflect new parameter names: README.md examples updated with new parameter conventions, tutorial documentation fully synchronized with API changes, and all code examples verified for consistency.

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
- **Every built-in model can now be fit with MCMC.** You can run MCMC fits against any of the built-in jet and medium models, not just a subset. See the documentation for details.
  - All built-in jet and medium models are now supported in MCMC fitting, including reverse shock, inverse Compton, magnetar injection, and etc.

#### ► **Adaptive Mesh Generation**
- **Fits use an angular grid tailored to your jet and viewing angle.** Grid resolution is placed where it matters, so runs are more efficient without extra tuning from you.
  - New adaptive angular grid generation based on jet properties and viewing angles. Grid points are distributed according to Doppler boosting factors for optimal efficiency.
- **Reverse-shock fits converge faster.** Reverse-shock runs reach convergence roughly 5x sooner.
  - ~5x faster convergence for reverse shocks.

#### ► **New Jet Models**
- **New StepPowerLawJet model for a uniform core with power-law wings.** Use it when your jet has a sharp core-to-wing transition.
  - StepPowerLawJet: uniform core with a sharp transition to power-law wings for realistic jet structures.
- **TwoComponentJet now supports fully independent narrow and wide components.** You can set separate energy and Lorentz factor profiles for each component.
  - Enhanced TwoComponentJet: separate narrow and wide components with independent energy and Lorentz factor profiles.
- **PowerLawJet now takes separate angular indices for energy and Lorentz factor.** Control the two angular dependencies independently instead of with one shared index.
  - Improved PowerLawJet: split power-law indices for energy (`k_e`) and Lorentz factor (`k_g`) angular dependence.

### **Performance & Computational Improvements**

#### ► **Shock Physics Enhancements**
- **Some internal variable names changed.** If you referenced these names directly, update them.
  - Major refactoring for clarity: `EAT_fwd` → `t_obs_fwd`, `Gamma_rel` → `Gamma_th`.
- **Reverse-shock dynamics rewritten into one unified model.** Shock crossing and post-crossing evolution now use a single, consistent treatment.
  - Major code refactoring for reverse shock dynamics: a unified model for shock crossing and post-crossing evolution.
- **More accurate reverse-shock emission during crossing — expect weaker reverse shocks.** If you fit reverse-shock data, re-evaluate your models, since predicted reverse-shock emission is now weaker than before.
  - Tracks the internal energy evolution during shock crossing for improved accuracy with more accurate shock heating and adiabatic cooling. **Note: The enhanced adiabatic cooling and detailed shock heating treatment leads to even weaker reverse shock emission compared to previous versions.**

#### ► **Numerical & Memory Optimizations**
- **Lower memory use.** Large runs need less memory.
  - Reduced memory footprint through optimized array operations and memory access patterns.
- **Repeated calculations are cached.** Frequently used grid values are reused instead of recomputed.
  - Enhanced caching strategies for frequently used calculations via grid pre-computation.

###  **Enhanced Python Interface**

#### ► **Data Management Improvements**
- **MultiBandData now handles light curves and spectra together with flexible weighting.** Manage both data types through one unified interface.
  - MultiBandData redesign: unified handling of light curves and spectra with flexible weighting.
- **New methods to compute flux at specific time-frequency pairs.** Get flux for exactly the epochs and frequencies you specify.
  - New Series Calculations methods for calculating flux at specific time-frequency pairs.
- **Large multi-wavelength datasets use memory more efficiently.** Big datasets take less memory to hold.
  - Optimized data structures for large multi-wavelength datasets.

###  **Development & Build System**

#### ► **Code Quality Enhancements**
- **Stronger compile-time checks under the hood.** Improves internal type safety and memory management; no action needed.
  - Enhanced template system: improved compile-time type checking and memory management.
- **Expanded API documentation.** Parameter descriptions are now more detailed.
  - Comprehensive API documentation with detailed parameter descriptions.
- **New MCMC tutorials with real-data fitting.** Updated examples walk through fitting real data.
  - New MCMC tutorials with real data fitting demonstrations.

### **API Changes & Migration**

#### ► **Parameter Interface Updates**
- **PowerLawJet: replace the single `k` with `k_e` and `k_g`.** Update your calls to pass separate energy and Lorentz factor indices.
  - Split single `k` parameter into separate `k_e` (energy) and `k_g` (Lorentz factor) indices for more flexible modeling.
- **TwoComponentJet parameter names are standardized.** Use the new names for the core and wide components in your calls.
  - Parameter names standardized to `(theta_c, E_iso, Gamma0, theta_w, E_iso_w, Gamma0_w)` for core and wide components.
- **New StepPowerLawJet model with its own parameter set.** Instantiate it with the listed core and power-law wing parameters.
  - New jet model with parameters `(theta_c, E_iso, Gamma0, E_iso_w, Gamma0_w, k_e, k_g)` for core and power-law wing components.
- **Wind medium can now include a stratified ISM/inner region.** Pass the optional parameters to model a stratified medium.
  - Wind extended with optional `n_ism` and `n_0` parameters for stratified medium modeling: `Wind(A_star, n_ism=0, n_0=inf)`.
- **Medium no longer takes a separate mass argument.** Drop the `mass` argument from your `Medium` calls.
  - Simplified from `Medium(rho, mass)` to `Medium(rho)` — removed separate mass parameter.
- **Model flux methods changed — update these calls.** One method was removed, one added for exposure-time averaging, and the `details` signature changed.
  - Removed `specific_flux_sorted_series()`; added `specific_flux_series_with_expo(t, nu, expo_time, num_points=10)` for exposure time averaging; changed `details(t_obs)` to `details(t_min, t_max)`.
- **Default resolution changed.** Fits use a new default resolution for a better performance/accuracy balance; results may shift slightly unless you set it explicitly.
  - Default resolution changed from `(0.3, 3.0, 5.0)` to `(0.1, 0.25, 10)`.
- **MCMC parameter set restructured across all jet, medium, and reverse-shock options.** Review your parameter lists: several were added, reverse-shock parameters were introduced, and `k_jet` was removed.
  - Major restructuring with new parameters for all jet types, medium configurations, and reverse shock physics. Added: `theta_v, n_ism, n0, A_star, k_e, k_g, duration, E_iso_w, Gamma0_w, theta_w, L0, t0, q`. Added reverse shock parameters: `p_r, eps_e_r, eps_B_r, xi_e_r`. Removed: `k_jet` (replaced by `k_e, k_g`).
- **MCMC data can now be weighted.** Pass the optional `weights` to `MultiBandData` for both light curves and spectra.
  - Enhanced `MultiBandData` with optional `weights` parameter for both light curves and spectra.
- **MCMC models now expose a single flux method.** Replace `light_curves()` and `spectra()` calls with the unified `specific_flux()`.
  - Replaced separate `light_curves()` and `spectra()` methods with unified `specific_flux()` method.

---

**Migration Notes**: This release includes some API changes. Most existing code will work with minimal modifications. See the documentation for detailed migration guidance.

**⚠️ Important Physics Changes**: With the improved reverse shock dynamics, we now find that the reverse shock emission is even weaker than what we reported in our previous code paper, which itself was already weaker than the analytical scalings. This further reduction is due to enhanced adiabatic cooling and a more detailed treatment of shock heating. Users fitting reverse shock data may therefore need to re-evaluate their models and parameters.

## [v0.2.8] - 2025-08-15

### Improved

- **Faster inverse Compton.** Inverse Compton scattering calculations run about 10x faster.
  - Major performance optimization (~10x speedup) for inverse Compton scattering calculations.
- **Photon parameter renamed to `I_nu_max`.** If your code uses `P_nu_max` in the photon interface, rename it to `I_nu_max`.
  - Changed parameter name from `P_nu_max` to `I_nu_max` for better consistency in photon interface.
- **Steadier spectra near `nu_max`.** Spectrum smoothing at `nu_max` is more numerically stable.
  - Enhanced spectrum smoothing at `nu_max` for better numerical stability.

### Documentation

- **Docs for simulation quantities evolution.** Updated documentation for detailed simulation quantities evolution.
- **More examples and API docs.** Enhanced examples and API documentation.

## [v0.2.7] - 2025-07-31

### Added

- **Inspect internal quantities during a run.** New Python interface lets you check how internal simulation quantities evolve under various reference frames.
  - New Python interface to check the evolution of internal simulation quantities under various reference frames.
- **Docs for simulation quantities evolution.** Comprehensive documentation for detailed simulation quantities evolution.

### Fixed

- **Faster runs from correct jet edge detection.** A misconfiguration made the computational domain unnecessarily large for every jet profile; runs are now sized correctly.
  - Fixed jet edge detection function that was mistakenly set to 90 degrees, making the computational domain unnecessarily large regardless of the jet profile.

### Documentation

- **Docs for simulation quantities evolution.** Updated documentation for detailed simulation quantities evolution.
- **Examples for tracking shocks and microphysics.** Enhanced examples showing how to track shock dynamics and microphysical parameters.

## [v0.2.6] - 2025-07-19

### Added

- **Built-in two-component jet.** Added support for built-in two-component jet configurations.

## [v0.2.5] - 2025-07-19

### Improved

- **Better jet edge detection for custom profiles.** Improved jet edge detection algorithm for user-defined jet profiles.

### Documentation

- **Updated README.** Updated README with latest features.
- **More examples.** Enhanced documentation with additional examples.
- **Updated API docs.** Updated API documentation.

## [v0.2.4] - 2025-06-28

### Added

- **Magnetar spin-down docs.** Added comprehensive documentation for magnetar spin-down energy injection.

### Documentation

- **Docs for magnetar energy injection.** Updated documentation for magnetar energy injection features.
- **More usage examples in README.** Enhanced README with additional usage examples.
- **Improved API reference.** Improved API reference documentation.

## [v0.2.3] - 2025-06-27

### Changed

- **Breaking: consistent units for custom medium/jet.** If you define your own medium or jet at the Python level, check your unit conventions against the new system.
  - API Breaking Change: modified Python-level user-defined medium/jet unit system for consistency.
- **Jet duration parameter renamed to `duration`.** Rename `T0` to `duration` in your jet setup.
  - Changed jet duration parameter name from `T0` (confusing to observer frame) to `duration`.

### Improved

- **Sorted flux density series.** Sorted series flux density for better data organization.

## [v0.2.2] - 2025-06-23

### Changed

- **New default resolution for faster runs.** Changed default resolution settings when unspecified for better performance.

### Improved

- **Smoother reverse shock.** Enhanced reverse shock smoothing algorithms.

## [v0.2.1] - 2025-06-22

### Improved

- **Better reverse shock modeling.** Reverse shock smoothing, post-crossing calibration, and single-electron peak power calibration are all improved.
  - Enhanced reverse shock smoothing; improved post-crossing reverse shock calibration; better single electron peak power calibration.

## [v0.1.9] - 2025-06-19

### Changed

- **Python 3.7 support dropped.** Removed support for Python 3.7; the minimum requirement is now Python 3.8+.

### Fixed

- **Reverse shock corrections.** Significant corrections to reverse shock calculations.

### Documentation

- **Updated docs and examples.** Updated documentation and examples.
- **Expanded API reference.** Enhanced API reference.

## [v0.1.8] - 2025-06-09

### Improved

- **Reverse shock refinements.** Refined reverse shock calculations and algorithms.
- **Code cleanup.** General code cleanup and optimization.

## [v0.1.7] - 2025-05-24

### Fixed

- **macOS build fixes.** Fixed macOS-specific bugs.
- **Deep Newtonian regime refinement.** Improved accuracy of calculations in the deep Newtonian regime.
  - Shocked electron calculations moved from energy space to momentum space; enhanced self-absorption refinement; interface improvements.

### Improved

- **Code quality updates.** General updates and optimizations.

## [v0.1.6] - 2025-05-15

### Added

- **Energy injection Python bindings.** Added the missing Python bindings header for energy injection functionality.

### Improved

- **Faster magnetar integration.** Enhanced C-level magnetar binding for better performance.

## [v0.1.5] - 2025-05-14

### Added

- **Inverse Compton enhancements.** Updated inverse Compton calculations and added a magnetar injection Python interface.
  - Updated IC spectrum calculations; added magnetar injection Python interface; enhanced IC code with cleanup and optimization.

### Changed

- **Python version support updated.** Removed Python 3.6 support and added Python 3.7 support.

### Fixed

- **Unit handling fixes.** Corrected unit handling in various calculations.
- **Build fixes.** Various build fixes and improvements.

### Documentation

- **Updated documentation examples.** Updated examples in documentation.
- **README updates.** Enhanced README with new features.
- **Improved API docs.** Improved API documentation.

## [v0.1.4] - 2025-05-11

### Added

- **Improved build system.** Improved CMake configuration and cross-platform support.

### Fixed

- **Windows build fixes.** Resolved various Windows build issues.
- **External dependency fixes.** Fixed external library integration.

## [v0.1.3] - 2025-05-07

### Improved

- **Build configuration.** Enhanced CMake build configuration.
- **Documentation.** Updated README and logo integration.

### Fixed

- **Cross-platform build fixes.** Resolved various build issues across platforms.

## [v0.1.2] - 2025-05-06

### Added

- **Documentation system.** Added a comprehensive documentation system.
- **Logo and branding.** Added project logo and branding.

### Improved

- **Code structure.** Enhanced code organization and commenting.
- **Build workflow.** Improved build workflow and packaging.

## [v0.1.1] - 2025-05-06

### Fixed

- **Build packaging fixes.** Fixed wheel building and packaging issues.
- **Installation instructions.** Updated README with installation instructions.

### Improved

- **README layout.** Enhanced README layout and organization.

## [v0.1.0] - 2025-05-05

### Added

- **Initial release.** First public release of VegasAfterglow.
- **Core features.** GRB afterglow modeling with forward and reverse shock dynamics, radiation, structured jets, and MCMC fitting.
  - High-performance C++ framework with Python interface; comprehensive GRB afterglow modeling; forward and reverse shock dynamics; synchrotron and inverse Compton radiation; structured jet configurations; MCMC parameter fitting capabilities.
- **Radiation mechanisms.** Synchrotron, inverse Compton, and SSC emission.
  - Synchrotron emission with self-absorption; inverse Compton scattering with Klein-Nishina corrections; Synchrotron Self-Compton (SSC).
- **Physical models.** Choice of ambient medium types, jet structures, shock regimes, and injection.
  - Multiple ambient medium types (ISM, wind, user-defined); various jet structures (top-hat, power-law, Gaussian, user-defined); relativistic and non-relativistic shock evolution; energy and mass injection.
- **Fast light curves.** Ultra-fast light curve computation on millisecond timescales.
- **Cross-platform support.** Linux, macOS, and Windows compatibility.
- **Python interface.** User-friendly Python bindings for easy integration.

### Documentation

- **API documentation.** Comprehensive API documentation.
- **Installation guides.** Installation guides for Python and C++.
- **Examples and tutorials.** Usage examples and tutorials.
- **Quick start guide.** Quick start guide with Jupyter notebooks.

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
