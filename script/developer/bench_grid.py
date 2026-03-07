"""Grid convergence and resolution benchmark tool.

Usage:
    python bench_grid.py                # quick accuracy check (default)
    python bench_grid.py --theta        # standard theta grid sweep
    python bench_grid.py --extended     # extended parameter space
    python bench_grid.py --phi          # phi grid sweep
    python bench_grid.py --resolutions  # resolution parameter sweep
    python bench_grid.py --combined     # full 5D combined sweep
"""

import argparse
import csv
import itertools
import time
from collections import defaultdict

import numpy as np

from VegasAfterglow import (
    GaussianJet, ISM, Model, Observer, PowerLawJet, Radiation,
    TophatJet, TwoComponentJet, Wind,
)

# ── Shared constants ───────────────────────────────────────────────────────────
THETA_C = 0.1
E_ISO = 1e52
GAMMA0 = 300
LUMI_DIST = 1e28
Z = 1.0
STANDARD_RAD = dict(eps_e=0.1, eps_B=0.01, p=2.3, xi_e=1.0)

TIMES = np.logspace(0, 8, 100)
BANDS = {"Radio": 1e9, "Optical": 4.84e14, "X-ray": 1e18}
SSC_BANDS = {**BANDS, "TeV": 2.4e26}
NUS_GRID = np.array([6e9, 4.56e14, 2.4e17])
RES_DEFAULT = (0.1, 0.5, 10)


# ── Helpers ────────────────────────────────────────────────────────────────────
def make_jet(jet_type, theta_c=THETA_C, gamma0=GAMMA0):
    if jet_type == "tophat":
        return TophatJet(theta_c=theta_c, E_iso=E_ISO, Gamma0=gamma0)
    elif jet_type == "gaussian":
        return GaussianJet(theta_c=theta_c, E_iso=E_ISO, Gamma0=gamma0)
    elif jet_type == "powerlaw":
        return PowerLawJet(theta_c=theta_c, E_iso=E_ISO, Gamma0=gamma0, k_e=2, k_g=2)
    raise ValueError(f"Unknown jet type: {jet_type}")


def max_rel_error(test, ref):
    """Max relative error with slope filtering to skip steep transitions."""
    valid = (ref > 0) & np.isfinite(ref) & (test > 0) & np.isfinite(test)
    if not np.any(valid):
        return np.nan
    log_ref = np.full_like(TIMES, np.nan)
    log_ref[ref > 0] = np.log10(ref[ref > 0])
    slope = np.gradient(log_ref, np.log10(TIMES))
    valid = valid & (np.abs(slope) <= 4)
    if np.any(valid):
        valid = valid & (ref > np.max(ref[valid]) * 1e-6)
    if not np.any(valid):
        return np.nan
    return float(np.max(np.abs(test[valid] - ref[valid]) / ref[valid]))


def simple_max_rel_error(test, ref):
    mask = np.abs(ref) > 1e-5 * np.max(np.abs(ref))
    if not np.any(mask):
        return 0.0
    return float(np.max(np.abs(test[mask] - ref[mask]) / np.abs(ref[mask])))


def build_model(case, res, **extra):
    jet = case["jet_factory"]()
    med = case["medium_factory"]()
    obs = Observer(lumi_dist=LUMI_DIST, z=Z, theta_obs=case["theta_obs"])
    fwd_rad = Radiation(**case["fwd_kw"])
    rvs_rad = Radiation(**case["rvs_kw"]) if case["rvs_kw"] else None
    return Model(jet, med, obs, fwd_rad, rvs_rad=rvs_rad, resolutions=res, **extra)


def extract_components(flux, components):
    accessors = {
        "fwd.sync": lambda f: np.asarray(f.fwd.sync),
        "fwd.ssc":  lambda f: np.asarray(f.fwd.ssc),
        "rvs.sync": lambda f: np.asarray(f.rvs.sync),
    }
    return {c: accessors[c](flux) for c in components if c in accessors}


def compute_refs(cases, res):
    refs = {}
    for i, case in enumerate(cases):
        t0 = time.time()
        model = build_model(case, res)
        ref_data = {}
        for band_name, nu in case["bands"].items():
            flux = model.flux_density(TIMES, np.full_like(TIMES, nu))
            ref_data[band_name] = extract_components(flux, case["components"])
        refs[i] = ref_data
        print(f"  [{i+1}/{len(cases)}] {case['label']} ({time.time()-t0:.1f}s)")
    return refs


# ── Case definitions ───────────────────────────────────────────────────────────
def make_theta_cases():
    """Standard cases: Gamma=300, all jets, ISM/wind, sync/SSC/rvs."""
    cases = []
    jet_factories = {
        "tophat":   (THETA_C, lambda: TophatJet(theta_c=THETA_C, E_iso=E_ISO, Gamma0=GAMMA0)),
        "gaussian": (THETA_C, lambda: GaussianJet(theta_c=THETA_C, E_iso=E_ISO, Gamma0=GAMMA0)),
        "powerlaw": (THETA_C, lambda: PowerLawJet(theta_c=THETA_C, E_iso=E_ISO, Gamma0=GAMMA0, k_e=2, k_g=2)),
        "two_comp": (0.05,    lambda: TwoComponentJet(
            theta_c=0.05, E_iso=E_ISO, Gamma0=GAMMA0, theta_w=0.3, E_iso_w=1e50, Gamma0_w=50)),
    }
    for jname, (tc, jf) in jet_factories.items():
        for ratio in [0, 2, 4]:
            cases.append(dict(
                label=f"{jname}/ISM/sync/{ratio}tc",
                jet_factory=jf, medium_factory=lambda: ISM(n_ism=1.0),
                fwd_kw={**STANDARD_RAD, "ssc": False}, rvs_kw=None,
                theta_obs=ratio * tc, components=["fwd.sync"], bands=BANDS,
            ))
    for ratio in [0, 2]:
        cases.append(dict(
            label=f"gaussian/wind/sync/{ratio}tc",
            jet_factory=lambda: GaussianJet(theta_c=THETA_C, E_iso=E_ISO, Gamma0=GAMMA0),
            medium_factory=lambda: Wind(A_star=0.1),
            fwd_kw={**STANDARD_RAD, "ssc": False}, rvs_kw=None,
            theta_obs=ratio * THETA_C, components=["fwd.sync"], bands=BANDS,
        ))
    for ratio in [0, 2]:
        cases.append(dict(
            label=f"gaussian/ISM/ssc/{ratio}tc",
            jet_factory=lambda: GaussianJet(theta_c=THETA_C, E_iso=E_ISO, Gamma0=GAMMA0),
            medium_factory=lambda: ISM(n_ism=1.0),
            fwd_kw={**STANDARD_RAD, "ssc": True}, rvs_kw=None,
            theta_obs=ratio * THETA_C, components=["fwd.sync", "fwd.ssc"], bands=SSC_BANDS,
        ))
    rvs_kw = {**STANDARD_RAD, "ssc": False, "kn": False}
    for ratio in [0, 2]:
        cases.append(dict(
            label=f"tophat/ISM/rvs/{ratio}tc",
            jet_factory=lambda: TophatJet(theta_c=THETA_C, E_iso=E_ISO, Gamma0=GAMMA0, duration=1),
            medium_factory=lambda: ISM(n_ism=1.0),
            fwd_kw={**STANDARD_RAD, "ssc": False}, rvs_kw=rvs_kw,
            theta_obs=ratio * THETA_C, components=["fwd.sync", "rvs.sync"], bands=BANDS,
        ))
    return cases


def make_extended_cases():
    """Extreme parameter space: high/low Gamma, wide/narrow theta_c."""
    cases = []
    rvs_kw = {**STANDARD_RAD, "ssc": False, "kn": False}

    def _sync(label, jf, theta_obs):
        return dict(label=label, jet_factory=jf, medium_factory=lambda: ISM(n_ism=1.0),
                    fwd_kw={**STANDARD_RAD, "ssc": False}, rvs_kw=None,
                    theta_obs=theta_obs, components=["fwd.sync"], bands=BANDS)

    def _rvs(label, jf, theta_obs):
        return dict(label=label, jet_factory=jf, medium_factory=lambda: ISM(n_ism=1.0),
                    fwd_kw={**STANDARD_RAD, "ssc": False}, rvs_kw=rvs_kw,
                    theta_obs=theta_obs, components=["fwd.sync", "rvs.sync"], bands=BANDS)

    # High Gamma
    G, TC = 1000, 0.1
    for ratio in [0, 2]:
        cases.append(_sync(f"tophat/highG/sync/{ratio}tc",   lambda G=G, TC=TC: TophatJet(theta_c=TC, E_iso=E_ISO, Gamma0=G),   ratio * TC))
        cases.append(_sync(f"gaussian/highG/sync/{ratio}tc", lambda G=G, TC=TC: GaussianJet(theta_c=TC, E_iso=E_ISO, Gamma0=G), ratio * TC))
        cases.append(_rvs(f"tophat/highG/rvs/{ratio}tc",    lambda G=G, TC=TC: TophatJet(theta_c=TC, E_iso=E_ISO, Gamma0=G, duration=1), ratio * TC))

    # Low Gamma
    G, TC = 5, 0.5
    cases.append(_sync("tophat/lowG/sync/0tc",   lambda G=G, TC=TC: TophatJet(theta_c=TC, E_iso=E_ISO, Gamma0=G),   0))
    cases.append(_sync("gaussian/lowG/sync/0tc", lambda G=G, TC=TC: GaussianJet(theta_c=TC, E_iso=E_ISO, Gamma0=G), 0))
    cases.append(_rvs("tophat/lowG/rvs/0tc",    lambda G=G, TC=TC: TophatJet(theta_c=TC, E_iso=E_ISO, Gamma0=G, duration=1), 0))

    # Wide theta_c
    TC, G = 0.5, 300
    for ratio in [0, 2]:
        cases.append(_sync(f"tophat/wide/sync/{ratio}tc",   lambda G=G, TC=TC: TophatJet(theta_c=TC, E_iso=E_ISO, Gamma0=G),   ratio * TC))
        cases.append(_sync(f"gaussian/wide/sync/{ratio}tc", lambda G=G, TC=TC: GaussianJet(theta_c=TC, E_iso=E_ISO, Gamma0=G), ratio * TC))

    # Narrow theta_c
    TC, G = 0.01, 1000
    cases.append(_sync("tophat/narrow/sync/0tc",   lambda G=G, TC=TC: TophatJet(theta_c=TC, E_iso=E_ISO, Gamma0=G),   0))
    cases.append(_sync("gaussian/narrow/sync/0tc", lambda G=G, TC=TC: GaussianJet(theta_c=TC, E_iso=E_ISO, Gamma0=G), 0))

    return cases


def make_phi_cases():
    """Off-axis cases for phi grid sweep."""
    cases = []
    jet_factories = {
        "tophat":   (THETA_C, lambda: TophatJet(theta_c=THETA_C, E_iso=E_ISO, Gamma0=GAMMA0)),
        "gaussian": (THETA_C, lambda: GaussianJet(theta_c=THETA_C, E_iso=E_ISO, Gamma0=GAMMA0)),
        "powerlaw": (THETA_C, lambda: PowerLawJet(theta_c=THETA_C, E_iso=E_ISO, Gamma0=GAMMA0, k_e=2, k_g=2)),
        "two_comp": (0.05,    lambda: TwoComponentJet(
            theta_c=0.05, E_iso=E_ISO, Gamma0=GAMMA0, theta_w=0.3, E_iso_w=1e50, Gamma0_w=50)),
    }
    for jname, (tc, jf) in jet_factories.items():
        for ratio in [2, 4]:
            cases.append(dict(
                label=f"{jname}/ISM/sync/{ratio}θc",
                jet_factory=jf, medium_factory=lambda: ISM(n_ism=1.0),
                fwd_kw={**STANDARD_RAD, "ssc": False}, rvs_kw=None,
                theta_obs=ratio * tc, components=["fwd.sync"], bands=BANDS,
            ))
    cases.append(dict(
        label="gaussian/wind/sync/2θc",
        jet_factory=lambda: GaussianJet(theta_c=THETA_C, E_iso=E_ISO, Gamma0=GAMMA0),
        medium_factory=lambda: Wind(A_star=0.1),
        fwd_kw={**STANDARD_RAD, "ssc": False}, rvs_kw=None,
        theta_obs=2 * THETA_C, components=["fwd.sync"], bands=BANDS,
    ))
    rvs_kw = {**STANDARD_RAD, "ssc": False, "kn": False}
    cases.append(dict(
        label="tophat/ISM/rvs/2θc",
        jet_factory=lambda: TophatJet(theta_c=THETA_C, E_iso=E_ISO, Gamma0=GAMMA0, duration=1),
        medium_factory=lambda: ISM(n_ism=1.0),
        fwd_kw={**STANDARD_RAD, "ssc": False}, rvs_kw=rvs_kw,
        theta_obs=2 * THETA_C, components=["fwd.sync", "rvs.sync"], bands=BANDS,
    ))
    return cases


# ── Runners ────────────────────────────────────────────────────────────────────
def run_quick():
    """Quick accuracy check: default resolution vs high-res reference."""
    RES_TEST = (0.1, 0.5, 5)
    RES_REF  = (0.3, 3.0, 20)
    jet_types = ["tophat", "gaussian", "powerlaw"]
    obs_factors = [0, 1, 2, 3]

    print("Computing references...")
    refs = {}
    for jt in jet_types:
        for of in obs_factors:
            t0 = time.time()
            m = Model(make_jet(jt), ISM(n_ism=1.0),
                      Observer(lumi_dist=LUMI_DIST, z=Z, theta_obs=of * THETA_C),
                      Radiation(eps_e=0.1, eps_B=0.01, p=2.2), resolutions=RES_REF)
            refs[(jt, of)] = np.asarray(m.flux_density_grid(TIMES, NUS_GRID).total)
            print(f"  ref {jt}/θ={of}θc: {time.time()-t0:.1f}s")

    print(f"\n{'jet':>10s} {'θ_obs':>6s} | {'radio':>8s} {'optical':>8s} {'X-ray':>8s} | {'worst':>8s}")
    print("-" * 60)
    worst_all = 0
    for jt in jet_types:
        for of in obs_factors:
            m = Model(make_jet(jt), ISM(n_ism=1.0),
                      Observer(lumi_dist=LUMI_DIST, z=Z, theta_obs=of * THETA_C),
                      Radiation(eps_e=0.1, eps_B=0.01, p=2.2), resolutions=RES_TEST)
            test = np.asarray(m.flux_density_grid(TIMES, NUS_GRID).total)
            ref = refs[(jt, of)]
            errs = [simple_max_rel_error(test[i], ref[i]) for i in range(3)]
            worst = max(errs)
            worst_all = max(worst_all, worst)
            print(f"{jt:>10s} {of:>5d}θc | {errs[0]:>7.1%} {errs[1]:>7.1%} {errs[2]:>7.1%} | {worst:>7.1%}")
    print(f"\nOverall worst: {worst_all:.1%}")


def run_theta_sweep(cases, sweep_values, check_val=32):
    """Sweep min_theta_num and report per-case convergence."""
    print(f"\n{len(cases)} cases × {len(sweep_values)} min_theta_num values\n")
    print("Computing references...")
    refs = compute_refs(cases, RES_DEFAULT)

    label_w = max(len(f"{c['label']}/{comp}") for c in cases for comp in c["components"]) + 2
    col_w   = max(8, max(len(f"mtn={m}") for m in sweep_values) + 2)
    header  = f"\n{'':>{label_w}s} " + "".join(f"{'mtn='+str(m):>{col_w}s}" for m in sweep_values) + f" | {'n_t':>4s}"
    print(header)
    print("-" * len(header))

    worst = {m: 0.0 for m in sweep_values}
    fails = []
    for ci, case in enumerate(cases):
        for comp in case["components"]:
            row = f"{case['label']}/{comp}"
            print(f"{row:>{label_w}s} ", end="")
            for mtn in sweep_values:
                model = build_model(case, RES_DEFAULT)
                w = 0.0
                for band_name, nu in case["bands"].items():
                    flux = model.flux_density(TIMES, np.full_like(TIMES, nu))
                    tc = extract_components(flux, [comp])
                    if comp in tc and comp in refs[ci][band_name]:
                        e = max_rel_error(tc[comp], refs[ci][band_name][comp])
                        if not np.isnan(e):
                            w = max(w, e)
                worst[mtn] = max(worst[mtn], w)
                if mtn == check_val and w > 0.05:
                    fails.append((row, w))
                print(f"{w:>{col_w-1}.1%} ", end="")
            n_t = len(build_model(case, RES_DEFAULT).details(TIMES[0], TIMES[-1]).theta)
            print(f"| {n_t:>4d}")

    print("-" * len(header))
    print(f"{'WORST ACROSS ALL':>{label_w}s} " + "".join(f"{worst[m]:>{col_w-1}.1%} " for m in sweep_values))
    if fails:
        print(f"\nFAILED at mtn={check_val}: {len(fails)} cases")
        for label, err in fails:
            print(f"  {label}: {err:.1%}")
    else:
        print(f"\nAll pass at mtn={check_val}")
    for mtn in sweep_values:
        if worst[mtn] <= 0.05:
            print(f"Smallest passing min_theta_num: {mtn}")
            break


def run_phi_sweep(cases, phi_resols, phi_ref=1.0, theta_resol=0.5, t_resol=10):
    """Sweep phi_resol values."""
    print(f"\n{len(cases)} cases × {len(phi_resols)} phi_resol values\n")
    print("Computing references...")
    refs = compute_refs(cases, (phi_ref, theta_resol, t_resol))

    label_w = 38
    col_w   = 10
    header  = f"\n{'':>{label_w}s} " + "".join(f"{'φ='+f'{p:.2f}':>{col_w}s}" for p in phi_resols) + f" | {'n_φ':>4s}"
    print(header)
    print("-" * len(header))

    worst = {p: 0.0 for p in phi_resols}
    for ci, case in enumerate(cases):
        for comp in case["components"]:
            row = f"{case['label']}/{comp}"
            print(f"{row:>{label_w}s} ", end="")
            for pr in phi_resols:
                model = build_model(case, (pr, theta_resol, t_resol))
                w = 0.0
                for band_name, nu in case["bands"].items():
                    flux = model.flux_density(TIMES, np.full_like(TIMES, nu))
                    tc = extract_components(flux, [comp])
                    if comp in tc and comp in refs[ci][band_name]:
                        e = max_rel_error(tc[comp], refs[ci][band_name][comp])
                        if not np.isnan(e):
                            w = max(w, e)
                worst[pr] = max(worst[pr], w)
                print(f"{w:>{col_w-1}.1%} ", end="")
            n_phi = len(build_model(case, (0.15, theta_resol, t_resol)).details(TIMES[0], TIMES[-1]).phi)
            print(f"| {n_phi:>4d}")

    print("-" * len(header))
    print(f"{'WORST ACROSS ALL':>{label_w}s} " + "".join(f"{worst[p]:>{col_w-1}.1%} " for p in phi_resols))
    for pr in phi_resols:
        if worst[pr] <= 0.05:
            print(f"Smallest phi_resol with <5%: {pr}")
            break


def _run_model_grid(jt, of, rvs, res, **extra):
    theta_obs = of * THETA_C
    rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
    m = Model(make_jet(jt), ISM(n_ism=1.0),
              Observer(lumi_dist=LUMI_DIST, z=Z, theta_obs=theta_obs),
              rad, rvs_rad=rad if rvs else None, resolutions=res, **extra)
    f = m.flux_density_grid(TIMES, NUS_GRID)
    return np.asarray(f.fwd.sync), (np.asarray(f.rvs.sync) if rvs else None)


def run_resolution_sweep():
    """Sweep (phi, theta, t) resolutions, output CSV + summary."""
    PHI_RESOLS   = [0.08, 0.1, 0.12, 0.15, 0.2]
    THETA_RESOLS = [0.3, 0.4, 0.5, 0.6, 0.8]
    T_RESOLS     = [6, 8, 10, 12, 15]
    REF_RES      = (0.3, 1.5, 25)
    jet_types    = ["tophat", "gaussian", "powerlaw"]
    obs_factors  = [0, 1, 2, 3]
    cases        = list(itertools.product(jet_types, obs_factors, [False, True]))

    print(f"Computing {len(cases)} references...")
    refs = {}
    for i, (jt, of, rvs) in enumerate(cases):
        t0 = time.time()
        refs[(jt, of, rvs)] = _run_model_grid(jt, of, rvs, REF_RES)
        print(f"  [{i+1}/{len(cases)}] {jt}/{'RS' if rvs else 'FS'}/{of}tc ({time.time()-t0:.1f}s)")

    combos = list(itertools.product(PHI_RESOLS, THETA_RESOLS, T_RESOLS))
    print(f"\nSweeping {len(combos)} combos × {len(cases)} cases = {len(combos)*len(cases)} runs...")

    results = []
    for ci, (pr, tr, ttr) in enumerate(combos):
        t0 = time.time()
        for jt, of, rvs in cases:
            fwd_test, rvs_test = _run_model_grid(jt, of, rvs, (pr, tr, ttr))
            fwd_ref, rvs_ref = refs[(jt, of, rvs)]
            fwd_err = max(simple_max_rel_error(fwd_test[i], fwd_ref[i]) for i in range(3))
            rvs_err = max(simple_max_rel_error(rvs_test[i], rvs_ref[i]) for i in range(3)) if rvs else None
            worst = max(fwd_err, rvs_err) if rvs else fwd_err
            results.append({"phi_r": pr, "theta_r": tr, "t_r": ttr, "jet": jt, "obs": of,
                            "shock": "RS" if rvs else "FS", "fwd_err": fwd_err,
                            "rvs_err": rvs_err, "worst": worst})
        if (ci + 1) % 10 == 0:
            print(f"  [{ci+1}/{len(combos)}] ({pr},{tr},{ttr}) {time.time()-t0:.1f}s")

    csv_path = "script/developer/bench_resolutions_results.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader()
        w.writerows(results)
    print(f"\nSaved {csv_path}")

    for shock in ["FS", "RS"]:
        print(f"\n--- {shock} (worst across jets/angles) ---")
        print(f"{'phi':>5} {'theta':>6} {'t':>4} | {'worst':>7}")
        print("-" * 35)
        groups = defaultdict(float)
        for r in results:
            if r["shock"] == shock:
                key = (r["phi_r"], r["theta_r"], r["t_r"])
                groups[key] = max(groups[key], r["worst"])
        for key in sorted(groups, key=groups.get)[:20]:
            print(f"{key[0]:>5.2f} {key[1]:>6.2f} {key[2]:>4} | {groups[key]:>6.1%}")


def run_combined_sweep():
    """Full 5D sweep: (min_theta, fwd_ratio, phi, theta, t). Pareto frontier of error vs cost."""
    MIN_THETA_NUMS = [24, 32, 40, 48]
    FWD_RATIOS     = [0.3, 0.4, 0.5, 0.6]
    PHI_RESOLS     = [0.1, 0.15, 0.2]
    THETA_RESOLS   = [0.3, 0.5, 0.7, 1.0]
    T_RESOLS       = [6, 10, 15]
    REF_RES        = (0.3, 1.5, 25)
    jet_types      = ["tophat", "gaussian", "powerlaw"]
    obs_factors    = [0, 1, 2, 3]
    cases          = list(itertools.product(jet_types, obs_factors, [False, True]))

    print(f"Computing {len(cases)} references...")
    refs = {}
    for i, (jt, of, rvs) in enumerate(cases):
        t0 = time.time()
        refs[(jt, of, rvs)] = _run_model_grid(jt, of, rvs, REF_RES, min_theta_num=80, fwd_ratio=0.4)
        print(f"  [{i+1}/{len(cases)}] {jt}/{'RS' if rvs else 'FS'}/{of}tc ({time.time()-t0:.1f}s)")

    combos = list(itertools.product(MIN_THETA_NUMS, FWD_RATIOS, PHI_RESOLS, THETA_RESOLS, T_RESOLS))
    print(f"\n{len(combos)} combos × {len(cases)} cases = {len(combos)*len(cases)} runs...")

    results = []
    for ci, (mtn, fr, pr, tr, ttr) in enumerate(combos):
        t0 = time.time()
        for jt, of, rvs in cases:
            fwd_test, rvs_test = _run_model_grid(jt, of, rvs, (pr, tr, ttr), min_theta_num=mtn, fwd_ratio=fr)
            fwd_ref, rvs_ref = refs[(jt, of, rvs)]
            fwd_err = max(simple_max_rel_error(fwd_test[i], fwd_ref[i]) for i in range(3))
            rvs_err = max(simple_max_rel_error(rvs_test[i], rvs_ref[i]) for i in range(3)) if rvs else None
            worst = max(fwd_err, rvs_err) if rvs else fwd_err
            results.append({"mtn": mtn, "fr": fr, "phi_r": pr, "theta_r": tr, "t_r": ttr,
                            "jet": jt, "obs": of, "shock": "RS" if rvs else "FS",
                            "fwd_err": fwd_err, "rvs_err": rvs_err, "worst": worst})
        if (ci + 1) % 20 == 0 or ci == 0:
            print(f"  [{ci+1}/{len(combos)}] ({mtn},{fr},{pr},{tr},{ttr}) {time.time()-t0:.1f}s")

    csv_path = "script/developer/bench_combined_results.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader()
        w.writerows(results)
    print(f"\nSaved {csv_path}")

    agg = defaultdict(lambda: {"fs": 0.0, "rs": 0.0})
    for r in results:
        key = (r["mtn"], r["fr"], r["phi_r"], r["theta_r"], r["t_r"])
        k = "rs" if r["shock"] == "RS" else "fs"
        agg[key][k] = max(agg[key][k], r["worst"])

    def cost(key):
        mtn, fr, pr, tr, ttr = key
        return 360 * pr * (mtn + 90 * tr) * ttr * 5

    ranked = sorted(agg.items(), key=lambda x: max(x[1]["fs"], x[1]["rs"]))
    print(f"\nTOP 30 by max(FS,RS) error:")
    print(f"{'mtn':>5} {'fr':>5} {'phi':>5} {'tht':>5} {'t':>4} | {'FS':>7} {'RS':>7} | {'cost':>8}")
    print("-" * 65)
    for key, v in ranked[:30]:
        mtn, fr, pr, tr, ttr = key
        print(f"{mtn:>5} {fr:>5.2f} {pr:>5.2f} {tr:>5.2f} {ttr:>4} | {v['fs']:>6.1%} {v['rs']:>6.1%} | {cost(key):>8.0f}")

    items = sorted([(key, max(v["fs"], v["rs"]), cost(key)) for key, v in agg.items()], key=lambda x: x[2])
    pareto, best = [], float("inf")
    for key, err, c in items:
        if err < best:
            best = err
            pareto.append((key, err, c))
    print(f"\nPARETO FRONTIER:")
    print(f"{'mtn':>5} {'fr':>5} {'phi':>5} {'tht':>5} {'t':>4} | {'err':>7} | {'cost':>8}")
    for key, err, c in pareto:
        mtn, fr, pr, tr, ttr = key
        print(f"{mtn:>5} {fr:>5.2f} {pr:>5.2f} {tr:>5.2f} {ttr:>4} | {err:>6.1%} | {c:>8.0f}")


# ── Entry point ────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Grid convergence benchmark")
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--theta",       action="store_true", help="Standard theta grid sweep")
    g.add_argument("--extended",    action="store_true", help="Extended parameter space")
    g.add_argument("--phi",         action="store_true", help="Phi grid sweep")
    g.add_argument("--resolutions", action="store_true", help="Resolution parameter sweep")
    g.add_argument("--combined",    action="store_true", help="Full 5D combined sweep")
    args = parser.parse_args()

    if args.theta:
        run_theta_sweep(make_theta_cases(), [8, 12, 16, 20, 24, 28, 32, 40, 48])
    elif args.extended:
        run_theta_sweep(make_extended_cases(), [8, 16, 24, 32, 40, 48])
    elif args.phi:
        run_phi_sweep(make_phi_cases(), [0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20, 0.30])
    elif args.resolutions:
        run_resolution_sweep()
    elif args.combined:
        run_combined_sweep()
    else:
        run_quick()


if __name__ == "__main__":
    main()
