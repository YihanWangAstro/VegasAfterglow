"""End-to-end closure relations: temporal and spectral power-law indices of the
computed light curves against standard afterglow theory (Sari, Piran & Narayan
1998; Granot & Sari 2002).

Tolerances are calibrated against the code's actual output (smooth spectral
shapes curve the local slope by ~0.04-0.10 relative to the asymptotic theory
values); a violation beyond them means the physics changed, not the numerics.
"""
import numpy as np
import pytest

from VegasAfterglow import ISM, Magnetar, Model, Observer, Radiation, TophatJet, Wind

pytestmark = pytest.mark.physics


def slope(x, y):
    return np.polyfit(np.log10(x), np.log10(np.asarray(y)), 1)[0]


def make(medium, p=2.5, eps_B=1e-3, theta_c=0.3, jet_kw=None, **rad_kw):
    jet_args = {"theta_c": theta_c, "E_iso": 1e53, "Gamma0": 300, **(jet_kw or {})}
    return Model(
        jet=TophatJet(**jet_args),
        medium=medium,
        observer=Observer(lumi_dist=3e28, z=0.5, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=eps_B, p=p, **rad_kw),
    )


T_MID = np.logspace(4.0, 5.5, 24)
NU_OPT = np.full_like(T_MID, 1e15)


# ---------------- temporal indices (alpha: F ~ t^-alpha) ----------------

def test_ism_temporal_index_mid_band():
    """ISM light curve at nu_m < nu < nu_c decays with the Granot & Sari slope alpha = 3(p-1)/4 = 1.125 within a calibrated 0.08."""
    # nu_m < nu < nu_c, ISM: alpha = 3(p-1)/4 = 1.125 (measured 1.168)
    alpha = -slope(T_MID, make(ISM(n_ism=1.0)).flux_density(T_MID, NU_OPT).total)
    assert abs(alpha - 1.125) < 0.08


def test_ism_temporal_index_scales_with_p():
    """Changing the electron index to p=2.2 shifts the ISM mid-band decay to alpha = 3(p-1)/4 = 0.9, so the closure relation tracks p."""
    # closure relation alpha(p): p=2.2 -> 3(p-1)/4 = 0.9
    alpha = -slope(T_MID, make(ISM(n_ism=1.0), p=2.2).flux_density(T_MID, NU_OPT).total)
    assert abs(alpha - 0.9) < 0.1


def test_ism_temporal_index_above_cooling():
    """ISM light curve above nu_c decays with alpha = (3p-2)/4 = 1.375 within a calibrated 0.15 that absorbs smooth-spectrum curvature."""
    # nu > nu_c: alpha = (3p-2)/4 = 1.375 (measured 1.475 -- smooth-spectrum curvature)
    m = make(ISM(n_ism=1.0), eps_B=0.1)
    alpha = -slope(T_MID, m.flux_density(T_MID, np.full_like(T_MID, 1e19)).total)
    assert abs(alpha - 1.375) < 0.15


def test_wind_temporal_index_mid_band():
    """Stellar-wind-medium light curve at nu_m < nu < nu_c decays with alpha = (3p-1)/4 = 1.625 within a calibrated 0.08."""
    # Wind: alpha = (3p-1)/4 = 1.625 (measured 1.620)
    alpha = -slope(T_MID, make(Wind(A_star=0.1)).flux_density(T_MID, NU_OPT).total)
    assert abs(alpha - 1.625) < 0.08


# ---------------- spectral indices (beta: F ~ nu^-beta) ----------------

def test_spectral_index_mid_band():
    """Spectral slope between nu_m and nu_c matches beta = (p-1)/2 = 0.75 within a calibrated 0.08."""
    # nu_m < nu < nu_c: beta = (p-1)/2 = 0.75 (measured 0.792)
    m = make(ISM(n_ism=1.0))
    nus = np.logspace(14, 16, 10)
    beta = -slope(nus, np.asarray(m.flux_density_grid(np.array([3e4]), nus).total)[:, 0])
    assert abs(beta - 0.75) < 0.08


def test_spectral_index_above_cooling():
    """Spectral slope above nu_c matches beta = p/2 = 1.25 within 0.1."""
    # nu > nu_c: beta = p/2 = 1.25 (measured 1.249)
    m = make(ISM(n_ism=1.0), eps_B=0.1)
    nus = np.logspace(18, 20, 8)
    beta = -slope(nus, np.asarray(m.flux_density_grid(np.array([3e4]), nus).total)[:, 0])
    assert abs(beta - 1.25) < 0.1


def test_spectrum_rises_below_peak():
    """The maximum local spectral slope in the nu_a-to-past-nu_m window must lie
    between +0.15 and +0.45 (rising segment near the nu^{1/3} asymptote; measured
    ~+0.34), and the local slope at the low-frequency end of the window must
    exceed that at the high-frequency end as the rise flattens across nu_m."""
    m = make(ISM(n_ism=1.0))
    nus = np.logspace(10.4, 13.5, 30)  # window between nu_a and past nu_m
    F = np.asarray(m.flux_density_grid(np.array([1e5]), nus).total)[:, 0]
    local = np.gradient(np.log10(F), np.log10(nus))
    assert 0.15 < np.max(local) < 0.45  # rising segment near nu^{1/3}
    assert local[0] > local[-1]         # turns over across nu_m


# ---------------- light-curve morphology ----------------

def test_jet_break_steepens_light_curve():
    """A narrow jet's post-jet-break temporal decay index exceeds the pre-break index by more than 0.7."""
    # narrow jet: post-break index exceeds pre-break by >~ 3/4 (measured +1.26)
    m = make(ISM(n_ism=1.0), theta_c=0.05)
    t_pre = np.logspace(3.3, 3.8, 10)
    t_post = np.logspace(6.0, 6.7, 10)
    a_pre = -slope(t_pre, m.flux_density(t_pre, np.full(10, 1e15)).total)
    a_post = -slope(t_post, m.flux_density(t_post, np.full(10, 1e15)).total)
    assert a_post - a_pre > 0.7


def test_magnetar_injection_flattens_decay():
    """Magnetar spin-down energy injection flattens the optical decay index by at least 0.1 relative to the no-injection model."""
    # measured: alpha 1.094 (no injection) -> 0.904 (L0=1e48, t0=1e4)
    t = np.logspace(3.5, 4.5, 12)
    nu = np.full_like(t, 1e15)
    m0 = make(ISM(n_ism=1.0), jet_kw={"E_iso": 1e52})
    mm = make(ISM(n_ism=1.0),
              jet_kw={"E_iso": 1e52, "magnetar": Magnetar(L0=1e48, t0=1e4, q=2)})
    a0 = -slope(t, m0.flux_density(t, nu).total)
    am = -slope(t, mm.flux_density(t, nu).total)
    assert am < a0 - 0.1


def test_thick_shell_reverse_shock_peaks_at_crossing():
    """In the thick-shell regime the reverse-shock synchrotron light curve peaks within a factor of ~3 of the shell-crossing time T*(1+z)."""
    # RS light curve peaks at ~T*(1+z) for thick shells (measured ratio 1.06)
    T_dur, z = 1000.0, 0.5
    m = Model(
        jet=TophatJet(theta_c=0.3, E_iso=1e53, Gamma0=100, duration=T_dur),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=3e28, z=z, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1e-3, p=2.5),
        rvs_rad=Radiation(eps_e=0.1, eps_B=1e-2, p=2.5),
    )
    t = np.logspace(1, 6, 60)
    rvs = np.asarray(m.flux_density(t, np.full_like(t, 1e14)).rvs.sync)
    t_peak = t[int(np.argmax(rvs))]
    assert 0.3 < t_peak / (T_dur * (1 + z)) < 3.0


def test_ssc_fraction_grows_with_eps_e_over_eps_B():
    """Lowering eps_B from 1e-2 to 1e-4 raises the SSC-to-synchrotron flux ratio more than tenfold (Compton Y ~ sqrt(eps_e/eps_B)), with SSC dominating at 1e24 Hz in both cases."""
    # Compton Y ~ sqrt(eps_e/eps_B) in the Y>>1 limit: lowering eps_B must
    # raise the SSC-to-synchrotron ratio
    t = np.logspace(4, 5, 8)
    nu = np.full_like(t, 1e24)  # SSC-dominated band (measured: 1.8e6 vs 220)

    def ssc_ratio(eps_B):
        f = make(ISM(n_ism=1.0), eps_B=eps_B, ssc=True).flux_density(t, nu)
        return np.max(np.asarray(f.fwd.ssc)) / np.max(np.asarray(f.fwd.sync))

    lo, hi = ssc_ratio(1e-2), ssc_ratio(1e-4)
    assert hi > 10 * lo
    assert lo > 1  # SSC dominates synchrotron at 1e24 Hz in both configs


def test_off_axis_dimmer_early_brighter_never():
    """Relativistic beaming makes the off-axis (theta_obs > theta_c) flux strictly lower than on-axis at every sampled early time."""
    # off-axis observer sees less flux at early times (relativistic beaming)
    t = np.logspace(3, 4, 8)
    nu = np.full_like(t, 1e15)
    on = np.asarray(make(ISM(n_ism=1.0), theta_c=0.1).flux_density(t, nu).total)
    m_off = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=3e28, z=0.5, theta_obs=0.4),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1e-3, p=2.5),
    )
    off = np.asarray(m_off.flux_density(t, nu).total)
    assert np.all(off < on)
