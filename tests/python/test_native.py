"""NativeFunc argument binding + the numba GIL-free path (if numba is installed)."""
import pytest

from VegasAfterglow import NativeFunc


class _StubCfunc:
    """Minimal stand-in for a numba cfunc: .address, ._pyfunc, callable."""

    def __init__(self, pyfunc, address=0):
        self.address = address
        self._pyfunc = pyfunc

    def __call__(self, *args):
        return self._pyfunc(*args)


def _profile(phi, theta, E_iso, theta_c):
    return E_iso * (1.0 - theta / theta_c)


# ----------------- pure-Python binding logic -----------------

def test_partitions_runtime_and_bound_params():
    """Binding trailing kwargs splits the cfunc signature so params holds the bound values in signature order while n_args still counts all four arguments."""
    nf = NativeFunc(_StubCfunc(_profile), E_iso=1e52, theta_c=0.1)
    assert nf.params == [1e52, 0.1]
    assert nf.n_args == 4


def test_call_fallback_applies_bound_params():
    """Calling a NativeFunc through the pure-Python fallback appends the bound params after the runtime args, reproducing the profile value E_iso*(1 - theta/theta_c)."""
    nf = NativeFunc(_StubCfunc(_profile), E_iso=2.0, theta_c=0.5)
    assert nf(0.0, 0.25) == pytest.approx(2.0 * (1.0 - 0.5))


def test_unknown_kwarg_raises_TypeError():
    """Binding a kwarg that is not in the cfunc signature raises a TypeError naming it an unexpected keyword."""
    with pytest.raises(TypeError, match="unexpected keyword"):
        NativeFunc(_StubCfunc(_profile), not_a_param=1.0)


def test_runtime_arg_after_bound_param_raises():
    """A signature where an unbound runtime argument follows a bound parameter is rejected with a ValueError saying it cannot appear after the bound one."""
    def bad_order(phi, E_iso, theta):  # runtime 'theta' after bound 'E_iso'
        return E_iso

    with pytest.raises(ValueError, match="cannot appear after"):
        NativeFunc(_StubCfunc(bad_order), E_iso=1e52)


# ----------------- numba integration (skipped when numba is absent) -----------------

def test_gil_free_ejecta_flux():
    """gil_free-compiled top-hat E_iso and Gamma0 profiles fed to Ejecta yield a flux_density light curve at 1e14 Hz that is finite and strictly positive everywhere."""
    pytest.importorskip("numba")
    import numpy as np

    from VegasAfterglow import Ejecta, ISM, Model, Observer, Radiation, gil_free

    @gil_free
    def top_hat_energy(phi, theta, E_iso, theta_c):
        return E_iso if theta < theta_c else 0.0

    @gil_free
    def top_hat_gamma(phi, theta, Gamma0, theta_c):
        return Gamma0 if theta < theta_c else 1.0

    jet = Ejecta(
        E_iso=top_hat_energy(E_iso=1e52, theta_c=0.1),
        Gamma0=top_hat_gamma(Gamma0=300.0, theta_c=0.1),
    )
    model = Model(
        jet=jet,
        medium=ISM(n_ism=1),
        observer=Observer(lumi_dist=1e28, z=0.1, theta_obs=0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=0.01, p=2.3),
    )
    t = np.logspace(2, 5, 8)
    nu = np.full_like(t, 1e14)
    flux = model.flux_density(t, nu)
    assert np.all(np.isfinite(flux.total))
    assert np.all(flux.total > 0)
