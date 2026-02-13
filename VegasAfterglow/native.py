import inspect


class NativeFunc:
    """Wrap a numba @cfunc with bound parameters for GIL-free C++ execution.

    Runtime arguments (phi, theta for BinaryFunc; phi, theta, r for TernaryFunc)
    are determined automatically: any cfunc parameter NOT provided as a keyword
    argument is treated as a runtime argument. Bound parameters must come after
    all runtime arguments in the cfunc signature.

    Example (BinaryFunc — custom jet)::

        @cfunc(float64(float64, float64, float64, float64))
        def gaussian_energy(phi, theta, E_iso, theta_c):
            return E_iso * math.exp(-(theta / theta_c) ** 2)

        nf = NativeFunc(gaussian_energy, E_iso=1e52, theta_c=0.1)

    Example (TernaryFunc — custom medium)::

        @cfunc(float64(float64, float64, float64, float64))
        def wind_rho(phi, theta, r, A_star):
            return A_star * 5e11 * mp / (r * r)

        nf = NativeFunc(wind_rho, A_star=1.0)
    """

    def __init__(self, cfunc_obj, **kwargs):
        self.address = cfunc_obj.address
        pyfunc = cfunc_obj._pyfunc
        all_names = list(inspect.signature(pyfunc).parameters.keys())

        # Reject kwargs that don't match any parameter name
        unknown = set(kwargs) - set(all_names)
        if unknown:
            raise TypeError(
                f"{pyfunc.__name__}() got unexpected keyword arguments: "
                + ", ".join(sorted(unknown))
            )

        # Partition into runtime args (not in kwargs) and bound params (in kwargs)
        runtime_names = [n for n in all_names if n not in kwargs]
        bound_names = [n for n in all_names if n in kwargs]

        # Validate: bound params must come after all runtime args in the signature
        if bound_names:
            first_bound = all_names.index(bound_names[0])
            for name in runtime_names:
                if all_names.index(name) > first_bound:
                    raise ValueError(
                        f"Runtime arg '{name}' cannot appear after "
                        f"bound param '{bound_names[0]}' in cfunc signature"
                    )

        self.n_args = len(all_names)
        self.params = [float(kwargs[name]) for name in bound_names]
        self._cfunc = cfunc_obj
        self._runtime_names = runtime_names

    def __call__(self, *args):
        """Fallback: call through Python (for testing without C++)."""
        return self._cfunc(*args, *self.params)


def gil_free(fn):
    """Decorator: compile a Python function with numba for GIL-free C++ execution.

    The decorated function becomes a factory — call it with keyword arguments
    to bind parameters and get a NativeFunc object.

    Example (custom jet)::

        @gil_free
        def gaussian_energy(phi, theta, E_iso, theta_c):
            return E_iso * math.exp(-0.5 * (theta / theta_c) ** 2)

        jet = Ejecta(
            E_iso=gaussian_energy(E_iso=1e52, theta_c=0.1),
            Gamma0=...,
        )

    Example (custom medium)::

        @gil_free
        def wind_rho(phi, theta, r, A_star):
            return A_star * 5e11 * 1.67e-24 / (r * r)

        medium = Medium(rho=wind_rho(A_star=1.0))

    Requires numba: ``pip install numba``
    """
    from numba import cfunc, float64

    n_params = len(inspect.signature(fn).parameters)
    compiled = cfunc(float64(*([float64] * n_params)))(fn)

    def factory(**kwargs):
        return NativeFunc(compiled, **kwargs)

    factory.__name__ = fn.__name__
    factory.__qualname__ = fn.__qualname__
    factory.__doc__ = fn.__doc__
    factory._cfunc = compiled
    return factory
