"""Parameter validation and sampler-input construction called from :meth:`Fitter.fit`.

* ``validate_parameters`` — cross-check ``ParamDef`` names against the chosen
  jet/medium/toggles.
* ``build_sampler_params`` — sampler-space labels + bounds + ``ndim``.
* ``build_prior_dict`` — bilby ``PriorDict`` (user priors override Uniform).
* ``generate_initial_positions`` — walker seed positions for emcee.

Factory-introspection helpers (``_factory_attrs_ast`` and
``collect_factory_accesses``) walk the user's custom jet/medium/extinction
factories to catch typos and dead parameters before the sampler starts.
"""

import ast
import inspect
import textwrap
from typing import List, Optional, Sequence, Tuple

import bilby
import numpy as np

from ..types import ModelParams, ParamDef, Scale
from .config import JET_RULES, MEDIUM_RULES, TOGGLE_RULES
from .utils import _get_latex_label, _param_label


def validate_parameters(fitter, param_defs: Sequence[ParamDef]) -> None:
    """Cross-check ``param_defs`` against the configured jet/medium/toggles."""
    # Skip validation when using custom factories (user manages their own params).
    if fitter._custom_jet or fitter._custom_medium:
        return

    param_names = {pd.name for pd in param_defs}
    missing, incompatible = [], []

    def add_violations(required: set, forbidden: set, context: str):
        missing.extend(f"{p} ({context})" for p in required - param_names)
        incompatible.extend(
            f"{p} (not used with {context})" for p in forbidden & param_names
        )

    for rules, config_attr, label in [
        (MEDIUM_RULES, fitter.medium, "medium"),
        (JET_RULES, fitter.jet, "jet"),
    ]:
        if config_attr in rules:
            required, forbidden = rules[config_attr]
            add_violations(required, forbidden, f"{config_attr} {label}")

    for toggle, (required_on, forbidden_off) in TOGGLE_RULES.items():
        enabled = toggle == "forward_shock" or getattr(fitter, toggle, False)
        if enabled:
            add_violations(required_on, set(), toggle.replace("_", " "))
        else:
            incompatible.extend(
                f"{p} ({toggle} disabled)" for p in forbidden_off & param_names
            )

    if missing or incompatible:
        msg = "Parameter validation failed:\n"
        if missing:
            msg += "Missing:\n  - " + "\n  - ".join(missing) + "\n"
        if incompatible:
            msg += "Incompatible:\n  - " + "\n  - ".join(incompatible) + "\n"
        msg += f"\nConfig: medium='{fitter.medium}', jet='{fitter.jet}', "
        msg += f"rvs_shock={fitter.rvs_shock}, magnetar={fitter.magnetar}"
        raise ValueError(msg)


def _factory_attrs_ast(factory, arg_index: int = 0):
    """Parse a factory's source and return every ``<arg>.X`` reference.

    Walks the full AST including nested function definitions and lambdas, so
    factories that return closures capturing ``params`` (the common pattern
    for ``Ejecta``-based jets) are handled correctly. ``arg_index`` selects
    which positional arg to inspect: 0 for jet/medium ``factory(params)``
    factories, 1 for extinction ``f(lam_cm, params)`` callables.

    Returns ``None`` when the source cannot be retrieved or parsed (e.g.,
    REPL-defined lambdas, bound C-extension callables); the caller should
    treat this as "introspection failed" and fall back to lenient validation
    rather than blocking the user.
    """
    try:
        src = inspect.getsource(factory)
    except (OSError, TypeError):
        return None
    src = textwrap.dedent(src)
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return None
    outermost = next(
        (
            n
            for n in ast.walk(tree)
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda))
        ),
        None,
    )
    if outermost is None or len(outermost.args.args) <= arg_index:
        return None
    arg_name = outermost.args.args[arg_index].arg
    return {
        node.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Name)
        and node.value.id == arg_name
    }


def collect_factory_accesses(fitter):
    """Discover which ``params`` attributes the custom factories read.

    Returns ``(accessed, ast_complete)`` where ``accessed`` is the union of
    names found across all configured custom factories and ``ast_complete``
    is False if any factory's source could not be parsed (in which case the
    caller should be lenient about ParamDefs that don't appear to match
    anything).

    Returns ``(set(), True)`` when no custom factories are configured.
    """
    if not (fitter._custom_jet or fitter._custom_medium or fitter._custom_extinction):
        return set(), True

    accessed: set = set()
    ast_complete = True

    for active, factory, arg_idx in (
        (fitter._custom_jet, fitter.jet_factory, 0),
        (fitter._custom_medium, fitter.medium_factory, 0),
        (fitter._custom_extinction, fitter._ext_law, 1),
    ):
        if not active:
            continue
        names = _factory_attrs_ast(factory, arg_idx)
        if names is None:
            ast_complete = False
        else:
            accessed |= names

    return accessed, ast_complete


def build_sampler_params(
    fitter, param_defs: List[ParamDef]
) -> Tuple[Tuple[str, ...], np.ndarray, np.ndarray, int]:
    """Build parameter labels and bounds for the sampler.

    Validates ``ParamDef`` names against the standard ``ModelParams`` fields
    and, when custom factories are configured, against the set of attributes
    those factories read at dry-run time. Catches:

    * ``ParamDef('typo')`` not used by any factory or standard model.
    * Custom factory reading ``params.X`` with no matching ``ParamDef`` for
      ``X`` (which would silently fail at fit time with -inf likelihood).
    """
    standard_fields = set(vars(ModelParams()).keys())
    custom_active = (
        fitter._custom_jet or fitter._custom_medium or fitter._custom_extinction
    )

    if custom_active:
        accessed, ast_complete = collect_factory_accesses(fitter)
        paramdef_names = {pd.name for pd in param_defs}

        for pd in param_defs:
            if pd.name in standard_fields or pd.name in accessed:
                continue
            if not ast_complete:
                continue
            raise AttributeError(
                f"ParamDef('{pd.name}') is neither a standard ModelParams "
                f"field nor read by any custom factory; likely a typo or "
                f"dead parameter. Custom factories read: {sorted(accessed)}"
            )

        if ast_complete:
            missing = accessed - standard_fields - paramdef_names
            if missing:
                raise AttributeError(
                    f"Custom factory reads {sorted(missing)} from `params` "
                    f"but no matching ParamDef is declared and these are "
                    f"not on ModelParams; the fit would silently return "
                    f"-inf likelihood. Add them as ParamDefs (Scale.fixed "
                    f"if not fitted), or use ``getattr(params, X, default)`` "
                    f"if the dynamic access is intentional."
                )
    else:
        for pd in param_defs:
            if pd.name not in standard_fields:
                raise AttributeError(f"'{pd.name}' is not a valid MCMC parameter")

    labels, lowers, uppers = zip(
        *(
            (
                _param_label(pd),
                np.log10(pd.lower) if pd.scale is Scale.log else pd.lower,
                np.log10(pd.upper) if pd.scale is Scale.log else pd.upper,
            )
            for pd in param_defs
            if pd.scale is not Scale.fixed
        )
    )
    return labels, np.array(lowers), np.array(uppers), len(labels)


def build_prior_dict(
    labels: Tuple[str, ...],
    pl: np.ndarray,
    pu: np.ndarray,
    defs: List[ParamDef],
    user_priors: Optional[dict],
):
    """Build bilby PriorDict: user priors where provided, Uniform for the rest."""
    label_to_def = {_param_label(pd): pd for pd in defs if pd.scale is not Scale.fixed}

    priors_dict = {}
    for i, name in enumerate(labels):
        if user_priors and name in user_priors:
            priors_dict[name] = user_priors[name]
        else:
            priors_dict[name] = bilby.core.prior.Uniform(
                pl[i], pu[i], name, _get_latex_label(label_to_def[name])
            )
    return bilby.core.prior.PriorDict(priors_dict)


def generate_initial_positions(
    param_defs: List[ParamDef],
    lower_bounds: np.ndarray,
    upper_bounds: np.ndarray,
    nwalkers: int,
    ndim: int,
) -> np.ndarray:
    """Generate initial walker positions centered around parameter initial values."""
    initial_vals = []
    idx = 0
    for pd in param_defs:
        if pd.scale is Scale.fixed:
            continue
        if pd.initial is not None:
            val = np.log10(pd.initial) if pd.scale is Scale.log else pd.initial
        else:
            val = 0.5 * (lower_bounds[idx] + upper_bounds[idx])
        initial_vals.append(val)
        idx += 1

    initial_vals = np.array(initial_vals)
    spread = 0.1 * (upper_bounds - lower_bounds)
    pos0 = initial_vals + spread * np.random.randn(nwalkers, ndim)

    eps = 1e-6 * (upper_bounds - lower_bounds)
    pos0 = np.clip(pos0, lower_bounds + eps, upper_bounds - eps)
    return pos0
