"""Save/load implementations for :class:`Fitter`.

The on-disk format is a bilby-native HDF5 / JSON file -- compatible with
``bilby.read_in_result`` for inspection-only access -- with a
``meta_data["vegasafterglow"]`` payload that holds the snapshot needed to
reconstruct the original ``Fitter`` end-to-end.

These are module-level functions; the corresponding ``Fitter.save`` and
``Fitter.load`` methods delegate here so the class itself stays focused on
fitting.
"""

import json

import numpy as np

from ..types import FitResult, ParamDef, Scale
from .config import _BandObs
from .utils import _build_transformer


def save_fitter(fitter, path) -> None:
    """Persist ``fitter`` to a bilby-native HDF5 / JSON file at ``path``.

    See :py:meth:`Fitter.save` for the full docstring.
    """
    if fitter.result is None:
        raise RuntimeError(
            "Fitter.save(...) requires a completed fit; call .fit(...) first."
        )

    r = fitter.result
    latex = list(r.latex_labels) if r.latex_labels else None

    # ParamDef + snapshot go through JSON strings: bilby's HDF5 writer rejects
    # nested list-of-dicts in meta_data, but a flat string field round-trips
    # cleanly.
    defs_json = (
        json.dumps(
            [
                {
                    "name": pd.name,
                    "lower": float(pd.lower),
                    "upper": float(pd.upper),
                    "scale": pd.scale.value,
                    "initial": None if pd.initial is None else float(pd.initial),
                }
                for pd in fitter._param_defs
            ]
        )
        if fitter._param_defs is not None
        else None
    )
    vg_meta = {
        "top_k_params": r.top_k_params,
        "top_k_log_probs": r.top_k_log_probs,
        "latex_labels": latex,
        "samples_shape": list(r.samples.shape),
        "param_defs_json": defs_json,
        "fitter_config_json": json.dumps(_snapshot_config(fitter)),
    }

    if r.bilby_result is not None:
        br = r.bilby_result
        br.meta_data = dict(br.meta_data or {})
        br.meta_data["vegasafterglow"] = vg_meta
        br.save_to_file(filename=str(path))
        return

    # Emcee path: synthesize a minimal bilby Result. The placeholder priors are
    # required by bilby's HDF5 reader but are not consulted on reload.
    import bilby
    import pandas as pd

    flat = r.samples.reshape(-1, r.samples.shape[-1])
    posterior = pd.DataFrame(flat, columns=list(r.labels))
    posterior["log_likelihood"] = r.log_probs.reshape(-1)
    priors = bilby.core.prior.PriorDict(
        {label: bilby.core.prior.Uniform(0, 1, label) for label in r.labels}
    )
    br = bilby.core.result.Result(
        label="vegasafterglow",
        outdir=".",
        sampler="emcee",
        posterior=posterior,
        search_parameter_keys=list(r.labels),
        parameter_labels=latex or list(r.labels),
        priors=priors,
        meta_data={"vegasafterglow": vg_meta},
    )
    br.save_to_file(filename=str(path))


def _snapshot_config(fitter) -> dict:
    """Capture constructor args + added observation data as a JSON-friendly dict.

    Custom (callable) jet / medium / extinction can't round-trip through
    serialisation and are recorded as ``None`` with ``had_custom_X`` flags
    set; the user must pass the same callable to ``Fitter.load(path, X=...)``
    in that case.
    """

    def _arr(x):
        return np.asarray(x).tolist()

    return {
        "version": 1,
        "constructor": {
            "z": float(fitter.z),
            "lumi_dist": float(fitter.lumi_dist),
            "jet": fitter.jet if isinstance(fitter.jet, str) else None,
            "medium": fitter.medium if isinstance(fitter.medium, str) else None,
            "fwd_ssc": bool(fitter.fwd_ssc),
            "rvs_ssc": bool(fitter.rvs_ssc),
            "rvs_shock": bool(fitter.rvs_shock),
            "kn": bool(fitter.kn),
            "cmb_cooling": bool(fitter.cmb_cooling),
            "magnetar": bool(fitter.magnetar),
            "rtol": float(fitter.rtol),
            "resolution": [
                float(fitter.phi_resol),
                float(fitter.theta_resol),
                float(fitter.t_resol),
            ],
            "extinction": fitter.extinction
            if isinstance(fitter.extinction, str)
            else None,
        },
        "had_custom_jet": bool(fitter._custom_jet),
        "had_custom_medium": bool(fitter._custom_medium),
        "had_custom_extinction": bool(fitter._custom_extinction),
        "point_data": [
            {
                "t": _arr(t),
                "nu": _arr(nu),
                "flux": _arr(flux),
                "err": _arr(err),
                "weights": _arr(w),
                "label": lbl,
            }
            for t, nu, flux, err, w, lbl in zip(
                fitter._point_t,
                fitter._point_nu,
                fitter._point_flux,
                fitter._point_err,
                fitter._point_weights,
                fitter._point_labels,
            )
        ],
        "band_data": [
            {
                "nu_min": float(b.nu_min),
                "nu_max": float(b.nu_max),
                "num_points": int(b.num_points),
                "t": _arr(b.t),
                "flux": _arr(b.flux),
                "err": _arr(b.err),
                "weights": _arr(b.weights),
                "name": b.name,
            }
            for b in fitter._band_obs
        ],
    }


def load_fitter(path, *, jet=None, medium=None, extinction=None):
    """Reload a saved fit as a fully-configured ``Fitter``.

    See :py:meth:`Fitter.load` for the full docstring.
    """
    import bilby

    from .fitter import Fitter  # local import to avoid circular dependency

    br = bilby.read_in_result(filename=str(path))
    vg = (br.meta_data or {}).get("vegasafterglow") or {}
    fitter_json = vg.get("fitter_config_json")
    if not fitter_json:
        raise ValueError(
            "Saved file does not include a Fitter snapshot (written by "
            "VegasAfterglow < v2.0.5). For inspection-only access to "
            "samples / corner / summary, use `bilby.read_in_result(path)` "
            "directly."
        )
    cfg = json.loads(fitter_json) if isinstance(fitter_json, str) else fitter_json
    cc = cfg["constructor"]

    def _resolve(name, original_str, override, had_custom):
        if had_custom:
            if override is None:
                raise ValueError(
                    f"Saved fit used a custom {name} (callable); pass it as "
                    f"`Fitter.load(path, {name}=...)`."
                )
            return override
        return override if override is not None else original_str

    fitter = Fitter(
        z=cc["z"],
        lumi_dist=cc["lumi_dist"],
        jet=_resolve("jet", cc["jet"], jet, cfg.get("had_custom_jet", False)),
        medium=_resolve(
            "medium", cc["medium"], medium, cfg.get("had_custom_medium", False)
        ),
        fwd_ssc=cc["fwd_ssc"],
        rvs_ssc=cc["rvs_ssc"],
        rvs_shock=cc["rvs_shock"],
        kn=cc["kn"],
        cmb_cooling=cc["cmb_cooling"],
        magnetar=cc["magnetar"],
        rtol=cc["rtol"],
        resolution=tuple(cc["resolution"]),
        extinction=_resolve(
            "extinction",
            cc["extinction"],
            extinction,
            cfg.get("had_custom_extinction", False),
        ),
    )

    for pt in cfg.get("point_data", []):
        fitter._point_t.append(np.asarray(pt["t"], dtype=np.float64))
        fitter._point_nu.append(np.asarray(pt["nu"], dtype=np.float64))
        fitter._point_flux.append(np.asarray(pt["flux"], dtype=np.float64))
        fitter._point_err.append(np.asarray(pt["err"], dtype=np.float64))
        fitter._point_weights.append(np.asarray(pt["weights"], dtype=np.float64))
        fitter._point_labels.append(pt.get("label"))
    for bd in cfg.get("band_data", []):
        fitter._band_obs.append(
            _BandObs(
                nu_min=bd["nu_min"],
                nu_max=bd["nu_max"],
                num_points=bd["num_points"],
                t=np.asarray(bd["t"], dtype=np.float64),
                flux=np.asarray(bd["flux"], dtype=np.float64),
                err=np.asarray(bd["err"], dtype=np.float64),
                weights=np.asarray(bd["weights"], dtype=np.float64),
                name=bd.get("name"),
            )
        )

    defs_json = vg.get("param_defs_json")
    if defs_json:
        defs_meta = json.loads(defs_json) if isinstance(defs_json, str) else defs_json
        param_defs = [
            ParamDef(
                name=str(d["name"]),
                lower=float(d["lower"]),
                upper=float(d["upper"]),
                scale=Scale(str(d["scale"])),
                initial=(None if d.get("initial") is None else float(d["initial"])),
            )
            for d in defs_meta
        ]
        fitter._param_defs = param_defs
        fitter._to_params = _build_transformer(param_defs)

    # Reconstruct the FitResult from the posterior + cached top-K metadata.
    labels = tuple(br.search_parameter_keys)
    ndim = len(labels)
    flat_samples = br.posterior[list(labels)].values
    flat_log_probs = br.posterior["log_likelihood"].values
    # HDF5 round-tripping can return ``samples_shape`` as a numpy array; coerce
    # to a tuple of ints for a clean reshape() call.
    shape = vg.get("samples_shape")
    if shape is not None and len(shape) == 3 and int(shape[-1]) == ndim:
        shape_tuple = tuple(int(s) for s in shape)
        samples = flat_samples.reshape(shape_tuple)
        log_probs = flat_log_probs.reshape(shape_tuple[0], shape_tuple[1])
    else:
        samples = flat_samples.reshape(-1, 1, ndim)
        log_probs = flat_log_probs.reshape(-1, 1)
    latex = vg.get("latex_labels") or br.parameter_labels

    fitter.result = FitResult(
        samples=samples,
        log_probs=log_probs,
        labels=labels,
        latex_labels=tuple(latex) if latex else None,
        top_k_params=vg.get("top_k_params"),
        top_k_log_probs=vg.get("top_k_log_probs"),
        bilby_result=br,
    )
    return fitter
