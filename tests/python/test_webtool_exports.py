import json

import numpy as np

from webtool.backend.app.exports import export_json, export_sed_json, export_skymap_json


def test_lightcurve_json_export_serializes_numpy_arrays():
    data = {
        "times": np.array([1.0, 10.0]),
        "frequencies": np.array([1e9]),
        "pt_components": [("total", np.array([[1.5, 2.5]]))],
        "band_data": [],
    }

    exported = export_json(data, "cgs", "s")
    parsed = json.loads(exported)

    assert parsed["units"] == {"time": "s", "flux_density": "cgs"}
    assert parsed["times"] == [1.0, 10.0]
    assert parsed["frequencies_Hz"] == [1e9]
    assert parsed["flux_density"]["1 GHz"]["total"] == [1.5, 2.5]


def test_spectrum_json_export_serializes_numpy_arrays():
    data = {
        "t_snapshots": np.array([1e3, 1e4]),
        "frequencies": np.array([1e9, 2e9]),
        "components": [("total", np.array([[1.0, 2.0], [3.0, 4.0]]))],
    }

    exported = export_sed_json(data, "cgs", "Hz")
    parsed = json.loads(exported)

    assert parsed["frequencies"] == [1e9, 2e9]
    assert parsed["t_snapshots_s"] == [1e3, 1e4]
    assert parsed["flux_density"]["16.7 min"]["total"] == [1.0, 3.0]
    assert parsed["flux_density"]["2.78 hr"]["total"] == [2.0, 4.0]


def test_skymap_json_export_serializes_numpy_arrays():
    data = {
        "t_obs_array": np.array([1e4, 1e5]),
        "nu_obs": 1e9,
        "fov_uas": 500.0,
        "extent_uas": [-250.0, 250.0, -250.0, 250.0],
        "images": np.array([[[1.0, 2.0], [3.0, 4.0]]]),
    }

    exported = export_skymap_json(data)
    parsed = json.loads(exported)

    assert parsed["t_obs_s"] == [1e4, 1e5]
    assert parsed["nu_obs_Hz"] == 1e9
    assert parsed["images"] == [[[1.0, 2.0], [3.0, 4.0]]]
