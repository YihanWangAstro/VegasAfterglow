from __future__ import annotations

import re
from typing import Any

from ..helpers import parse_entry, z_from_lumi_dist_mpc
from ..schemas import ObservationGroup, SharedParams

_FREQ_SPLIT = re.compile(r",(?![^\[]*\])")
_OBS_SPLIT = re.compile(r"[,\s\t]+")


def shared_to_physics(shared: SharedParams) -> dict[str, Any]:
    # Exclude UI-only fields; add derived d_L_cm and z in place of d_L_mpc.
    physics = shared.model_dump(exclude={"flux_unit", "time_unit", "d_L_mpc", "num_t"})
    physics["d_L_cm"] = shared.d_L_mpc * 3.0856775814913673e24
    physics["z"] = z_from_lumi_dist_mpc(shared.d_L_mpc)
    return physics


def parse_frequency_input(frequencies_input: str) -> tuple[list[float], list[list[Any]], list[str]]:
    frequencies: list[float] = []
    bands: list[list[Any]] = []
    warnings: list[str] = []
    tokens = [tok.strip() for tok in _FREQ_SPLIT.split(frequencies_input) if tok.strip()]

    for token in tokens:
        try:
            entry = parse_entry(token)
            if isinstance(entry, tuple):
                bands.append([float(entry[0]), float(entry[1]), str(entry[2])])
            else:
                frequencies.append(float(entry))
        except Exception:
            warnings.append(f"Unknown frequency or filter: '{token}'")

    return sorted(frequencies), bands, warnings


def parse_observation_rows(groups: list[ObservationGroup], is_lc: bool) -> tuple:
    rows = []
    x_default = "day" if is_lc else "Hz"
    for group in groups:
        if not group.visible:
            continue
        label = group.legend or "data"
        x_unit = group.x_unit or x_default
        y_unit = group.y_unit
        for line in group.text.strip().splitlines():
            parts = _OBS_SPLIT.split(line.strip())
            if len(parts) < 2:
                continue
            try:
                x_val = float(parts[0])
                y_val = float(parts[1])
            except (ValueError, IndexError):
                continue
            try:
                err_val = float(parts[2]) if len(parts) > 2 else 0.0
            except ValueError:
                err_val = 0.0
            rows.append((label, x_val, x_unit, y_val, err_val, y_unit))
    return tuple(rows)
