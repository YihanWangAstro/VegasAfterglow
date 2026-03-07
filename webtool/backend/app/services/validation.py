from __future__ import annotations

from fastapi import HTTPException

from ..constants import FLUX_SCALES, INSTRUMENTS, TIME_SCALES
from ..schemas import SharedParams


def validate_shared(shared: SharedParams) -> None:
    if shared.flux_unit not in list(FLUX_SCALES.keys()) + ["AB mag"]:
        raise HTTPException(status_code=400, detail=f"Unsupported flux_unit: {shared.flux_unit}")
    if shared.time_unit not in TIME_SCALES:
        raise HTTPException(status_code=400, detail=f"Unsupported time_unit: {shared.time_unit}")


def validate_instruments(selected_instruments: list[str]) -> list[str]:
    unknown = [name for name in selected_instruments if name not in INSTRUMENTS]
    if unknown:
        raise HTTPException(status_code=400, detail=f"Unknown instruments: {', '.join(unknown)}")
    return selected_instruments


def resolve_export_kinds(requested: list[str], defaults: set[str], allowed: set[str]) -> set[str]:
    if not requested:
        return defaults
    normalized = {kind.strip().lower() for kind in requested if kind.strip()}
    unknown = normalized - allowed
    if unknown:
        raise HTTPException(status_code=400, detail=f"Unknown export kind(s): {', '.join(sorted(unknown))}")
    return normalized
