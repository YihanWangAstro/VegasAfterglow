import { useCallback, useRef } from "react";
import type { AxisRanges, AxisSignatures, Mode } from "../lib/types";

export function useAxisPersistence() {
  const axisRangesRef = useRef<Record<Mode, AxisRanges>>({
    lightcurve: { xaxis: null, yaxis: null, yaxis2: null },
    spectrum: { xaxis: null, yaxis: null, yaxis2: null },
    skymap: { xaxis: null, yaxis: null, yaxis2: null },
  });

  const axisSignaturesRef = useRef<Record<Mode, AxisSignatures>>({
    lightcurve: { yaxis: null, yaxis2: null },
    spectrum: { yaxis: null, yaxis2: null },
    skymap: { yaxis: null, yaxis2: null },
  });

  const clearSavedXAxis = useCallback((mode: Mode) => {
    const current = axisRangesRef.current[mode];
    axisRangesRef.current[mode] = { ...current, xaxis: null };
  }, []);

  const clearSavedYAxes = useCallback((mode: Mode) => {
    const current = axisRangesRef.current[mode];
    axisRangesRef.current[mode] = { ...current, yaxis: null, yaxis2: null };
  }, []);

  return {
    axisRangesRef,
    axisSignaturesRef,
    clearSavedXAxis,
    clearSavedYAxes,
  };
}
