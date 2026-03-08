import {
  useCallback,
  useDeferredValue,
  useEffect,
  useMemo,
  type Dispatch,
  type MutableRefObject,
  type ReactNode,
  type SetStateAction,
} from "react";
import { ALL_AXES } from "../lib/constants";
import type {
  AxisName,
  AxisRanges,
  AxisSignatures,
  BookmarkEntry,
  LcPlotData,
  Mode,
  RunResponse,
  SedPlotData,
  SharedParams,
  SkymapPlotData,
} from "../lib/types";
import { formatParamValueNode } from "../lib/utils/math";
import {
  applySkymapAxisPolicy,
  axisSignature,
  legendFontSizeForWidth,
  parseAxisRange,
  parseRelayoutRange,
  rangesEqual,
  shouldPreserveAxisRangesOnAutorange,
} from "../lib/utils/plot";
import { buildLcFigure, buildSedFigure, buildSkymapFigure } from "../lib/utils/plot-builders";

type UseFigurePresentationArgs = {
  mode: Mode;
  result: RunResponse | null;
  resultMode: Mode | null;
  compareEnabled: boolean;
  compareResult: RunResponse | null;
  compareResultMode: Mode | null;
  selectedCompareBookmark: BookmarkEntry | null;
  plotWidthPx: number;
  zoomRevision: number;
  shared: SharedParams;
  axisRangesRef: MutableRefObject<Record<Mode, AxisRanges>>;
  axisSignaturesRef: MutableRefObject<Record<Mode, AxisSignatures>>;
  setZoomRevision: Dispatch<SetStateAction<number>>;
};

export function useFigurePresentation({
  mode,
  result,
  resultMode,
  compareEnabled,
  compareResult,
  compareResultMode,
  selectedCompareBookmark,
  plotWidthPx,
  zoomRevision,
  shared,
  axisRangesRef,
  axisSignaturesRef,
  setZoomRevision,
}: UseFigurePresentationArgs) {
  const warnings = resultMode === mode ? result?.meta?.warnings ?? [] : [];

  const handlePlotRelayout = useCallback(
    (event: Record<string, unknown>) => {
      const current = axisRangesRef.current[mode];
      const next: AxisRanges = { ...current };
      let changed = false;
      const preserveOnAutorange = shouldPreserveAxisRangesOnAutorange(mode);

      for (const axis of ALL_AXES) {
        if (event[`${axis}.autorange`] === true && next[axis] !== null) {
          if (preserveOnAutorange) {
            continue;
          }
          next[axis] = null;
          changed = true;
        }
        const range = parseRelayoutRange(event, axis);
        if (!range) continue;
        const prev = next[axis];
        if (!rangesEqual(prev, range)) {
          next[axis] = range;
          changed = true;
        }
      }

      if (changed) {
        axisRangesRef.current[mode] = next;
        setZoomRevision((value) => value + 1);
      }
    },
    [axisRangesRef, mode, setZoomRevision],
  );

  // Build figure client-side from plot_data for all modes.
  const plotData = result?.plot_data;
  const builtFigure = useMemo(() => {
    if (resultMode !== mode) return null;
    if (!plotData) return null;
    if (mode === "lightcurve") return buildLcFigure(plotData as LcPlotData);
    if (mode === "spectrum") return buildSedFigure(plotData as SedPlotData);
    if (mode === "skymap") return buildSkymapFigure(plotData as SkymapPlotData);
    return null;
  }, [plotData, resultMode, mode]);

  const comparePlotData = compareResult?.plot_data;
  const compareBuiltFigure = useMemo(() => {
    if (!compareEnabled || compareResultMode !== mode) return null;
    if (!comparePlotData) return null;
    if (mode === "lightcurve") return buildLcFigure(comparePlotData as LcPlotData);
    if (mode === "spectrum") return buildSedFigure(comparePlotData as SedPlotData);
    return null;
  }, [compareEnabled, comparePlotData, compareResultMode, mode]);

  useEffect(() => {
    if (resultMode !== mode) return;
    const layout = builtFigure?.layout;
    if (!layout) return;

    const current = axisRangesRef.current[mode];
    const next: AxisRanges = { ...current };
    let changed = false;

    for (const axis of ALL_AXES) {
      if (next[axis] !== null) continue;
      const axisObj = layout[axis];
      if (!axisObj || typeof axisObj !== "object") continue;
      const range = parseAxisRange((axisObj as { range?: unknown }).range);
      if (!range) continue;
      next[axis] = range;
      changed = true;
    }

    if (changed) {
      axisRangesRef.current[mode] = next;
      setZoomRevision((value) => value + 1);
    }
  }, [axisRangesRef, builtFigure, mode, resultMode, setZoomRevision]);

  useEffect(() => {
    if (resultMode !== mode) return;
    const layout = builtFigure?.layout;
    if (!layout) return;

    const nextYSig = axisSignature(layout, "yaxis");
    const nextY2Sig = axisSignature(layout, "yaxis2");
    const prevSig = axisSignaturesRef.current[mode];
    const clearY = prevSig.yaxis !== null && prevSig.yaxis !== nextYSig;
    const clearY2 = prevSig.yaxis2 !== null && prevSig.yaxis2 !== nextY2Sig;

    axisSignaturesRef.current[mode] = { yaxis: nextYSig, yaxis2: nextY2Sig };
    if (!clearY && !clearY2) return;

    const current = axisRangesRef.current[mode];
    const next: AxisRanges = { ...current };
    let changed = false;

    if (clearY && next.yaxis !== null) {
      next.yaxis = null;
      changed = true;
    }
    if (clearY2 && next.yaxis2 !== null) {
      next.yaxis2 = null;
      changed = true;
    }

    if (changed) {
      axisRangesRef.current[mode] = next;
      setZoomRevision((value) => value + 1);
    }
  }, [axisRangesRef, axisSignaturesRef, builtFigure, mode, resultMode, setZoomRevision]);

  // baseFigure: compare overlay, legend font, skymap layout, uirevision.
  // Runs only when computation results or compare settings change, NOT on every zoom.
  const baseFigure = useMemo(() => {
    if (resultMode !== mode) return null;
    const figure = builtFigure;
    if (!figure?.data) return null;

    const rawData = Array.isArray(figure.data) ? [...figure.data] : [];
    let mergedData = rawData;

    if (mode !== "skymap" && compareEnabled && compareResultMode === mode && selectedCompareBookmark) {
      const compareData = compareBuiltFigure?.data;
      if (Array.isArray(compareData) && compareData.length > 0) {
        const compareLabel = selectedCompareBookmark.name;
        const overlayData = compareData.map((trace, index) => {
          if (!trace || typeof trace !== "object") return trace;
          const nextTrace = { ...(trace as Record<string, unknown>) };
          const name = typeof nextTrace.name === "string" ? nextTrace.name : `trace ${index + 1}`;
          nextTrace.name = `${name} (${compareLabel})`;
          nextTrace.legendgroup =
            typeof nextTrace.legendgroup === "string" ? `cmp-${nextTrace.legendgroup}` : `cmp-${index + 1}`;
          nextTrace.opacity = typeof nextTrace.opacity === "number" ? Math.min(0.9, nextTrace.opacity) : 0.85;
          if (typeof nextTrace.fill === "string" && nextTrace.fill !== "none") {
            nextTrace.fill = "none";
          }

          const traceType = typeof nextTrace.type === "string" ? nextTrace.type : "";
          if (traceType === "" || traceType === "scatter" || traceType === "scattergl") {
            const lineObj =
              nextTrace.line && typeof nextTrace.line === "object"
                ? { ...(nextTrace.line as Record<string, unknown>) }
                : {};
            lineObj.dash = "dot";
            nextTrace.line = lineObj;
          }
          return nextTrace;
        });
        mergedData = [...rawData, ...overlayData];
      }
    }

    const layout = { ...(figure.layout ?? {}) } as Record<string, unknown>;

    if (mode === "skymap") {
      delete layout.updatemenus;
      applySkymapAxisPolicy(layout);
    }

    layout.uirevision = `${mode}-zoom-state`;
    return { ...figure, data: mergedData, layout };
  }, [builtFigure, compareBuiltFigure, compareEnabled, compareResultMode, mode, resultMode, selectedCompareBookmark]);

  // displayFigure: cheap zoom application — runs on every zoom/pan interaction.
  const displayFigure = useMemo(() => {
    if (!baseFigure) return null;

    const zoom = axisRangesRef.current[mode];
    const axes: AxisName[] = ["xaxis", "yaxis", "yaxis2"];
    const layout = { ...baseFigure.layout } as Record<string, unknown>;

    for (const axis of axes) {
      const saved = zoom[axis];
      if (!saved) continue;
      const axisObj = layout[axis];
      const nextAxis = axisObj && typeof axisObj === "object" ? { ...(axisObj as Record<string, unknown>) } : {};
      nextAxis.range = [saved[0], saved[1]];
      nextAxis.autorange = false;
      layout[axis] = nextAxis;
    }

    if (mode === "skymap") {
      applySkymapAxisPolicy(layout);
    } else {
      const legendFontSize = legendFontSizeForWidth(plotWidthPx);
      const legendObj =
        layout.legend && typeof layout.legend === "object" ? { ...(layout.legend as Record<string, unknown>) } : {};
      const legendFontObj =
        legendObj.font && typeof legendObj.font === "object" ? { ...(legendObj.font as Record<string, unknown>) } : {};
      legendFontObj.size = legendFontSize;
      legendObj.font = legendFontObj;
      layout.legend = legendObj;
    }

    return { ...baseFigure, layout };
  }, [axisRangesRef, baseFigure, mode, plotWidthPx, zoomRevision]);

  const deferredFigure = useDeferredValue(displayFigure);

  const figureCaption = useMemo(() => {
    const title =
      mode === "lightcurve"
        ? "Multi-band GRB afterglow light curves"
        : mode === "spectrum"
          ? "GRB afterglow spectra"
          : "GRB afterglow sky images";
    const mediumLabel =
      shared.medium_type === "ISM"
        ? "uniform ISM"
        : shared.medium_type === "Wind"
          ? "stellar-wind medium"
          : "wind-bubble medium";
    const mediumParam =
      shared.medium_type === "ISM" ? (
        <>
          n<sub>ism</sub>={formatParamValueNode(shared.n_ism)} cm<sup>-3</sup>
        </>
      ) : shared.medium_type === "Wind" ? (
        <>
          A<sub>*</sub>={formatParamValueNode(shared.A_star)}
        </>
      ) : (
        <>
          A<sub>*</sub>={formatParamValueNode(shared.A_star)}, n<sub>floor</sub>={formatParamValueNode(shared.n_ism)} cm<sup>-3</sup>
        </>
      );
    const reverseShockTerms: ReactNode[] = [];
    if (shared.enable_rvs) {
      reverseShockTerms.push(
        <>
          jet duration={formatParamValueNode(shared.duration)} s
        </>,
      );
      reverseShockTerms.push(
        <>
          ε<sub>e,r</sub>={formatParamValueNode(shared.eps_e_r)}
        </>,
      );
      reverseShockTerms.push(
        <>
          ε<sub>B,r</sub>={formatParamValueNode(shared.eps_B_r)}
        </>,
      );
      reverseShockTerms.push(
        <>
          p<sub>r</sub>={formatParamValueNode(shared.p_r)}
        </>,
      );
      reverseShockTerms.push(
        <>
          ξ<sub>e,r</sub>={formatParamValueNode(shared.xi_e_r)}
        </>,
      );
      if (shared.rvs_ssc) reverseShockTerms.push("SSC_r enabled");
      if (shared.rvs_kn) reverseShockTerms.push("KN_r enabled");
    }

    return (
      <>
        {title} generated with VegasAfterglow for a {shared.jet_type.toLowerCase()} jet in a {mediumLabel}, with E
        <sub>iso</sub>={formatParamValueNode(shared.E_iso)} erg, Γ<sub>0</sub>={formatParamValueNode(shared.Gamma0)}, θ
        <sub>c</sub>={formatParamValueNode(shared.theta_c)} rad, {mediumParam}, ε<sub>e</sub>={formatParamValueNode(shared.eps_e)}, ε
        <sub>B</sub>={formatParamValueNode(shared.eps_B)}, p={formatParamValueNode(shared.p)}, ξ<sub>e</sub>=
        {formatParamValueNode(shared.xi_e)}, d<sub>L</sub>={formatParamValueNode(shared.d_L_mpc)} Mpc, and θ<sub>obs</sub>=
        {formatParamValueNode(shared.theta_obs)} rad
        {shared.enable_rvs ? (
          <>
            ; reverse shock with{" "}
            {reverseShockTerms.map((item, idx) => (
              <span key={`rvs-node-${idx}`}>
                {idx > 0 ? ", " : ""}
                {item}
              </span>
            ))}
          </>
        ) : null}
        .
      </>
    );
  }, [mode, shared]);

  return {
    warnings,
    handlePlotRelayout,
    displayFigure,
    deferredFigure,
    figureCaption,
  };
}
