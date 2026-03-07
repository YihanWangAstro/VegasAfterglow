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
  Mode,
  RunResponse,
  SharedParams,
} from "../lib/types";
import { formatParamValueNode } from "../lib/utils/math";
import {
  axisSignature,
  legendFontSizeForWidth,
  parseAxisRange,
  parseRelayoutRange,
  rangesEqual,
  remapScientificHoverTemplate,
} from "../lib/utils/plot";

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

      for (const axis of ALL_AXES) {
        if (event[`${axis}.autorange`] === true && next[axis] !== null) {
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

  useEffect(() => {
    if (resultMode !== mode) return;
    const layout = result?.figure?.layout;
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
  }, [axisRangesRef, mode, result, resultMode, setZoomRevision]);

  useEffect(() => {
    if (resultMode !== mode) return;
    const layout = result?.figure?.layout;
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
  }, [axisRangesRef, axisSignaturesRef, mode, result, resultMode, setZoomRevision]);

  // baseFigure: expensive work — data remap, compare overlay, legend, skymap layout.
  // Runs only when computation results or compare settings change, NOT on every zoom.
  const baseFigure = useMemo(() => {
    if (resultMode !== mode) return null;
    const figure = result?.figure;
    if (!figure?.data) return null;

    const rawData = Array.isArray(figure.data) ? [...figure.data] : [];
    const baseData =
      mode === "skymap"
        ? rawData
        : rawData.map((trace) => {
            if (!trace || typeof trace !== "object") return trace;
            return remapScientificHoverTemplate(trace as Record<string, unknown>);
          });
    let mergedData = baseData;

    if (mode !== "skymap" && compareEnabled && compareResultMode === mode && selectedCompareBookmark) {
      const compareData = compareResult?.figure?.data;
      if (Array.isArray(compareData) && compareData.length > 0) {
        const compareLabel = selectedCompareBookmark.name;
        const overlayData = compareData.map((trace, index) => {
          if (!trace || typeof trace !== "object") return trace;
          const sourceTrace = remapScientificHoverTemplate(trace as Record<string, unknown>);
          const nextTrace = { ...sourceTrace };
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
        mergedData = [...baseData, ...overlayData];
      }
    }

    const layout = { ...(figure.layout ?? {}) } as Record<string, unknown>;

    if (mode === "skymap") {
      delete layout.updatemenus;
      const xAxisObj =
        layout.xaxis && typeof layout.xaxis === "object" ? { ...(layout.xaxis as Record<string, unknown>) } : {};
      const yAxisObj =
        layout.yaxis && typeof layout.yaxis === "object" ? { ...(layout.yaxis as Record<string, unknown>) } : {};
      xAxisObj.constrain = "domain";
      yAxisObj.constrain = "domain";
      yAxisObj.scaleanchor = "x";
      yAxisObj.scaleratio = 1;
      layout.xaxis = xAxisObj;
      layout.yaxis = yAxisObj;
    }

    if (mode !== "skymap") {
      const legendFontSize = legendFontSizeForWidth(plotWidthPx);
      const legendObj =
        layout.legend && typeof layout.legend === "object" ? { ...(layout.legend as Record<string, unknown>) } : {};
      const legendFontObj =
        legendObj.font && typeof legendObj.font === "object" ? { ...(legendObj.font as Record<string, unknown>) } : {};
      legendFontObj.size = legendFontSize;
      legendObj.font = legendFontObj;
      layout.legend = legendObj;
    }

    layout.uirevision = `${mode}-zoom-state`;
    return { ...figure, data: mergedData, layout };
  }, [compareEnabled, compareResult, compareResultMode, mode, plotWidthPx, result, resultMode, selectedCompareBookmark]);

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
      const xAxisObj =
        layout.xaxis && typeof layout.xaxis === "object" ? { ...(layout.xaxis as Record<string, unknown>) } : {};
      const yAxisObj =
        layout.yaxis && typeof layout.yaxis === "object" ? { ...(layout.yaxis as Record<string, unknown>) } : {};
      const xRange = parseAxisRange((xAxisObj as { range?: unknown }).range);
      const yRange = parseAxisRange((yAxisObj as { range?: unknown }).range);

      if (xRange && yRange) {
        const xCenter = 0.5 * (xRange[0] + xRange[1]);
        const yCenter = 0.5 * (yRange[0] + yRange[1]);
        const span = Math.max(Math.abs(xRange[1] - xRange[0]), Math.abs(yRange[1] - yRange[0]));
        const half = span * 0.5;
        xAxisObj.range = [xCenter - half, xCenter + half];
        yAxisObj.range = [yCenter - half, yCenter + half];
        xAxisObj.autorange = false;
        yAxisObj.autorange = false;
      }

      layout.xaxis = xAxisObj;
      layout.yaxis = yAxisObj;
    }

    return { ...baseFigure, layout };
  }, [axisRangesRef, baseFigure, mode, zoomRevision]);

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
