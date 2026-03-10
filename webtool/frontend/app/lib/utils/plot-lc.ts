/**
 * Light curve Plotly figure builder.
 */

import type { LcPlotData } from "../types";
import {
  AXIS_COMMON,
  AXIS_TITLE_FONT,
  COMP_DASHES,
  COMP_LABELS,
  COMP_ORDER,
  DEFAULT_INST_COLOR,
  FBAND_TITLE,
  FLUX_AXIS_LABELS,
  FLUX_LATEX,
  FLUX_SCALES_CGS,
  HOVERLABEL_COMMON,
  LAYOUT_COMMON,
  LEGEND_COMMON,
  TIME_AXIS_LABELS,
  TIME_SCALES_S,
  axisTitle,
  type PlotlyFigure,
} from "./plot-constants";
import {
  formatBandLabel,
  formatFreqLabel,
  formatShiftSuffix,
  freqColors,
  resolveInstruments,
  resolveObsStyle,
  tab10Colors,
} from "./plot-formatters";
import {
  arrayMinMax,
  buildYAxisConfig,
  cgsToAbMag,
  magYRange,
  pushAll,
  smartYlim,
} from "./plot-helpers";

export type LcDisplayOptions = { fluxUnit: string; timeUnit: string; instruments: string[] };

export function buildLcFigure(pd: LcPlotData, opts: LcDisplayOptions, obsShifts?: Map<string, number>): PlotlyFigure {
  const { t_min_s, t_max_s, times_s, pt, bands, obs = [] } = pd;
  const instruments = resolveInstruments(opts.instruments);
  const flux_unit = opts.fluxUnit;
  const time_unit = opts.timeUnit;
  const isABmag = flux_unit === "AB mag";
  const tScale = TIME_SCALES_S[time_unit] ?? 86400;
  const fluxScale = FLUX_SCALES_CGS[flux_unit] ?? 1e-26;
  const tDisp = times_s.map((t) => t / tScale);
  const hasPt = pt != null && pt.freq_hz.length > 0;
  const hasBands = bands.length > 0;
  const hasFbandInst = instruments.some((i) => i.kind === "Fband");
  const useSecondary = hasPt && (hasBands || hasFbandInst);

  // Compute colors for all frequencies (pt + band centroids) together.
  const allFreqs = [...(pt?.freq_hz ?? []), ...bands.map((b) => b.nu_cen)];
  const allColors = freqColors(allFreqs);
  const ptColors = allColors.slice(0, pt?.freq_hz.length ?? 0);
  const bandColors = allColors.slice(pt?.freq_hz.length ?? 0);

  // Compute labels once for reuse in traces and obs matching.
  const ptLabels = pt?.freq_hz.map(formatFreqLabel) ?? [];
  const bandLabels = bands.map((b) => formatBandLabel(b.nu_min, b.nu_max, b.name));

  const traces: Record<string, unknown>[] = [];
  const tUnit = TIME_AXIS_LABELS[time_unit] ?? time_unit;
  const xHover = `t=%{x} ${tUnit}`;

  // PT traces (per frequency)
  if (hasPt && pt) {
    const ptOrderedComps = COMP_ORDER.filter((n) => n in pt.components);
    for (let i = 0; i < pt.freq_hz.length; i++) {
      const label = ptLabels[i];
      const color = ptColors[i];
      const shift = obsShifts?.get(label) ?? 1;
      const shiftSuffix = formatShiftSuffix(shift);
      for (const compName of ptOrderedComps) {
        const isTotal = compName === "total";
        const baseTraceName = isTotal ? label : `${label} [${COMP_LABELS[compName] ?? compName}]`;
        const traceName = `${baseTraceName}${shiftSuffix}`;
        const fluxCgs = pt.components[compName][i];
        const yDisp = isABmag ? fluxCgs.map((f) => cgsToAbMag(f * shift)) : fluxCgs.map((f) => (f * shift) / fluxScale);
        const yUnit = isABmag ? "mag" : FLUX_AXIS_LABELS[flux_unit];
        const hoverY = isABmag ? "mag=%{y:.2f}" : `F\u03bd=%{y} ${yUnit}`;
        traces.push({
          type: tDisp.length >= 400 ? "scattergl" : "scatter",
          x: tDisp,
          y: yDisp,
          mode: "lines",
          name: traceName,
          line: {
            color,
            width: isTotal ? 1.2 : 0.9,
            dash: COMP_DASHES[compName] ?? "solid",
          },
          opacity: isTotal ? 1.0 : 0.75,
          legendgroup: traceName,
          hoverlabel: HOVERLABEL_COMMON,
          hovertemplate: `${xHover}<br>${hoverY}<extra>${traceName}</extra>`,
          ...(useSecondary ? { yaxis: "y" } : {}),
        });
      }
    }
  }

  // Band traces
  if (hasBands) {
    for (let bIdx = 0; bIdx < bands.length; bIdx++) {
      const band = bands[bIdx];
      const bandLabel = bandLabels[bIdx];
      const shift = obsShifts?.get(bandLabel) ?? 1;
      const shiftSuffix = formatShiftSuffix(shift);
      const orderedComps = COMP_ORDER.filter((n) => n in band.components);
      for (const compName of orderedComps) {
        const isTotal = compName === "total";
        const baseTraceName = isTotal
          ? bandLabel
          : `${bandLabel} [${COMP_LABELS[compName] ?? compName}]`;
        const traceName = `${baseTraceName}${shiftSuffix}`;
        const fluxCgs = band.components[compName];
        const yDisp = shift !== 1 ? fluxCgs.map((f: number) => f * shift) : fluxCgs;
        traces.push({
          type: tDisp.length >= 400 ? "scattergl" : "scatter",
          x: tDisp,
          y: yDisp,
          mode: "lines",
          name: traceName,
          line: {
            color: bandColors[bIdx],
            width: isTotal ? 1.2 : 0.9,
            dash: COMP_DASHES[compName] ?? "solid",
          },
          opacity: isTotal ? 1.0 : 0.75,
          legendgroup: traceName,
          hoverlabel: HOVERLABEL_COMMON,
          hovertemplate: `${xHover}<br>F=%{y} erg/cm\u00b2/s<extra>${traceName}</extra>`,
          ...(useSecondary ? { yaxis: "y2" } : hasPt ? { yaxis: "y" } : {}),
        });
      }
    }
  }

  // Instrument sensitivity traces
  const instPalette = tab10Colors(instruments.length);
  const tLo = t_min_s / tScale;
  const tHi = t_max_s / tScale;
  const logXRange = Math.log10(tHi / tLo);
  type Arrow = { xArrow: number; yVal: number; instColor: string; isFnu: boolean; instYAxis: Record<string, unknown> };
  const pendingArrows: Arrow[] = [];
  for (let ii = 0; ii < instruments.length; ii++) {
    const inst = instruments[ii];
    const isFnu = inst.kind === "Fnu";
    const yVal = isFnu ? (isABmag ? NaN : inst.sensitivity / fluxScale) : inst.sensitivity;
    if (!isFinite(yVal)) continue;
    const instColor = instPalette[ii] ?? DEFAULT_INST_COLOR;
    const instYAxis = useSecondary
      ? { yaxis: isFnu ? "y" : "y2" }
      : hasPt
        ? { yaxis: "y" }
        : {};
    traces.push({
      type: "scatter",
      x: [tLo, tHi],
      y: [yVal, yVal],
      mode: "lines",
      name: inst.name,
      legendgroup: "inst",
      legendgrouptitle: { text: "Instruments" },
      showlegend: true,
      line: { color: instColor, width: 1, dash: "dot" },
      hovertemplate: isFnu
        ? `${inst.name}<br>F\u03bd=%{y}<extra></extra>`
        : `${inst.name}<br>F=%{y} erg/cm\u00b2/s<extra></extra>`,
      ...instYAxis,
    });
    const xArrow = isFnu ? tLo * Math.pow(10, logXRange * 0.02) : tHi / Math.pow(10, logXRange * 0.02);
    pendingArrows.push({ xArrow, yVal, instColor, isFnu, instYAxis });
  }

  // Obs traces -- resolve colors from curve labels/colors
  const curveLabels = [...ptLabels, ...bandLabels];
  const curveColors = [...ptColors, ...bandColors];
  const allFnuYs: number[] = [];
  let obsUnmatched = 0;
  for (const group of obs) {
    const { label, fnu, fband } = group;
    const shift = obsShifts?.get(label) ?? 1;
    const shiftSuffix = formatShiftSuffix(shift);
    const displayLabel = `${label}${shiftSuffix}`;
    const { color, showLegend } = resolveObsStyle(label, curveLabels, curveColors, obsUnmatched);
    if (showLegend) obsUnmatched++;
    if (fnu.length > 0) {
      const xs = fnu.map((r) => r[0] / tScale);
      let ys: number[], errs: number[], hoverY: string;
      if (isABmag) {
        ys = fnu.map((r) => cgsToAbMag(r[1] * shift));
        errs = fnu.map((r) => (r[1] > 0 ? (r[2] / r[1]) * 2.5 / Math.LN10 : 0));
        hoverY = "mag=%{y:.2f}";
      } else {
        ys = fnu.map((r) => (r[1] * shift) / fluxScale);
        errs = fnu.map((r) => (r[2] * shift) / fluxScale);
        hoverY = `F\u03bd=%{y} ${FLUX_AXIS_LABELS[flux_unit]}`;
      }
      pushAll(allFnuYs, ys.filter((v) => isFinite(v)));
      traces.push({
        type: "scatter",
        x: xs,
        y: ys,
        mode: "markers",
        name: displayLabel,
        legendgroup: `obs_${label}`,
        showlegend: showLegend,
        marker: { color, size: 8, symbol: "circle", line: { color: "black", width: 0.5 } },
        error_y: { type: "data", array: errs, visible: true, color, thickness: 1.0, width: 3 },
        hovertemplate: `${xHover}<br>${hoverY}<extra>${displayLabel}</extra>`,
        ...(useSecondary ? { yaxis: "y" } : {}),
      });
    }
    if (fband.length > 0 && (useSecondary || (!hasPt && !hasBands))) {
      const xs = fband.map((r) => r[0] / tScale);
      const ys = fband.map((r) => r[1] * shift);
      const errs = fband.map((r) => r[2] * shift);
      traces.push({
        type: "scatter",
        x: xs,
        y: ys,
        mode: "markers",
        name: displayLabel,
        legendgroup: `obs_${label}`,
        showlegend: showLegend,
        marker: { color, size: 8, symbol: "diamond", line: { color: "black", width: 0.5 } },
        error_y: { type: "data", array: errs, visible: true, color, thickness: 1.0, width: 3 },
        hovertemplate: `${xHover}<br>F=%{y} erg/cm\u00b2/s<extra>${displayLabel}</extra>`,
        yaxis: "y2",
      });
    }
  }

  // Compute y-axis ranges
  const xLo = Math.log10(t_min_s / tScale);
  const xHi = Math.log10(t_max_s / tScale);

  let yRangePt: [number, number] | null | "reversed" = null;
  if (hasPt && pt) {
    if (isABmag) {
      yRangePt = magYRange(pt.components);
    } else {
      const allFlux: number[][] = Object.values(pt.components).flatMap((byFreq) =>
        byFreq.map((ts) => ts.map((f) => f / fluxScale)),
      );
      const lim = smartYlim(allFlux);
      yRangePt = lim ? [Math.log10(lim[0]), Math.log10(lim[1])] : null;
    }
  }

  // Expand y range to include obs data
  if (!isABmag && allFnuYs.length > 0 && yRangePt && yRangePt !== "reversed") {
    const posObs = allFnuYs.filter((v) => v > 0);
    if (posObs.length > 0) {
      const [obsMin, obsMax] = arrayMinMax(posObs);
      const newLo = Math.min(yRangePt[0], Math.log10(obsMin * 0.3));
      const newHi = Math.max(yRangePt[1], Math.log10(obsMax * 3));
      yRangePt = [newLo, newHi];
    }
  }

  let yRangeBd: [number, number] | null = null;
  if (hasBands) {
    const allBandFlux: number[][] = bands.flatMap((b) => Object.values(b.components));
    const lim = smartYlim(allBandFlux);
    yRangeBd = lim ? [Math.log10(lim[0]), Math.log10(lim[1])] : null;
  }

  // If no band data but Fband instrument limits exist, set a reasonable y2 range
  if (!yRangeBd) {
    const fbandSens = instruments
      .filter((inst) => inst.kind === "Fband" && isFinite(inst.sensitivity) && inst.sensitivity > 0)
      .map((inst) => Math.log10(inst.sensitivity));
    if (fbandSens.length > 0) {
      const [lo, hi] = arrayMinMax(fbandSens);
      yRangeBd = [lo - 2, hi + 2];
    }
  }

  // Build layout
  const xAxis: Record<string, unknown> = {
    type: "log",
    title: axisTitle(`$t_{\\rm obs} - t_0\\;(\\mathrm{${tUnit}})$`),
    range: [xLo, xHi],
    ...AXIS_COMMON,
  };

  const yAxisPt = buildYAxisConfig(
    isABmag,
    isABmag ? "AB mag" : `$F_\\nu\\;(\\mathrm{${FLUX_LATEX[flux_unit] ?? flux_unit}})$`,
    yRangePt,
  );

  // Add arrow traces now that y-ranges are known
  for (const a of pendingArrows) {
    const yRange = a.isFnu ? yRangePt : yRangeBd;
    const span = yRange && yRange !== "reversed" ? yRange[1] - yRange[0] : 4;
    const arrowTop = a.yVal * Math.pow(10, span * 0.02);
    traces.push({
      type: "scatter",
      x: [a.xArrow, a.xArrow],
      y: [a.yVal, arrowTop],
      mode: "lines+markers",
      line: { color: a.instColor, width: 1 },
      marker: { symbol: ["circle", "triangle-up"], size: [0, 7], color: a.instColor },
      showlegend: false,
      hoverinfo: "skip",
      legendgroup: "inst",
      ...a.instYAxis,
    });
  }

  const layout: Record<string, unknown> = {
    xaxis: xAxis,
    legend: LEGEND_COMMON,
    ...LAYOUT_COMMON,
  };

  if (useSecondary) {
    layout.yaxis = yAxisPt;
    const y2Range = yRangeBd ? { range: yRangeBd } : {};
    layout.yaxis2 = {
      type: "log",
      ...y2Range,
      overlaying: "y",
      side: "right",
      ...AXIS_COMMON,
      title: { text: FBAND_TITLE, font: AXIS_TITLE_FONT, standoff: 8 },
      showgrid: false,
    };
    layout.margin = { l: 65, r: 75, t: 15, b: 55 };
  } else if (hasBands && !hasPt) {
    layout.yaxis = {
      type: "log",
      title: axisTitle(FBAND_TITLE),
      ...(yRangeBd ? { range: yRangeBd } : {}),
      ...AXIS_COMMON,
    };
  } else {
    layout.yaxis = yAxisPt;
  }

  return { data: traces, layout };
}
