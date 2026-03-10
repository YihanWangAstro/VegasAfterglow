/**
 * Spectrum / SED Plotly figure builder.
 */

import type { SedPlotData } from "../types";
import {
  AXIS_COMMON,
  COMP_DASHES,
  COMP_LABELS,
  COMP_ORDER,
  DEFAULT_INST_COLOR,
  FLUX_AXIS_LABELS,
  FLUX_LATEX,
  FLUX_SCALES_CGS,
  FREQ_DISP_SCALES,
  HOVERLABEL_COMMON,
  LAYOUT_COMMON,
  LEGEND_COMMON,
  axisTitle,
  type PlotlyFigure,
} from "./plot-constants";
import {
  formatShiftSuffix,
  formatTimeLabel,
  resolveInstruments,
  resolveObsStyle,
  tab10Colors,
  timeColors,
} from "./plot-formatters";
import {
  buildYAxisConfig,
  cgsToAbMag,
  magYRange,
  smartYlim,
} from "./plot-helpers";

export type SedDisplayOptions = { fluxUnit: string; freqUnit: string; nufnu: boolean; instruments: string[] };

export function buildSedFigure(pd: SedPlotData, opts: SedDisplayOptions, obsShifts?: Map<string, number>): PlotlyFigure {
  const {
    freq_hz,
    t_snapshots_s,
    components,
    obs = [],
  } = pd;
  const instruments = resolveInstruments(opts.instruments);
  const flux_unit = opts.fluxUnit;
  const freq_unit = opts.freqUnit;
  const nufnu = opts.nufnu;
  const t_labels = t_snapshots_s.map(formatTimeLabel);
  const isABmag = flux_unit === "AB mag";
  const actualNufnu = nufnu && !isABmag;
  const fluxScale = FLUX_SCALES_CGS[flux_unit] ?? 1e-26;
  const freqScale = FREQ_DISP_SCALES[freq_unit] ?? 1;
  const nuDisp = freq_hz.map((nu) => nu / freqScale);

  const traces: Record<string, unknown>[] = [];
  const xHover = `\u03bd=%{x} ${freq_unit}`;

  const orderedComps = COMP_ORDER.filter((n) => n in components);
  const t_colors = timeColors(t_snapshots_s.length);

  for (let j = 0; j < t_snapshots_s.length; j++) {
    const label = t_labels[j];
    const color = t_colors[j];
    const shift = obsShifts?.get(label) ?? 1;
    const shiftSuffix = formatShiftSuffix(shift);
    for (const compName of orderedComps) {
      const isTotal = compName === "total";
      const baseTraceName = isTotal ? label : `${label} [${COMP_LABELS[compName] ?? compName}]`;
      const traceName = `${baseTraceName}${shiftSuffix}`;
      const fluxCgs = components[compName][j]; // [nu_idx], CGS
      let yDisp: number[], hoverY: string;
      if (isABmag) {
        yDisp = fluxCgs.map((f) => cgsToAbMag(f * shift));
        hoverY = "mag=%{y:.2f}";
      } else if (actualNufnu) {
        yDisp = freq_hz.map((nu, k) => nu * fluxCgs[k] * shift);
        hoverY = "\u03bdF\u03bd=%{y} erg/cm\u00b2/s";
      } else {
        yDisp = fluxCgs.map((f) => (f * shift) / fluxScale);
        hoverY = `F\u03bd=%{y} ${FLUX_AXIS_LABELS[flux_unit]}`;
      }
      traces.push({
        type: "scatter",
        x: nuDisp,
        y: yDisp,
        mode: "lines",
        name: traceName,
        line: {
          color,
          width: isTotal ? 1.0 : 0.9,
          dash: COMP_DASHES[compName] ?? "solid",
        },
        opacity: isTotal ? 1.0 : 0.75,
        legendgroup: traceName,
        hoverlabel: HOVERLABEL_COMMON,
        hovertemplate: `${xHover}<br>${hoverY}<extra>${traceName}</extra>`,
      });
    }
  }

  // Instrument traces for spectrum
  const sedInstPalette = tab10Colors(instruments.length);
  for (let ii = 0; ii < instruments.length; ii++) {
    const inst = instruments[ii];
    const isFnu = inst.kind === "Fnu";
    const xLo = inst.nu_min / freqScale;
    const xHi = inst.nu_max / freqScale;
    const yVal = isFnu ? (isABmag ? NaN : inst.sensitivity / fluxScale) : inst.sensitivity;
    if (!isFinite(yVal)) continue;
    const instColor = sedInstPalette[ii] ?? DEFAULT_INST_COLOR;
    traces.push({
      type: "scatter",
      x: [xLo, xHi],
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
    });
    const tickTop = isABmag ? yVal - 0.15 : yVal * Math.SQRT2;
    const tickBottom = isABmag ? yVal + 0.15 : yVal / Math.SQRT2;
    for (const xEnd of [xLo, xHi]) {
      traces.push({
        type: "scatter",
        x: [xEnd, xEnd],
        y: [tickTop, tickBottom],
        mode: "lines",
        line: { color: instColor, width: 1 },
        showlegend: false,
        hoverinfo: "skip",
        legendgroup: "inst",
      });
    }
  }

  // Obs traces for spectrum -- resolve colors from time snapshot labels/colors
  let sedObsUnmatched = 0;
  for (const group of obs) {
    const { label, fnu } = group;
    const shift = obsShifts?.get(label) ?? 1;
    const shiftSuffix = formatShiftSuffix(shift);
    const displayLabel = `${label}${shiftSuffix}`;
    const { color, showLegend } = resolveObsStyle(label, t_labels, t_colors, sedObsUnmatched);
    if (showLegend) sedObsUnmatched++;
    if (fnu.length > 0) {
      const xs = fnu.map((r) => r[0] / freqScale);
      let ys: number[], errs: number[], hoverY: string;
      if (isABmag) {
        ys = fnu.map((r) => cgsToAbMag(r[1] * shift));
        errs = fnu.map((r) => (r[1] > 0 ? (r[2] / r[1]) * 2.5 / Math.LN10 : 0));
        hoverY = "mag=%{y:.2f}";
      } else if (actualNufnu) {
        ys = fnu.map((r) => r[0] * r[1] * shift); // nu * Fnu in CGS
        errs = fnu.map((r) => r[0] * r[2] * shift);
        hoverY = "\u03bdF\u03bd=%{y} erg/cm\u00b2/s";
      } else {
        ys = fnu.map((r) => (r[1] * shift) / fluxScale);
        errs = fnu.map((r) => (r[2] * shift) / fluxScale);
        hoverY = `F\u03bd=%{y} ${FLUX_AXIS_LABELS[flux_unit]}`;
      }
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
      });
    }
  }

  // Compute y range
  const nuMin = freq_hz[0];
  const nuMax = freq_hz[freq_hz.length - 1];
  const xLo = Math.log10(nuMin / freqScale);
  const xHi = Math.log10(nuMax / freqScale);

  let yRange: [number, number] | null | "reversed" = null;
  if (isABmag) {
    yRange = magYRange(components);
  } else {
    const allFlux: number[][] = actualNufnu
      ? Object.values(components).flatMap((byT) =>
          byT.map((snapshot) => freq_hz.map((nu, k) => nu * snapshot[k])),
        )
      : Object.values(components).flatMap((byT) =>
          byT.map((snapshot) => snapshot.map((f) => f / fluxScale)),
        );
    const lim = smartYlim(allFlux);
    yRange = lim ? [Math.log10(lim[0]), Math.log10(lim[1])] : null;
  }

  const yLabel = isABmag
    ? "AB mag"
    : actualNufnu
      ? "$\\nu F_\\nu\\;(\\mathrm{erg\\,cm^{-2}\\,s^{-1}})$"
      : `$F_\\nu\\;(\\mathrm{${FLUX_LATEX[flux_unit] ?? flux_unit}})$`;

  const layout: Record<string, unknown> = {
    xaxis: {
      type: "log",
      title: axisTitle(`$\\nu\\;(\\mathrm{${freq_unit}})$`),
      range: [xLo, xHi],
      ...AXIS_COMMON,
    },
    yaxis: buildYAxisConfig(isABmag, yLabel, yRange),
    legend: LEGEND_COMMON,
    ...LAYOUT_COMMON,
  };

  return { data: traces, layout };
}
