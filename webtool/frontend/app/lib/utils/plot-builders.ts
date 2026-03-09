/**
 * Client-side Plotly figure builders for LC, SED/spectrum, and skymap plots.
 * Receives raw arrays from the backend and builds Plotly-compatible
 * trace + layout objects.
 */

import { INSTRUMENT_CATALOG } from "../constants";
import type { LcPlotData, SedPlotData, SkymapPlotData } from "../types";

function shiftSuperscript(n: number): string {
  const map: Record<string, string> = { "-": "\u207B", "0": "\u2070", "1": "\u00B9", "2": "\u00B2", "3": "\u00B3", "4": "\u2074", "5": "\u2075", "6": "\u2076", "7": "\u2077", "8": "\u2078", "9": "\u2079" };
  return String(n).split("").map((c) => map[c] ?? c).join("");
}

function formatShiftSuffix(shift: number): string {
  if (shift === 1) return "";
  const exp = Math.floor(Math.log10(shift));
  const leading = Math.round(shift / Math.pow(10, exp));
  let base: string;
  if (exp === 0) base = String(leading);
  else if (leading === 1) base = `10${shiftSuperscript(exp)}`;
  else base = `${leading}\u00b710${shiftSuperscript(exp)}`;
  return ` \u00d7 ${base}`;
}

export type PlotlyFigure = {
  data: Record<string, unknown>[];
  layout: Record<string, unknown>;
  frames?: unknown[];
};

// ---------------------------------------------------------------------------
// Unit constants
// ---------------------------------------------------------------------------

export const TIME_SCALES_S: Record<string, number> = { s: 1, day: 86400, hr: 3600, min: 60 };
const TIME_AXIS_LABELS: Record<string, string> = { s: "s", day: "day", hr: "hr", min: "min" };

// Scale: display = cgs / scale.  e.g. 1 mJy = 1e-26 CGS, so cgs / 1e-26 = mJy
export const FLUX_SCALES_CGS: Record<string, number> = {
  mJy: 1e-26,
  Jy: 1e-23,
  uJy: 1e-29,
  cgs: 1,
};

// Plain-text flux labels for hover templates
const FLUX_AXIS_LABELS: Record<string, string> = {
  mJy: "mJy",
  Jy: "Jy",
  uJy: "\u03bcJy",
  cgs: "erg cm\u207b\u00b2 s\u207b\u00b9 Hz\u207b\u00b9",
  "AB mag": "AB mag",
};

// LaTeX flux unit labels for axis titles
const FLUX_LATEX: Record<string, string> = {
  mJy: "mJy",
  Jy: "Jy",
  uJy: "\\mu Jy",
  cgs: "erg\\,cm^{-2}\\,s^{-1}\\,Hz^{-1}",
};

// Physical constants used for label formatting and frequency scaling
const KEV_HZ = 2.417989242084918e17;
const C_ANGSTROM_PER_S = 2.99792458e18;

// Freq scales: display = hz / scale
export const FREQ_DISP_SCALES: Record<string, number> = {
  Hz: 1,
  GHz: 1e9,
  keV: KEV_HZ,
  MeV: KEV_HZ * 1e3,
};

const EV_HZ = KEV_HZ / 1e3;

/** Format number to 3 significant figures, stripping trailing zeros (matches Python :.3g). */
function sig3(val: number): string {
  return parseFloat(val.toPrecision(3)).toString();
}

/** Return the broad-band category for a frequency. */
function broadBand(nu: number): string {
  if (nu < 1e12) return "Radio";
  if (nu < 4e14) return "IR";
  if (nu < 7.5e14) return "Optical";
  const E_eV = nu / EV_HZ;
  if (E_eV < 100) return "UV";
  if (E_eV < 1e5) return "X-ray";
  if (E_eV < 1e9) return "\u03b3-ray";
  if (E_eV < 1e12) return "GeV";
  return "TeV";
}

/** Format a frequency as energy/wavelength/frequency string (no category prefix). */
function formatEnergy(nu: number): string {
  const E_keV = nu / KEV_HZ;
  const lam_nm = C_ANGSTROM_PER_S / nu / 10;
  const v = (val: number, unit: string) => `${sig3(val)} ${unit}`;

  if (E_keV >= 1e9) return v(E_keV / 1e9, "TeV");
  if (E_keV >= 1e6) return v(E_keV / 1e6, "GeV");
  if (E_keV >= 1e3) return v(E_keV / 1e3, "MeV");
  if (E_keV >= 0.1) return v(E_keV, "keV");
  if (lam_nm > 100 && lam_nm < 10000) return v(lam_nm, "nm");
  if (nu >= 1e12) return v(nu / 1e12, "THz");
  if (nu >= 1e9) return `${sig3(nu / 1e9)} GHz`;
  if (nu >= 1e6) return `${Math.round(nu / 1e6)} MHz`;
  return `${Math.round(nu)} Hz`;
}

/** Format a frequency with broad-band category prefix, e.g. "Radio (1 GHz)". */
export function formatFreqLabel(nu: number): string {
  return `${broadBand(nu)} (${formatEnergy(nu)})`;
}

/** Format a frequency band label. */
export function formatBandLabel(nuMin: number, nuMax: number, name: string | null): string {
  const range = `${formatEnergy(nuMin)}\u2013${formatEnergy(nuMax)}`;
  return name ? `${name} (${range})` : range;
}

/** Format observer time for display labels. */
export function formatTimeLabel(tSec: number): string {
  if (tSec >= 365.25 * 86400) return `${sig3(tSec / (365.25 * 86400))} yr`;
  if (tSec >= 86400) return `${sig3(tSec / 86400)} day`;
  if (tSec >= 3600) return `${sig3(tSec / 3600)} hr`;
  if (tSec >= 60) return `${sig3(tSec / 60)} min`;
  return `${sig3(tSec)} s`;
}

const COMP_DASHES: Record<string, string> = {
  total: "solid",
  fwd_sync: "dash",
  fwd_ssc: "dot",
  rvs_sync: "dashdot",
  rvs_ssc: "longdashdot",
};

const COMP_LABELS: Record<string, string> = {
  total: "total",
  fwd_sync: "fwd syn",
  fwd_ssc: "fwd SSC",
  rvs_sync: "rvs syn",
  rvs_ssc: "rvs SSC",
};

export const COMP_ORDER = ["total", "fwd_sync", "fwd_ssc", "rvs_sync", "rvs_ssc"];

const FBAND_TITLE = "$F\\;(\\mathrm{erg\\,cm^{-2}\\,s^{-1}})$";

// Discrete qualitative palette ordered warm → cool (radio → X-ray/gamma).
// Colors are assigned by frequency rank so closely-spaced frequencies stay distinct.
const FREQ_PALETTE = [
  "#E03530", "#E8872E", "#D4A017", "#8C564B", "#BCBD22",
  "#2AB07E", "#17BECF", "#2878B5", "#1D3557", "#7B3FA0",
  "#E377C2", "#555555",
];

const TIME_PALETTE = ["#E03530", "#E8872E", "#2AB07E", "#2878B5", "#7B3FA0", "#D4A017"];

// Unmatched obs colors — intentionally distinct from FREQ_PALETTE to avoid confusion.
const OBS_COLORS = [
  "#000000", "#FF1493", "#00CED1", "#FF8C00", "#4169E1",
  "#32CD32", "#DC143C", "#8B008B", "#008080", "#A0522D",
];

/** Assign colors from FREQ_PALETTE by frequency rank (warm→cool).
 *  Spreads evenly across the full palette so even a few frequencies use the full warm→cool range. */
function freqColors(freqHz: number[]): string[] {
  const n = freqHz.length;
  if (n === 0) return [];
  const order = freqHz
    .map((nu, i) => ({ nu, i }))
    .sort((a, b) => a.nu - b.nu);
  const colors = new Array<string>(n);
  const P = FREQ_PALETTE.length;
  const step = n <= P ? (P - 1) / Math.max(1, n - 1) : 1;
  for (let rank = 0; rank < n; rank++) {
    const idx = n <= P ? Math.round(rank * step) : rank % P;
    colors[order[rank].i] = FREQ_PALETTE[idx];
  }
  return colors;
}

/** Assign colors from TIME_PALETTE by index. */
function timeColors(n: number): string[] {
  return Array.from({ length: n }, (_, i) => TIME_PALETTE[i % TIME_PALETTE.length]);
}

type ResolvedInstrument = { name: string; nu_min: number; nu_max: number; sensitivity: number; kind: "Fnu" | "Fband" };

function resolveInstruments(names: string[]): ResolvedInstrument[] {
  const result: ResolvedInstrument[] = [];
  for (const name of names) {
    const entry = INSTRUMENT_CATALOG[name];
    if (entry) result.push({ name, ...entry });
  }
  return result;
}

/** Pick n colors cycling through the matplotlib tab10 colormap. */
const TAB10 = [
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
];
function tab10Colors(n: number): string[] {
  return Array.from({ length: n }, (_, i) => TAB10[i % TAB10.length]);
}
const DEFAULT_INST_COLOR = "#757575";

/**
 * Resolve color and showlegend for each obs group by matching its label
 * against the curve/snapshot labels. Matched groups adopt the curve color
 * and hide their legend entry; unmatched groups get a default palette color.
 */
function resolveObsStyle(
  obsLabel: string,
  curveLabels: string[],
  curveColors: string[],
  unmatchedIndex: number,
): { color: string; showLegend: boolean } {
  const idx = curveLabels.indexOf(obsLabel);
  if (idx >= 0) {
    return { color: curveColors[idx], showLegend: false };
  }
  return { color: OBS_COLORS[unmatchedIndex % OBS_COLORS.length], showLegend: true };
}
const AXIS_TITLE_FONT = { size: 13, color: "#000" };
function axisTitle(text: string) { return { text, font: AXIS_TITLE_FONT }; }

const AXIS_COMMON: Record<string, unknown> = {
  showline: true,
  linewidth: 0.8,
  linecolor: "#000",
  mirror: true,
  ticks: "inside",
  ticklen: 5,
  tickwidth: 0.8,
  tickcolor: "#000",
  minor: {
    ticks: "inside",
    ticklen: 2.5,
    showgrid: true,
    gridcolor: "rgba(0,0,0,0.06)",
    griddash: "dot",
    gridwidth: 0.3,
  },
  showgrid: true,
  gridcolor: "rgba(0,0,0,0.10)",
  griddash: "dot",
  gridwidth: 0.3,
  tickfont: { size: 11, color: "#000" },
  exponentformat: "power",
};

const LEGEND_COMMON: Record<string, unknown> = {
  orientation: "v",
  yanchor: "bottom",
  y: 0.02,
  xanchor: "left",
  x: 0.02,
  font: { size: 12, color: "#333" },
  tracegroupgap: 1,
  borderwidth: 0,
  bgcolor: "rgba(255,255,255,0.75)",
};

const LAYOUT_COMMON: Record<string, unknown> = {
  hovermode: "closest",
  autosize: true,
  margin: { l: 65, r: 20, t: 15, b: 55 },
  plot_bgcolor: "#ffffff",
  paper_bgcolor: "#ffffff",
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function cgsToAbMag(f: number): number {
  if (f <= 0) return NaN;
  return -2.5 * Math.log10(f) - 48.6;
}

/**
 * Collect all positive values from 2D arrays, compute smart log y-axis limits.
 * Returns null if no positive values found.
 */
function smartYlim(arrays: number[][]): [number, number] | null {
  let yMin = Infinity;
  let yMax = -Infinity;
  for (const arr of arrays) {
    for (const v of arr) {
      if (v > 0) {
        if (v < yMin) yMin = v;
        if (v > yMax) yMax = v;
      }
    }
  }
  if (!isFinite(yMin) || !isFinite(yMax)) return null;
  let yTop = yMax * 10;
  let yBot = yMin / 10;
  if (yTop / yBot > 1e12) {
    yBot = yTop * 1e-12;
  }
  return [yBot, yTop];
}

/** Compute AB mag y-axis range from nested component arrays (Record<name, number[][]>). */
function magYRange(components: Record<string, number[][]>): [number, number] | "reversed" {
  let magMin = Infinity;
  let magMax = -Infinity;
  for (const compFlux of Object.values(components)) {
    for (const arr of compFlux) {
      for (const f of arr) {
        const m = cgsToAbMag(f);
        if (isFinite(m)) {
          if (m < magMin) magMin = m;
          if (m > magMax) magMax = m;
        }
      }
    }
  }
  return isFinite(magMin) ? [magMax + 1, magMin - 1] : "reversed";
}

/** Build y-axis layout object for AB mag or log scale. */
function buildYAxisConfig(
  isABmag: boolean,
  label: string,
  range: [number, number] | null | "reversed",
): Record<string, unknown> {
  if (isABmag) {
    return {
      type: "linear",
      title: axisTitle(label),
      ...(range === "reversed"
        ? { autorange: "reversed" }
        : range
          ? { range }
          : { autorange: "reversed" }),
      ...AXIS_COMMON,
    };
  }
  return {
    type: "log",
    title: axisTitle(label),
    ...(range && range !== "reversed" ? { range } : {}),
    ...AXIS_COMMON,
  };
}

/** Loop-based min/max for an array (avoids stack overflow from Math.min/max(...arr)). */
function arrayMinMax(arr: number[]): [number, number] {
  let lo = Infinity;
  let hi = -Infinity;
  for (const v of arr) {
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  return [lo, hi];
}

/** Append values to target array (avoids stack overflow from push(...spread)). */
function pushAll(target: number[], source: number[]) {
  for (const v of source) target.push(v);
}

// ---------------------------------------------------------------------------
// LC figure builder
// ---------------------------------------------------------------------------

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
          hoverlabel: { bgcolor: "white", font: { color: "black" }, bordercolor: "#ccc" },
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
          hoverlabel: { bgcolor: "white", font: { color: "black" }, bordercolor: "#ccc" },
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

  // Obs traces — resolve colors from curve labels/colors
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
        marker: { color, size: 6, symbol: "circle" },
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
        marker: { color, size: 6, symbol: "diamond" },
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

// ---------------------------------------------------------------------------
// Spectrum / SED figure builder
// ---------------------------------------------------------------------------

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
          width: isTotal ? 1.2 : 0.9,
          dash: COMP_DASHES[compName] ?? "solid",
        },
        opacity: isTotal ? 1.0 : 0.75,
        legendgroup: traceName,
        hoverlabel: { bgcolor: "white", font: { color: "black" }, bordercolor: "#ccc" },
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

  // Obs traces for spectrum — resolve colors from time snapshot labels/colors
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
        marker: { color, size: 6, symbol: "circle" },
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

// ---------------------------------------------------------------------------
// Skymap figure builder
// ---------------------------------------------------------------------------

function decodeFloat32Frame(b64: string, ny: number, nx: number): number[][] {
  const binary = atob(b64);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  const floats = new Float32Array(bytes.buffer);
  const rows: number[][] = [];
  for (let r = 0; r < ny; r++) {
    const row = new Array<number>(nx);
    for (let c = 0; c < nx; c++) {
      const v = floats[r * nx + c];
      row[c] = isFinite(v) ? v : NaN;
    }
    rows.push(row);
  }
  return rows;
}

/** Decode all base64 frames upfront (expensive). Cache this result. */
export function decodeSkymapFrames(pd: SkymapPlotData): number[][][] {
  const { nx, ny, frames_b64f32 } = pd;
  return frames_b64f32.map((b64) => decodeFloat32Frame(b64, ny, nx));
}

// Skymap axes: inherit common styling but omit minor ticks, tickfont, and exponentformat.
const { minor: _, tickfont: _tf, exponentformat: _ef, ...SKYMAP_AXIS_STYLE } = AXIS_COMMON;

// Log10 offsets to convert from CGS specific intensity (erg/cm²/s/Hz/sr) to display units.
// Display value = log10(I_cgs) + offset.
const SKY_UNIT_LOG_OFFSETS: Record<string, number> = {
  "cgs": 0,
  "mJy/sr": 26,      // 1 mJy = 1e-26 erg/cm²/s/Hz
  "Jy/sr": 23,       // 1 Jy  = 1e-23
  "MJy/sr": 17,      // 1 MJy = 1e-17
  "μJy/μas²": Math.log10(1e29 / (206265e6) ** 2),  // ≈ 6.371
};

const SKY_UNIT_LABELS: Record<string, string> = {
  "cgs": "erg cm<sup>-2</sup> s<sup>-1</sup> Hz<sup>-1</sup> sr<sup>-1</sup>",
  "mJy/sr": "mJy sr<sup>-1</sup>",
  "Jy/sr": "Jy sr<sup>-1</sup>",
  "MJy/sr": "MJy sr<sup>-1</sup>",
  "μJy/μas²": "\u03bcJy \u03bcas<sup>-2</sup>",
};

const SKY_FOV_SCALES: Record<string, number> = { "μas": 1, "mas": 1e3, "arcsec": 1e6 };
const SKY_FOV_LABELS: Record<string, string> = { "μas": "\u03bcas", "mas": "mas", "arcsec": "arcsec" };

type SkymapOptions = {
  intensityUnit?: string;
  fovUnit?: string;
};

function offsetFrame(frame: number[][], offset: number): number[][] {
  if (offset === 0) return frame;
  return frame.map((row) => row.map((v) => Number.isFinite(v) ? v + offset : NaN));
}

export function buildSkymapFigure(pd: SkymapPlotData, decodedFrames: number[][][], opts?: SkymapOptions): PlotlyFigure {
  const { extent_uas, t_obs_s, nu_obs_hz, dx, dy, x0, y0, z_min, z_max } = pd;
  const t_labels = t_obs_s.map(formatTimeLabel);
  const nu_label = formatFreqLabel(nu_obs_hz);

  const fovUnit = opts?.fovUnit ?? "μas";
  const fovScale = SKY_FOV_SCALES[fovUnit] ?? 1;
  const fovLabel = SKY_FOV_LABELS[fovUnit] ?? fovUnit;
  const [x_min, x_max, y_min, y_max] = extent_uas.map((v) => v / fovScale) as [number, number, number, number];
  const dxScaled = dx / fovScale;
  const dyScaled = dy / fovScale;
  const x0Scaled = x0 / fovScale;
  const y0Scaled = y0 / fovScale;

  const unit = opts?.intensityUnit ?? "cgs";
  const logOffset = SKY_UNIT_LOG_OFFSETS[unit] ?? 0;
  const unitLabel = SKY_UNIT_LABELS[unit] ?? unit;
  const colorbarTitle = `log<sub>10</sub> I (${unitLabel})`;

  const zMinDisp = z_min + logOffset;
  const zMaxDisp = z_max + logOffset;

  const offsetFrames = decodedFrames.map((f) => offsetFrame(f, logOffset));

  const data: Record<string, unknown>[] = [{
    type: "heatmap",
    z: offsetFrames[0], x0: x0Scaled, dx: dxScaled, y0: y0Scaled, dy: dyScaled,
    colorscale: "Electric",
    zmin: zMinDisp, zmax: zMaxDisp,
    zauto: false,
    colorbar: {
      title: { text: colorbarTitle, side: "right", font: { size: 12 } },
      xref: "paper",
      x: 1.01,
      xanchor: "left",
      len: 0.9,
      thickness: 28,
      tickfont: { size: 12 },
      dtick: 1,
    },
    hovertemplate: `\u0394x=%{x} ${fovLabel}<br>\u0394y=%{y} ${fovLabel}<br>log\u2081\u2080 I=%{z:.3g}<extra></extra>`,
  }];

  const isAnimated = decodedFrames.length > 1;

  const frames: Record<string, unknown>[] = offsetFrames.map((z, i) => ({
    name: `frame_${i}`,
    data: [{ z }],
    traces: [0],
  }));

  const layout: Record<string, unknown> = {
    title: isAnimated ? `\u03bd = ${nu_label}` : `t = ${t_labels[0]}, \u03bd = ${nu_label}`,
    xaxis: { title: axisTitle(`\u0394x (${fovLabel})`), range: [x_min, x_max], autorange: false, ...SKYMAP_AXIS_STYLE },
    yaxis: {
      title: axisTitle(`\u0394y (${fovLabel})`), range: [y_min, y_max], autorange: false,
      ...SKYMAP_AXIS_STYLE,
    },
    template: "none",
    plot_bgcolor: "#ffffff",
    paper_bgcolor: "#ffffff",
    margin: decodedFrames.length > 1 ? { l: 60, r: 120, t: 50, b: 130 } : { l: 60, r: 120, t: 50, b: 60 },
  };

  if (decodedFrames.length > 1) {
    const steps = t_labels.map((label, i) => ({
      label,
      method: "animate",
      args: [[`frame_${i}`], { mode: "immediate", frame: { duration: 0, redraw: true }, transition: { duration: 0 } }],
    }));
    layout.sliders = [{
      active: 0, x: 0.13, y: -0.08, len: 0.8, pad: { t: 8, b: 0 },
      currentvalue: { prefix: "t = ", font: { size: 13 }, xanchor: "left" }, steps,
    }];
  }

  return { data, layout, frames: decodedFrames.length > 1 ? frames : [] };
}
