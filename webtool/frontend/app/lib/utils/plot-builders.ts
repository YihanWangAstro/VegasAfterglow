/**
 * Client-side Plotly figure builders for LC, SED/spectrum, and skymap plots.
 * Receives raw arrays from the backend and builds Plotly-compatible
 * trace + layout objects.
 */

import type { LcPlotData, SedPlotData, SkymapPlotData } from "../types";

export type PlotlyFigure = {
  data: Record<string, unknown>[];
  layout: Record<string, unknown>;
  frames?: unknown[];
};

// ---------------------------------------------------------------------------
// Unit constants
// ---------------------------------------------------------------------------

const TIME_SCALES_S: Record<string, number> = { s: 1, day: 86400, hr: 3600, min: 60 };
const TIME_AXIS_LABELS: Record<string, string> = { s: "s", day: "day", hr: "hr", min: "min" };

// Scale: display = cgs / scale.  e.g. 1 mJy = 1e-26 CGS, so cgs / 1e-26 = mJy
const FLUX_SCALES_CGS: Record<string, number> = {
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
const FREQ_DISP_SCALES: Record<string, number> = {
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
  const lo = formatFreqLabel(nuMin);
  const hi = formatFreqLabel(nuMax);
  const range = `${lo}-${hi}`;
  return name ? `${name}(${range})` : range;
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
  fwd_ssc: "fwd IC",
  rvs_sync: "rvs syn",
  rvs_ssc: "rvs IC",
};

const COMP_ORDER = ["total", "fwd_sync", "fwd_ssc", "rvs_sync", "rvs_ssc"];

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

const INSTRUMENT_COLORS: Record<string, string> = {
  "VLA": "#E63946", "ALMA": "#457B9D", "MeerKAT": "#2A9D8F", "ngVLA": "#E9C46A",
  "Rubin/LSST": "#F4A261", "JWST": "#264653", "WFST": "#8B5E3C", "SVOM/VT": "#D4A017",
  "Swift/XRT": "#A8DADC", "Chandra": "#1D3557", "EP/WXT": "#FF6B6B", "EP/FXT": "#4ECDC4",
  "SVOM/MXT": "#45B7D1", "SVOM/ECLAIRs": "#96CEB4", "Swift/BAT": "#7A5C61",
  "Fermi/GBM": "#BC6C25", "Fermi/LAT": "#6A0572", "CTA": "#C44536",
};
const DEFAULT_INST_COLOR = "#888888";

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

export function buildLcFigure(pd: LcPlotData): PlotlyFigure {
  const { flux_unit, time_unit, t_min_s, t_max_s, times_s, pt, bands, obs, instruments } = pd;
  const isABmag = flux_unit === "AB mag";
  const tScale = TIME_SCALES_S[time_unit] ?? 86400;
  const fluxScale = FLUX_SCALES_CGS[flux_unit] ?? 1e-26;
  const tDisp = times_s.map((t) => t / tScale);
  const hasPt = pt != null && pt.freq_hz.length > 0;
  const hasBands = bands.length > 0;
  const useSecondary = hasPt && hasBands;

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
      for (const compName of ptOrderedComps) {
        const isTotal = compName === "total";
        const traceName = isTotal ? label : `${label} (${COMP_LABELS[compName] ?? compName})`;
        const fluxCgs = pt.components[compName][i];
        const yDisp = isABmag ? fluxCgs.map(cgsToAbMag) : fluxCgs.map((f) => f / fluxScale);
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
      const orderedComps = COMP_ORDER.filter((n) => n in band.components);
      for (const compName of orderedComps) {
        const isTotal = compName === "total";
        const traceName = isTotal
          ? bandLabel
          : `${bandLabel} (${COMP_LABELS[compName] ?? compName})`;
        const fluxCgs = band.components[compName];
        traces.push({
          type: tDisp.length >= 400 ? "scattergl" : "scatter",
          x: tDisp,
          y: fluxCgs,
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
  for (const inst of instruments) {
    const isFnu = inst.kind === "Fnu";
    const tLo = inst.t_lo_s / tScale;
    const tHi = inst.t_hi_s / tScale;
    const yVal = isFnu ? (isABmag ? NaN : inst.sensitivity / fluxScale) : inst.sensitivity;
    if (!isFinite(yVal)) continue;
    traces.push({
      type: "scatter",
      x: [tLo, tHi],
      y: [yVal, yVal],
      mode: "lines",
      name: inst.name,
      legendgroup: "inst",
      legendgrouptitle: { text: "Instruments" },
      showlegend: true,
      line: { color: INSTRUMENT_COLORS[inst.name] ?? DEFAULT_INST_COLOR, width: 1, dash: isFnu ? "dash" : "solid" },
      hovertemplate: isFnu
        ? `${inst.name}<br>F\u03bd=%{y}<extra></extra>`
        : `${inst.name}<br>F=%{y} erg/cm\u00b2/s<extra></extra>`,
      ...(useSecondary
        ? { yaxis: isFnu ? "y" : "y2" }
        : hasPt
          ? { yaxis: "y" }
          : {}),
    });
  }

  // Obs traces — resolve colors from curve labels/colors
  const curveLabels = [...ptLabels, ...bandLabels];
  const curveColors = [...ptColors, ...bandColors];
  const allFnuYs: number[] = [];
  let obsUnmatched = 0;
  for (const group of obs) {
    const { label, fnu, fband } = group;
    const { color, showLegend } = resolveObsStyle(label, curveLabels, curveColors, obsUnmatched);
    if (showLegend) obsUnmatched++;
    if (fnu.length > 0) {
      const xs = fnu.map((r) => r[0] / tScale);
      let ys: number[], errs: number[], hoverY: string;
      if (isABmag) {
        ys = fnu.map((r) => cgsToAbMag(r[1]));
        errs = fnu.map((r) => (r[1] > 0 ? (r[2] / r[1]) * 2.5 / Math.LN10 : 0));
        hoverY = "mag=%{y:.2f}";
      } else {
        ys = fnu.map((r) => r[1] / fluxScale);
        errs = fnu.map((r) => r[2] / fluxScale);
        hoverY = `F\u03bd=%{y} ${FLUX_AXIS_LABELS[flux_unit]}`;
      }
      pushAll(allFnuYs, ys.filter((v) => isFinite(v)));
      traces.push({
        type: "scatter",
        x: xs,
        y: ys,
        mode: "markers",
        name: label,
        legendgroup: `obs_${label}`,
        showlegend: showLegend,
        marker: { color, size: 6, symbol: "circle" },
        error_y: { type: "data", array: errs, visible: true, color, thickness: 1.0, width: 3 },
        hovertemplate: `${xHover}<br>${hoverY}<extra>${label}</extra>`,
        ...(useSecondary ? { yaxis: "y" } : {}),
      });
    }
    if (fband.length > 0 && (useSecondary || (!hasPt && !hasBands))) {
      const xs = fband.map((r) => r[0] / tScale);
      const ys = fband.map((r) => r[1]);
      const errs = fband.map((r) => r[2]);
      traces.push({
        type: "scatter",
        x: xs,
        y: ys,
        mode: "markers",
        name: label,
        legendgroup: `obs_${label}`,
        showlegend: showLegend,
        marker: { color, size: 6, symbol: "diamond" },
        error_y: { type: "data", array: errs, visible: true, color, thickness: 1.0, width: 3 },
        hovertemplate: `${xHover}<br>F=%{y} erg/cm\u00b2/s<extra>${label}</extra>`,
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

  // Build layout
  const xAxis: Record<string, unknown> = {
    type: "log",
    title: axisTitle(`$t_{\\rm obs}\\;(\\mathrm{${tUnit}})$`),
    range: [xLo, xHi],
    ...AXIS_COMMON,
  };

  const yAxisPt = buildYAxisConfig(
    isABmag,
    isABmag ? "AB mag" : `$F_\\nu\\;(\\mathrm{${FLUX_LATEX[flux_unit] ?? flux_unit}})$`,
    yRangePt,
  );

  const layout: Record<string, unknown> = {
    xaxis: xAxis,
    legend: LEGEND_COMMON,
    ...LAYOUT_COMMON,
  };

  if (useSecondary) {
    layout.yaxis = yAxisPt;
    layout.yaxis2 = {
      type: "log",
      ...(yRangeBd ? { range: yRangeBd } : {}),
      overlaying: "y",
      side: "right",
      ...AXIS_COMMON,
      title: axisTitle(FBAND_TITLE),
      showgrid: false,
    };
    layout.margin = { l: 65, r: 65, t: 15, b: 55 };
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

export function buildSedFigure(pd: SedPlotData): PlotlyFigure {
  const {
    flux_unit,
    freq_unit,
    nufnu,
    freq_hz,
    t_snapshots_s,
    components,
    obs,
    instruments,
  } = pd;
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
    for (const compName of orderedComps) {
      const isTotal = compName === "total";
      const traceName = isTotal ? label : `${label} (${COMP_LABELS[compName] ?? compName})`;
      const fluxCgs = components[compName][j]; // [nu_idx], CGS
      let yDisp: number[], hoverY: string;
      if (isABmag) {
        yDisp = fluxCgs.map(cgsToAbMag);
        hoverY = "mag=%{y:.2f}";
      } else if (actualNufnu) {
        yDisp = freq_hz.map((nu, k) => nu * fluxCgs[k]);
        hoverY = "\u03bdF\u03bd=%{y} erg/cm\u00b2/s";
      } else {
        yDisp = fluxCgs.map((f) => f / fluxScale);
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
  for (const inst of instruments) {
    const isFnu = inst.kind === "Fnu";
    const xLo = inst.nu_min / freqScale;
    const xHi = inst.nu_max / freqScale;
    const yVal = isFnu ? (isABmag ? NaN : inst.sensitivity / fluxScale) : inst.sensitivity;
    if (!isFinite(yVal)) continue;
    traces.push({
      type: "scatter",
      x: [xLo, xHi],
      y: [yVal, yVal],
      mode: "lines",
      name: inst.name,
      legendgroup: "inst",
      legendgrouptitle: { text: "Instruments" },
      showlegend: true,
      line: { color: INSTRUMENT_COLORS[inst.name] ?? DEFAULT_INST_COLOR, width: 1, dash: isFnu ? "dash" : "solid" },
      hovertemplate: isFnu
        ? `${inst.name}<br>F\u03bd=%{y}<extra></extra>`
        : `${inst.name}<br>F=%{y} erg/cm\u00b2/s<extra></extra>`,
    });
  }

  // Obs traces for spectrum — resolve colors from time snapshot labels/colors
  let sedObsUnmatched = 0;
  for (const group of obs) {
    const { label, fnu } = group;
    const { color, showLegend } = resolveObsStyle(label, t_labels, t_colors, sedObsUnmatched);
    if (showLegend) sedObsUnmatched++;
    if (fnu.length > 0) {
      const xs = fnu.map((r) => r[0] / freqScale);
      let ys: number[], errs: number[], hoverY: string;
      if (isABmag) {
        ys = fnu.map((r) => cgsToAbMag(r[1]));
        errs = fnu.map((r) => (r[1] > 0 ? (r[2] / r[1]) * 2.5 / Math.LN10 : 0));
        hoverY = "mag=%{y:.2f}";
      } else if (actualNufnu) {
        ys = fnu.map((r) => r[0] * r[1]); // nu * Fnu in CGS
        errs = fnu.map((r) => r[0] * r[2]);
        hoverY = "\u03bdF\u03bd=%{y} erg/cm\u00b2/s";
      } else {
        ys = fnu.map((r) => r[1] / fluxScale);
        errs = fnu.map((r) => r[2] / fluxScale);
        hoverY = `F\u03bd=%{y} ${FLUX_AXIS_LABELS[flux_unit]}`;
      }
      traces.push({
        type: "scatter",
        x: xs,
        y: ys,
        mode: "markers",
        name: label,
        legendgroup: `obs_${label}`,
        showlegend: showLegend,
        marker: { color, size: 6, symbol: "circle" },
        error_y: { type: "data", array: errs, visible: true, color, thickness: 1.0, width: 3 },
        hovertemplate: `${xHover}<br>${hoverY}<extra>${label}</extra>`,
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

// Skymap axes: inherit common styling but omit minor ticks, tickfont, and exponentformat.
const { minor: _, tickfont: _tf, exponentformat: _ef, ...SKYMAP_AXIS_STYLE } = AXIS_COMMON;

export function buildSkymapFigure(pd: SkymapPlotData): PlotlyFigure {
  const { nx, ny, frames_b64f32, extent_uas, t_obs_s, nu_obs_hz, dx, dy, x0, y0, z_min, z_max } = pd;
  const t_labels = t_obs_s.map(formatTimeLabel);
  const nu_label = formatFreqLabel(nu_obs_hz);
  const [x_min, x_max, y_min, y_max] = extent_uas;

  const z0 = decodeFloat32Frame(frames_b64f32[0], ny, nx);
  const colorbarTitle = "log<sub>10</sub> I (erg cm<sup>-2</sup> s<sup>-1</sup> Hz<sup>-1</sup> sr<sup>-1</sup>)";

  const data: Record<string, unknown>[] = [{
    type: "heatmap",
    z: z0, x0, dx, y0, dy,
    colorscale: "Electric",
    zmin: z_min, zmax: z_max,
    zauto: false,
    colorbar: {
      title: { text: colorbarTitle, side: "right", font: { size: 12 } },
      xref: "paper",
      x: 1.01,
      xanchor: "left",
      len: 0.85,
      thickness: 28,
      tickfont: { size: 12 },
    },
    hovertemplate: "\u0394x=%{x} \u03bcas<br>\u0394y=%{y} \u03bcas<br>log\u2081\u2080 I=%{z:.3g}<extra></extra>",
  }];

  const isAnimated = frames_b64f32.length > 1;

  const frames: Record<string, unknown>[] = frames_b64f32.map((b64, i) => ({
    name: `frame_${i}`,
    data: [{ z: i === 0 ? z0 : decodeFloat32Frame(b64, ny, nx) }],
    traces: [0],
  }));

  const layout: Record<string, unknown> = {
    title: isAnimated ? `\u03bd = ${nu_label}` : `t = ${t_labels[0]}, \u03bd = ${nu_label}`,
    xaxis: { title: axisTitle("\u0394x (\u03bcas)"), range: [x_min, x_max], autorange: false, ...SKYMAP_AXIS_STYLE },
    yaxis: {
      title: axisTitle("\u0394y (\u03bcas)"), range: [y_min, y_max], autorange: false,
      ...SKYMAP_AXIS_STYLE,
    },
    template: "none",
    plot_bgcolor: "#ffffff",
    paper_bgcolor: "#ffffff",
    margin: frames_b64f32.length > 1 ? { l: 60, r: 120, t: 50, b: 130 } : { l: 60, r: 120, t: 50, b: 60 },
  };

  if (frames_b64f32.length > 1) {
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

  return { data, layout, frames: frames_b64f32.length > 1 ? frames : [] };
}
