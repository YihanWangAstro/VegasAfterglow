/**
 * Shared constants for Plotly figure builders: unit scales, palettes,
 * axis/layout/legend presets, and component style maps.
 */

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

export type PlotlyFigure = {
  data: Record<string, unknown>[];
  layout: Record<string, unknown>;
  frames?: unknown[];
};

// ---------------------------------------------------------------------------
// Unit constants
// ---------------------------------------------------------------------------

/** Physical constants used for label formatting and frequency scaling. */
export const KEV_HZ = 2.417989242084918e17;
export const C_ANGSTROM_PER_S = 2.99792458e18;
export const EV_HZ = KEV_HZ / 1e3;

export const TIME_SCALES_S: Record<string, number> = { s: 1, day: 86400, hr: 3600, min: 60 };
export const TIME_AXIS_LABELS: Record<string, string> = { s: "s", day: "day", hr: "hr", min: "min" };

// Scale: display = cgs / scale.  e.g. 1 mJy = 1e-26 CGS, so cgs / 1e-26 = mJy
export const FLUX_SCALES_CGS: Record<string, number> = {
  mJy: 1e-26,
  Jy: 1e-23,
  uJy: 1e-29,
  cgs: 1,
};

// Plain-text flux labels for hover templates
export const FLUX_AXIS_LABELS: Record<string, string> = {
  mJy: "mJy",
  Jy: "Jy",
  uJy: "\u03bcJy",
  cgs: "erg cm\u207b\u00b2 s\u207b\u00b9 Hz\u207b\u00b9",
  "AB mag": "AB mag",
};

// LaTeX flux unit labels for axis titles
export const FLUX_LATEX: Record<string, string> = {
  mJy: "mJy",
  Jy: "Jy",
  uJy: "\\mu Jy",
  cgs: "erg\\,cm^{-2}\\,s^{-1}\\,Hz^{-1}",
};

// Freq scales: display = hz / scale
export const FREQ_DISP_SCALES: Record<string, number> = {
  Hz: 1,
  GHz: 1e9,
  keV: KEV_HZ,
  MeV: KEV_HZ * 1e3,
};

// ---------------------------------------------------------------------------
// Component style maps
// ---------------------------------------------------------------------------

export const COMP_DASHES: Record<string, string> = {
  total: "solid",
  fwd_sync: "dash",
  fwd_ssc: "dot",
  rvs_sync: "dashdot",
  rvs_ssc: "longdashdot",
};

export const COMP_LABELS: Record<string, string> = {
  total: "total",
  fwd_sync: "fwd syn",
  fwd_ssc: "fwd SSC",
  rvs_sync: "rvs syn",
  rvs_ssc: "rvs SSC",
};

export const COMP_ORDER = ["total", "fwd_sync", "fwd_ssc", "rvs_sync", "rvs_ssc"];

export const FBAND_TITLE = "$F\\;(\\mathrm{erg\\,cm^{-2}\\,s^{-1}})$";

// ---------------------------------------------------------------------------
// Color palettes
// ---------------------------------------------------------------------------

// Discrete qualitative palette ordered warm -> cool (radio -> X-ray/gamma).
// Colors are assigned by frequency rank so closely-spaced frequencies stay distinct.
export const FREQ_PALETTE = [
  "#E03530", "#E8872E", "#D4A017", "#8C564B", "#BCBD22",
  "#2AB07E", "#17BECF", "#2878B5", "#1D3557", "#7B3FA0",
  "#E377C2", "#555555",
];

export const TIME_PALETTE = ["#E03530", "#E8872E", "#2AB07E", "#2878B5", "#7B3FA0", "#D4A017"];

// Unmatched obs colors -- intentionally distinct from FREQ_PALETTE to avoid confusion.
export const OBS_COLORS = [
  "#000000", "#FF1493", "#00CED1", "#FF8C00", "#4169E1",
  "#32CD32", "#DC143C", "#8B008B", "#008080", "#A0522D",
];

export const TAB10 = [
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
];

export const DEFAULT_INST_COLOR = "#757575";

// ---------------------------------------------------------------------------
// Axis / layout / legend presets
// ---------------------------------------------------------------------------

export const AXIS_TITLE_FONT = { size: 13, color: "#000" };

export function axisTitle(text: string) { return { text, font: AXIS_TITLE_FONT }; }

export const AXIS_COMMON: Record<string, unknown> = {
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

export const LEGEND_COMMON: Record<string, unknown> = {
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

export const HOVERLABEL_COMMON = { bgcolor: "white", font: { color: "black" }, bordercolor: "#ccc" };

export const LAYOUT_COMMON: Record<string, unknown> = {
  hovermode: "closest",
  autosize: true,
  margin: { l: 65, r: 20, t: 15, b: 55 },
  plot_bgcolor: "#ffffff",
  paper_bgcolor: "#ffffff",
};
