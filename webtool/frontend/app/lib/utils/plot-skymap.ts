/**
 * Skymap Plotly figure builder.
 */

import type { SkymapPlotData } from "../types";
import {
  AXIS_COMMON,
  axisTitle,
  type PlotlyFigure,
} from "./plot-constants";
import { formatFreqLabel, formatTimeLabel } from "./plot-formatters";

// Skymap axes: inherit common styling but omit minor ticks, tickfont, and exponentformat.
const { minor: _, tickfont: _tf, exponentformat: _ef, ...SKYMAP_AXIS_STYLE } = AXIS_COMMON;

// Log10 offsets to convert from CGS specific intensity (erg/cm^2/s/Hz/sr) to display units.
// Display value = log10(I_cgs) + offset.
const SKY_UNIT_LOG_OFFSETS: Record<string, number> = {
  "cgs": 0,
  "mJy/sr": 26,      // 1 mJy = 1e-26 erg/cm^2/s/Hz
  "Jy/sr": 23,       // 1 Jy  = 1e-23
  "MJy/sr": 17,      // 1 MJy = 1e-17
  "μJy/μas²": Math.log10(1e29 / (206265e6) ** 2),  // ~= 6.371
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

export type SkymapOptions = {
  intensityUnit?: string;
  fovUnit?: string;
};

// ---------------------------------------------------------------------------
// Frame decoding
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

/** Decode the base64 frame. Cache this result. */
export function decodeSkymapFrame(pd: SkymapPlotData): number[][] {
  return decodeFloat32Frame(pd.frame_b64f32, pd.ny, pd.nx);
}

function offsetFrame(frame: number[][], offset: number): number[][] {
  if (offset === 0) return frame;
  return frame.map((row) => row.map((v) => Number.isFinite(v) ? v + offset : NaN));
}

// ---------------------------------------------------------------------------
// Figure builder
// ---------------------------------------------------------------------------

export function buildSkymapFigure(pd: SkymapPlotData, decodedFrame: number[][], opts?: SkymapOptions): PlotlyFigure {
  const { extent_uas, t_obs_s, nu_obs_hz, dx, dy, x0, y0, z_min, z_max } = pd;
  const t_label = formatTimeLabel(t_obs_s);
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

  const z = offsetFrame(decodedFrame, logOffset);

  const data: Record<string, unknown>[] = [{
    type: "heatmap",
    z, x0: x0Scaled, dx: dxScaled, y0: y0Scaled, dy: dyScaled,
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

  const layout: Record<string, unknown> = {
    title: `t = ${t_label}, \u03bd = ${nu_label}`,
    xaxis: { title: axisTitle(`\u0394x (${fovLabel})`), range: [x_min, x_max], autorange: false, ...SKYMAP_AXIS_STYLE },
    yaxis: {
      title: axisTitle(`\u0394y (${fovLabel})`), range: [y_min, y_max], autorange: false,
      ...SKYMAP_AXIS_STYLE,
    },
    template: "none",
    plot_bgcolor: "#ffffff",
    paper_bgcolor: "#ffffff",
    margin: { l: 60, r: 120, t: 50, b: 60 },
  };

  return { data, layout, frames: [] };
}
