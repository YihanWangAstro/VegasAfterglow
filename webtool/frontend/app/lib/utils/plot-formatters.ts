/**
 * Formatting functions for frequency labels, time labels, color assignment,
 * instrument resolution, and observation style matching.
 */

import { INSTRUMENT_CATALOG } from "../constants";
import {
  C_ANGSTROM_PER_S,
  DEFAULT_INST_COLOR,
  EV_HZ,
  FREQ_PALETTE,
  KEV_HZ,
  OBS_COLORS,
  TAB10,
  TIME_PALETTE,
} from "./plot-constants";

// ---------------------------------------------------------------------------
// Superscript / shift formatting
// ---------------------------------------------------------------------------

function shiftSuperscript(n: number): string {
  const map: Record<string, string> = { "-": "\u207B", "0": "\u2070", "1": "\u00B9", "2": "\u00B2", "3": "\u00B3", "4": "\u2074", "5": "\u2075", "6": "\u2076", "7": "\u2077", "8": "\u2078", "9": "\u2079" };
  return String(n).split("").map((c) => map[c] ?? c).join("");
}

export function formatShiftSuffix(shift: number): string {
  if (shift === 1) return "";
  const exp = Math.floor(Math.log10(shift));
  const leading = Math.round(shift / Math.pow(10, exp));
  let base: string;
  if (exp === 0) base = String(leading);
  else if (leading === 1) base = `10${shiftSuperscript(exp)}`;
  else base = `${leading}\u00b710${shiftSuperscript(exp)}`;
  return ` \u00d7 ${base}`;
}

// ---------------------------------------------------------------------------
// Number formatting
// ---------------------------------------------------------------------------

/** Format number to 3 significant figures, stripping trailing zeros (matches Python :.3g). */
export function sig3(val: number): string {
  return parseFloat(val.toPrecision(3)).toString();
}

// ---------------------------------------------------------------------------
// Frequency / energy / time formatting
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// Color assignment
// ---------------------------------------------------------------------------

/** Assign colors from FREQ_PALETTE by frequency rank (warm->cool).
 *  Spreads evenly across the full palette so even a few frequencies use the full warm->cool range. */
export function freqColors(freqHz: number[]): string[] {
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
export function timeColors(n: number): string[] {
  return Array.from({ length: n }, (_, i) => TIME_PALETTE[i % TIME_PALETTE.length]);
}

/** Pick n colors cycling through the matplotlib tab10 colormap. */
export function tab10Colors(n: number): string[] {
  return Array.from({ length: n }, (_, i) => TAB10[i % TAB10.length]);
}

// ---------------------------------------------------------------------------
// Instrument resolution
// ---------------------------------------------------------------------------

export type ResolvedInstrument = { name: string; nu_min: number; nu_max: number; sensitivity: number; kind: "Fnu" | "Fband" };

export function resolveInstruments(names: string[]): ResolvedInstrument[] {
  const result: ResolvedInstrument[] = [];
  for (const name of names) {
    const entry = INSTRUMENT_CATALOG[name];
    if (entry) result.push({ name, ...entry });
  }
  return result;
}

// ---------------------------------------------------------------------------
// Observation style resolution
// ---------------------------------------------------------------------------

/**
 * Resolve color and showlegend for each obs group by matching its label
 * against the curve/snapshot labels. Matched groups adopt the curve color
 * and hide their legend entry; unmatched groups get a default palette color.
 */
export function resolveObsStyle(
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
