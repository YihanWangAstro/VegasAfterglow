/**
 * Barrel re-export for plot builder modules.
 *
 * All plotting logic has been split into focused modules:
 * - plot-constants: unit scales, palettes, axis/layout presets
 * - plot-formatters: frequency/time/band label formatting, color assignment
 * - plot-helpers: shared math helpers (y-range, AB mag, min/max)
 * - plot-lc: light curve figure builder
 * - plot-sed: spectrum/SED figure builder
 * - plot-skymap: skymap figure builder
 */

// Re-export everything that was previously exported from this file,
// so existing imports continue to work.

export type { PlotlyFigure } from "./plot-constants";
export {
  COMP_ORDER,
  FLUX_SCALES_CGS,
  FREQ_DISP_SCALES,
  TIME_SCALES_S,
} from "./plot-constants";

export {
  formatBandLabel,
  formatFreqLabel,
  formatTimeLabel,
} from "./plot-formatters";

export { type LcDisplayOptions, buildLcFigure } from "./plot-lc";
export { type SedDisplayOptions, buildSedFigure } from "./plot-sed";
export { type SkymapOptions, buildSkymapFigure, decodeSkymapFrame } from "./plot-skymap";
