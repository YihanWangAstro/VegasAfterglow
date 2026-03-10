/**
 * Shared helper functions used by multiple plot builders (LC, SED).
 */

import {
  AXIS_COMMON,
  axisTitle,
} from "./plot-constants";

/** Convert CGS flux to AB magnitude. */
export function cgsToAbMag(f: number): number {
  if (f <= 0) return NaN;
  return -2.5 * Math.log10(f) - 48.6;
}

/**
 * Collect all positive values from 2D arrays, compute smart log y-axis limits.
 * Returns null if no positive values found.
 */
export function smartYlim(arrays: number[][]): [number, number] | null {
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
export function magYRange(components: Record<string, number[][]>): [number, number] | "reversed" {
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
export function buildYAxisConfig(
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
export function arrayMinMax(arr: number[]): [number, number] {
  let lo = Infinity;
  let hi = -Infinity;
  for (const v of arr) {
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  return [lo, hi];
}

/** Append values to target array (avoids stack overflow from push(...spread)). */
export function pushAll(target: number[], source: number[]) {
  for (const v of source) target.push(v);
}
