import { AXIS_EPS } from "../constants";
import type { AxisName, AxisRange, Mode } from "../types";

export function parseAxisRange(value: unknown): AxisRange | null {
  if (!Array.isArray(value) || value.length !== 2) return null;
  const lo = Number(value[0]);
  const hi = Number(value[1]);
  if (!Number.isFinite(lo) || !Number.isFinite(hi)) return null;
  return [lo, hi];
}

export function parseRelayoutRange(event: Record<string, unknown>, axis: AxisName): AxisRange | null {
  const direct = parseAxisRange(event[`${axis}.range`]);
  if (direct) return direct;
  const lo = Number(event[`${axis}.range[0]`]);
  const hi = Number(event[`${axis}.range[1]`]);
  if (!Number.isFinite(lo) || !Number.isFinite(hi)) return null;
  return [lo, hi];
}

export function rangesEqual(a: AxisRange | null, b: AxisRange | null): boolean {
  if (a === b) return true;
  if (!a || !b) return false;
  return Math.abs(a[0] - b[0]) <= AXIS_EPS && Math.abs(a[1] - b[1]) <= AXIS_EPS;
}

export function axisTitleText(axisObj: Record<string, unknown>): string {
  const title = axisObj.title;
  if (typeof title === "string") return title;
  if (title && typeof title === "object") {
    const text = (title as { text?: unknown }).text;
    if (typeof text === "string") return text;
  }
  return "";
}

export function axisSignature(layout: Record<string, unknown>, axis: AxisName): string | null {
  const axisObj = layout[axis];
  if (!axisObj || typeof axisObj !== "object") return null;
  const axisRecord = axisObj as Record<string, unknown>;
  const axisType = typeof axisRecord.type === "string" ? axisRecord.type : "";
  const title = axisTitleText(axisRecord);
  return `${axisType}|${title}`;
}

function asAxisRecord(layout: Record<string, unknown>, axis: "xaxis" | "yaxis"): Record<string, unknown> {
  const axisObj = layout[axis];
  if (!axisObj || typeof axisObj !== "object") return {};
  return { ...(axisObj as Record<string, unknown>) };
}

export function shouldPreserveAxisRangesOnAutorange(mode: Mode): boolean {
  return mode === "skymap";
}

export function applySkymapAxisPolicy(layout: Record<string, unknown>): void {
  const xAxisObj = asAxisRecord(layout, "xaxis");
  const yAxisObj = asAxisRecord(layout, "yaxis");

  const xRange = parseAxisRange((xAxisObj as { range?: unknown }).range);
  const yRange = parseAxisRange((yAxisObj as { range?: unknown }).range);
  if (xRange && yRange) {
    const xCenter = 0.5 * (xRange[0] + xRange[1]);
    const yCenter = 0.5 * (yRange[0] + yRange[1]);
    const span = Math.max(Math.abs(xRange[1] - xRange[0]), Math.abs(yRange[1] - yRange[0]));
    const half = span * 0.5;
    xAxisObj.range = [xCenter - half, xCenter + half];
    yAxisObj.range = [yCenter - half, yCenter + half];
  }

  xAxisObj.autorange = false;
  yAxisObj.autorange = false;
  xAxisObj.constrain = "domain";
  yAxisObj.constrain = "domain";
  yAxisObj.scaleanchor = "x";
  yAxisObj.scaleratio = 1;

  layout.xaxis = xAxisObj;
  layout.yaxis = yAxisObj;
}

export function legendFontSizeForWidth(plotWidthPx: number): number {
  if (plotWidthPx <= 300) return 6;
  if (plotWidthPx <= 380) return 7;
  if (plotWidthPx <= 480) return 8;
  if (plotWidthPx <= 640) return 9;
  if (plotWidthPx <= 860) return 10;
  return 11;
}
