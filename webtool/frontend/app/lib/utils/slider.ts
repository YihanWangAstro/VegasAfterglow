import type { CSSProperties, KeyboardEvent as ReactKeyboardEvent } from "react";
import { clampRangeValue } from "./math";

type RangeKeyAction = "dec" | "inc" | "min" | "max" | null;

function rangeKeyAction(event: ReactKeyboardEvent<HTMLInputElement>): RangeKeyAction {
  const key = event.key;
  const code = event.code;
  if (key === "ArrowLeft" || key === "ArrowDown") return "dec";
  if (key === "ArrowRight" || key === "ArrowUp") return "inc";
  if (key === "Home") return "min";
  if (key === "End") return "max";
  if (code === "ArrowLeft" || code === "ArrowDown") return "dec";
  if (code === "ArrowRight" || code === "ArrowUp") return "inc";
  if (code === "Home") return "min";
  if (code === "End") return "max";
  return null;
}

export function stepRangeValue(
  event: ReactKeyboardEvent<HTMLInputElement>,
  current: number,
  min: number,
  max: number,
  step: number,
  apply: (next: number) => void,
): void {
  let next: number | null = null;
  switch (rangeKeyAction(event)) {
    case "dec":
      next = current - step;
      break;
    case "inc":
      next = current + step;
      break;
    case "min":
      next = min;
      break;
    case "max":
      next = max;
      break;
    default:
      return;
  }

  event.preventDefault();
  const snapped = clampRangeValue(Number(next.toFixed(10)), min, max);
  apply(snapped);
}

export function sliderFillStyle(value: number, min: number, max: number): CSSProperties {
  const span = max - min;
  const raw = span > 0 ? ((value - min) / span) * 100 : 0;
  const bounded = Math.max(0, Math.min(100, raw));
  return { "--fill-percent": `${bounded}%` } as CSSProperties;
}
