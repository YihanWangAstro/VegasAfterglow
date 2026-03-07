import type { MutableRefObject } from "react";

export function clearTimeoutRef(timerRef: MutableRefObject<number | null>): void {
  if (timerRef.current === null) return;
  window.clearTimeout(timerRef.current);
  timerRef.current = null;
}

export function isAbortLikeError(err: unknown, controller: AbortController): boolean {
  return controller.signal.aborted || (err instanceof DOMException && err.name === "AbortError");
}
