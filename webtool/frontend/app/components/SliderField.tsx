import {
  useCallback,
  useEffect,
  useMemo,
  useRef,
  useState,
  type CSSProperties,
  type KeyboardEvent as ReactKeyboardEvent,
  type PointerEvent as ReactPointerEvent,
  type ReactNode,
} from "react";
import { SLIDER_COMMIT_INTERVAL_MS } from "../lib/constants";
import { clearTimeoutRef } from "../lib/utils/async";
import { safeLog10 } from "../lib/utils/math";
import { sliderFillStyle, stepRangeValue } from "../lib/utils/slider";

type SliderDraftHandlers = {
  draft: number;
  fillStyle: CSSProperties;
  handleInput: (raw: string) => void;
  handleKeyDown: (event: ReactKeyboardEvent<HTMLInputElement>) => void;
  handleStart: (event: ReactPointerEvent<HTMLInputElement>) => void;
  handleEnd: () => void;
};

function useSliderDraftState({
  value,
  min,
  max,
  step,
  onCommit,
}: {
  value: number;
  min: number;
  max: number;
  step: number;
  onCommit: (next: number) => void;
}): SliderDraftHandlers {
  const [draft, setDraft] = useState(value);
  const draggingRef = useRef(false);
  const timerRef = useRef<number | null>(null);
  const pendingRef = useRef<number | null>(null);
  const fillStyle = useMemo(() => sliderFillStyle(draft, min, max), [draft, min, max]);

  useEffect(() => {
    if (!draggingRef.current) {
      setDraft(value);
    }
  }, [value]);

  const flushCommit = useCallback(
    (next: number) => {
      clearTimeoutRef(timerRef);
      pendingRef.current = null;
      onCommit(next);
    },
    [onCommit],
  );

  const scheduleCommit = useCallback(
    (next: number) => {
      pendingRef.current = next;
      if (timerRef.current !== null) return;
      timerRef.current = window.setTimeout(() => {
        timerRef.current = null;
        const pending = pendingRef.current;
        pendingRef.current = null;
        if (pending !== null) {
          onCommit(pending);
        }
      }, SLIDER_COMMIT_INTERVAL_MS);
    },
    [onCommit],
  );

  useEffect(() => {
    return () => clearTimeoutRef(timerRef);
  }, []);

  const handleInput = useCallback(
    (raw: string) => {
      const next = Number(raw);
      setDraft(next);
      scheduleCommit(next);
    },
    [scheduleCommit],
  );

  const handleKeyDown = useCallback(
    (event: ReactKeyboardEvent<HTMLInputElement>) => {
      stepRangeValue(event, draft, min, max, step, (next) => {
        setDraft(next);
        flushCommit(next);
      });
    },
    [draft, flushCommit, max, min, step],
  );

  const handleStart = useCallback((event: ReactPointerEvent<HTMLInputElement>) => {
    draggingRef.current = true;
    event.currentTarget.focus();
  }, []);

  const handleEnd = useCallback(() => {
    draggingRef.current = false;
    flushCommit(draft);
  }, [draft, flushCommit]);

  return {
    draft,
    fillStyle,
    handleInput,
    handleKeyDown,
    handleStart,
    handleEnd,
  };
}

export type LogSliderProps = {
  label: ReactNode;
  minExp: number;
  maxExp: number;
  step: number;
  value: number;
  defaultExp: number;
  precision?: number;
  onChange: (next: number) => void;
};

export function LogSliderField({ label, minExp, maxExp, step, value, defaultExp, precision = 1, onChange }: LogSliderProps) {
  const exp = safeLog10(value, defaultExp);
  const {
    draft: draftExp,
    fillStyle,
    handleInput,
    handleKeyDown,
    handleStart,
    handleEnd,
  } = useSliderDraftState({
    value: exp,
    min: minExp,
    max: maxExp,
    step,
    onCommit: (nextExp) => onChange(10 ** nextExp),
  });

  return (
    <label className="sb-field sb-slider">
      <span className="sb-label">{label}</span>
      <div className="sb-slider-track">
        <input
          type="range"
          min={minExp}
          max={maxExp}
          step={step}
          value={draftExp}
          style={fillStyle}
          onPointerDown={handleStart}
          onPointerUp={handleEnd}
          onPointerCancel={handleEnd}
          onInput={(e) => handleInput((e.target as HTMLInputElement).value)}
          onChange={(e) => handleInput((e.target as HTMLInputElement).value)}
          onKeyDown={handleKeyDown}
          aria-keyshortcuts="ArrowLeft ArrowRight ArrowUp ArrowDown Home End"
        />
        <span className="sb-value">{draftExp.toFixed(precision)}</span>
      </div>
    </label>
  );
}

export type SliderProps = {
  label: ReactNode;
  min: number;
  max: number;
  step: number;
  value: number;
  decimals?: number;
  formatValue?: (value: number) => string;
  onChange: (next: number) => void;
};

export function SliderField({ label, min, max, step, value, decimals = 2, formatValue, onChange }: SliderProps) {
  const { draft, fillStyle, handleInput, handleKeyDown, handleStart, handleEnd } = useSliderDraftState({
    value,
    min,
    max,
    step,
    onCommit: onChange,
  });

  return (
    <label className="sb-field sb-slider">
      <span className="sb-label">{label}</span>
      <div className="sb-slider-track">
        <input
          type="range"
          min={min}
          max={max}
          step={step}
          value={draft}
          style={fillStyle}
          onPointerDown={handleStart}
          onPointerUp={handleEnd}
          onPointerCancel={handleEnd}
          onInput={(e) => handleInput((e.target as HTMLInputElement).value)}
          onChange={(e) => handleInput((e.target as HTMLInputElement).value)}
          onKeyDown={handleKeyDown}
          aria-keyshortcuts="ArrowLeft ArrowRight ArrowUp ArrowDown Home End"
        />
        <span className="sb-value">{formatValue ? formatValue(draft) : draft.toFixed(decimals)}</span>
      </div>
    </label>
  );
}
