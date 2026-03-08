import { useCallback, useEffect, useRef, useState, type Dispatch, type SetStateAction, type TransitionStartFunction } from "react";
import { AUTO_RUN_DEBOUNCE_IDLE_MS, COLD_START_HINT_MS } from "../lib/constants";
import type { ComputationSpec, Mode, RunResponse } from "../lib/types";
import { clearTimeoutRef, isAbortLikeError } from "../lib/utils/async";
import { isEmptyPrimaryInputSpec } from "../lib/utils/snapshot";

type Args = {
  bootReady: boolean;
  computationSpec: ComputationSpec;
  compareSpec: ComputationSpec | null;
  compareEnabled: boolean;
  mode: Mode;
  sliderInteracting: boolean;
  apiCandidates: string[];
  postFigureRequest: (spec: ComputationSpec, controller: AbortController) => Promise<RunResponse>;
  startFigureTransition: TransitionStartFunction;
  setResult: Dispatch<SetStateAction<RunResponse | null>>;
  setResultMode: Dispatch<SetStateAction<Mode | null>>;
  setError: Dispatch<SetStateAction<string>>;
  clearCompareOverlay: (clearError?: boolean) => void;
  setCompareRunning: Dispatch<SetStateAction<boolean>>;
  setCompareError: Dispatch<SetStateAction<string>>;
  setCompareResult: Dispatch<SetStateAction<RunResponse | null>>;
  setCompareResultMode: Dispatch<SetStateAction<Mode | null>>;
};

export function useComputationEffects({
  bootReady,
  computationSpec,
  compareSpec,
  compareEnabled,
  mode,
  sliderInteracting,
  apiCandidates,
  postFigureRequest,
  startFigureTransition,
  setResult,
  setResultMode,
  setError,
  clearCompareOverlay,
  setCompareRunning,
  setCompareError,
  setCompareResult,
  setCompareResultMode,
}: Args) {
  const [showRunning, setShowRunning] = useState(false);
  const [showColdStartHint, setShowColdStartHint] = useState(false);

  const requestSeqRef = useRef(0);
  const runTimerRef = useRef<number | null>(null);
  const activeRequestRef = useRef<AbortController | null>(null);
  const coldStartTimerRef = useRef<number | null>(null);
  const activeRunRequestSeqRef = useRef<number | null>(null);

  const compareRequestSeqRef = useRef(0);
  const compareTimerRef = useRef<number | null>(null);
  const compareActiveRequestRef = useRef<AbortController | null>(null);
  const requestInFlightRef = useRef(false);
  const pendingSpecRef = useRef<ComputationSpec | null>(null);
  const [running, setRunning] = useState(false);

  const runComputation = useCallback(
    async (spec: ComputationSpec) => {
      const requestMode = spec.endpoint;
      const controller = new AbortController();
      activeRequestRef.current = controller;
      const requestSeq = ++requestSeqRef.current;
      activeRunRequestSeqRef.current = requestSeq;

      setError("");
      setShowColdStartHint(false);
      clearTimeoutRef(coldStartTimerRef);
      coldStartTimerRef.current = window.setTimeout(() => {
        if (activeRunRequestSeqRef.current === requestSeq) {
          setShowColdStartHint(true);
        }
      }, COLD_START_HINT_MS);

      try {
        const data = await postFigureRequest(spec, controller);
        if (requestSeq !== requestSeqRef.current) return;

        startFigureTransition(() => {
          setResult(data);
          setResultMode(requestMode);
        });
      } catch (err) {
        if (isAbortLikeError(err, controller)) return;
        if (requestSeq !== requestSeqRef.current) return;

        const message = err instanceof Error ? err.message : "Unknown error";
        const lower = message.toLowerCase();
        const isNetworkFailure = lower.includes("load failed") || lower.includes("failed to fetch");
        setError(isNetworkFailure ? `Cannot reach backend. Checked ${apiCandidates.join(" or ")}.` : message);
      } finally {
        if (activeRequestRef.current === controller) {
          activeRequestRef.current = null;
        }
        if (activeRunRequestSeqRef.current === requestSeq) {
          activeRunRequestSeqRef.current = null;
          clearTimeoutRef(coldStartTimerRef);
          setShowColdStartHint(false);
        }
      }
    },
    [apiCandidates, postFigureRequest, setError, setResult, setResultMode, startFigureTransition],
  );

  const runCompareComputation = useCallback(
    async (spec: ComputationSpec) => {
      compareActiveRequestRef.current?.abort();
      const controller = new AbortController();
      compareActiveRequestRef.current = controller;
      const requestSeq = ++compareRequestSeqRef.current;
      setCompareRunning(true);
      setCompareError("");

      try {
        const data = await postFigureRequest(spec, controller);
        if (requestSeq !== compareRequestSeqRef.current) return;
        setCompareResult(data);
        setCompareResultMode(spec.endpoint);
      } catch (err) {
        if (isAbortLikeError(err, controller)) return;
        if (requestSeq !== compareRequestSeqRef.current) return;
        const message = err instanceof Error ? err.message : "Compare request failed";
        setCompareError(message);
        setCompareResult(null);
        setCompareResultMode(null);
      } finally {
        if (compareActiveRequestRef.current === controller) {
          compareActiveRequestRef.current = null;
        }
        if (requestSeq === compareRequestSeqRef.current) {
          setCompareRunning(false);
        }
      }
    },
    [postFigureRequest, setCompareError, setCompareResult, setCompareResultMode, setCompareRunning],
  );

  const cancelCompareRun = useCallback(() => {
    compareActiveRequestRef.current?.abort();
    setCompareRunning(false);
  }, [setCompareRunning]);

  const flushQueuedComputation = useCallback(async () => {
    if (requestInFlightRef.current) return;
    if (!pendingSpecRef.current) return;

    requestInFlightRef.current = true;
    setRunning(true);
    try {
      let nextSpec: ComputationSpec | null = pendingSpecRef.current;
      while (nextSpec) {
        pendingSpecRef.current = null;
        await runComputation(nextSpec);
        nextSpec = pendingSpecRef.current;
      }
    } finally {
      requestInFlightRef.current = false;
      setRunning(false);
    }
  }, [runComputation]);

  const enqueueComputation = useCallback((spec: ComputationSpec) => {
    pendingSpecRef.current = spec;
  }, []);

  const clearQueuedComputation = useCallback(() => {
    pendingSpecRef.current = null;
  }, []);

  useEffect(() => {
    if (!bootReady) return;
    if (isEmptyPrimaryInputSpec(computationSpec)) {
      clearQueuedComputation();
      clearTimeoutRef(runTimerRef);
      activeRequestRef.current?.abort();
      setError("");
      startFigureTransition(() => {
        setResult(null);
        setResultMode(null);
      });
      return;
    }
    enqueueComputation(computationSpec);
    clearTimeoutRef(runTimerRef);
    if (sliderInteracting) {
      // Fire immediately — SLIDER_COMMIT_INTERVAL_MS already throttles commit rate
      void flushQueuedComputation();
    } else {
      runTimerRef.current = window.setTimeout(() => {
        runTimerRef.current = null;
        void flushQueuedComputation();
      }, AUTO_RUN_DEBOUNCE_IDLE_MS);
    }
    return () => clearTimeoutRef(runTimerRef);
  }, [
    bootReady,
    clearQueuedComputation,
    computationSpec,
    enqueueComputation,
    flushQueuedComputation,
    setError,
    setResult,
    setResultMode,
    sliderInteracting,
    startFigureTransition,
  ]);

  useEffect(() => {
    clearTimeoutRef(compareTimerRef);

    if (!bootReady || !compareSpec) {
      cancelCompareRun();
      if (!compareEnabled || mode === "skymap") {
        clearCompareOverlay(true);
      }
      return;
    }

    if (isEmptyPrimaryInputSpec(compareSpec)) {
      cancelCompareRun();
      clearCompareOverlay(true);
      return;
    }

    const debounceMs = sliderInteracting ? 0 : AUTO_RUN_DEBOUNCE_IDLE_MS;
    compareTimerRef.current = window.setTimeout(() => {
      compareTimerRef.current = null;
      void runCompareComputation(compareSpec);
    }, debounceMs);

    return () => clearTimeoutRef(compareTimerRef);
  }, [
    bootReady,
    cancelCompareRun,
    clearCompareOverlay,
    compareEnabled,
    compareSpec,
    mode,
    runCompareComputation,
    sliderInteracting,
  ]);

  useEffect(() => {
    return () => {
      clearTimeoutRef(runTimerRef);
      clearTimeoutRef(coldStartTimerRef);
      clearQueuedComputation();
      activeRequestRef.current?.abort();
      clearTimeoutRef(compareTimerRef);
      compareActiveRequestRef.current?.abort();
    };
  }, [clearQueuedComputation]);

  useEffect(() => {
    if (!running) {
      setShowRunning(false);
      return;
    }
    const timer = window.setTimeout(() => {
      setShowRunning(true);
    }, 120);
    return () => window.clearTimeout(timer);
  }, [running]);

  return { showRunning, showColdStartHint };
}
