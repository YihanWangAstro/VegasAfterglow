"use client";

import dynamic from "next/dynamic";
import { memo, useCallback, useEffect, useRef, useState } from "react";
import { PLOT_CONFIG } from "../lib/constants";
import type { Mode } from "../lib/types";
import type { PlotlyFigure } from "../lib/utils/plot-builders";

const Plot = dynamic(() => import("./PlotlyCustom").then((m) => ({ default: m.default })), { ssr: false });

export const PlotFigure = memo(function PlotFigure({
  figure,
  mode,
  onRelayout,
}: {
  figure: PlotlyFigure;
  mode: Mode;
  onRelayout: (event: Record<string, unknown>) => void;
}) {
  const graphRef = useRef<any>(null);
  const [graphRevision, setGraphRevision] = useState(0);

  const bindGraphRef = useCallback((_nextFigure: unknown, graphDiv: unknown) => {
    const nextGraph = graphDiv as any;
    if (graphRef.current !== nextGraph) {
      graphRef.current = nextGraph;
      setGraphRevision((v) => v + 1);
    }
  }, []);

  useEffect(() => {
    const hasFrames = Array.isArray(figure.frames) && figure.frames.length > 1;
    if (mode !== "skymap" || !hasFrames) return;

    const graphDiv = graphRef.current as any;
    if (!graphDiv?.animate) return;
    const frameNames = (figure.frames ?? [])
      .map((frame) => (frame && typeof frame === "object" ? (frame as { name?: unknown }).name : undefined))
      .filter((name): name is string => typeof name === "string" && name.length > 0);
    if (frameNames.length === 0) return;

    let cancelled = false;
    let starter: number | null = null;
    let frameIndex = 0;
    const playArgs: Record<string, unknown> = {
      frame: { duration: 300, redraw: true },
      transition: { duration: 0 },
      mode: "immediate",
    };
    const stopArgs = { frame: { duration: 0, redraw: false }, mode: "immediate" };
    const stopAnimationSafely = () => {
      try {
        const maybePromise = graphDiv.animate?.([null], stopArgs);
        if (maybePromise && typeof maybePromise.then === "function") {
          void maybePromise.catch(() => undefined);
        }
      } catch {
        // Ignore transient Plotly stop failures.
      }
    };

    const playSequence = async () => {
      if (cancelled) return;

      const nextFrame = frameNames[frameIndex % frameNames.length];
      frameIndex += 1;
      try {
        await graphDiv.animate([nextFrame], playArgs);
        if (cancelled) return;
        starter = window.setTimeout(() => {
          void playSequence();
        }, 0);
      } catch {
        // Transient Plotly frame registration race; retry shortly.
        if (cancelled) return;
        starter = window.setTimeout(() => {
          void playSequence();
        }, 60);
      }
    };

    starter = window.setTimeout(() => {
      void playSequence();
    }, 80);

    return () => {
      cancelled = true;
      if (starter !== null) {
        window.clearTimeout(starter);
      }
      stopAnimationSafely();
    };
  }, [figure.frames, graphRevision, mode]);

  return (
    <div className={`plot-wrap${mode === "skymap" ? " plot-wrap-square" : ""}`}>
      <Plot
        data={figure.data as any}
        layout={(figure.layout ?? {}) as any}
        frames={(figure.frames ?? []) as any}
        config={PLOT_CONFIG as any}
        onRelayout={onRelayout as any}
        onInitialized={bindGraphRef}
        onUpdate={bindGraphRef}
        style={{ width: "100%", height: "100%" }}
        useResizeHandler
      />
    </div>
  );
});
