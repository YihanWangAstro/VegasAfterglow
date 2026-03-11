"use client";

import dynamic from "next/dynamic";
import { memo, useCallback, useRef } from "react";
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

  const bindGraphRef = useCallback((_nextFigure: unknown, graphDiv: unknown) => {
    graphRef.current = graphDiv as any;
  }, []);

  return (
    <div className={`plot-wrap${mode === "skymap" ? " plot-wrap-square" : ""}`}>
      <Plot
        data={figure.data as any}
        layout={(figure.layout ?? {}) as any}
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
