import type { ReactNode } from "react";
import { FREQ_HELP_TEXT } from "../lib/constants";
import type { ModeLogSliderSpec } from "../lib/types";
import { HelpHint } from "./HelpHint";
import { SliderField } from "./SliderField";

type Props = {
  lcFreqDraft: string;
  setLcFreqDraft: (value: string) => void;
  commitLcFreq: () => void;
  renderDistanceObserverControls: () => ReactNode;
  renderModeLogSliderRow: (specs: ModeLogSliderSpec[], className?: string) => ReactNode;
  lcTMin: number;
  setLcTMin: (value: number) => void;
  lcTMax: number;
  setLcTMax: (value: number) => void;
  numT: number;
  setNumT: (value: number) => void;
  renderPlotUnitControls: (timeDisabled: boolean, showTime?: boolean) => ReactNode;
  renderSharedControls: (options: { hideDistanceObserver?: boolean; hidePlotUnits?: boolean }) => ReactNode;
};

export function ModePanelLightcurve({
  lcFreqDraft,
  setLcFreqDraft,
  commitLcFreq,
  renderDistanceObserverControls,
  renderModeLogSliderRow,
  lcTMin,
  setLcTMin,
  lcTMax,
  setLcTMax,
  numT,
  setNumT,
  renderPlotUnitControls,
  renderSharedControls,
}: Props) {
  return (
    <>
      <div className="sb-stack">
        <label className="sb-field">
          <span className="sb-label sb-label-with-help">
            <span>Frequencies</span>
            <HelpHint text={FREQ_HELP_TEXT} ariaLabel="Frequency input help" />
          </span>
          <input
            value={lcFreqDraft}
            onChange={(e) => setLcFreqDraft(e.target.value)}
            onBlur={commitLcFreq}
            onKeyDown={(e) => { if (e.key === "Enter") commitLcFreq(); }}
            placeholder="e.g. 1e9, R, 1keV, XRT, [0.3keV,10keV]"
          />
        </label>
        {renderDistanceObserverControls()}
      </div>

      <div className="sb-stack">
        {renderModeLogSliderRow([
          {
            key: "lc-t-min",
            label: (
              <>
                log10(t<sub>min</sub> (s))
              </>
            ),
            minExp: -1,
            maxExp: 6,
            step: 0.05,
            value: lcTMin,
            defaultExp: 0,
            onChange: setLcTMin,
          },
          {
            key: "lc-t-max",
            label: (
              <>
                log10(t<sub>max</sub> (s))
              </>
            ),
            minExp: 3,
            maxExp: 10,
            step: 0.05,
            value: lcTMax,
            defaultExp: 8,
            onChange: setLcTMax,
          },
        ])}
        <SliderField
          label={<>t points</>}
          min={50}
          max={300}
          step={10}
          decimals={0}
          value={numT}
          onChange={(v) => setNumT(Math.round(v))}
        />
        {renderPlotUnitControls(false)}
      </div>

      {renderSharedControls({ hideDistanceObserver: true, hidePlotUnits: true })}
    </>
  );
}
