import type { ReactNode } from "react";
import { FREQ_UNIT_OPTIONS } from "../lib/constants";
import type { ModeLogSliderSpec } from "../lib/types";
import { SliderField } from "./SliderField";

type Props = {
  sedTimesDraft: string;
  setSedTimesDraft: (value: string) => void;
  commitSedTimes: () => void;
  renderDistanceObserverControls: () => ReactNode;
  renderModeLogSliderRow: (specs: ModeLogSliderSpec[], className?: string) => ReactNode;
  sedNuMin: number;
  setSedNuMin: (value: number) => void;
  sedNuMax: number;
  setSedNuMax: (value: number) => void;
  sedNumNu: number;
  setSedNumNu: (value: number) => void;
  sedNuFNu: boolean;
  setSedNuFNu: (value: boolean) => void;
  fluxIsAbMag: boolean;
  sedFreqUnit: string;
  setSedFreqUnit: (value: string) => void;
  renderStringOptions: (options: readonly string[]) => ReactNode;
  renderFluxUnitControl: () => ReactNode;
  renderSharedControls: (options: { hideDistanceObserver?: boolean; hidePlotUnits?: boolean }) => ReactNode;
};

export function ModePanelSpectrum({
  sedTimesDraft,
  setSedTimesDraft,
  commitSedTimes,
  renderDistanceObserverControls,
  renderModeLogSliderRow,
  sedNuMin,
  setSedNuMin,
  sedNuMax,
  setSedNuMax,
  sedNumNu,
  setSedNumNu,
  sedNuFNu,
  setSedNuFNu,
  fluxIsAbMag,
  sedFreqUnit,
  setSedFreqUnit,
  renderStringOptions,
  renderFluxUnitControl,
  renderSharedControls,
}: Props) {
  return (
    <>
      <div className="sb-stack">
        <label className="sb-field">
          <span className="sb-label">
            t<sub>obs</sub> (s)
          </span>
          <input
            value={sedTimesDraft}
            onChange={(e) => setSedTimesDraft(e.target.value)}
            onBlur={commitSedTimes}
            placeholder="e.g. 1e3, 1e4, 1e5"
          />
        </label>
        {renderDistanceObserverControls()}
      </div>

      <div className="sb-stack">
        {renderModeLogSliderRow([
          {
            key: "sed-nu-min",
            label: (
              <>
                log10(ν<sub>min</sub> (Hz))
              </>
            ),
            minExp: 6,
            maxExp: 20,
            step: 0.05,
            value: sedNuMin,
            defaultExp: 8,
            onChange: setSedNuMin,
          },
          {
            key: "sed-nu-max",
            label: (
              <>
                log10(ν<sub>max</sub> (Hz))
              </>
            ),
            minExp: 10,
            maxExp: 35,
            step: 0.05,
            value: sedNuMax,
            defaultExp: 20,
            onChange: setSedNuMax,
          },
        ])}
        <div className="sb-row-2 sb-row-fit-right">
          <SliderField
            label={<>ν points</>}
            min={50}
            max={500}
            step={10}
            decimals={0}
            value={sedNumNu}
            onChange={(v) => setSedNumNu(Math.round(v))}
          />
          <label className="sb-checkbox-inline sb-nufnu-inline">
            <input type="checkbox" checked={sedNuFNu} disabled={fluxIsAbMag} onChange={(e) => setSedNuFNu(e.target.checked)} />
            ν F<sub>ν</sub>
          </label>
        </div>
        <div className="sb-row-2">
          <label className="sb-field">
            <span className="sb-label">ν unit</span>
            <select value={sedFreqUnit} onChange={(e) => setSedFreqUnit(e.target.value)}>
              {renderStringOptions(FREQ_UNIT_OPTIONS)}
            </select>
          </label>
          {renderFluxUnitControl()}
        </div>
        {fluxIsAbMag ? <p className="sb-muted">ν F<sub>ν</sub> is disabled for AB mag.</p> : null}
      </div>

      {renderSharedControls({ hideDistanceObserver: true, hidePlotUnits: true })}
    </>
  );
}
