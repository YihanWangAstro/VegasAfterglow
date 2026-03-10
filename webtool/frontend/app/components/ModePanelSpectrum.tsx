import { useCallback, type ReactNode } from "react";
import { FREQ_UNIT_OPTIONS } from "../lib/constants";
import { useParameters } from "../lib/ParameterContext";
import type { ModeLogSliderSpec } from "../lib/types";
import { SliderField } from "./SliderField";

function renderStringOptions(options: readonly string[]): ReactNode[] {
  return options.map((item) => (
    <option key={item} value={item}>
      {item}
    </option>
  ));
}

type Props = {
  renderDistanceObserverControls: () => ReactNode;
  renderModeLogSliderRow: (specs: ModeLogSliderSpec[], className?: string) => ReactNode;
  renderFluxUnitControl: () => ReactNode;
  renderSharedControls: (options: { hideDistanceObserver?: boolean; hidePlotUnits?: boolean }) => ReactNode;
};

export function ModePanelSpectrum({
  renderDistanceObserverControls,
  renderModeLogSliderRow,
  renderFluxUnitControl,
  renderSharedControls,
}: Props) {
  const { state, actions } = useParameters();
  const { sedTimesDraft, sedNuMin, sedNuMax, sedNumNu, sedNuFNu, sedFreqUnit, shared } = state;
  const fluxIsAbMag = shared.flux_unit === "AB mag";

  const commitSedTimes = useCallback(() => actions.setSedTimes(sedTimesDraft), [sedTimesDraft, actions.setSedTimes]);

  return (
    <>
      <div className="sb-stack">
        <label className="sb-field">
          <span className="sb-label">
            t<sub>obs</sub> (s)
          </span>
          <input
            value={sedTimesDraft}
            onChange={(e) => actions.setSedTimesDraft(e.target.value)}
            onBlur={commitSedTimes}
            onKeyDown={(e) => { if (e.key === "Enter") commitSedTimes(); }}
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
            onChange: actions.setSedNuMin,
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
            onChange: actions.setSedNuMax,
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
            onChange={(v) => actions.setSedNumNu(Math.round(v))}
          />
          <label className="sb-checkbox-inline sb-nufnu-inline">
            <input type="checkbox" checked={sedNuFNu} disabled={fluxIsAbMag} onChange={(e) => actions.setSedNuFNu(e.target.checked)} />
            ν F<sub>ν</sub>
          </label>
        </div>
        <div className="sb-row-2">
          <label className="sb-field">
            <span className="sb-label">ν unit</span>
            <select value={sedFreqUnit} onChange={(e) => actions.setSedFreqUnit(e.target.value)}>
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
