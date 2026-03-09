import type { ReactNode } from "react";
import { FREQ_HELP_TEXT, SKY_FOV_UNIT_OPTIONS, SKY_INTENSITY_UNIT_OPTIONS } from "../lib/constants";
import { HelpHint } from "./HelpHint";
import { LogSliderField } from "./SliderField";

type Props = {
  skyNpixel: number;
  setSkyNpixel: (value: number) => void;
  skyTObs: number;
  setSkyTObs: (value: number) => void;
  skyNuInputDraft: string;
  setSkyNuInputDraft: (value: string) => void;
  commitSkyNuInput: () => void;
  skyFov: number;
  setSkyFov: (value: number) => void;
  skyFovUnit: string;
  setSkyFovUnit: (value: string) => void;
  skyIntensityUnit: string;
  setSkyIntensityUnit: (value: string) => void;
  skyPixelOptions: readonly number[];
  renderDistanceObserverControls: () => ReactNode;
  renderSharedControls: (options?: { hideDistanceObserver?: boolean; hidePlotUnits?: boolean }) => ReactNode;
};

export function ModePanelSkymap({
  skyNpixel,
  setSkyNpixel,
  skyTObs,
  setSkyTObs,
  skyNuInputDraft,
  setSkyNuInputDraft,
  commitSkyNuInput,
  skyFov,
  setSkyFov,
  skyFovUnit,
  setSkyFovUnit,
  skyIntensityUnit,
  setSkyIntensityUnit,
  skyPixelOptions,
  renderDistanceObserverControls,
  renderSharedControls,
}: Props) {
  return (
    <>
      <div className="sb-stack">
        <label className="sb-field">
          <span className="sb-label sb-label-with-help">
            <span>
              ν<sub>obs</sub>
            </span>
            <HelpHint text={FREQ_HELP_TEXT} ariaLabel="Sky frequency input help" />
          </span>
          <input
            value={skyNuInputDraft}
            onChange={(e) => setSkyNuInputDraft(e.target.value)}
            onBlur={commitSkyNuInput}
            onKeyDown={(e) => { if (e.key === "Enter") commitSkyNuInput(); }}
            placeholder="e.g. 1GHz"
          />
        </label>
        {renderDistanceObserverControls()}
        <LogSliderField
          label={
            <>
              log10(t<sub>obs</sub> (s))
            </>
          }
          minExp={2}
          maxExp={9}
          step={0.05}
          value={skyTObs}
          defaultExp={6}
          onChange={setSkyTObs}
        />
      </div>

      <div className="sb-stack">

        <div className="sb-row-2 sb-row-gap-top">
          <LogSliderField
            label="log10(FOV (μas))"
            minExp={1}
            maxExp={5}
            step={0.05}
            value={skyFov}
            defaultExp={Math.log10(500)}
            onChange={setSkyFov}
          />
          <label className="sb-field">
            <span className="sb-label">Pixels</span>
            <select value={String(skyNpixel)} onChange={(e) => setSkyNpixel(Number(e.target.value))}>
              {skyPixelOptions.map((v) => (
                <option key={v} value={v}>
                  {v}
                </option>
              ))}
            </select>
          </label>
        </div>

        <div className="sb-row-2 sb-row-gap-top">
          <label className="sb-field">
            <span className="sb-label">FOV unit</span>
            <select value={skyFovUnit} onChange={(e) => setSkyFovUnit(e.target.value)}>
              {SKY_FOV_UNIT_OPTIONS.map((u) => (
                <option key={u} value={u}>
                  {u}
                </option>
              ))}
            </select>
          </label>
          <label className="sb-field">
            <span className="sb-label">Intensity unit</span>
            <select value={skyIntensityUnit} onChange={(e) => setSkyIntensityUnit(e.target.value)}>
              {SKY_INTENSITY_UNIT_OPTIONS.map((u) => (
                <option key={u} value={u}>
                  {u}
                </option>
              ))}
            </select>
          </label>
        </div>
      </div>

      {renderSharedControls({ hideDistanceObserver: true, hidePlotUnits: true })}
    </>
  );
}
