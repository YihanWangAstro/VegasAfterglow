import type { ReactNode } from "react";
import { FREQ_HELP_TEXT, SKY_MAX_N_FRAMES, SKY_MAX_PIXEL_ANIMATE } from "../lib/constants";
import type { ModeLogSliderSpec } from "../lib/types";
import { HelpHint } from "./HelpHint";
import { LogSliderField, SliderField } from "./SliderField";

type Props = {
  skyAnimate: boolean;
  setSkyAnimate: (value: boolean) => void;
  skyNpixel: number;
  setSkyNpixel: (value: number) => void;
  renderModeLogSliderRow: (specs: ModeLogSliderSpec[], className?: string) => ReactNode;
  skyTMin: number;
  setSkyTMin: (value: number) => void;
  skyTMax: number;
  setSkyTMax: (value: number) => void;
  skyNFrames: number;
  setSkyNFrames: (value: number) => void;
  skyTObs: number;
  setSkyTObs: (value: number) => void;
  skyNuInput: string;
  setSkyNuInput: (value: string) => void;
  skyFov: number;
  setSkyFov: (value: number) => void;
  skyPixelOptions: readonly number[];
  renderSharedControls: (options?: { hideDistanceObserver?: boolean; hidePlotUnits?: boolean }) => ReactNode;
};

export function ModePanelSkymap({
  skyAnimate,
  setSkyAnimate,
  skyNpixel,
  setSkyNpixel,
  renderModeLogSliderRow,
  skyTMin,
  setSkyTMin,
  skyTMax,
  setSkyTMax,
  skyNFrames,
  setSkyNFrames,
  skyTObs,
  setSkyTObs,
  skyNuInput,
  setSkyNuInput,
  skyFov,
  setSkyFov,
  skyPixelOptions,
  renderSharedControls,
}: Props) {
  return (
    <>
      <label className="sb-checkbox-inline">
        <input
          type="checkbox"
          checked={skyAnimate}
          onChange={(e) => {
            const nextAnimate = e.target.checked;
            setSkyAnimate(nextAnimate);
            if (nextAnimate && skyNpixel > SKY_MAX_PIXEL_ANIMATE) {
              setSkyNpixel(SKY_MAX_PIXEL_ANIMATE);
            }
          }}
        />
        Animate
      </label>

      {skyAnimate ? (
        <>
          {renderModeLogSliderRow(
            [
              {
                key: "sky-t-min",
                label: (
                  <>
                    log10(t<sub>min</sub> (s))
                  </>
                ),
                minExp: 2,
                maxExp: 9,
                step: 0.05,
                value: skyTMin,
                defaultExp: 4,
                onChange: setSkyTMin,
              },
              {
                key: "sky-t-max",
                label: (
                  <>
                    log10(t<sub>max</sub> (s))
                  </>
                ),
                minExp: 2,
                maxExp: 9,
                step: 0.05,
                value: skyTMax,
                defaultExp: 7,
                onChange: setSkyTMax,
              },
            ],
            "sb-row-2",
          )}
          <SliderField
            label="Frames"
            min={3}
            max={SKY_MAX_N_FRAMES}
            step={1}
            decimals={0}
            value={skyNFrames}
            onChange={(v) => setSkyNFrames(Math.max(3, Math.min(SKY_MAX_N_FRAMES, Math.round(v))))}
          />
        </>
      ) : (
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
      )}

      <label className="sb-field">
        <span className="sb-label sb-label-with-help">
          <span>
            ν<sub>obs</sub>
          </span>
          <HelpHint text={FREQ_HELP_TEXT} ariaLabel="Sky frequency input help" />
        </span>
        <input value={skyNuInput} onChange={(e) => setSkyNuInput(e.target.value)} placeholder="e.g. 1e9, 1GHz, 1keV" />
      </label>

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

      {renderSharedControls()}
    </>
  );
}
