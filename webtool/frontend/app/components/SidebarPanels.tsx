import type { ReactNode } from "react";
import { useParameters } from "../lib/ParameterContext";
import { PLOT_FLUX_UNIT_OPTIONS, SKY_MAX_PIXEL_STATIC, SKY_PIXEL_OPTIONS, TIME_UNIT_OPTIONS } from "../lib/constants";
import type {
  BookmarkEntry,
  InstrumentGroup,
  ModeLogSliderSpec,
  SelectOption,
  SharedParams,
} from "../lib/types";
import { BookmarkPanel } from "./BookmarkPanel";
import { DistanceControls } from "./DistanceControls";
import { ModePanelLightcurve } from "./ModePanelLightcurve";
import { ModePanelSkymap } from "./ModePanelSkymap";
import { ModePanelSpectrum } from "./ModePanelSpectrum";
import { SharedControls } from "./SharedControls";
import { LogSliderField } from "./SliderField";

type BookmarkSection = {
  bookmarkNameDraft: string;
  setBookmarkNameDraft: (value: string) => void;
  saveCurrentBookmark: () => void;
  setupStatusText: string;
  bookmarks: BookmarkEntry[];
  loadBookmarkById: (id: string) => void;
  copyShareLinkForSnapshot: (snapshot: BookmarkEntry["snapshot"]) => Promise<void>;
  removeBookmarkById: (id: string) => void;
};

type CompareSection = {
  compareEnabled: boolean;
  setCompareEnabled: (value: boolean) => void;
  compareBookmarkId: string;
  setCompareBookmarkId: (value: string) => void;
  compareRunning: boolean;
  compareError: string;
};

type Props = {
  instrumentGroups: InstrumentGroup[];
  setError: (value: string) => void;
  lcCurveOptions?: SelectOption[];
  sedCurveOptions?: SelectOption[];
  bookmark: BookmarkSection;
  compare: CompareSection;
};

const skyPixelOptions = SKY_PIXEL_OPTIONS.filter((v) => v <= SKY_MAX_PIXEL_STATIC);

function renderStringOptions(options: readonly string[]): ReactNode[] {
  return options.map((item) => (
    <option key={item} value={item}>
      {item}
    </option>
  ));
}

export function SidebarPanels({
  instrumentGroups,
  setError,
  lcCurveOptions,
  sedCurveOptions,
  bookmark,
  compare,
}: Props) {
  const { state, actions, setSharedField, setDistanceMpc, setDistanceRedshift } = useParameters();
  const { mode, shared } = state;

  function renderDistanceObserverControls() {
    return (
      <DistanceControls
        distanceDriver={state.distanceDriver}
        distanceLinked={state.distanceLinked}
        shared={shared}
        setDistanceLinked={actions.setDistanceLinked}
        setDistanceDriver={actions.setDistanceDriver}
        setDistanceMpc={setDistanceMpc}
        setDistanceRedshift={setDistanceRedshift}
        setSharedField={setSharedField}
      />
    );
  }

  function renderFluxUnitControl() {
    return (
      <label className="sb-field">
        <span className="sb-label">Flux</span>
        <select value={shared.flux_unit} onChange={(e) => setSharedField("flux_unit", e.target.value as SharedParams["flux_unit"])}>
          {renderStringOptions(PLOT_FLUX_UNIT_OPTIONS)}
        </select>
      </label>
    );
  }

  function renderPlotUnitControls(timeDisabled: boolean, showTime = true) {
    if (!showTime) return renderFluxUnitControl();

    return (
      <div className="sb-row-2">
        {renderFluxUnitControl()}
        <label className="sb-field">
          <span className="sb-label">Time</span>
          <select value={shared.time_unit} disabled={timeDisabled} onChange={(e) => setSharedField("time_unit", e.target.value as SharedParams["time_unit"])}>
            {renderStringOptions(TIME_UNIT_OPTIONS)}
          </select>
        </label>
      </div>
    );
  }

  function renderModeLogSlider(spec: ModeLogSliderSpec) {
    return (
      <LogSliderField
        key={spec.key}
        label={spec.label}
        minExp={spec.minExp}
        maxExp={spec.maxExp}
        step={spec.step}
        value={spec.value}
        defaultExp={spec.defaultExp}
        onChange={spec.onChange}
      />
    );
  }

  function renderModeLogSliderRow(specs: ModeLogSliderSpec[], className = "sb-row-2 sb-row-gap-top") {
    return <div className={className}>{specs.map((spec) => renderModeLogSlider(spec))}</div>;
  }

  function renderSharedControls(options?: { hideDistanceObserver?: boolean; hidePlotUnits?: boolean }) {
    return (
      <SharedControls
        instrumentGroups={instrumentGroups}
        setError={setError}
        lcCurveOptions={lcCurveOptions}
        sedCurveOptions={sedCurveOptions}
        hideDistanceObserver={options?.hideDistanceObserver}
        hidePlotUnits={options?.hidePlotUnits}
        distanceControls={renderDistanceObserverControls()}
        plotUnitControls={renderPlotUnitControls(mode !== "lightcurve", mode === "lightcurve")}
      />
    );
  }

  return (
    <>
      <BookmarkPanel
        mode={mode}
        {...bookmark}
        {...compare}
      />

      {mode === "lightcurve" ? (
        <ModePanelLightcurve
          renderDistanceObserverControls={renderDistanceObserverControls}
          renderModeLogSliderRow={renderModeLogSliderRow}
          renderPlotUnitControls={renderPlotUnitControls}
          renderSharedControls={renderSharedControls}
        />
      ) : null}

      {mode === "spectrum" ? (
        <ModePanelSpectrum
          renderDistanceObserverControls={renderDistanceObserverControls}
          renderModeLogSliderRow={renderModeLogSliderRow}
          renderFluxUnitControl={renderFluxUnitControl}
          renderSharedControls={renderSharedControls}
        />
      ) : null}

      {mode === "skymap" ? (
        <ModePanelSkymap
          skyPixelOptions={skyPixelOptions}
          renderDistanceObserverControls={renderDistanceObserverControls}
          renderSharedControls={renderSharedControls}
        />
      ) : null}
    </>
  );
}
