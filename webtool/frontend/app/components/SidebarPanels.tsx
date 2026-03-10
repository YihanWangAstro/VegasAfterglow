import type { ReactNode } from "react";
import type { UpdateValue } from "../hooks/useParameterState";
import { PLOT_FLUX_UNIT_OPTIONS, SKY_MAX_PIXEL_STATIC, SKY_PIXEL_OPTIONS, TIME_UNIT_OPTIONS } from "../lib/constants";
import type {
  BookmarkEntry,
  DistanceDriver,
  InstrumentGroup,
  Mode,
  ModeLogSliderSpec,
  ObservationGroup,
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

type DistanceSection = {
  distanceDriver: DistanceDriver;
  distanceLinked: boolean;
  setDistanceLinked: (next: boolean) => void;
  setDistanceDriver: (next: DistanceDriver) => void;
  setDistanceMpc: (next: number) => void;
  setDistanceRedshift: (next: number) => void;
  setSharedField: <K extends keyof SharedParams>(key: K, value: SharedParams[K]) => void;
};

type SharedControlsSection = {
  instrumentGroups: InstrumentGroup[];
  lcInstruments: string[];
  sedInstruments: string[];
  setLcInstruments: (value: UpdateValue<string[]>) => void;
  setSedInstruments: (value: UpdateValue<string[]>) => void;
  lcObsGroups: ObservationGroup[];
  sedObsGroups: ObservationGroup[];
  activeLcObsTab: number;
  activeSedObsTab: number;
  setLcObsGroups: (value: UpdateValue<ObservationGroup[]>) => void;
  setSedObsGroups: (value: UpdateValue<ObservationGroup[]>) => void;
  setActiveLcObsTab: (value: UpdateValue<number>) => void;
  setActiveSedObsTab: (value: UpdateValue<number>) => void;
  setError: (value: string) => void;
  lcCurveOptions?: SelectOption[];
  sedCurveOptions?: SelectOption[];
};

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

type LightcurveSection = {
  lcFreqDraft: string;
  setLcFreqDraft: (value: string) => void;
  commitLcFreq: () => void;
  lcTMin: number;
  setLcTMin: (value: number) => void;
  lcTMax: number;
  setLcTMax: (value: number) => void;
};

type SpectrumSection = {
  sedTimesDraft: string;
  setSedTimesDraft: (value: string) => void;
  commitSedTimes: () => void;
  sedNuMin: number;
  setSedNuMin: (value: number) => void;
  sedNuMax: number;
  setSedNuMax: (value: number) => void;
  sedNumNu: number;
  setSedNumNu: (value: number) => void;
  sedNuFNu: boolean;
  setSedNuFNu: (value: boolean) => void;
  sedFreqUnit: string;
  setSedFreqUnit: (value: string) => void;
};

type SkymapSection = {
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
};

type Props = {
  mode: Mode;
  shared: SharedParams;
  distance: DistanceSection;
  sharedControls: SharedControlsSection;
  bookmark: BookmarkSection;
  compare: CompareSection;
  lightcurve: LightcurveSection;
  spectrum: SpectrumSection;
  skymap: SkymapSection;
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
  mode,
  shared,
  distance,
  sharedControls,
  bookmark,
  compare,
  lightcurve,
  spectrum,
  skymap,
}: Props) {
  function renderDistanceObserverControls() {
    return (
      <DistanceControls
        distanceDriver={distance.distanceDriver}
        distanceLinked={distance.distanceLinked}
        shared={shared}
        setDistanceLinked={distance.setDistanceLinked}
        setDistanceDriver={distance.setDistanceDriver}
        setDistanceMpc={distance.setDistanceMpc}
        setDistanceRedshift={distance.setDistanceRedshift}
        setSharedField={distance.setSharedField}
      />
    );
  }

  function renderFluxUnitControl() {
    return (
      <label className="sb-field">
        <span className="sb-label">Flux</span>
        <select value={shared.flux_unit} onChange={(e) => distance.setSharedField("flux_unit", e.target.value as SharedParams["flux_unit"])}>
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
          <select value={shared.time_unit} disabled={timeDisabled} onChange={(e) => distance.setSharedField("time_unit", e.target.value as SharedParams["time_unit"])}>
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
        mode={mode}
        shared={shared}
        setSharedField={distance.setSharedField}
        hideDistanceObserver={options?.hideDistanceObserver}
        hidePlotUnits={options?.hidePlotUnits}
        distanceControls={renderDistanceObserverControls()}
        plotUnitControls={renderPlotUnitControls(mode !== "lightcurve", mode === "lightcurve")}
        {...sharedControls}
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
          {...lightcurve}
          renderDistanceObserverControls={renderDistanceObserverControls}
          renderModeLogSliderRow={renderModeLogSliderRow}
          numT={shared.num_t}
          setNumT={(value) => distance.setSharedField("num_t", value)}
          renderPlotUnitControls={renderPlotUnitControls}
          renderSharedControls={renderSharedControls}
        />
      ) : null}

      {mode === "spectrum" ? (
        <ModePanelSpectrum
          {...spectrum}
          renderDistanceObserverControls={renderDistanceObserverControls}
          renderModeLogSliderRow={renderModeLogSliderRow}
          fluxIsAbMag={shared.flux_unit === "AB mag"}
          renderStringOptions={renderStringOptions}
          renderFluxUnitControl={renderFluxUnitControl}
          renderSharedControls={renderSharedControls}
        />
      ) : null}

      {mode === "skymap" ? (
        <ModePanelSkymap
          {...skymap}
          skyPixelOptions={skyPixelOptions}
          renderDistanceObserverControls={renderDistanceObserverControls}
          renderSharedControls={renderSharedControls}
        />
      ) : null}
    </>
  );
}
