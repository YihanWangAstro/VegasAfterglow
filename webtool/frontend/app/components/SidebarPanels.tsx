import type { ReactNode } from "react";
import type { UpdateValue } from "../hooks/useParameterState";
import { PLOT_FLUX_UNIT_OPTIONS, SKY_MAX_PIXEL_ANIMATE, SKY_MAX_PIXEL_STATIC, SKY_PIXEL_OPTIONS, TIME_UNIT_OPTIONS } from "../lib/constants";
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
  lcFreq: string;
  setLcFreq: (value: string) => void;
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
  skyAnimate: boolean;
  setSkyAnimate: (value: boolean) => void;
  skyNpixel: number;
  setSkyNpixel: (value: number) => void;
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
  const { distanceDriver, distanceLinked, setDistanceLinked, setDistanceDriver, setDistanceMpc, setDistanceRedshift, setSharedField } =
    distance;
  const {
    instrumentGroups,
    lcInstruments,
    sedInstruments,
    setLcInstruments,
    setSedInstruments,
    lcObsGroups,
    sedObsGroups,
    activeLcObsTab,
    activeSedObsTab,
    setLcObsGroups,
    setSedObsGroups,
    setActiveLcObsTab,
    setActiveSedObsTab,
    setError,
    lcCurveOptions,
    sedCurveOptions,
  } = sharedControls;
  const {
    bookmarkNameDraft,
    setBookmarkNameDraft,
    saveCurrentBookmark,
    setupStatusText,
    bookmarks,
    loadBookmarkById,
    copyShareLinkForSnapshot,
    removeBookmarkById,
  } = bookmark;
  const { compareEnabled, setCompareEnabled, compareBookmarkId, setCompareBookmarkId, compareRunning, compareError } = compare;
  const { lcFreq, setLcFreq, lcTMin, setLcTMin, lcTMax, setLcTMax } = lightcurve;
  const { sedTimesDraft, setSedTimesDraft, commitSedTimes, sedNuMin, setSedNuMin, sedNuMax, setSedNuMax, sedNumNu, setSedNumNu, sedNuFNu, setSedNuFNu, sedFreqUnit, setSedFreqUnit } =
    spectrum;
  const { skyAnimate, setSkyAnimate, skyNpixel, setSkyNpixel, skyTMin, setSkyTMin, skyTMax, setSkyTMax, skyNFrames, setSkyNFrames, skyTObs, setSkyTObs, skyNuInput, setSkyNuInput, skyFov, setSkyFov } =
    skymap;

  function renderDistanceObserverControls() {
    return (
      <DistanceControls
        distanceDriver={distanceDriver}
        distanceLinked={distanceLinked}
        shared={shared}
        setDistanceLinked={setDistanceLinked}
        setDistanceDriver={setDistanceDriver}
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
        mode={mode}
        shared={shared}
        setSharedField={setSharedField}
        hideDistanceObserver={options?.hideDistanceObserver}
        hidePlotUnits={options?.hidePlotUnits}
        distanceControls={renderDistanceObserverControls()}
        plotUnitControls={renderPlotUnitControls(mode !== "lightcurve", mode === "lightcurve")}
        instrumentGroups={instrumentGroups}
        lcInstruments={lcInstruments}
        sedInstruments={sedInstruments}
        setLcInstruments={setLcInstruments}
        setSedInstruments={setSedInstruments}
        lcObsGroups={lcObsGroups}
        sedObsGroups={sedObsGroups}
        activeLcObsTab={activeLcObsTab}
        activeSedObsTab={activeSedObsTab}
        setLcObsGroups={setLcObsGroups}
        setSedObsGroups={setSedObsGroups}
        setActiveLcObsTab={setActiveLcObsTab}
        setActiveSedObsTab={setActiveSedObsTab}
        setError={setError}
        lcCurveOptions={lcCurveOptions}
        sedCurveOptions={sedCurveOptions}
      />
    );
  }

  const skyPixelMax = skyAnimate ? SKY_MAX_PIXEL_ANIMATE : SKY_MAX_PIXEL_STATIC;
  const skyPixelOptions = SKY_PIXEL_OPTIONS.filter((value) => value <= skyPixelMax);

  return (
    <>
      <BookmarkPanel
        mode={mode}
        bookmarkNameDraft={bookmarkNameDraft}
        setBookmarkNameDraft={setBookmarkNameDraft}
        saveCurrentBookmark={saveCurrentBookmark}
        setupStatusText={setupStatusText}
        bookmarks={bookmarks}
        loadBookmarkById={loadBookmarkById}
        copyShareLinkForSnapshot={copyShareLinkForSnapshot}
        removeBookmarkById={removeBookmarkById}
        compareEnabled={compareEnabled}
        setCompareEnabled={setCompareEnabled}
        compareBookmarkId={compareBookmarkId}
        setCompareBookmarkId={setCompareBookmarkId}
        compareRunning={compareRunning}
        compareError={compareError}
      />

      {mode === "lightcurve" ? (
        <ModePanelLightcurve
          lcFreq={lcFreq}
          setLcFreq={setLcFreq}
          renderDistanceObserverControls={renderDistanceObserverControls}
          renderModeLogSliderRow={renderModeLogSliderRow}
          lcTMin={lcTMin}
          setLcTMin={setLcTMin}
          lcTMax={lcTMax}
          setLcTMax={setLcTMax}
          numT={shared.num_t}
          setNumT={(value) => setSharedField("num_t", value)}
          renderPlotUnitControls={renderPlotUnitControls}
          renderSharedControls={renderSharedControls}
        />
      ) : null}

      {mode === "spectrum" ? (
        <ModePanelSpectrum
          sedTimesDraft={sedTimesDraft}
          setSedTimesDraft={setSedTimesDraft}
          commitSedTimes={commitSedTimes}
          renderDistanceObserverControls={renderDistanceObserverControls}
          renderModeLogSliderRow={renderModeLogSliderRow}
          sedNuMin={sedNuMin}
          setSedNuMin={setSedNuMin}
          sedNuMax={sedNuMax}
          setSedNuMax={setSedNuMax}
          sedNumNu={sedNumNu}
          setSedNumNu={setSedNumNu}
          sedNuFNu={sedNuFNu}
          setSedNuFNu={setSedNuFNu}
          fluxIsAbMag={shared.flux_unit === "AB mag"}
          sedFreqUnit={sedFreqUnit}
          setSedFreqUnit={setSedFreqUnit}
          renderStringOptions={renderStringOptions}
          renderFluxUnitControl={renderFluxUnitControl}
          renderSharedControls={renderSharedControls}
        />
      ) : null}

      {mode === "skymap" ? (
        <ModePanelSkymap
          skyAnimate={skyAnimate}
          setSkyAnimate={setSkyAnimate}
          skyNpixel={skyNpixel}
          setSkyNpixel={setSkyNpixel}
          renderModeLogSliderRow={renderModeLogSliderRow}
          skyTMin={skyTMin}
          setSkyTMin={setSkyTMin}
          skyTMax={skyTMax}
          setSkyTMax={setSkyTMax}
          skyNFrames={skyNFrames}
          setSkyNFrames={setSkyNFrames}
          skyTObs={skyTObs}
          setSkyTObs={setSkyTObs}
          skyNuInput={skyNuInput}
          setSkyNuInput={setSkyNuInput}
          skyFov={skyFov}
          setSkyFov={setSkyFov}
          skyPixelOptions={skyPixelOptions}
          renderSharedControls={renderSharedControls}
        />
      ) : null}
    </>
  );
}
