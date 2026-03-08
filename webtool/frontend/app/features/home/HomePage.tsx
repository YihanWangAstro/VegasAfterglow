"use client";

import {
  useCallback,
  useEffect,
  useMemo,
  useRef,
  useState,
  useTransition,
} from "react";
import { PlotFigure } from "../../components/PlotFigure";
import { SidebarFooter } from "../../components/SidebarFooter";
import { SidebarPanels } from "../../components/SidebarPanels";
import { ServerStatus } from "../../components/ServerStatus";
import {
  BIBTEX_TEXT,
  DEFAULT_APP_VERSION,
  MODE_OPTIONS,
  SKY_MAX_PIXEL_ANIMATE,
  SPECTRUM_TEXT_COMMIT_DEBOUNCE_MS,
} from "../../lib/constants";
import type {
  ComputationSpec,
  InstrumentGroup,
  LcPlotData,
  Mode,
  RunResponse,
  SedPlotData,
  SelectOption,
  UiSnapshot,
} from "../../lib/types";
import {
  buildComputationSpecFromSnapshot,
  cloneSnapshot,
} from "../../lib/utils/snapshot";
import { formatFreqLabel, formatBandLabel, formatTimeLabel } from "../../lib/utils/plot-builders";
import { useApiStatus } from "../../hooks/useApiStatus";
import { useAxisPersistence } from "../../hooks/useAxisPersistence";
import { useBookmarks } from "../../hooks/useBookmarks";
import { useBootstrapState } from "../../hooks/useBootstrapState";
import { useComputationEffects } from "../../hooks/useComputationEffects";
import { useComputationApi } from "../../hooks/useComputationApi";
import { useDistanceLinking } from "../../hooks/useDistanceLinking";
import { useDownloadExports } from "../../hooks/useDownloadExports";
import { useFigurePresentation } from "../../hooks/useFigurePresentation";
import { useHomeSyncEffects } from "../../hooks/useHomeSyncEffects";
import { usePlotWidthObserver, useViewportLayout } from "../../hooks/useLayoutEffects";
import { useParameterState } from "../../hooks/useParameterState";
import { useShareActions } from "../../hooks/useShareActions";
import { useSliderInteracting } from "../../hooks/useSliderInteracting";

export default function HomePage() {
  const {
    apiCandidates,
    fetchFromApi,
    activeApiStatusRow,
    otherApiStatusRows,
    probeOtherServersOnce,
  } = useApiStatus();

  const { parameterState, actions } = useParameterState();
  const {
    mode,
    distanceLinked,
    distanceDriver,
    shared,
    lcFreq,
    lcTMin,
    lcTMax,
    lcInstruments,
    lcObsGroups,
    activeLcObsTab,
    sedTimes,
    sedTimesDraft,
    sedNuMin,
    sedNuMax,
    sedNumNu,
    sedFreqUnit,
    sedNuFNu,
    sedInstruments,
    sedObsGroups,
    activeSedObsTab,
    skyAnimate,
    skyTObs,
    skyTMin,
    skyTMax,
    skyNFrames,
    skyNuInput,
    skyFov,
    skyNpixel,
  } = parameterState;

  const {
    setMode,
    setDistanceLinked,
    setDistanceDriver,
    setShared,
    setLcFreq,
    setLcTMin,
    setLcTMax,
    setLcInstruments,
    setLcObsGroups,
    setActiveLcObsTab,
    setSedTimes,
    setSedTimesDraft,
    setSedNuMin,
    setSedNuMax,
    setSedNumNu,
    setSedFreqUnit,
    setSedNuFNu,
    setSedInstruments,
    setSedObsGroups,
    setActiveSedObsTab,
    setSkyAnimate,
    setSkyTObs,
    setSkyTMin,
    setSkyTMax,
    setSkyNFrames,
    setSkyNuInput,
    setSkyFov,
    setSkyNpixel,
  } = actions;

  const [sidebarOpen, setSidebarOpen] = useState(false);

  const [error, setError] = useState<string>("");
  const [result, setResult] = useState<RunResponse | null>(null);
  const [resultMode, setResultMode] = useState<Mode | null>(null);
  const [compareResult, setCompareResult] = useState<RunResponse | null>(null);
  const [compareResultMode, setCompareResultMode] = useState<Mode | null>(null);
  const [compareRunning, setCompareRunning] = useState(false);
  const [compareError, setCompareError] = useState("");
  const [plotWidthPx, setPlotWidthPx] = useState(960);
  const [isFigurePending, startFigureTransition] = useTransition();
  const [appVersion, setAppVersion] = useState(DEFAULT_APP_VERSION);
  const [citeLinkText, setCiteLinkText] = useState("cite");
  const [zoomRevision, setZoomRevision] = useState(0);
  const workspaceRef = useRef<HTMLElement | null>(null);
  const { axisRangesRef, axisSignaturesRef, clearSavedXAxis: clearSavedXAxisBase, clearSavedYAxes: clearSavedYAxesBase } =
    useAxisPersistence();
  const sliderInteracting = useSliderInteracting();

  const [instrumentGroups, setInstrumentGroups] = useState<InstrumentGroup[]>([]);
  const {
    bookmarks,
    setBookmarks,
    bookmarkNameDraft,
    setBookmarkNameDraft,
    setupStatusText,
    setSetupStatusText,
    saveSnapshotAsBookmark,
    removeBookmarkById,
    buildShareLink,
  } = useBookmarks();
  const [compareEnabled, setCompareEnabled] = useState(false);
  const [compareBookmarkId, setCompareBookmarkId] = useState("");

  const { setSharedField, setDistanceMpc, setDistanceRedshift } = useDistanceLinking({
    distanceLinked,
    distanceDriver,
    shared,
    setShared,
  });

  const currentSnapshot = useMemo<UiSnapshot>(
    () => ({
      mode,
      distance_linked: distanceLinked,
      distance_driver: distanceDriver,
      shared: { ...shared },
      lightcurve: {
        frequencies_input: lcFreq,
        t_min: lcTMin,
        t_max: lcTMax,
        selected_instruments: [...lcInstruments],
        observation_groups: lcObsGroups.map((group) => ({ ...group })),
      },
      spectrum: {
        t_snapshots_input: sedTimes,
        nu_min: sedNuMin,
        nu_max: sedNuMax,
        num_nu: sedNumNu,
        freq_unit: sedFreqUnit,
        show_nufnu: sedNuFNu,
        selected_instruments: [...sedInstruments],
        observation_groups: sedObsGroups.map((group) => ({ ...group })),
      },
      skymap: {
        animate: skyAnimate,
        t_obs: skyTObs,
        t_min: skyTMin,
        t_max: skyTMax,
        n_frames: skyNFrames,
        nu_input: skyNuInput,
        fov: skyFov,
        npixel: skyNpixel,
      },
    }),
    [parameterState],
  );

  const applySnapshot = useCallback((snapshot: UiSnapshot) => {
    actions.applySnapshot(snapshot);
    setError("");
    setCompareError("");
  }, [actions.applySnapshot]);

  const { copyShareLinkForSnapshot, copyBibtex } = useShareActions({
    buildShareLink,
    setSetupStatusText,
    setCiteLinkText,
    bibtexText: BIBTEX_TEXT,
  });

  const saveCurrentBookmark = useCallback(() => {
    saveSnapshotAsBookmark(currentSnapshot, bookmarkNameDraft);
  }, [bookmarkNameDraft, currentSnapshot, saveSnapshotAsBookmark]);

  const clearCompareOverlay = useCallback((clearError = true) => {
    setCompareResult(null);
    setCompareResultMode(null);
    if (clearError) {
      setCompareError("");
    }
  }, []);

  const loadBookmarkById = useCallback(
    (bookmarkId: string) => {
      const target = bookmarks.find((bookmark) => bookmark.id === bookmarkId);
      if (!target) return;
      try {
        const snapshot = cloneSnapshot(target.snapshot);
        if (!snapshot || typeof snapshot !== "object") throw new Error("Invalid snapshot");
        if (!snapshot.lightcurve || !snapshot.spectrum || !snapshot.skymap) throw new Error("Invalid snapshot");
        applySnapshot(snapshot);
      } catch {
        setError("Bookmark is incompatible with current app version.");
      }
    },
    [applySnapshot, bookmarks],
  );

  const selectedCompareBookmark = useMemo(
    () => bookmarks.find((bookmark) => bookmark.id === compareBookmarkId) ?? null,
    [bookmarks, compareBookmarkId],
  );

  const compareSpec = useMemo<ComputationSpec | null>(() => {
    if (!compareEnabled || mode === "skymap" || !selectedCompareBookmark) return null;
    return buildComputationSpecFromSnapshot(selectedCompareBookmark.snapshot, mode);
  }, [compareEnabled, mode, selectedCompareBookmark]);

  const { bootReady, urlStateReady } = useBootstrapState({
    fetchFromApi,
    actions,
    setInstrumentGroups,
    setAppVersion,
    bookmarks,
    setBookmarks,
    lcObsGroups,
    sedObsGroups,
    applySnapshot,
    setError,
  });

  const { warnings, handlePlotRelayout, displayFigure, deferredFigure, figureCaption } = useFigurePresentation({
    mode,
    result,
    resultMode,
    compareEnabled,
    compareResult,
    compareResultMode,
    selectedCompareBookmark,
    plotWidthPx,
    zoomRevision,
    shared,
    lcObsGroups,
    sedObsGroups,
    axisRangesRef,
    axisSignaturesRef,
    setZoomRevision,
  });
  usePlotWidthObserver({ workspaceRef, setPlotWidthPx });

  const clearSavedXAxis = useCallback(
    (targetMode: Mode) => {
      const before = axisRangesRef.current[targetMode].xaxis;
      if (before === null) return;
      clearSavedXAxisBase(targetMode);
      if (mode === targetMode) {
        setZoomRevision((value) => value + 1);
      }
    },
    [axisRangesRef, clearSavedXAxisBase, mode],
  );

  const clearSavedYAxes = useCallback(
    (targetMode: Mode) => {
      const current = axisRangesRef.current[targetMode];
      if (current.yaxis === null && current.yaxis2 === null) return;
      clearSavedYAxesBase(targetMode);
      if (mode === targetMode) {
        setZoomRevision((value) => value + 1);
      }
    },
    [axisRangesRef, clearSavedYAxesBase, mode],
  );
  const { fullComputationSpec, computationSpec, postRunRequest, postFigureRequest } = useComputationApi({
    parameterState,
    sliderInteracting,
    fetchFromApi,
  });

  const { showRunning, showColdStartHint } = useComputationEffects({
    bootReady,
    computationSpec,
    compareSpec,
    compareEnabled,
    mode,
    sliderInteracting,
    apiCandidates,
    postFigureRequest,
    startFigureTransition,
    setResult,
    setResultMode,
    setError,
    clearCompareOverlay,
    setCompareRunning,
    setCompareError,
    setCompareResult,
    setCompareResultMode,
  });

  useHomeSyncEffects({
    bookmarksLength: bookmarks.length,
    hasCompareBookmarkId: (id) => bookmarks.some((bookmark) => bookmark.id === id),
    firstBookmarkId: bookmarks[0]?.id ?? "",
    setCompareBookmarkId,
    setCompareEnabled,
    lcObsGroupsLength: lcObsGroups.length,
    sedObsGroupsLength: sedObsGroups.length,
    setActiveLcObsTab,
    setActiveSedObsTab,
    urlStateReady,
    sliderInteracting,
    buildShareLink,
    currentSnapshot,
    sedTimesDraft,
    sedTimes,
    setSedTimes,
    spectrumTextCommitDebounceMs: SPECTRUM_TEXT_COMMIT_DEBOUNCE_MS,
    skyAnimate,
    skyNpixel,
    setSkyNpixel,
    skyMaxPixelAnimate: SKY_MAX_PIXEL_ANIMATE,
    clearSavedXAxis,
    clearSavedYAxes,
    lcTMin,
    lcTMax,
    sedNuMin,
    sedNuMax,
    sedNuFNu,
    sedFreqUnit,
    shared,
    skyFov,
    mode,
    axisRangesRef,
    setZoomRevision,
  });

  useViewportLayout(sidebarOpen);

  const { downloading, downloadCurrentExport } = useDownloadExports({
    mode,
    fullComputationSpec,
    postRunRequest,
    exportsCache: result?.exports,
    setResult,
    setError,
  });

  const resultPlotData = result?.plot_data;
  const lcCurveOptions = useMemo<SelectOption[] | undefined>(() => {
    if (!resultPlotData) return undefined;
    const pd = resultPlotData as LcPlotData;
    const ptOpts = pd?.pt
      ? pd.pt.freq_hz.map((nu) => ({ label: formatFreqLabel(nu), value: String(nu) }))
      : [];
    const bandOpts = (pd?.bands ?? []).map((b) => {
      const label = formatBandLabel(b.nu_min, b.nu_max, b.name);
      return { label, value: `band:${label}` };
    });
    return [...ptOpts, ...bandOpts];
  }, [resultPlotData]);

  const sedCurveOptions = useMemo<SelectOption[] | undefined>(() => {
    if (!resultPlotData) return undefined;
    const pd = resultPlotData as SedPlotData;
    if (!pd?.t_snapshots_s) return [];
    return pd.t_snapshots_s.map((t) => ({ label: formatTimeLabel(t), value: String(t) }));
  }, [resultPlotData]);

  const sidebarPanelsProps = {
    mode,
    shared,
    distance: {
      distanceDriver,
      distanceLinked,
      setDistanceLinked,
      setDistanceDriver,
      setDistanceMpc,
      setDistanceRedshift,
      setSharedField,
    },
    sharedControls: {
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
    },
    bookmark: {
      bookmarkNameDraft,
      setBookmarkNameDraft,
      saveCurrentBookmark,
      setupStatusText,
      bookmarks,
      loadBookmarkById,
      copyShareLinkForSnapshot,
      removeBookmarkById,
    },
    compare: {
      compareEnabled,
      setCompareEnabled,
      compareBookmarkId,
      setCompareBookmarkId,
      compareRunning,
      compareError,
    },
    lightcurve: { lcFreq, setLcFreq, lcTMin, setLcTMin, lcTMax, setLcTMax },
    spectrum: {
      sedTimesDraft,
      setSedTimesDraft,
      commitSedTimes: () => setSedTimes(sedTimesDraft),
      sedNuMin,
      setSedNuMin,
      sedNuMax,
      setSedNuMax,
      sedNumNu,
      setSedNumNu,
      sedNuFNu,
      setSedNuFNu,
      sedFreqUnit,
      setSedFreqUnit,
    },
    skymap: {
      skyAnimate,
      setSkyAnimate,
      skyNpixel,
      setSkyNpixel,
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
    },
  };

  const emptyStateMessage =
    mode === "lightcurve" && lcFreq.trim().length === 0
      ? "Enter at least one frequency, filter, or band to render a light curve."
      : mode === "spectrum" && sedTimes.trim().length === 0
        ? "Enter at least one observation time to render a spectrum."
        : "Adjust parameters in the sidebar to update the plot.";

  return (
    <main className={`app-shell ${sidebarOpen ? "sidebar-open" : ""}`}>
      {!sidebarOpen ? (
        <button className="sidebar-toggle" onClick={() => setSidebarOpen(true)} aria-label="Open controls">
          Controls
        </button>
      ) : null}

      <aside className={`sidebar ${sidebarOpen ? "open" : ""}`}>
        <div className="sidebar-header">
          <a
            className="sidebar-logo-link"
            href="https://github.com/YihanWangAstro/VegasAfterglow"
            target="_blank"
            rel="noreferrer"
          >
            <img src="/logo-horizontal.svg" alt="VegasAfterglow" />
          </a>
          <button className="sidebar-close" onClick={() => setSidebarOpen(false)}>
            Close
          </button>
        </div>
        <ServerStatus
          activeApiStatusRow={activeApiStatusRow}
          otherApiStatusRows={otherApiStatusRows}
          probeOtherServersOnce={probeOtherServersOnce}
        />
        <div className="mode-block">
          <span className="sb-label sb-section-label">Mode</span>
          <div className="mode-radios">
            {MODE_OPTIONS.map((item) => (
              <label key={item.value} className="mode-radio">
                <input
                  type="radio"
                  name="plot-mode"
                  checked={mode === item.value}
                  onChange={() => setMode(item.value)}
                />
                <span>{item.label}</span>
              </label>
            ))}
          </div>
        </div>

        <SidebarPanels {...sidebarPanelsProps} />

        <SidebarFooter
          mode={mode}
          downloading={downloading}
          hasFigureData={Boolean(displayFigure?.data)}
          appVersion={appVersion}
          citeLinkText={citeLinkText}
          onDownload={(kind) => void downloadCurrentExport(kind)}
          onCopyBibtex={() => void copyBibtex()}
        />

        {showRunning || isFigurePending ? <div className="status">{`Running ${mode} ...`}</div> : null}
        {error ? <div className="error">{error}</div> : null}
      </aside>

      {sidebarOpen ? (
        <button className="sidebar-backdrop" onClick={() => setSidebarOpen(false)} aria-label="Close sidebar" />
      ) : null}

      <section ref={workspaceRef} className="workspace">
        {showColdStartHint ? (
          <div className="startup-hint" role="status" aria-live="polite">
            Processing request... this can take longer for complex parameters or a waking backend.
          </div>
        ) : null}
        {warnings.length > 0 ? (
          <div className="warn-box">
            {warnings.map((warning, idx) => (
              <p key={idx}>{warning}</p>
            ))}
          </div>
        ) : null}

        {deferredFigure?.data ? (
          <div className={`figure-stack${mode === "skymap" ? " figure-stack-square" : ""}`}>
            <PlotFigure figure={deferredFigure} mode={mode} onRelayout={handlePlotRelayout} />
            <p className="figure-caption">{figureCaption}</p>
          </div>
        ) : (
          <div className="workspace-empty">
            <p className="muted">{emptyStateMessage}</p>
            {error ? <p className="workspace-inline-error">{error}</p> : null}
          </div>
        )}

      </section>
    </main>
  );
}
