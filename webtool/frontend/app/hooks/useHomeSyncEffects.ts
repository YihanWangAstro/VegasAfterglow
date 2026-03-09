import { useEffect, useRef, type Dispatch, type MutableRefObject, type SetStateAction } from "react";
import type { AxisRanges, Mode, SharedParams, UiSnapshot } from "../lib/types";
import { clampTabIndex } from "../lib/utils/snapshot";
import { clearTimeoutRef } from "../lib/utils/async";

type Args = {
  bookmarksLength: number;
  hasCompareBookmarkId: (id: string) => boolean;
  firstBookmarkId: string;
  setCompareBookmarkId: Dispatch<SetStateAction<string>>;
  setCompareEnabled: Dispatch<SetStateAction<boolean>>;

  lcObsGroupsLength: number;
  sedObsGroupsLength: number;
  setActiveLcObsTab: (value: number | ((prev: number) => number)) => void;
  setActiveSedObsTab: (value: number | ((prev: number) => number)) => void;

  urlStateReady: boolean;
  sliderInteracting: boolean;
  buildShareLink: (snapshot: UiSnapshot) => string | null;
  currentSnapshot: UiSnapshot;

  clearSavedXAxis: (mode: Mode) => void;
  clearSavedYAxes: (mode: Mode) => void;
  lcTMin: number;
  lcTMax: number;
  sedNuMin: number;
  sedNuMax: number;
  sedNuFNu: boolean;
  sedFreqUnit: string;
  shared: SharedParams;
  skyFov: number;
  skyNpixel: number;
  mode: Mode;
  axisRangesRef: MutableRefObject<Record<Mode, AxisRanges>>;
  setZoomRevision: Dispatch<SetStateAction<number>>;
};

export function useHomeSyncEffects({
  bookmarksLength,
  hasCompareBookmarkId,
  firstBookmarkId,
  setCompareBookmarkId,
  setCompareEnabled,
  lcObsGroupsLength,
  sedObsGroupsLength,
  setActiveLcObsTab,
  setActiveSedObsTab,
  urlStateReady,
  sliderInteracting,
  buildShareLink,
  currentSnapshot,
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
  skyNpixel,
  mode,
  axisRangesRef,
  setZoomRevision,
}: Args) {
  const urlSyncTimerRef = useRef<number | null>(null);

  useEffect(() => {
    if (bookmarksLength === 0) {
      setCompareBookmarkId("");
      setCompareEnabled(false);
      return;
    }
    setCompareBookmarkId((prev) => (hasCompareBookmarkId(prev) ? prev : firstBookmarkId));
  }, [bookmarksLength, firstBookmarkId, hasCompareBookmarkId, setCompareBookmarkId, setCompareEnabled]);

  useEffect(() => {
    setActiveLcObsTab((prev) => clampTabIndex(prev, lcObsGroupsLength));
    setActiveSedObsTab((prev) => clampTabIndex(prev, sedObsGroupsLength));
  }, [lcObsGroupsLength, sedObsGroupsLength, setActiveLcObsTab, setActiveSedObsTab]);

  useEffect(() => {
    if (!urlStateReady || typeof window === "undefined") return;
    if (sliderInteracting) return;
    clearTimeoutRef(urlSyncTimerRef);
    const link = buildShareLink(currentSnapshot);
    if (!link) return;
    urlSyncTimerRef.current = window.setTimeout(() => {
      urlSyncTimerRef.current = null;
      const current = window.location.href;
      if (current === link) return;
      window.history.replaceState({}, "", link);
    }, 120);

    return () => {
      clearTimeoutRef(urlSyncTimerRef);
    };
  }, [buildShareLink, currentSnapshot, sliderInteracting, urlStateReady]);

  useEffect(() => {
    clearSavedXAxis("lightcurve");
    clearSavedXAxis("spectrum");
    clearSavedYAxes("spectrum");
    clearSavedXAxis("skymap");
    clearSavedYAxes("skymap");
  }, [
    clearSavedXAxis,
    clearSavedYAxes,
    lcTMax,
    lcTMin,
    sedNuFNu,
    sedNuMax,
    sedNuMin,
    sedFreqUnit,
    shared.flux_unit,
    shared.time_unit,
    skyFov,
    skyNpixel,
  ]);

  useEffect(() => {
    axisRangesRef.current[mode] = { xaxis: null, yaxis: null, yaxis2: null };
    setZoomRevision((value) => value + 1);
  }, [axisRangesRef, mode, setZoomRevision]);

  useEffect(() => {
    return () => {
      clearTimeoutRef(urlSyncTimerRef);
    };
  }, []);
}
