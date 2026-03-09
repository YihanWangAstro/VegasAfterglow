import { useEffect, useRef, useState, type Dispatch, type SetStateAction } from "react";
import {
  BOOKMARKS_STORAGE_KEY,
  FALLBACK_SHARED,
  INSTRUMENT_CATALOG,
  MAX_BOOKMARKS,
  OBS_STORAGE_LC_KEY,
  OBS_STORAGE_SED_KEY,
  SKY_MAX_PIXEL_STATIC,
  URL_STATE_PARAM,
} from "../lib/constants";
import type { BookmarkEntry, DefaultsResponse, InstrumentGroup, ObservationGroup, UiSnapshot } from "../lib/types";
import { parseStoredObsGroups } from "../lib/utils/obs";
import {
  decodeUrlState,
  groupInstrumentsByType,
  parseStoredBookmarks,
  setIfArray,
  setIfBoolean,
  setIfNumber,
  setIfString,
} from "../lib/utils/snapshot";
import { readLocalStorageItem, writeLocalStorageJson } from "../lib/utils/storage";
import type { ParameterActions } from "./useParameterState";

type Args = {
  fetchFromApi: (path: string, init?: RequestInit) => Promise<Response>;
  actions: ParameterActions;
  setInstrumentGroups: Dispatch<SetStateAction<InstrumentGroup[]>>;
  setAppVersion: Dispatch<SetStateAction<string>>;
  bookmarks: BookmarkEntry[];
  setBookmarks: Dispatch<SetStateAction<BookmarkEntry[]>>;
  lcObsGroups: ObservationGroup[];
  sedObsGroups: ObservationGroup[];
  applySnapshot: (snapshot: UiSnapshot) => void;
  setError: (message: string) => void;
};

export function useBootstrapState({
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
}: Args) {
  const [bootReady, setBootReady] = useState(false);
  const [obsStorageReady, setObsStorageReady] = useState(false);
  const [urlStateReady, setUrlStateReady] = useState(false);
  const urlStateAppliedRef = useRef(false);

  useEffect(() => {
    async function boot() {
      try {
        const [defaultsRes, optionsRes] = await Promise.all([
          fetchFromApi("defaults", { cache: "no-store" }),
          fetchFromApi("options", { cache: "no-store" }),
        ]);

        if (defaultsRes.ok) {
          const data = (await defaultsRes.json()) as DefaultsResponse;

          if (data.shared) {
            actions.setShared({ ...FALLBACK_SHARED, ...data.shared });
          }

          if (data.lightcurve) {
            setIfString(data.lightcurve.frequencies_input, actions.setLcFreq);
            setIfNumber(data.lightcurve.t_min, actions.setLcTMin);
            setIfNumber(data.lightcurve.t_max, actions.setLcTMax);
            setIfArray<ObservationGroup>(data.lightcurve.observation_groups, actions.setLcObsGroups);
          }

          if (data.spectrum) {
            setIfString(data.spectrum.t_snapshots_input, (value) => {
              actions.setSedTimes(value);
              actions.setSedTimesDraft(value);
            });
            setIfNumber(data.spectrum.nu_min, actions.setSedNuMin);
            setIfNumber(data.spectrum.nu_max, actions.setSedNuMax);
            setIfNumber(data.spectrum.num_nu, actions.setSedNumNu);
            setIfString(data.spectrum.freq_unit, actions.setSedFreqUnit);
            setIfBoolean(data.spectrum.show_nufnu, actions.setSedNuFNu);
            setIfArray<ObservationGroup>(data.spectrum.observation_groups, actions.setSedObsGroups);
          }

          if (data.skymap) {
            setIfNumber(data.skymap.t_obs, actions.setSkyTObs);
            setIfString(data.skymap.nu_input, actions.setSkyNuInput);
            setIfNumber(data.skymap.fov, actions.setSkyFov);
            setIfNumber(data.skymap.npixel, (value) => {
              actions.setSkyNpixel(Math.max(64, Math.min(SKY_MAX_PIXEL_STATIC, Math.round(value))));
            });
          }
        }

        // Instrument catalog is defined locally — group from local constant.
        const catalogNames = Object.keys(INSTRUMENT_CATALOG);
        setInstrumentGroups(groupInstrumentsByType(catalogNames));
        const allowed = new Set(catalogNames);
        actions.setLcInstruments((prev) => prev.filter((name) => allowed.has(name)));
        actions.setSedInstruments((prev) => prev.filter((name) => allowed.has(name)));

        if (optionsRes.ok) {
          const options = (await optionsRes.json()) as Record<string, unknown>;
          if (typeof options.version === "string" && options.version.trim()) {
            setAppVersion(options.version.trim());
          }
        }
      } catch {
        // Keep defaults if backend/options unavailable.
      } finally {
        setBootReady(true);
      }
    }

    boot();
  }, [fetchFromApi]);

  useEffect(() => {
    if (!bootReady) return;
    const storedLc = parseStoredObsGroups(readLocalStorageItem(OBS_STORAGE_LC_KEY), true);
    if (storedLc !== null) {
      actions.setLcObsGroups(storedLc);
    }
    const storedSed = parseStoredObsGroups(readLocalStorageItem(OBS_STORAGE_SED_KEY), false);
    if (storedSed !== null) {
      actions.setSedObsGroups(storedSed);
    }
    setObsStorageReady(true);
  }, [actions, bootReady]);

  useEffect(() => {
    if (!obsStorageReady) return;
    writeLocalStorageJson(OBS_STORAGE_LC_KEY, lcObsGroups);
    writeLocalStorageJson(OBS_STORAGE_SED_KEY, sedObsGroups);
  }, [lcObsGroups, obsStorageReady, sedObsGroups]);

  useEffect(() => {
    const stored = parseStoredBookmarks(readLocalStorageItem(BOOKMARKS_STORAGE_KEY));
    if (stored !== null) {
      setBookmarks(stored.slice(0, MAX_BOOKMARKS));
    }
  }, [setBookmarks]);

  useEffect(() => {
    writeLocalStorageJson(BOOKMARKS_STORAGE_KEY, bookmarks);
  }, [bookmarks]);

  useEffect(() => {
    if (!bootReady || typeof window === "undefined") return;
    if (urlStateAppliedRef.current) return;
    urlStateAppliedRef.current = true;

    const raw = new URL(window.location.href).searchParams.get(URL_STATE_PARAM);
    if (!raw) {
      setUrlStateReady(true);
      return;
    }
    const snapshot = decodeUrlState(raw);
    if (!snapshot) {
      setError("URL state is invalid or incompatible.");
      setUrlStateReady(true);
      return;
    }
    applySnapshot(snapshot);
    setUrlStateReady(true);
  }, [applySnapshot, bootReady, setError]);

  return { bootReady, urlStateReady };
}
