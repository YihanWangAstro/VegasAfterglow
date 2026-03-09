import { useCallback, useMemo, useReducer } from "react";
import { FALLBACK_SHARED } from "../lib/constants";
import type { DistanceDriver, Mode, ObservationGroup, SharedParams, UiSnapshot } from "../lib/types";

export type UpdateValue<T> = T | ((prev: T) => T);

export const INITIAL_PARAMETER_STATE = {
  mode: "lightcurve" as Mode,
  distanceLinked: true,
  distanceDriver: "dL" as DistanceDriver,
  shared: FALLBACK_SHARED as SharedParams,
  lcFreq: "1e9, R, 1keV",
  lcFreqDraft: "1e9, R, 1keV",
  lcTMin: 1,
  lcTMax: 1e8,
  lcInstruments: [] as string[],
  lcObsGroups: [] as ObservationGroup[],
  activeLcObsTab: 0,
  sedTimes: "1e3, 1e4, 1e5, 1e6",
  sedTimesDraft: "1e3, 1e4, 1e5, 1e6",
  sedNuMin: 1e8,
  sedNuMax: 1e20,
  sedNumNu: 200,
  sedFreqUnit: "Hz",
  sedNuFNu: false,
  sedInstruments: [] as string[],
  sedObsGroups: [] as ObservationGroup[],
  activeSedObsTab: 0,
  skyTObs: 1e6,
  skyNuInput: "1e9",
  skyNuInputDraft: "1e9",
  skyFov: 500,
  skyFovUnit: "μas" as string,
  skyNpixel: 256,
  skyIntensityUnit: "cgs" as string,
};

export type ParameterState = typeof INITIAL_PARAMETER_STATE;

type ParameterAction =
  | {
      type: "set";
      key: keyof ParameterState;
      value: unknown;
    }
  | {
      type: "apply_snapshot";
      snapshot: UiSnapshot;
    };

function resolveUpdateValue<T>(prev: T, value: unknown): T {
  if (typeof value === "function") {
    return (value as (prev: T) => T)(prev);
  }
  return value as T;
}

export function parameterReducer(state: ParameterState, action: ParameterAction): ParameterState {
  if (action.type === "set") {
    const key = action.key;
    const prevValue = state[key] as unknown;
    const nextValue = resolveUpdateValue(prevValue, action.value);
    if (Object.is(prevValue, nextValue)) return state;
    return { ...state, [key]: nextValue } as ParameterState;
  }

  const { snapshot } = action;
  return {
    ...state,
    mode: snapshot.mode,
    distanceLinked: snapshot.distance_linked,
    distanceDriver: snapshot.distance_driver,
    shared: { ...FALLBACK_SHARED, ...snapshot.shared },
    lcFreq: snapshot.lightcurve.frequencies_input,
    lcFreqDraft: snapshot.lightcurve.frequencies_input,
    lcTMin: snapshot.lightcurve.t_min,
    lcTMax: snapshot.lightcurve.t_max,
    lcInstruments: [...snapshot.lightcurve.selected_instruments],
    lcObsGroups: snapshot.lightcurve.observation_groups.map((group) => ({ ...group })),
    activeLcObsTab: 0,
    sedTimes: snapshot.spectrum.t_snapshots_input,
    sedTimesDraft: snapshot.spectrum.t_snapshots_input,
    sedNuMin: snapshot.spectrum.nu_min,
    sedNuMax: snapshot.spectrum.nu_max,
    sedNumNu: snapshot.spectrum.num_nu,
    sedFreqUnit: snapshot.spectrum.freq_unit,
    sedNuFNu: snapshot.spectrum.show_nufnu,
    sedInstruments: [...snapshot.spectrum.selected_instruments],
    sedObsGroups: snapshot.spectrum.observation_groups.map((group) => ({ ...group })),
    activeSedObsTab: 0,
    skyTObs: snapshot.skymap.t_obs,
    skyNuInput: snapshot.skymap.nu_input,
    skyNuInputDraft: snapshot.skymap.nu_input,
    skyFov: snapshot.skymap.fov,
    skyFovUnit: snapshot.skymap.fov_unit ?? "μas",
    skyNpixel: snapshot.skymap.npixel,
    skyIntensityUnit: snapshot.skymap.intensity_unit ?? "cgs",
  };
}

export function useParameterState() {
  const [parameterState, dispatchParameter] = useReducer(parameterReducer, INITIAL_PARAMETER_STATE);

  const setParam = useCallback(<K extends keyof ParameterState>(key: K, value: UpdateValue<ParameterState[K]>): void => {
    dispatchParameter({ type: "set", key, value });
  }, []);

  const actions = useMemo(() => ({
    setMode: (value: UpdateValue<Mode>) => setParam("mode", value),
    setDistanceLinked: (value: UpdateValue<boolean>) => setParam("distanceLinked", value),
    setDistanceDriver: (value: UpdateValue<DistanceDriver>) => setParam("distanceDriver", value),
    setShared: (value: UpdateValue<SharedParams>) => setParam("shared", value),
    setLcFreq: (value: UpdateValue<string>) => setParam("lcFreq", value),
    setLcFreqDraft: (value: UpdateValue<string>) => setParam("lcFreqDraft", value),
    setLcTMin: (value: UpdateValue<number>) => setParam("lcTMin", value),
    setLcTMax: (value: UpdateValue<number>) => setParam("lcTMax", value),
    setLcInstruments: (value: UpdateValue<string[]>) => setParam("lcInstruments", value),
    setLcObsGroups: (value: UpdateValue<ObservationGroup[]>) => setParam("lcObsGroups", value),
    setActiveLcObsTab: (value: UpdateValue<number>) => setParam("activeLcObsTab", value),
    setSedTimes: (value: UpdateValue<string>) => setParam("sedTimes", value),
    setSedTimesDraft: (value: UpdateValue<string>) => setParam("sedTimesDraft", value),
    setSedNuMin: (value: UpdateValue<number>) => setParam("sedNuMin", value),
    setSedNuMax: (value: UpdateValue<number>) => setParam("sedNuMax", value),
    setSedNumNu: (value: UpdateValue<number>) => setParam("sedNumNu", value),
    setSedFreqUnit: (value: UpdateValue<string>) => setParam("sedFreqUnit", value),
    setSedNuFNu: (value: UpdateValue<boolean>) => setParam("sedNuFNu", value),
    setSedInstruments: (value: UpdateValue<string[]>) => setParam("sedInstruments", value),
    setSedObsGroups: (value: UpdateValue<ObservationGroup[]>) => setParam("sedObsGroups", value),
    setActiveSedObsTab: (value: UpdateValue<number>) => setParam("activeSedObsTab", value),
    setSkyTObs: (value: UpdateValue<number>) => setParam("skyTObs", value),
    setSkyNuInput: (value: UpdateValue<string>) => setParam("skyNuInput", value),
    setSkyNuInputDraft: (value: UpdateValue<string>) => setParam("skyNuInputDraft", value),
    setSkyFov: (value: UpdateValue<number>) => setParam("skyFov", value),
    setSkyFovUnit: (value: UpdateValue<string>) => setParam("skyFovUnit", value),
    setSkyNpixel: (value: UpdateValue<number>) => setParam("skyNpixel", value),
    setSkyIntensityUnit: (value: UpdateValue<string>) => setParam("skyIntensityUnit", value),
    applySnapshot: (snapshot: UiSnapshot) => dispatchParameter({ type: "apply_snapshot", snapshot }),
  }), [setParam]);

  return { parameterState, actions };
}

export type ParameterActions = ReturnType<typeof useParameterState>["actions"];
