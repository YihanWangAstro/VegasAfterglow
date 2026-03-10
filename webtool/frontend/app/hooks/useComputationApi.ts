import { useCallback, useMemo } from "react";
import { ENABLE_INTERACTIVE_DOWNSAMPLE } from "../lib/constants";
import type { ComputationSpec, RunResponse, SharedParams } from "../lib/types";
import { buildInteractiveSpec, normalizeShared } from "../lib/utils/snapshot";
import type { ParameterState } from "./useParameterState";

/** Strip presentation-only fields from shared params so they don't affect computation identity. */
function stripPresentationParams(shared: SharedParams): Omit<SharedParams, "flux_unit" | "time_unit"> {
  const { flux_unit: _, time_unit: _t, ...physics } = shared;
  return physics;
}

type Args = {
  parameterState: ParameterState;
  sliderInteracting: boolean;
  fetchFromApi: (path: string, init?: RequestInit) => Promise<Response>;
};

export function useComputationApi({ parameterState, sliderInteracting, fetchFromApi }: Args) {
  // Stable identity key for physics params — changes only when a physics-relevant
  // parameter changes, not when presentation-only params (flux_unit, time_unit) change.
  const { mode, shared, lcFreq, lcTMin, lcTMax, sedTimes, sedNuMin, sedNuMax, sedNumNu, skyTObs, skyNuInput, skyFov, skyNpixel } = parameterState;
  const physicsKey = useMemo(
    () => JSON.stringify(stripPresentationParams(normalizeShared(shared, mode))),
    [shared, mode],
  );

  const fullComputationSpec = useMemo<ComputationSpec>(() => {
    const normalized = JSON.parse(physicsKey);

    if (mode === "lightcurve") {
      return {
        endpoint: "lightcurve",
        payload: { shared: normalized, frequencies_input: lcFreq, t_min: lcTMin, t_max: lcTMax },
      };
    }

    if (mode === "spectrum") {
      return {
        endpoint: "spectrum",
        payload: { shared: normalized, t_snapshots_input: sedTimes, nu_min: sedNuMin, nu_max: sedNuMax, num_nu: sedNumNu },
      };
    }

    return {
      endpoint: "skymap",
      payload: { shared: normalized, animate: false, t_obs: skyTObs, t_min: skyTObs, t_max: skyTObs, n_frames: 1, nu_input: skyNuInput, fov: skyFov, npixel: skyNpixel },
    };
  }, [physicsKey, mode, lcFreq, lcTMin, lcTMax, sedTimes, sedNuMin, sedNuMax, sedNumNu, skyTObs, skyNuInput, skyFov, skyNpixel]);

  const computationSpec = useMemo<ComputationSpec>(
    () => buildInteractiveSpec(fullComputationSpec, ENABLE_INTERACTIVE_DOWNSAMPLE && sliderInteracting),
    [fullComputationSpec, sliderInteracting],
  );

  const postRunRequest = useCallback(
    async (
      spec: ComputationSpec,
      extraPayload: Record<string, unknown>,
      signal?: AbortSignal,
    ): Promise<RunResponse> => {
      const response = await fetchFromApi(spec.endpoint, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          ...spec.payload,
          ...extraPayload,
        }),
        signal,
      });
      if (!response.ok) {
        const details = await response.text();
        throw new Error(details || `HTTP ${response.status}`);
      }
      const data = (await response.json()) as RunResponse;
      return data;
    },
    [fetchFromApi],
  );

  const postFigureRequest = useCallback(
    async (spec: ComputationSpec, controller: AbortController): Promise<RunResponse> => {
      return postRunRequest(spec, {}, controller.signal);
    },
    [postRunRequest],
  );

  return { computationSpec, postFigureRequest };
}
