import { useCallback, useMemo } from "react";
import { ENABLE_INTERACTIVE_DOWNSAMPLE } from "../lib/constants";
import { decodePlotlyFigureBinaryInPlace } from "../lib/plotly/decode";
import type { ComputationSpec, RunResponse } from "../lib/types";
import { compactObservationGroups } from "../lib/utils/obs";
import { buildInteractiveSpec, normalizeShared } from "../lib/utils/snapshot";
import type { ParameterState } from "./useParameterState";

type Args = {
  parameterState: ParameterState;
  sliderInteracting: boolean;
  fetchFromApi: (path: string, init?: RequestInit) => Promise<Response>;
};

export function useComputationApi({ parameterState, sliderInteracting, fetchFromApi }: Args) {
  const fullComputationSpec = useMemo<ComputationSpec>(() => {
    const {
      mode,
      shared,
      lcFreq,
      lcTMin,
      lcTMax,
      lcInstruments,
      lcObsGroups,
      sedTimes,
      sedNuMin,
      sedNuMax,
      sedNumNu,
      sedFreqUnit,
      sedNuFNu,
      sedInstruments,
      sedObsGroups,
      skyAnimate,
      skyTObs,
      skyTMin,
      skyTMax,
      skyNFrames,
      skyNuInput,
      skyFov,
      skyNpixel,
    } = parameterState;

    const normalized = normalizeShared(shared, mode);
    const lcObsPayload = compactObservationGroups(lcObsGroups);
    const sedObsPayload = compactObservationGroups(sedObsGroups);

    if (mode === "lightcurve") {
      return {
        endpoint: "lightcurve",
        payload: {
          shared: normalized,
          frequencies_input: lcFreq,
          t_min: lcTMin,
          t_max: lcTMax,
          selected_instruments: lcInstruments,
          observation_groups: lcObsPayload,
        },
      };
    }

    if (mode === "spectrum") {
      return {
        endpoint: "spectrum",
        payload: {
          shared: normalized,
          t_snapshots_input: sedTimes,
          nu_min: sedNuMin,
          nu_max: sedNuMax,
          num_nu: sedNumNu,
          freq_unit: sedFreqUnit,
          show_nufnu: sedNuFNu,
          selected_instruments: sedInstruments,
          observation_groups: sedObsPayload,
        },
      };
    }

    return {
      endpoint: "skymap",
      payload: {
        shared: normalized,
        animate: skyAnimate,
        t_obs: skyTObs,
        t_min: skyTMin,
        t_max: skyTMax,
        n_frames: skyNFrames,
        nu_input: skyNuInput,
        fov: skyFov,
        npixel: skyNpixel,
      },
    };
  }, [parameterState]);

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
      if (data.figure) {
        decodePlotlyFigureBinaryInPlace(data.figure);
      }
      return data;
    },
    [fetchFromApi],
  );

  const postFigureRequest = useCallback(
    async (spec: ComputationSpec, controller: AbortController): Promise<RunResponse> => {
      return postRunRequest(
        spec,
        {
          include_figure: true,
          include_exports: true,
        },
        controller.signal,
      );
    },
    [postRunRequest],
  );

  return { fullComputationSpec, computationSpec, postRunRequest, postFigureRequest };
}
