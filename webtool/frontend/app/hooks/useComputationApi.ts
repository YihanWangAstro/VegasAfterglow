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
  // Physics-only spec — excludes presentation params (flux_unit, time_unit,
  // freq_unit, show_nufnu) so changing display units doesn't trigger a backend
  // re-computation.  The frontend handles all unit conversion locally.
  const fullComputationSpec = useMemo<ComputationSpec>(() => {
    const {
      mode,
      shared,
      lcFreq,
      lcTMin,
      lcTMax,
      sedTimes,
      sedNuMin,
      sedNuMax,
      sedNumNu,
      skyTObs,
      skyNuInput,
      skyFov,
      skyNpixel,
    } = parameterState;

    const normalized = stripPresentationParams(normalizeShared(shared, mode));

    if (mode === "lightcurve") {
      return {
        endpoint: "lightcurve",
        payload: {
          shared: normalized,
          frequencies_input: lcFreq,
          t_min: lcTMin,
          t_max: lcTMax,
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
        },
      };
    }

    return {
      endpoint: "skymap",
      payload: {
        shared: normalized,
        animate: false,
        t_obs: skyTObs,
        t_min: skyTObs,
        t_max: skyTObs,
        n_frames: 1,
        nu_input: skyNuInput,
        fov: skyFov,
        npixel: skyNpixel,
      },
    };
  }, [
    parameterState.mode,
    parameterState.shared.E_iso, parameterState.shared.Gamma0, parameterState.shared.theta_c,
    parameterState.shared.theta_w, parameterState.shared.E_iso_w, parameterState.shared.Gamma0_w,
    parameterState.shared.k_e, parameterState.shared.k_g,
    parameterState.shared.jet_type, parameterState.shared.medium_type,
    parameterState.shared.n_ism, parameterState.shared.A_star,
    parameterState.shared.eps_e, parameterState.shared.eps_B, parameterState.shared.p, parameterState.shared.xi_e,
    parameterState.shared.d_L_mpc, parameterState.shared.z, parameterState.shared.theta_obs,
    parameterState.shared.enable_rvs, parameterState.shared.duration,
    parameterState.shared.eps_e_r, parameterState.shared.eps_B_r, parameterState.shared.p_r, parameterState.shared.xi_e_r,
    parameterState.shared.spreading,
    parameterState.shared.ssc, parameterState.shared.kn, parameterState.shared.rvs_ssc, parameterState.shared.rvs_kn,
    parameterState.shared.num_t, parameterState.shared.k_m,
    parameterState.shared.res_phi, parameterState.shared.res_theta, parameterState.shared.res_t,
    parameterState.lcFreq, parameterState.lcTMin, parameterState.lcTMax,
    parameterState.sedTimes, parameterState.sedNuMin, parameterState.sedNuMax, parameterState.sedNumNu,
    parameterState.skyTObs, parameterState.skyNuInput, parameterState.skyFov, parameterState.skyNpixel,
  ]);

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
