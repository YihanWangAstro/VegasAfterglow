import type { ReactNode } from "react";

export type Mode = "lightcurve" | "spectrum" | "skymap";

export type SelectOption = { label: string; value: string };

export type ComputationSpec = {
  endpoint: Mode;
  payload: Record<string, unknown>;
};

export type AxisRange = [number, number];
export type AxisName = "xaxis" | "yaxis" | "yaxis2";
export type AxisRanges = {
  xaxis: AxisRange | null;
  yaxis: AxisRange | null;
  yaxis2: AxisRange | null;
};
export type AxisSignatures = {
  yaxis: string | null;
  yaxis2: string | null;
};

export type ApiEndpoint = {
  key: string;
  healthTarget: string;
};

export type ApiHealthStatus = {
  up: boolean;
  latencyMs: number | null;
  statusCode: number | null;
  region: string | null;
  locationLabel: string | null;
};

export type ObsEntry = {
  label: string;
  fnu: [number, number, number][];   // [x_s, y_cgs, err_cgs]
  fband: [number, number, number][]; // [x_s, y_cgs, err_cgs]
};

export type LcPlotData = {
  flux_unit: string;
  time_unit: string;
  t_min_s: number;
  t_max_s: number;
  times_s: number[];
  pt: {
    freq_hz: number[];
    components: Record<string, number[][]>; // name -> [freq_idx][time_idx]
  } | null;
  bands: {
    nu_min: number;
    nu_max: number;
    name: string | null;
    nu_cen: number;
    components: Record<string, number[]>; // name -> [time_idx]
  }[];
  obs: ObsEntry[];
  instruments: {
    name: string;
    nu_min: number;
    nu_max: number;
    sensitivity: number;
    kind: string;
    t_lo_s: number;
    t_hi_s: number;
  }[];
};

export type SedPlotData = {
  flux_unit: string;
  freq_unit: string;
  nufnu: boolean;
  freq_hz: number[];
  t_snapshots_s: number[];
  components: Record<string, number[][]>; // name -> [t_idx][nu_idx]
  obs: ObsEntry[];
  instruments: {
    name: string;
    nu_min: number;
    nu_max: number;
    sensitivity: number;
    kind: string;
  }[];
};

export type SkymapPlotData = {
  nx: number;
  ny: number;
  frames_b64f32: string[];  // base64-encoded float32 arrays, shape [ny][nx] per frame
  extent_uas: [number, number, number, number];  // [x_min, x_max, y_min, y_max]
  t_obs_s: number[];
  nu_obs_hz: number;
  dx: number;
  dy: number;
  x0: number;
  y0: number;
  z_min: number;
  z_max: number;
};

export type RunResponse = {
  meta?: {
    compute_seconds?: number;
    warnings?: string[];
    [key: string]: unknown;
  };
  plot_data?: LcPlotData | SedPlotData | SkymapPlotData;
  exports?: Record<string, string>;
};

export type OptionsResponse = {
  instruments?: string[];
  instrument_groups?: { label?: string; items?: string[] }[];
  version?: string;
};

export type DownloadKind = "csv" | "json" | "gif";
export type DistanceDriver = "dL" | "z";

export type ObservationGroup = {
  legend: string;
  x_unit: string;
  y_unit: string;
  text: string;
  visible: boolean;
  freq?: string;
};

export type SharedParams = {
  d_L_mpc: number;
  z: number;
  theta_obs: number;
  flux_unit: string;
  time_unit: string;
  jet_type: string;
  theta_c: number;
  E_iso: number;
  Gamma0: number;
  spreading: boolean;
  duration: number;
  k_e: number;
  k_g: number;
  theta_w: number;
  E_iso_w: number;
  Gamma0_w: number;
  medium_type: string;
  n_ism: number;
  A_star: number;
  k_m: number;
  eps_e: number;
  eps_B: number;
  p: number;
  xi_e: number;
  ssc: boolean;
  kn: boolean;
  enable_rvs: boolean;
  eps_e_r: number;
  eps_B_r: number;
  p_r: number;
  xi_e_r: number;
  rvs_ssc: boolean;
  rvs_kn: boolean;
  num_t: number;
  res_phi: number;
  res_theta: number;
  res_t: number;
};

export type UiSnapshot = {
  mode: Mode;
  distance_linked: boolean;
  distance_driver: DistanceDriver;
  shared: SharedParams;
  lightcurve: {
    frequencies_input: string;
    t_min: number;
    t_max: number;
    selected_instruments: string[];
    observation_groups: ObservationGroup[];
  };
  spectrum: {
    t_snapshots_input: string;
    nu_min: number;
    nu_max: number;
    num_nu: number;
    freq_unit: string;
    show_nufnu: boolean;
    selected_instruments: string[];
    observation_groups: ObservationGroup[];
  };
  skymap: {
    animate: boolean;
    t_obs: number;
    t_min: number;
    t_max: number;
    n_frames: number;
    nu_input: string;
    fov: number;
    npixel: number;
  };
};

export type BookmarkEntry = {
  id: string;
  name: string;
  created_at: string;
  snapshot: UiSnapshot;
};

export type SharedLogSliderSpec = {
  kind: "log";
  key: keyof SharedParams;
  label: ReactNode;
  minExp: number;
  maxExp: number;
  step: number;
  defaultExp: number;
};

export type SharedLinearSliderSpec = {
  kind: "linear";
  key: keyof SharedParams;
  label: ReactNode;
  min: number;
  max: number;
  step: number;
  decimals?: number;
};

export type SharedSliderSpec = SharedLogSliderSpec | SharedLinearSliderSpec;

export type ModeLogSliderSpec = {
  key: string;
  label: ReactNode;
  minExp: number;
  maxExp: number;
  step: number;
  value: number;
  defaultExp: number;
  onChange: (next: number) => void;
};

export type InstrumentGroup = {
  label: string;
  items: string[];
};

export type DefaultsResponse = {
  shared?: Partial<SharedParams>;
  lightcurve?: {
    frequencies_input?: string;
    t_min?: number;
    t_max?: number;
    selected_instruments?: string[];
    observation_groups?: ObservationGroup[];
  };
  spectrum?: {
    t_snapshots_input?: string;
    nu_min?: number;
    nu_max?: number;
    num_nu?: number;
    freq_unit?: string;
    show_nufnu?: boolean;
    selected_instruments?: string[];
    observation_groups?: ObservationGroup[];
  };
  skymap?: {
    animate?: boolean;
    t_obs?: number;
    t_min?: number;
    t_max?: number;
    n_frames?: number;
    nu_input?: string;
    fov?: number;
    npixel?: number;
  };
};
