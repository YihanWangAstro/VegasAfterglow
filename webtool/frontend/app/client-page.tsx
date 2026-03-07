"use client";

import dynamic from "next/dynamic";
import {
  memo,
  useCallback,
  useDeferredValue,
  useEffect,
  useMemo,
  useRef,
  useState,
  useTransition,
  type CSSProperties,
  type KeyboardEvent as ReactKeyboardEvent,
  type PointerEvent as ReactPointerEvent,
  type ReactNode,
} from "react";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });
const PLOT_CONFIG = { responsive: true, displaylogo: false } as const;

type Mode = "lightcurve" | "spectrum" | "skymap";
type ComputationSpec = {
  endpoint: Mode;
  payload: Record<string, unknown>;
};
type AxisRange = [number, number];
type AxisName = "xaxis" | "yaxis" | "yaxis2";
type AxisRanges = {
  xaxis: AxisRange | null;
  yaxis: AxisRange | null;
  yaxis2: AxisRange | null;
};
type AxisSignatures = {
  yaxis: string | null;
  yaxis2: string | null;
};
type ApiEndpoint = {
  key: string;
  label: string;
  healthTarget: string;
};
type ApiHealthStatus = {
  up: boolean;
  latencyMs: number | null;
  statusCode: number | null;
  region: string | null;
  locationLabel: string | null;
  backend: string | null;
};

type RunResponse = {
  meta?: {
    compute_seconds?: number;
    warnings?: string[];
    [key: string]: unknown;
  };
  figure?: {
    data?: unknown[];
    layout?: Record<string, unknown>;
    frames?: unknown[];
  };
  exports?: Record<string, string>;
};

type OptionsResponse = {
  instruments?: string[];
  instrument_groups?: { label?: string; items?: string[] }[];
  version?: string;
};
type DownloadKind = "csv" | "json" | "gif";

type ObservationGroup = {
  legend: string;
  x_unit: string;
  y_unit: string;
  text: string;
  visible: boolean;
};
type DistanceDriver = "dL" | "z";

type UiSnapshot = {
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

type BookmarkEntry = {
  id: string;
  name: string;
  created_at: string;
  snapshot: UiSnapshot;
};

const FALLBACK_SHARED = {
  d_L_mpc: 100,
  z: 0.022,
  theta_obs: 0,
  flux_unit: "mJy",
  time_unit: "s",
  jet_type: "Top-hat",
  theta_c: 0.1,
  E_iso: 1e52,
  Gamma0: 300,
  spreading: false,
  duration: 1.0,
  k_e: 2.0,
  k_g: 2.0,
  theta_w: 0.3,
  E_iso_w: 1e51,
  Gamma0_w: 100,
  medium_type: "ISM",
  n_ism: 1.0,
  A_star: 0.1,
  k_m: 2.0,
  eps_e: 0.1,
  eps_B: 1e-3,
  p: 2.3,
  xi_e: 1.0,
  ssc: false,
  kn: false,
  enable_rvs: false,
  eps_e_r: 0.1,
  eps_B_r: 1e-3,
  p_r: 2.3,
  xi_e_r: 1.0,
  rvs_ssc: false,
  rvs_kn: false,
  num_t: 100,
  res_phi: 0.1,
  res_theta: 0.25,
  res_t: 10,
};

type SharedParams = typeof FALLBACK_SHARED;
type InstrumentGroup = {
  label: string;
  items: string[];
};

const Y_UNIT_OPTIONS = ["mJy", "Jy", "uJy", "cgs", "AB mag", "erg/cm²/s"];
const TIME_UNIT_OPTIONS = ["s", "day", "hr", "min"];
const FREQ_UNIT_OPTIONS = ["Hz", "GHz", "keV", "MeV"];
const SKY_DEFAULT_N_FRAMES = 15;
const SKY_MAX_N_FRAMES = 30;
const SKY_MAX_PIXEL_ANIMATE = 512;
const SKY_MAX_PIXEL_STATIC = 1024;
const SKY_PIXEL_OPTIONS = [64, 128, 256, 512, 1024] as const;
const ENABLE_INTERACTIVE_DOWNSAMPLE = false;
const INTERACTIVE_LIGHTCURVE_NUM_T_MAX = 120;
const INTERACTIVE_SPECTRUM_NUM_NU_MAX = 120;
const INTERACTIVE_SKY_PIXEL_MAX_STATIC = 256;
const INTERACTIVE_SKY_PIXEL_MAX_ANIMATE = 256;
const INTERACTIVE_SKY_FRAMES_MAX = 8;
const AUTO_RUN_DEBOUNCE_IDLE_MS = 10;
const AUTO_RUN_DEBOUNCE_SLIDING_MS = 20;
const SLIDER_COMMIT_INTERVAL_MS = 20;
const SPECTRUM_TEXT_COMMIT_DEBOUNCE_MS = 220;
const API_STATUS_REFRESH_DIRECT_MS = 2000;
const AXIS_EPS = 1e-6;
const H0_KM_S_MPC = 67.4;
const C_KM_S = 299792.458;
const MODE_OPTIONS: { value: Mode; label: string }[] = [
  { value: "lightcurve", label: "Light Curve" },
  { value: "spectrum", label: "Spectrum" },
  { value: "skymap", label: "Sky Image" },
];
const INSTRUMENT_TYPE_BY_NAME: Record<string, string> = {
  VLA: "Radio",
  ALMA: "Radio",
  MeerKAT: "Radio",
  ngVLA: "Radio",
  "Rubin/LSST": "Optical/IR",
  JWST: "Optical/IR",
  WFST: "Optical/IR",
  "SVOM/VT": "Optical/IR",
  "Swift/XRT": "X-ray",
  Chandra: "X-ray",
  "EP/WXT": "X-ray",
  "EP/FXT": "X-ray",
  "SVOM/MXT": "X-ray",
  "SVOM/ECLAIRs": "X-ray",
  "Swift/BAT": "Gamma-ray",
  "Fermi/GBM": "Gamma-ray",
  "Fermi/LAT": "Gamma-ray",
  CTA: "Gamma-ray",
};
const INSTRUMENT_GROUP_ORDER = ["Radio", "Optical/IR", "X-ray", "Gamma-ray", "Other"] as const;
const OBS_STORAGE_LC_KEY = "afterglow:webtool:obs:lightcurve";
const OBS_STORAGE_SED_KEY = "afterglow:webtool:obs:spectrum";
const BOOKMARKS_STORAGE_KEY = "afterglow:webtool:bookmarks:v1";
const MAX_BOOKMARKS = 30;
const URL_STATE_PARAM = "state";
const URL_STATE_VERSION = 1;
const URL_STATE_MAX_CHARS = 3500;
const FREQ_HELP_TEXT = [
  "Accepts Hz values (1e9), unit suffixes (1GHz, 1keV), filter names, instrument bands, or custom ranges ([0.3keV,10keV]).",
  "",
  "Filters: U B V R I J H Ks | u b v uvw1 uvm2 uvw2 | g r i z | w VT_B VT_R |",
  "F225W F275W F336W F438W F475W F555W F606W F625W F775W F814W F850LP F105W F110W F125W F140W F160W",
  "",
  "Bands: XRT BAT FXT WXT MXT ECLAIRs LAT GBM",
  "",
  "Units: Hz kHz MHz GHz | eV keV MeV GeV",
].join("\n");
const JET_TYPE_OPTIONS: SharedParams["jet_type"][] = ["Top-hat", "Gaussian", "Power-law", "Two-component"];
const MEDIUM_TYPE_OPTIONS: { value: SharedParams["medium_type"]; label: string }[] = [
  { value: "ISM", label: "ISM" },
  { value: "Wind", label: "Wind" },
  { value: "Wind bubble", label: "Wind+ISM" },
];
const DEFAULT_APP_VERSION = "2.0.1";
const BIBTEX_TEXT = `@ARTICLE{2026JHEAp..5000490W,
       author = {{Wang}, Yihan and {Chen}, Connery and {Zhang}, Bing},
        title = "{VegasAfterglow: A high-performance framework for gamma-ray burst afterglows}",
      journal = {Journal of High Energy Astrophysics},
         year = 2026,
        month = feb,
       volume = {50},
          eid = {100490},
        pages = {100490},
          doi = {10.1016/j.jheap.2025.100490},
archivePrefix = {arXiv},
       eprint = {2507.10829},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2026JHEAp..5000490W},
}

@book{Zhang2018,
  author    = {Zhang, Bing},
  title     = {{The Physics of Gamma-Ray Bursts}},
  publisher = {Cambridge University Press},
  year      = {2018},
  doi       = {10.1017/9781139226530}
}`;

const PlotFigure = memo(function PlotFigure({
  figure,
  mode,
  onRelayout,
}: {
  figure: NonNullable<RunResponse["figure"]>;
  mode: Mode;
  onRelayout: (event: Record<string, unknown>) => void;
}) {
  const graphRef = useRef<any>(null);
  const [graphRevision, setGraphRevision] = useState(0);

  const bindGraphRef = useCallback((_nextFigure: unknown, graphDiv: unknown) => {
    const nextGraph = graphDiv as any;
    if (graphRef.current !== nextGraph) {
      graphRef.current = nextGraph;
      setGraphRevision((v) => v + 1);
    }
  }, []);

  useEffect(() => {
    const hasFrames = Array.isArray(figure.frames) && figure.frames.length > 1;
    if (mode !== "skymap" || !hasFrames) return;

    const graphDiv = graphRef.current as any;
    if (!graphDiv?.animate) return;
    const frameNames = (figure.frames ?? [])
      .map((frame) => (frame && typeof frame === "object" ? (frame as { name?: unknown }).name : undefined))
      .filter((name): name is string => typeof name === "string" && name.length > 0);
    if (frameNames.length === 0) return;

    let cancelled = false;
    let starter: number | null = null;
    let frameIndex = 0;
    const playArgs: Record<string, unknown> = {
      frame: { duration: 300, redraw: true },
      transition: { duration: 0 },
      mode: "immediate",
    };
    const stopArgs = { frame: { duration: 0, redraw: false }, mode: "immediate" };
    const stopAnimationSafely = () => {
      try {
        const maybePromise = graphDiv.animate?.([null], stopArgs);
        if (maybePromise && typeof maybePromise.then === "function") {
          void maybePromise.catch(() => undefined);
        }
      } catch {
        // Ignore transient Plotly stop failures.
      }
    };

    const playSequence = async () => {
      if (cancelled) return;

      const nextFrame = frameNames[frameIndex % frameNames.length];
      frameIndex += 1;
      try {
        await graphDiv.animate([nextFrame], playArgs);
        if (cancelled) return;
        starter = window.setTimeout(() => {
          void playSequence();
        }, 0);
      } catch {
        // Transient Plotly frame registration race; retry shortly.
        if (cancelled) return;
        starter = window.setTimeout(() => {
          void playSequence();
        }, 60);
      }
    };

    // Give Plotly one tick to register frames before first animation call.
    starter = window.setTimeout(() => {
      void playSequence();
    }, 80);

    return () => {
      cancelled = true;
      if (starter !== null) {
        window.clearTimeout(starter);
      }
      stopAnimationSafely();
    };
  }, [figure.frames, graphRevision, mode]);

  return (
    <div className={`plot-wrap${mode === "skymap" ? " plot-wrap-square" : ""}`}>
      <Plot
        data={figure.data as any}
        layout={(figure.layout ?? {}) as any}
        frames={(figure.frames ?? []) as any}
        config={PLOT_CONFIG as any}
        onRelayout={(event) => onRelayout((event ?? {}) as Record<string, unknown>)}
        onInitialized={bindGraphRef}
        onUpdate={bindGraphRef}
        style={{ width: "100%", height: "100%" }}
        useResizeHandler
      />
    </div>
  );
});

function stripTrailingSlash(url: string): string {
  return url.replace(/\/+$/, "");
}

function apiServerFromTarget(target: string): string {
  try {
    return new URL(target, "http://localhost").host;
  } catch {
    return target;
  }
}

function backendFromHealthPayload(payload: unknown): string | null {
  if (!payload || typeof payload !== "object") return null;
  const record = payload as Record<string, unknown>;
  const revision = typeof record.revision === "string" ? record.revision.trim() : "";
  const serviceName = typeof record.service_name === "string" ? record.service_name.trim() : "";
  const rawInstance = typeof record.instance === "string" ? record.instance.trim() : "";
  const instance =
    rawInstance && rawInstance !== "localhost" && rawInstance !== "127.0.0.1" ? rawInstance : "";

  const parts: string[] = [];
  if (revision) {
    parts.push(revision);
  } else if (serviceName) {
    parts.push(serviceName);
  }
  if (instance) parts.push(instance);
  return parts.length > 0 ? parts.join(" | ") : null;
}

function regionFromHealthPayload(payload: unknown): string | null {
  if (!payload || typeof payload !== "object") return null;
  const record = payload as Record<string, unknown>;
  const region = typeof record.region === "string" ? record.region.trim() : "";
  return region || null;
}

function locationLabelFromHealthPayload(payload: unknown): string | null {
  if (!payload || typeof payload !== "object") return null;
  const record = payload as Record<string, unknown>;
  const label = typeof record.location_label === "string" ? record.location_label.trim() : "";
  return label || null;
}

function safeLog10(value: number, fallback: number): number {
  return Number.isFinite(value) && value > 0 ? Math.log10(value) : fallback;
}

function formatParamValueNode(value: number): ReactNode {
  if (!Number.isFinite(value)) return String(value);
  if (value === 0) return "0";
  const abs = Math.abs(value);
  const exponent = Math.floor(Math.log10(abs));
  const mantissa = value / Math.pow(10, exponent);
  const mantissaText = mantissa.toFixed(2).replace(/\.?0+$/, "");
  if (exponent === 0) return mantissaText;
  if (mantissaText === "1") {
    return (
      <>
        10<sup>{exponent}</sup>
      </>
    );
  }
  if (mantissaText === "-1") {
    return (
      <>
        -10<sup>{exponent}</sup>
      </>
    );
  }
  return (
    <>
      {mantissaText} x 10<sup>{exponent}</sup>
    </>
  );
}

function clampRangeValue(value: number, min: number, max: number): number {
  return Math.max(min, Math.min(max, value));
}

function formatPowerHoverValue(value: unknown): string {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return String(value ?? "");
  if (numeric === 0) return "0";
  const abs = Math.abs(numeric);
  const exponent = Math.floor(Math.log10(abs));
  const mantissa = numeric / Math.pow(10, exponent);
  const mantissaText = Number(mantissa.toPrecision(3)).toString();
  if (exponent === 0) return mantissaText;
  return `${mantissaText}\u00d710<sup>${exponent}</sup>`;
}

function remapScientificHoverTemplate(trace: Record<string, unknown>): Record<string, unknown> {
  const hovertemplate = typeof trace.hovertemplate === "string" ? trace.hovertemplate : "";
  if (!hovertemplate) return trace;

  let needsX = false;
  let needsY = false;
  let needsZ = false;
  let rewritten = hovertemplate;

  rewritten = rewritten.replace(/%\{x:[^}]*e\}/g, () => {
    needsX = true;
    return "%{customdata[0]}";
  });
  rewritten = rewritten.replace(/%\{y:[^}]*e\}/g, () => {
    needsY = true;
    return "%{customdata[1]}";
  });
  rewritten = rewritten.replace(/%\{z:[^}]*e\}/g, () => {
    needsZ = true;
    return "%{customdata[2]}";
  });

  if (!needsX && !needsY && !needsZ) return trace;

  const xSeries = Array.isArray(trace.x) ? trace.x : null;
  const ySeries = Array.isArray(trace.y) ? trace.y : null;
  const zSeries = Array.isArray(trace.z) ? trace.z : null;
  const zIs1D = zSeries ? zSeries.every((item) => !Array.isArray(item)) : false;

  const length = Math.max(
    needsX && xSeries ? xSeries.length : 0,
    needsY && ySeries ? ySeries.length : 0,
    needsZ && zIs1D && zSeries ? zSeries.length : 0,
  );
  if (length === 0) return trace;

  const customdata = Array.from({ length }, (_, idx) => [
    needsX && xSeries ? formatPowerHoverValue(xSeries[idx]) : "",
    needsY && ySeries ? formatPowerHoverValue(ySeries[idx]) : "",
    needsZ && zIs1D && zSeries ? formatPowerHoverValue(zSeries[idx]) : "",
  ]);

  return {
    ...trace,
    customdata,
    hovertemplate: rewritten,
  };
}

type RangeKeyAction = "dec" | "inc" | "min" | "max" | null;

function rangeKeyAction(event: ReactKeyboardEvent<HTMLInputElement>): RangeKeyAction {
  const key = event.key;
  const code = event.code;
  if (key === "ArrowLeft" || key === "ArrowDown" || key === "Left" || key === "Down") return "dec";
  if (key === "ArrowRight" || key === "ArrowUp" || key === "Right" || key === "Up") return "inc";
  if (key === "Home") return "min";
  if (key === "End") return "max";
  if (code === "ArrowLeft" || code === "ArrowDown") return "dec";
  if (code === "ArrowRight" || code === "ArrowUp") return "inc";
  if (code === "Home") return "min";
  if (code === "End") return "max";
  if (key === "UIKeyInputLeftArrow" || key === "UIKeyInputDownArrow") return "dec";
  if (key === "UIKeyInputRightArrow" || key === "UIKeyInputUpArrow") return "inc";

  const nativeEvent = event.nativeEvent as KeyboardEvent & { keyIdentifier?: string };
  const keyIdentifier = nativeEvent.keyIdentifier;
  if (keyIdentifier === "Left" || keyIdentifier === "Down") return "dec";
  if (keyIdentifier === "Right" || keyIdentifier === "Up") return "inc";
  if (keyIdentifier === "Home") return "min";
  if (keyIdentifier === "End") return "max";

  const keyCode = Number(nativeEvent.keyCode);
  if (keyCode === 37 || keyCode === 40 || keyCode === 63234 || keyCode === 63233) return "dec";
  if (keyCode === 39 || keyCode === 38 || keyCode === 63235 || keyCode === 63232) return "inc";
  if (keyCode === 36 || keyCode === 63273) return "min";
  if (keyCode === 35 || keyCode === 63275) return "max";
  return null;
}

function stepRangeValue(
  event: ReactKeyboardEvent<HTMLInputElement>,
  current: number,
  min: number,
  max: number,
  step: number,
  apply: (next: number) => void,
): void {
  let next: number | null = null;
  switch (rangeKeyAction(event)) {
    case "dec":
      next = current - step;
      break;
    case "inc":
      next = current + step;
      break;
    case "min":
      next = min;
      break;
    case "max":
      next = max;
      break;
    default:
      return;
  }

  event.preventDefault();
  const snapped = clampRangeValue(Number(next.toFixed(10)), min, max);
  apply(snapped);
}

function sliderFillStyle(value: number, min: number, max: number): CSSProperties {
  const span = max - min;
  const raw = span > 0 ? ((value - min) / span) * 100 : 0;
  const bounded = Math.max(0, Math.min(100, raw));
  return { "--fill-percent": `${bounded}%` } as CSSProperties;
}

function legendFontSizeForWidth(plotWidthPx: number): number {
  if (plotWidthPx <= 300) return 6;
  if (plotWidthPx <= 380) return 7;
  if (plotWidthPx <= 480) return 8;
  if (plotWidthPx <= 640) return 9;
  if (plotWidthPx <= 860) return 10;
  return 11;
}

function downloadText(content: string, filename: string, mime: string): void {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}

function downloadBase64(content: string, filename: string, mime: string): void {
  const raw = content.includes(",") ? content.split(",").pop() ?? "" : content;
  const binary = window.atob(raw);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i += 1) {
    bytes[i] = binary.charCodeAt(i);
  }
  const blob = new Blob([bytes], { type: mime });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}

function defaultObsGroup(isLc: boolean, idx: number): ObservationGroup {
  return {
    legend: `data ${idx}`,
    x_unit: isLc ? "day" : "Hz",
    y_unit: "mJy",
    text: "",
    visible: true,
  };
}

function parseStoredObsGroups(raw: string | null, isLc: boolean): ObservationGroup[] | null {
  if (!raw) return null;
  try {
    const parsed = JSON.parse(raw);
    if (!Array.isArray(parsed)) return null;

    const normalized = parsed
      .map((item, idx) => {
        if (!item || typeof item !== "object") return null;
        const base = defaultObsGroup(isLc, idx + 1);
        const group = item as Partial<ObservationGroup>;
        return {
          legend: typeof group.legend === "string" ? group.legend : base.legend,
          x_unit: typeof group.x_unit === "string" ? group.x_unit : base.x_unit,
          y_unit: typeof group.y_unit === "string" ? group.y_unit : base.y_unit,
          text: typeof group.text === "string" ? group.text : "",
          visible: typeof group.visible === "boolean" ? group.visible : true,
        };
      })
      .filter((group): group is ObservationGroup => group !== null);

    return normalized;
  } catch {
    return null;
  }
}

function compactObservationGroups(groups: ObservationGroup[]): ObservationGroup[] {
  return groups
    .filter((group) => group.visible && group.text.trim().length > 0)
    .map((group) => ({
      legend: group.legend.trim() || "data",
      x_unit: group.x_unit,
      y_unit: group.y_unit,
      text: group.text.trim(),
      visible: true,
    }));
}

function cloneSnapshot(snapshot: UiSnapshot): UiSnapshot {
  return JSON.parse(JSON.stringify(snapshot)) as UiSnapshot;
}

function asFiniteNumber(value: unknown, fallback: number): number {
  const next = Number(value);
  return Number.isFinite(next) ? next : fallback;
}

function clampNumber(value: number, min: number, max: number): number {
  return Math.max(min, Math.min(max, value));
}

function asString(value: unknown, fallback: string): string {
  return typeof value === "string" ? value : fallback;
}

function asBoolean(value: unknown, fallback: boolean): boolean {
  return typeof value === "boolean" ? value : fallback;
}

function luminosityDistanceMpcFromRedshift(z: number): number {
  const zSafe = Math.max(0, z);
  return (C_KM_S / H0_KM_S_MPC) * zSafe * (1 + zSafe);
}

function redshiftFromLuminosityDistanceMpc(dLmpc: number): number {
  const dSafe = Math.max(0, dLmpc);
  const x = (dSafe * H0_KM_S_MPC) / C_KM_S;
  return (-1 + Math.sqrt(1 + 4 * x)) / 2;
}

function normalizeDistanceMpc(value: number): number {
  if (!Number.isFinite(value)) return 0;
  return Math.max(0, Math.trunc(value));
}

function roundToSignificant(value: number, digits: number): number {
  if (!Number.isFinite(value) || value === 0) return value;
  const abs = Math.abs(value);
  const exponent = Math.floor(Math.log10(abs));
  const scale = Math.pow(10, exponent - digits + 1);
  return Math.round(value / scale) * scale;
}

function normalizeRedshift(value: number): number {
  if (!Number.isFinite(value)) return 0;
  return Math.max(0, roundToSignificant(value, 3));
}

function nearlyEqualRelative(a: number, b: number, relTol = 1e-6): boolean {
  return Math.abs(a - b) <= relTol * Math.max(1, Math.abs(a), Math.abs(b));
}

function asStringArray(value: unknown): string[] {
  if (!Array.isArray(value)) return [];
  return value.filter((item): item is string => typeof item === "string");
}

function normalizeSharedParams(value: unknown): SharedParams {
  const defaults = FALLBACK_SHARED;
  const next: SharedParams = { ...defaults };
  if (!value || typeof value !== "object") return next;
  const record = value as Record<string, unknown>;
  const writable = next as Record<string, unknown>;

  for (const key of Object.keys(defaults) as (keyof SharedParams)[]) {
    const defaultValue = defaults[key];
    const candidate = record[key as string];
    if (typeof defaultValue === "number") {
      const numeric = Number(candidate);
      if (Number.isFinite(numeric)) {
        writable[key as string] = numeric;
      }
      continue;
    }
    if (typeof defaultValue === "boolean") {
      if (typeof candidate === "boolean") {
        writable[key as string] = candidate;
      }
      continue;
    }
    if (typeof defaultValue === "string" && typeof candidate === "string") {
      writable[key as string] = candidate;
    }
  }

  return next;
}

function normalizeObsGroupArray(value: unknown, isLc: boolean): ObservationGroup[] {
  if (!Array.isArray(value)) return [];
  return value
    .map((item, idx) => {
      if (!item || typeof item !== "object") return null;
      const base = defaultObsGroup(isLc, idx + 1);
      const group = item as Partial<ObservationGroup>;
      return {
        legend: asString(group.legend, base.legend),
        x_unit: asString(group.x_unit, base.x_unit),
        y_unit: asString(group.y_unit, base.y_unit),
        text: asString(group.text, ""),
        visible: asBoolean(group.visible, true),
      } satisfies ObservationGroup;
    })
    .filter((group): group is ObservationGroup => group !== null);
}

function defaultUiSnapshot(): UiSnapshot {
  return {
    mode: "lightcurve",
    distance_linked: true,
    distance_driver: "dL",
    shared: { ...FALLBACK_SHARED },
    lightcurve: {
      frequencies_input: "1e9, R, 1keV",
      t_min: 1,
      t_max: 1e8,
      selected_instruments: [],
      observation_groups: [],
    },
    spectrum: {
      t_snapshots_input: "1e3, 1e4, 1e5, 1e6",
      nu_min: 1e8,
      nu_max: 1e20,
      num_nu: 200,
      freq_unit: "Hz",
      show_nufnu: false,
      selected_instruments: [],
      observation_groups: [],
    },
    skymap: {
      animate: false,
      t_obs: 1e6,
      t_min: 1e4,
      t_max: 1e7,
      n_frames: SKY_DEFAULT_N_FRAMES,
      nu_input: "1e9",
      fov: 500,
      npixel: 256,
    },
  };
}

function normalizeUiSnapshot(value: unknown): UiSnapshot | null {
  if (!value || typeof value !== "object") return null;
  const defaults = defaultUiSnapshot();
  const record = value as Record<string, unknown>;

  const modeCandidate = record.mode;
  const mode: Mode =
    modeCandidate === "lightcurve" || modeCandidate === "spectrum" || modeCandidate === "skymap"
      ? modeCandidate
      : defaults.mode;
  const distanceLinked = asBoolean(record.distance_linked, defaults.distance_linked);
  const distanceDriver =
    record.distance_driver === "dL" || record.distance_driver === "z"
      ? record.distance_driver
      : defaults.distance_driver;

  const lcRaw = record.lightcurve && typeof record.lightcurve === "object" ? (record.lightcurve as Record<string, unknown>) : {};
  const sedRaw = record.spectrum && typeof record.spectrum === "object" ? (record.spectrum as Record<string, unknown>) : {};
  const skyRaw = record.skymap && typeof record.skymap === "object" ? (record.skymap as Record<string, unknown>) : {};

  return {
    mode,
    distance_linked: distanceLinked,
    distance_driver: distanceDriver,
    shared: normalizeSharedParams(record.shared),
    lightcurve: {
      frequencies_input: asString(lcRaw.frequencies_input, defaults.lightcurve.frequencies_input),
      t_min: asFiniteNumber(lcRaw.t_min, defaults.lightcurve.t_min),
      t_max: asFiniteNumber(lcRaw.t_max, defaults.lightcurve.t_max),
      selected_instruments: asStringArray(lcRaw.selected_instruments),
      observation_groups: normalizeObsGroupArray(lcRaw.observation_groups, true),
    },
    spectrum: {
      t_snapshots_input: asString(sedRaw.t_snapshots_input, defaults.spectrum.t_snapshots_input),
      nu_min: asFiniteNumber(sedRaw.nu_min, defaults.spectrum.nu_min),
      nu_max: asFiniteNumber(sedRaw.nu_max, defaults.spectrum.nu_max),
      num_nu: clampNumber(asFiniteNumber(sedRaw.num_nu, defaults.spectrum.num_nu), 50, 2000),
      freq_unit: asString(sedRaw.freq_unit, defaults.spectrum.freq_unit),
      show_nufnu: asBoolean(sedRaw.show_nufnu, defaults.spectrum.show_nufnu),
      selected_instruments: asStringArray(sedRaw.selected_instruments),
      observation_groups: normalizeObsGroupArray(sedRaw.observation_groups, false),
    },
    skymap: {
      animate: asBoolean(skyRaw.animate, defaults.skymap.animate),
      t_obs: asFiniteNumber(skyRaw.t_obs, defaults.skymap.t_obs),
      t_min: asFiniteNumber(skyRaw.t_min, defaults.skymap.t_min),
      t_max: asFiniteNumber(skyRaw.t_max, defaults.skymap.t_max),
      n_frames: clampNumber(asFiniteNumber(skyRaw.n_frames, defaults.skymap.n_frames), 3, SKY_MAX_N_FRAMES),
      nu_input: asString(skyRaw.nu_input, defaults.skymap.nu_input),
      fov: asFiniteNumber(skyRaw.fov, defaults.skymap.fov),
      npixel: clampNumber(asFiniteNumber(skyRaw.npixel, defaults.skymap.npixel), 64, SKY_MAX_PIXEL_STATIC),
    },
  };
}

function makeShareableSnapshot(snapshot: UiSnapshot): UiSnapshot {
  const shared = cloneSnapshot(snapshot);
  // Keep URL state compact; uploaded observation tables stay local.
  shared.lightcurve.observation_groups = [];
  shared.spectrum.observation_groups = [];
  return shared;
}

function encodeUrlState(snapshot: UiSnapshot): string | null {
  try {
    const payload = JSON.stringify({ v: URL_STATE_VERSION, snapshot });
    const bytes = new TextEncoder().encode(payload);
    let binary = "";
    for (let i = 0; i < bytes.length; i += 1) {
      binary += String.fromCharCode(bytes[i]);
    }
    const encoded = window.btoa(binary).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/g, "");
    if (encoded.length > URL_STATE_MAX_CHARS) return null;
    return encoded;
  } catch {
    return null;
  }
}

function decodeUrlState(raw: string): UiSnapshot | null {
  try {
    const normalized = raw.replace(/-/g, "+").replace(/_/g, "/");
    const padded = normalized + "=".repeat((4 - (normalized.length % 4)) % 4);
    const binary = window.atob(padded);
    const bytes = Uint8Array.from(binary, (char) => char.charCodeAt(0));
    const decoded = new TextDecoder().decode(bytes);
    const parsed = JSON.parse(decoded);
    if (!parsed || typeof parsed !== "object") return null;
    const payload = parsed as { v?: unknown; snapshot?: unknown };
    if (payload.v !== undefined && payload.v !== URL_STATE_VERSION) return null;
    return normalizeUiSnapshot(payload.snapshot ?? parsed);
  } catch {
    return null;
  }
}

function parseStoredBookmarks(raw: string | null): BookmarkEntry[] | null {
  if (!raw) return null;
  try {
    const parsed = JSON.parse(raw);
    if (!Array.isArray(parsed)) return null;
    const normalized = parsed
      .map((item) => {
        if (!item || typeof item !== "object") return null;
        const candidate = item as Partial<BookmarkEntry>;
        if (typeof candidate.id !== "string" || typeof candidate.name !== "string") return null;
        if (!candidate.snapshot || typeof candidate.snapshot !== "object") return null;
        const snapshot = normalizeUiSnapshot(candidate.snapshot);
        if (!snapshot) return null;
        const createdAt = typeof candidate.created_at === "string" ? candidate.created_at : new Date().toISOString();
        return {
          id: candidate.id,
          name: candidate.name,
          created_at: createdAt,
          snapshot,
        } satisfies BookmarkEntry;
      })
      .filter((item): item is BookmarkEntry => item !== null);
    return normalized;
  } catch {
    return null;
  }
}

function buildComputationSpecFromSnapshot(snapshot: UiSnapshot, mode: Mode): ComputationSpec {
  const normalized = normalizeShared(snapshot.shared, mode);

  if (mode === "lightcurve") {
    return {
      endpoint: "lightcurve",
      payload: {
        shared: normalized,
        frequencies_input: snapshot.lightcurve.frequencies_input,
        t_min: snapshot.lightcurve.t_min,
        t_max: snapshot.lightcurve.t_max,
        selected_instruments: [],
        observation_groups: [],
      },
    };
  }

  if (mode === "spectrum") {
    return {
      endpoint: "spectrum",
      payload: {
        shared: normalized,
        t_snapshots_input: snapshot.spectrum.t_snapshots_input,
        nu_min: snapshot.spectrum.nu_min,
        nu_max: snapshot.spectrum.nu_max,
        num_nu: snapshot.spectrum.num_nu,
        freq_unit: snapshot.spectrum.freq_unit,
        show_nufnu: snapshot.spectrum.show_nufnu,
        selected_instruments: [],
        observation_groups: [],
      },
    };
  }

  return {
    endpoint: "skymap",
    payload: {
      shared: normalized,
      animate: snapshot.skymap.animate,
      t_obs: snapshot.skymap.t_obs,
      t_min: snapshot.skymap.t_min,
      t_max: snapshot.skymap.t_max,
      n_frames: snapshot.skymap.n_frames,
      nu_input: snapshot.skymap.nu_input,
      fov: snapshot.skymap.fov,
      npixel: snapshot.skymap.npixel,
    },
  };
}

function groupInstrumentsByType(instrumentNames: string[]): InstrumentGroup[] {
  const grouped = new Map<string, string[]>();
  for (const label of INSTRUMENT_GROUP_ORDER) {
    grouped.set(label, []);
  }

  for (const name of instrumentNames) {
    const label = INSTRUMENT_TYPE_BY_NAME[name] ?? "Other";
    grouped.get(label)?.push(name);
  }

  return INSTRUMENT_GROUP_ORDER.map((label) => ({ label, items: grouped.get(label) ?? [] })).filter(
    (group) => group.items.length > 0,
  );
}

function normalizeUploadCell(value: unknown): string {
  if (value === null || value === undefined) return "";
  if (typeof value === "number") {
    return Number.isFinite(value) ? String(value) : "";
  }
  return String(value).trim();
}

function normalizeObservationTokens(cells: unknown[]): string[] {
  const values = cells.map(normalizeUploadCell).filter((item) => item.length > 0);
  if (values.length < 2) return [];
  const x = values[0];
  const y = values[1];
  const err = values[2] ?? "";

  if (!Number.isFinite(Number(x)) || !Number.isFinite(Number(y))) {
    return [];
  }
  if (err.length > 0 && !Number.isFinite(Number(err))) {
    return [x, y];
  }
  return err.length > 0 ? [x, y, err] : [x, y];
}

function rowsFromDelimitedText(text: string): string[] {
  const lines = text.split(/\r?\n/);
  const rows: string[] = [];
  for (const line of lines) {
    const trimmed = line.trim();
    if (!trimmed || trimmed.startsWith("#")) continue;
    const tokens = normalizeObservationTokens(trimmed.split(/[,\t ]+/));
    if (tokens.length >= 2) {
      rows.push(tokens.join("  "));
    }
  }
  return rows;
}

function rowsFromMatrix(rows: unknown[][]): string[] {
  const out: string[] = [];
  for (const row of rows) {
    const tokens = normalizeObservationTokens(row);
    if (tokens.length >= 2) {
      out.push(tokens.join("  "));
    }
  }
  return out;
}

async function parseObservationUpload(file: File): Promise<string> {
  const lowered = file.name.toLowerCase();
  const isExcel = lowered.endsWith(".xlsx") || lowered.endsWith(".xls");
  let rows: string[] = [];
  if (isExcel) {
    const XLSX = await import("xlsx");
    const bytes = await file.arrayBuffer();
    const workbook = XLSX.read(bytes, { type: "array" });
    const firstSheet = workbook.SheetNames[0];
    if (firstSheet) {
      const sheet = workbook.Sheets[firstSheet];
      const matrix = XLSX.utils.sheet_to_json(sheet, { header: 1, blankrows: false, raw: true }) as unknown[][];
      rows = rowsFromMatrix(matrix);
    }
  } else {
    rows = rowsFromDelimitedText(await file.text());
  }

  if (rows.length === 0) {
    throw new Error("No valid rows found. Expected numeric columns: x, value, error(optional).");
  }
  return rows.join("\n");
}

function HelpHint({ text, ariaLabel }: { text: string; ariaLabel: string }) {
  return (
    <span className="sb-help" tabIndex={0} role="button" aria-label={ariaLabel}>
      ?
      <span className="sb-help-tip">{text}</span>
    </span>
  );
}

function parseAxisRange(value: unknown): AxisRange | null {
  if (!Array.isArray(value) || value.length !== 2) return null;
  const lo = Number(value[0]);
  const hi = Number(value[1]);
  if (!Number.isFinite(lo) || !Number.isFinite(hi)) return null;
  return [lo, hi];
}

function parseRelayoutRange(event: Record<string, unknown>, axis: AxisName): AxisRange | null {
  const direct = parseAxisRange(event[`${axis}.range`]);
  if (direct) return direct;
  const lo = Number(event[`${axis}.range[0]`]);
  const hi = Number(event[`${axis}.range[1]`]);
  if (!Number.isFinite(lo) || !Number.isFinite(hi)) return null;
  return [lo, hi];
}

function rangesEqual(a: AxisRange | null, b: AxisRange | null): boolean {
  if (a === b) return true;
  if (!a || !b) return false;
  return Math.abs(a[0] - b[0]) <= AXIS_EPS && Math.abs(a[1] - b[1]) <= AXIS_EPS;
}

function axisTitleText(axisObj: Record<string, unknown>): string {
  const title = axisObj.title;
  if (typeof title === "string") return title;
  if (title && typeof title === "object") {
    const text = (title as { text?: unknown }).text;
    if (typeof text === "string") return text;
  }
  return "";
}

function axisSignature(layout: Record<string, unknown>, axis: AxisName): string | null {
  const axisObj = layout[axis];
  if (!axisObj || typeof axisObj !== "object") return null;
  const axisRecord = axisObj as Record<string, unknown>;
  const axisType = typeof axisRecord.type === "string" ? axisRecord.type : "";
  const title = axisTitleText(axisRecord);
  return `${axisType}|${title}`;
}

function normalizeShared(shared: SharedParams, mode: Mode): SharedParams {
  const next = { ...shared };
  next.d_L_mpc = normalizeDistanceMpc(next.d_L_mpc);

  if (next.jet_type !== "Power-law") {
    next.k_e = 2.0;
    next.k_g = 2.0;
  }

  if (next.jet_type !== "Two-component") {
    next.theta_w = 0.3;
    next.E_iso_w = 1e51;
    next.Gamma0_w = 100;
  }

  if (next.medium_type === "ISM") {
    next.A_star = 0.1;
    next.k_m = 2.0;
  } else if (next.medium_type === "Wind") {
    next.n_ism = 1.0;
  }

  if (!next.enable_rvs) {
    next.duration = 1.0;
    next.eps_e_r = next.eps_e;
    next.eps_B_r = next.eps_B;
    next.p_r = next.p;
    next.xi_e_r = next.xi_e;
    next.rvs_ssc = false;
    next.rvs_kn = false;
  }

  if (mode === "skymap") {
    next.flux_unit = "cgs";
    next.time_unit = "s";
  }

  return next;
}

function clampInt(value: unknown, min: number, max: number): number {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return min;
  return Math.max(min, Math.min(max, Math.round(numeric)));
}

function buildInteractiveSpec(base: ComputationSpec, interactive: boolean): ComputationSpec {
  if (!interactive) return base;

  if (base.endpoint === "lightcurve") {
    const payload = { ...base.payload };
    const maybeShared = payload.shared;
    if (maybeShared && typeof maybeShared === "object") {
      const nextShared = { ...(maybeShared as Record<string, unknown>) };
      nextShared.num_t = clampInt(nextShared.num_t, 50, INTERACTIVE_LIGHTCURVE_NUM_T_MAX);
      payload.shared = nextShared;
    }
    return { ...base, payload };
  }

  if (base.endpoint === "spectrum") {
    const payload = { ...base.payload };
    payload.num_nu = clampInt(payload.num_nu, 50, INTERACTIVE_SPECTRUM_NUM_NU_MAX);
    return { ...base, payload };
  }

  const payload = { ...base.payload };
  const animate = Boolean(payload.animate);
  payload.npixel = clampInt(
    payload.npixel,
    64,
    animate ? INTERACTIVE_SKY_PIXEL_MAX_ANIMATE : INTERACTIVE_SKY_PIXEL_MAX_STATIC,
  );
  if (animate) {
    payload.n_frames = clampInt(payload.n_frames, 3, INTERACTIVE_SKY_FRAMES_MAX);
  }
  return { ...base, payload };
}

function isEmptyLightcurveFrequencySpec(spec: ComputationSpec): boolean {
  if (spec.endpoint !== "lightcurve") return false;
  const raw = spec.payload.frequencies_input;
  return typeof raw !== "string" || raw.trim().length === 0;
}

type LogSliderProps = {
  label: ReactNode;
  minExp: number;
  maxExp: number;
  step: number;
  value: number;
  defaultExp: number;
  precision?: number;
  onChange: (next: number) => void;
};

function LogSliderField({ label, minExp, maxExp, step, value, defaultExp, precision = 1, onChange }: LogSliderProps) {
  const exp = safeLog10(value, defaultExp);
  const [draftExp, setDraftExp] = useState(exp);
  const draggingRef = useRef(false);
  const timerRef = useRef<number | null>(null);
  const pendingRef = useRef<number | null>(null);
  const fillStyle = useMemo(() => sliderFillStyle(draftExp, minExp, maxExp), [draftExp, minExp, maxExp]);

  useEffect(() => {
    if (!draggingRef.current) {
      setDraftExp(exp);
    }
  }, [exp]);

  const flushCommit = useCallback(
    (nextExp: number) => {
      if (timerRef.current !== null) {
        window.clearTimeout(timerRef.current);
        timerRef.current = null;
      }
      pendingRef.current = null;
      onChange(10 ** nextExp);
    },
    [onChange],
  );

  const scheduleCommit = useCallback(
    (nextExp: number) => {
      pendingRef.current = nextExp;
      if (timerRef.current !== null) return;
      timerRef.current = window.setTimeout(() => {
        timerRef.current = null;
        const pending = pendingRef.current;
        pendingRef.current = null;
        if (pending !== null) {
          onChange(10 ** pending);
        }
      }, SLIDER_COMMIT_INTERVAL_MS);
    },
    [onChange],
  );

  useEffect(() => {
    return () => {
      if (timerRef.current !== null) {
        window.clearTimeout(timerRef.current);
      }
    };
  }, []);

  const handleInput = (raw: string) => {
    const nextExp = Number(raw);
    setDraftExp(nextExp);
    scheduleCommit(nextExp);
  };

  const handleKeyDown = (event: ReactKeyboardEvent<HTMLInputElement>) => {
    stepRangeValue(event, draftExp, minExp, maxExp, step, (nextExp) => {
      setDraftExp(nextExp);
      flushCommit(nextExp);
    });
  };

  const handleStart = (event: ReactPointerEvent<HTMLInputElement>) => {
    draggingRef.current = true;
    event.currentTarget.focus();
  };

  const handleEnd = () => {
    draggingRef.current = false;
    flushCommit(draftExp);
  };

  return (
    <label className="sb-field sb-slider">
      <span className="sb-label">{label}</span>
      <div className="sb-slider-track">
        <input
          type="range"
          min={minExp}
          max={maxExp}
          step={step}
          value={draftExp}
          style={fillStyle}
          onPointerDown={handleStart}
          onPointerUp={handleEnd}
          onPointerCancel={handleEnd}
          onInput={(e) => handleInput((e.target as HTMLInputElement).value)}
          onChange={(e) => handleInput((e.target as HTMLInputElement).value)}
          onKeyDown={handleKeyDown}
          aria-keyshortcuts="ArrowLeft ArrowRight ArrowUp ArrowDown Home End"
        />
        <span className="sb-value">{draftExp.toFixed(precision)}</span>
      </div>
    </label>
  );
}

type SliderProps = {
  label: ReactNode;
  min: number;
  max: number;
  step: number;
  value: number;
  decimals?: number;
  formatValue?: (value: number) => string;
  onChange: (next: number) => void;
};

function SliderField({ label, min, max, step, value, decimals = 2, formatValue, onChange }: SliderProps) {
  const [draft, setDraft] = useState(value);
  const draggingRef = useRef(false);
  const timerRef = useRef<number | null>(null);
  const pendingRef = useRef<number | null>(null);
  const fillStyle = useMemo(() => sliderFillStyle(draft, min, max), [draft, min, max]);

  useEffect(() => {
    if (!draggingRef.current) {
      setDraft(value);
    }
  }, [value]);

  const flushCommit = useCallback(
    (next: number) => {
      if (timerRef.current !== null) {
        window.clearTimeout(timerRef.current);
        timerRef.current = null;
      }
      pendingRef.current = null;
      onChange(next);
    },
    [onChange],
  );

  const scheduleCommit = useCallback(
    (next: number) => {
      pendingRef.current = next;
      if (timerRef.current !== null) return;
      timerRef.current = window.setTimeout(() => {
        timerRef.current = null;
        const pending = pendingRef.current;
        pendingRef.current = null;
        if (pending !== null) {
          onChange(pending);
        }
      }, SLIDER_COMMIT_INTERVAL_MS);
    },
    [onChange],
  );

  useEffect(() => {
    return () => {
      if (timerRef.current !== null) {
        window.clearTimeout(timerRef.current);
      }
    };
  }, []);

  const handleInput = (raw: string) => {
    const next = Number(raw);
    setDraft(next);
    scheduleCommit(next);
  };

  const handleKeyDown = (event: ReactKeyboardEvent<HTMLInputElement>) => {
    stepRangeValue(event, draft, min, max, step, (next) => {
      setDraft(next);
      flushCommit(next);
    });
  };

  const handleStart = (event: ReactPointerEvent<HTMLInputElement>) => {
    draggingRef.current = true;
    event.currentTarget.focus();
  };

  const handleEnd = () => {
    draggingRef.current = false;
    flushCommit(draft);
  };

  return (
    <label className="sb-field sb-slider">
      <span className="sb-label">{label}</span>
      <div className="sb-slider-track">
        <input
          type="range"
          min={min}
          max={max}
          step={step}
          value={draft}
          style={fillStyle}
          onPointerDown={handleStart}
          onPointerUp={handleEnd}
          onPointerCancel={handleEnd}
          onInput={(e) => handleInput((e.target as HTMLInputElement).value)}
          onChange={(e) => handleInput((e.target as HTMLInputElement).value)}
          onKeyDown={handleKeyDown}
          aria-keyshortcuts="ArrowLeft ArrowRight ArrowUp ArrowDown Home End"
        />
        <span className="sb-value">{formatValue ? formatValue(draft) : draft.toFixed(decimals)}</span>
      </div>
    </label>
  );
}

export default function HomePage() {
  const apiCandidates = useMemo(() => {
    const configured = process.env.NEXT_PUBLIC_API_URL;
    const host = typeof window !== "undefined" ? window.location.hostname : "";
    const fallback = host === "127.0.0.1" ? "http://127.0.0.1:8000" : "http://localhost:8000";
    const primary = stripTrailingSlash(configured ?? fallback);
    const candidates = [primary];
    if (primary.includes("localhost")) {
      candidates.push(primary.replace("localhost", "127.0.0.1"));
    }
    if (primary.includes("127.0.0.1")) {
      candidates.push(primary.replace("127.0.0.1", "localhost"));
    }
    return Array.from(new Set(candidates));
  }, []);
  const apiStatusExtraCandidates = useMemo(() => {
    const extraRaw = process.env.NEXT_PUBLIC_API_STATUS_URLS ?? "";
    return extraRaw
      .split(",")
      .map((item) => stripTrailingSlash(item.trim()))
      .filter((item) => item.length > 0);
  }, []);
  const apiStatusCandidates = useMemo(() => {
    // Always probe both direct backends and primary API domain. This lets us
    // map "(in use)" from global API to the concrete backend region.
    return Array.from(new Set([...apiStatusExtraCandidates, ...apiCandidates]));
  }, [apiCandidates, apiStatusExtraCandidates]);
  const apiStatusDisplaySet = useMemo(() => {
    return apiStatusExtraCandidates.length > 0 ? new Set(apiStatusExtraCandidates) : new Set(apiStatusCandidates);
  }, [apiStatusCandidates, apiStatusExtraCandidates]);

  const [mode, setMode] = useState<Mode>("lightcurve");
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const [shared, setShared] = useState<SharedParams>(FALLBACK_SHARED);
  const [distanceLinked, setDistanceLinked] = useState(true);
  const [distanceDriver, setDistanceDriver] = useState<DistanceDriver>("dL");

  const [error, setError] = useState<string>("");
  const [running, setRunning] = useState(false);
  const [showRunning, setShowRunning] = useState(false);
  const [result, setResult] = useState<RunResponse | null>(null);
  const [resultMode, setResultMode] = useState<Mode | null>(null);
  const [compareResult, setCompareResult] = useState<RunResponse | null>(null);
  const [compareResultMode, setCompareResultMode] = useState<Mode | null>(null);
  const [compareRunning, setCompareRunning] = useState(false);
  const [compareError, setCompareError] = useState("");
  const [plotWidthPx, setPlotWidthPx] = useState(960);
  const [isFigurePending, startFigureTransition] = useTransition();
  const [downloading, setDownloading] = useState<DownloadKind | null>(null);
  const [appVersion, setAppVersion] = useState(DEFAULT_APP_VERSION);
  const [citeLinkText, setCiteLinkText] = useState("cite");
  const [setupStatusText, setSetupStatusText] = useState("");
  const [bootReady, setBootReady] = useState(false);
  const [obsStorageReady, setObsStorageReady] = useState(false);
  const [urlStateReady, setUrlStateReady] = useState(false);
  const [zoomRevision, setZoomRevision] = useState(0);
  const requestSeqRef = useRef(0);
  const requestInFlightRef = useRef(false);
  const pendingSpecRef = useRef<ComputationSpec | null>(null);
  const runTimerRef = useRef<number | null>(null);
  const activeRequestRef = useRef<AbortController | null>(null);
  const urlSyncTimerRef = useRef<number | null>(null);
  const setupStatusTimerRef = useRef<number | null>(null);
  const compareRequestSeqRef = useRef(0);
  const compareTimerRef = useRef<number | null>(null);
  const compareActiveRequestRef = useRef<AbortController | null>(null);
  const [sliderInteracting, setSliderInteracting] = useState(false);
  const sliderInteractingRef = useRef(false);
  const workspaceRef = useRef<HTMLElement | null>(null);
  const axisRangesRef = useRef<Record<Mode, AxisRanges>>({
    lightcurve: { xaxis: null, yaxis: null, yaxis2: null },
    spectrum: { xaxis: null, yaxis: null, yaxis2: null },
    skymap: { xaxis: null, yaxis: null, yaxis2: null },
  });
  const axisSignaturesRef = useRef<Record<Mode, AxisSignatures>>({
    lightcurve: { yaxis: null, yaxis2: null },
    spectrum: { yaxis: null, yaxis2: null },
    skymap: { yaxis: null, yaxis2: null },
  });

  const [instrumentGroups, setInstrumentGroups] = useState<InstrumentGroup[]>([]);

  const [lcFreq, setLcFreq] = useState("1e9, R, 1keV");
  const [lcTMin, setLcTMin] = useState(1);
  const [lcTMax, setLcTMax] = useState(1e8);
  const [lcInstruments, setLcInstruments] = useState<string[]>([]);
  const [lcObsGroups, setLcObsGroups] = useState<ObservationGroup[]>([]);
  const [activeLcObsTab, setActiveLcObsTab] = useState(0);

  const [sedTimes, setSedTimes] = useState("1e3, 1e4, 1e5, 1e6");
  const [sedTimesDraft, setSedTimesDraft] = useState("1e3, 1e4, 1e5, 1e6");
  const [sedNuMin, setSedNuMin] = useState(1e8);
  const [sedNuMax, setSedNuMax] = useState(1e20);
  const [sedNumNu, setSedNumNu] = useState(200);
  const [sedFreqUnit, setSedFreqUnit] = useState("Hz");
  const [sedNuFNu, setSedNuFNu] = useState(false);
  const [sedInstruments, setSedInstruments] = useState<string[]>([]);
  const [sedObsGroups, setSedObsGroups] = useState<ObservationGroup[]>([]);
  const [activeSedObsTab, setActiveSedObsTab] = useState(0);

  const [skyAnimate, setSkyAnimate] = useState(false);
  const [skyTObs, setSkyTObs] = useState(1e6);
  const [skyTMin, setSkyTMin] = useState(1e4);
  const [skyTMax, setSkyTMax] = useState(1e7);
  const [skyNFrames, setSkyNFrames] = useState(SKY_DEFAULT_N_FRAMES);
  const [skyNuInput, setSkyNuInput] = useState("1e9");
  const [skyFov, setSkyFov] = useState(500);
  const [skyNpixel, setSkyNpixel] = useState(256);
  const [bookmarks, setBookmarks] = useState<BookmarkEntry[]>([]);
  const [bookmarkNameDraft, setBookmarkNameDraft] = useState("");
  const [compareEnabled, setCompareEnabled] = useState(false);
  const [compareBookmarkId, setCompareBookmarkId] = useState("");
  const [activeApiKey, setActiveApiKey] = useState<string>("");
  const [showAllServerStatus, setShowAllServerStatus] = useState(false);
  const [apiStatusByKey, setApiStatusByKey] = useState<Record<string, ApiHealthStatus>>({});
  const apiStatusByKeyRef = useRef<Record<string, ApiHealthStatus>>({});

  const setSharedField = useCallback(<K extends keyof SharedParams>(key: K, value: SharedParams[K]) => {
    setShared((prev) => ({ ...prev, [key]: value }));
  }, []);

  const setDistanceMpc = useCallback(
    (nextMpcRaw: number) => {
      if (!Number.isFinite(nextMpcRaw)) return;
      const nextMpc = normalizeDistanceMpc(nextMpcRaw);
      setShared((prev) => {
        const next: SharedParams = { ...prev, d_L_mpc: nextMpc };
        if (distanceLinked) {
          next.z = normalizeRedshift(redshiftFromLuminosityDistanceMpc(nextMpc));
        }
        return next;
      });
    },
    [distanceLinked],
  );

  const setDistanceRedshift = useCallback(
    (nextZRaw: number) => {
      if (!Number.isFinite(nextZRaw)) return;
      const nextZ = Math.max(0, nextZRaw);
      setShared((prev) => {
        const next: SharedParams = { ...prev, z: nextZ };
        if (distanceLinked) {
          next.d_L_mpc = normalizeDistanceMpc(luminosityDistanceMpcFromRedshift(nextZ));
        }
        return next;
      });
    },
    [distanceLinked],
  );

  useEffect(() => {
    if (!distanceLinked) return;
    setShared((prev) => {
      if (distanceDriver === "dL") {
        const nextMpc = normalizeDistanceMpc(prev.d_L_mpc);
        const nextZ = normalizeRedshift(redshiftFromLuminosityDistanceMpc(nextMpc));
        if (nextMpc !== prev.d_L_mpc) {
          return { ...prev, d_L_mpc: nextMpc, z: nextZ };
        }
        if (nearlyEqualRelative(prev.z, nextZ)) return prev;
        return { ...prev, z: nextZ };
      }
      const nextDL = normalizeDistanceMpc(luminosityDistanceMpcFromRedshift(prev.z));
      if (prev.d_L_mpc === nextDL) return prev;
      return { ...prev, d_L_mpc: nextDL };
    });
  }, [distanceDriver, distanceLinked, shared.d_L_mpc, shared.z]);

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
    [
      mode,
      shared,
      distanceDriver,
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
      distanceLinked,
    ],
  );

  const applySnapshot = useCallback((snapshot: UiSnapshot) => {
    setMode(snapshot.mode);
    setDistanceLinked(snapshot.distance_linked);
    setDistanceDriver(snapshot.distance_driver);
    setShared({ ...FALLBACK_SHARED, ...snapshot.shared });
    setLcFreq(snapshot.lightcurve.frequencies_input);
    setLcTMin(snapshot.lightcurve.t_min);
    setLcTMax(snapshot.lightcurve.t_max);
    setLcInstruments([...snapshot.lightcurve.selected_instruments]);
    setLcObsGroups(snapshot.lightcurve.observation_groups.map((group) => ({ ...group })));
    setActiveLcObsTab(0);

    setSedTimes(snapshot.spectrum.t_snapshots_input);
    setSedTimesDraft(snapshot.spectrum.t_snapshots_input);
    setSedNuMin(snapshot.spectrum.nu_min);
    setSedNuMax(snapshot.spectrum.nu_max);
    setSedNumNu(snapshot.spectrum.num_nu);
    setSedFreqUnit(snapshot.spectrum.freq_unit);
    setSedNuFNu(snapshot.spectrum.show_nufnu);
    setSedInstruments([...snapshot.spectrum.selected_instruments]);
    setSedObsGroups(snapshot.spectrum.observation_groups.map((group) => ({ ...group })));
    setActiveSedObsTab(0);

    setSkyAnimate(snapshot.skymap.animate);
    setSkyTObs(snapshot.skymap.t_obs);
    setSkyTMin(snapshot.skymap.t_min);
    setSkyTMax(snapshot.skymap.t_max);
    setSkyNFrames(snapshot.skymap.n_frames);
    setSkyNuInput(snapshot.skymap.nu_input);
    setSkyFov(snapshot.skymap.fov);
    setSkyNpixel(snapshot.skymap.npixel);
    setError("");
    setCompareError("");
  }, []);

  const copyTextToClipboard = useCallback(async (text: string): Promise<boolean> => {
    try {
      await navigator.clipboard.writeText(text);
      return true;
    } catch {
      try {
        const ta = document.createElement("textarea");
        ta.value = text;
        ta.style.position = "fixed";
        ta.style.left = "-9999px";
        document.body.appendChild(ta);
        ta.select();
        const copied = document.execCommand("copy");
        document.body.removeChild(ta);
        return copied;
      } catch {
        return false;
      }
    }
  }, []);

  const saveSnapshotAsBookmark = useCallback((snapshot: UiSnapshot, proposedName?: string) => {
    const nameHint = proposedName?.trim();
    setBookmarks((prev) => {
      const nextSetupIndex =
        prev.reduce((max, entry) => {
          const matched = /^setup\s+(\d+)$/i.exec(entry.name.trim());
          if (!matched) return max;
          const value = Number(matched[1]);
          return Number.isFinite(value) ? Math.max(max, value) : max;
        }, 0) + 1;
      const fallbackName = `Setup ${nextSetupIndex}`;
      const name = nameHint && nameHint.length > 0 ? nameHint : fallbackName;
      const entry: BookmarkEntry = {
        id: `bm-${Date.now()}-${Math.random().toString(36).slice(2, 8)}`,
        name,
        created_at: new Date().toISOString(),
        snapshot: cloneSnapshot(snapshot),
      };
      return [entry, ...prev].slice(0, MAX_BOOKMARKS);
    });
    setBookmarkNameDraft("");
  }, []);

  const showSetupStatus = useCallback((text: string, timeoutMs = 1500) => {
    setSetupStatusText(text);
    if (setupStatusTimerRef.current !== null) {
      window.clearTimeout(setupStatusTimerRef.current);
      setupStatusTimerRef.current = null;
    }
    setupStatusTimerRef.current = window.setTimeout(() => {
      setupStatusTimerRef.current = null;
      setSetupStatusText("");
    }, timeoutMs);
  }, []);

  const buildShareLink = useCallback((snapshot: UiSnapshot): string | null => {
    if (typeof window === "undefined") return null;
    const encoded = encodeUrlState(makeShareableSnapshot(snapshot));
    if (!encoded) return null;
    const url = new URL(window.location.href);
    url.searchParams.set(URL_STATE_PARAM, encoded);
    return url.toString();
  }, []);

  const copyShareLinkForSnapshot = useCallback(
    async (snapshot: UiSnapshot) => {
      const link = buildShareLink(snapshot);
      if (!link) {
        showSetupStatus("Link is too large to share.");
        return;
      }
      const copied = await copyTextToClipboard(link);
      showSetupStatus(copied ? "setup url link copied" : "Copy failed.");
    },
    [buildShareLink, copyTextToClipboard, showSetupStatus],
  );

  const saveCurrentBookmark = useCallback(() => {
    saveSnapshotAsBookmark(currentSnapshot, bookmarkNameDraft);
  }, [bookmarkNameDraft, currentSnapshot, saveSnapshotAsBookmark]);

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

  const removeBookmarkById = useCallback((bookmarkId: string) => {
    setBookmarks((prev) => prev.filter((bookmark) => bookmark.id !== bookmarkId));
  }, []);

  const selectedCompareBookmark = useMemo(
    () => bookmarks.find((bookmark) => bookmark.id === compareBookmarkId) ?? null,
    [bookmarks, compareBookmarkId],
  );

  const apiEndpoints = useMemo<ApiEndpoint[]>(
    () =>
      apiStatusCandidates.map((base) => ({
        key: base,
        label: apiServerFromTarget(base),
        healthTarget: `${base}/api/health`,
      })),
    [apiStatusCandidates],
  );

  const resolveApiKeyFromTarget = useCallback(
    (target: string): string | null => {
      for (const base of apiCandidates) {
        if (target.startsWith(`${base}/api/`)) return base;
      }
      return null;
    },
    [apiCandidates],
  );

  const noteApiTarget = useCallback((target: string) => {
    const key = resolveApiKeyFromTarget(target);
    if (!key) return;
    setActiveApiKey((prev) => (prev === key ? prev : key));
  }, [resolveApiKeyFromTarget]);

  const updateApiStatus = useCallback((key: string, next: ApiHealthStatus) => {
    setApiStatusByKey((prev) => {
      const current = prev[key];
      if (
        current &&
        current.up === next.up &&
        current.latencyMs === next.latencyMs &&
        current.statusCode === next.statusCode &&
        current.region === next.region &&
        current.locationLabel === next.locationLabel &&
        current.backend === next.backend
      ) {
        return prev;
      }
      return { ...prev, [key]: next };
    });
  }, []);

  useEffect(() => {
    apiStatusByKeyRef.current = apiStatusByKey;
  }, [apiStatusByKey]);

  const probeEndpoint = useCallback(
    async (endpoint: ApiEndpoint, cancelledRef?: { current: boolean }): Promise<void> => {
      const started = performance.now();
      try {
        const response = await fetch(endpoint.healthTarget, { cache: "no-store" });
        if (cancelledRef?.current) return;
        if (!response.ok) {
          updateApiStatus(endpoint.key, {
            up: false,
            latencyMs: null,
            statusCode: response.status,
            region: null,
            locationLabel: null,
            backend: null,
          });
          return;
        }
        let backend: string | null = null;
        let region: string | null = null;
        let locationLabel: string | null = null;
        try {
          const payload = (await response.clone().json()) as unknown;
          region = regionFromHealthPayload(payload);
          locationLabel = locationLabelFromHealthPayload(payload);
          backend = backendFromHealthPayload(payload);
        } catch {
          region = null;
          locationLabel = null;
          backend = null;
        }
        const elapsed = Math.max(1, Math.round(performance.now() - started));
        updateApiStatus(endpoint.key, {
          up: true,
          latencyMs: elapsed,
          statusCode: response.status,
          region,
          locationLabel,
          backend,
        });
      } catch {
        if (cancelledRef?.current) return;
        updateApiStatus(endpoint.key, {
          up: false,
          latencyMs: null,
          statusCode: null,
          region: null,
          locationLabel: null,
          backend: null,
        });
      }
    },
    [updateApiStatus],
  );

  useEffect(() => {
    setApiStatusByKey((prev) => {
      const allowed = new Set(apiEndpoints.map((endpoint) => endpoint.key));
      let changed = false;
      const next: Record<string, ApiHealthStatus> = {};
      for (const [key, value] of Object.entries(prev)) {
        if (!allowed.has(key)) {
          changed = true;
          continue;
        }
        next[key] = value;
      }
      return changed ? next : prev;
    });
    if (activeApiKey && !apiEndpoints.some((endpoint) => endpoint.key === activeApiKey)) {
      setActiveApiKey("");
    }
  }, [activeApiKey, apiEndpoints]);

  useEffect(() => {
    const cancelledRef = { current: false };
    let directTimer: number | null = null;

    const endpointByKey = (key: string): ApiEndpoint | null =>
      apiEndpoints.find((endpoint) => endpoint.key === key) ?? null;
    const isDisplayEndpoint = (key: string): boolean => apiStatusDisplaySet.has(key);
    const mapDisplayEndpointByRegion = (region: string | null): ApiEndpoint | null => {
      if (!region) return null;
      return (
        apiEndpoints.find((endpoint) => {
          if (!isDisplayEndpoint(endpoint.key)) return false;
          const endpointRegion = apiStatusByKeyRef.current[endpoint.key]?.region ?? null;
          return endpointRegion === region;
        }) ?? null
      );
    };
    const pickDisplayEndpoint = (): ApiEndpoint | null => {
      if (activeApiKey && isDisplayEndpoint(activeApiKey)) {
        const direct = endpointByKey(activeApiKey);
        if (direct) return direct;
      }
      if (activeApiKey) {
        const region = apiStatusByKeyRef.current[activeApiKey]?.region ?? null;
        const mapped = mapDisplayEndpointByRegion(region);
        if (mapped) return mapped;
      }
      return apiEndpoints.find((endpoint) => isDisplayEndpoint(endpoint.key)) ?? apiEndpoints[0] ?? null;
    };

    const pickActiveEndpoint = (): ApiEndpoint | null => {
      if (activeApiKey) {
        const source = endpointByKey(activeApiKey);
        if (source) return source;
      }
      return pickDisplayEndpoint();
    };

    const probeActive = async () => {
      const sourceEndpoint = pickActiveEndpoint();
      if (sourceEndpoint) {
        await probeEndpoint(sourceEndpoint, cancelledRef);
      }
      const displayEndpoint = pickDisplayEndpoint();
      if (displayEndpoint && (!sourceEndpoint || displayEndpoint.key !== sourceEndpoint.key)) {
        await probeEndpoint(displayEndpoint, cancelledRef);
      }
    };

    void probeActive();
    directTimer = window.setInterval(() => {
      void probeActive();
    }, API_STATUS_REFRESH_DIRECT_MS);

    return () => {
      cancelledRef.current = true;
      if (directTimer !== null) {
        window.clearInterval(directTimer);
      }
    };
  }, [activeApiKey, apiEndpoints, apiStatusDisplaySet, probeEndpoint]);

  useEffect(() => {
    const cancelledRef = { current: false };
    const targets = apiEndpoints.filter((endpoint) => apiStatusDisplaySet.has(endpoint.key));
    if (targets.length === 0) return;
    // Measure all visible backends once on load so dropdown rows are pre-populated.
    void Promise.all(targets.map((endpoint) => probeEndpoint(endpoint, cancelledRef)));
    return () => {
      cancelledRef.current = true;
    };
  }, [apiEndpoints, apiStatusDisplaySet, probeEndpoint]);

  const compareSpec = useMemo<ComputationSpec | null>(() => {
    if (!compareEnabled || mode === "skymap" || !selectedCompareBookmark) return null;
    return buildComputationSpecFromSnapshot(selectedCompareBookmark.snapshot, mode);
  }, [compareEnabled, mode, selectedCompareBookmark]);

  const fetchFromApi = useCallback(async (path: string, init?: RequestInit): Promise<Response> => {
    const cleanPath = path.replace(/^\/+/, "");
    let lastErr: unknown = null;

    const directTargets = apiCandidates.map((base) => `${base}/api/${cleanPath}`);
    for (const target of directTargets) {
      try {
        const response = await fetch(target, init);
        noteApiTarget(target);
        return response;
      } catch (err) {
        lastErr = err;
      }
    }
    throw lastErr instanceof Error ? lastErr : new Error("Load failed");
  }, [apiCandidates, noteApiTarget]);

  useEffect(() => {
    async function boot() {
      try {
        const [defaultsRes, optionsRes] = await Promise.all([
          fetchFromApi("defaults", { cache: "no-store" }),
          fetchFromApi("options", { cache: "no-store" }),
        ]);

        if (defaultsRes.ok) {
          const data = (await defaultsRes.json()) as {
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

          if (data.shared) {
            setShared({ ...FALLBACK_SHARED, ...data.shared });
          }

          if (data.lightcurve) {
            if (typeof data.lightcurve.frequencies_input === "string") setLcFreq(data.lightcurve.frequencies_input);
            if (typeof data.lightcurve.t_min === "number") setLcTMin(data.lightcurve.t_min);
            if (typeof data.lightcurve.t_max === "number") setLcTMax(data.lightcurve.t_max);
            if (Array.isArray(data.lightcurve.selected_instruments)) setLcInstruments(data.lightcurve.selected_instruments);
            if (Array.isArray(data.lightcurve.observation_groups)) setLcObsGroups(data.lightcurve.observation_groups);
          }

          if (data.spectrum) {
            if (typeof data.spectrum.t_snapshots_input === "string") {
              setSedTimes(data.spectrum.t_snapshots_input);
              setSedTimesDraft(data.spectrum.t_snapshots_input);
            }
            if (typeof data.spectrum.nu_min === "number") setSedNuMin(data.spectrum.nu_min);
            if (typeof data.spectrum.nu_max === "number") setSedNuMax(data.spectrum.nu_max);
            if (typeof data.spectrum.num_nu === "number") setSedNumNu(data.spectrum.num_nu);
            if (typeof data.spectrum.freq_unit === "string") setSedFreqUnit(data.spectrum.freq_unit);
            if (typeof data.spectrum.show_nufnu === "boolean") setSedNuFNu(data.spectrum.show_nufnu);
            if (Array.isArray(data.spectrum.selected_instruments)) setSedInstruments(data.spectrum.selected_instruments);
            if (Array.isArray(data.spectrum.observation_groups)) setSedObsGroups(data.spectrum.observation_groups);
          }

          if (data.skymap) {
            if (typeof data.skymap.animate === "boolean") setSkyAnimate(data.skymap.animate);
            if (typeof data.skymap.t_obs === "number") setSkyTObs(data.skymap.t_obs);
            if (typeof data.skymap.t_min === "number") setSkyTMin(data.skymap.t_min);
            if (typeof data.skymap.t_max === "number") setSkyTMax(data.skymap.t_max);
            if (typeof data.skymap.n_frames === "number") {
              setSkyNFrames(Math.max(3, Math.min(SKY_MAX_N_FRAMES, Math.round(data.skymap.n_frames))));
            }
            if (typeof data.skymap.nu_input === "string") setSkyNuInput(data.skymap.nu_input);
            if (typeof data.skymap.fov === "number") setSkyFov(data.skymap.fov);
            if (typeof data.skymap.npixel === "number") {
              setSkyNpixel(Math.max(64, Math.min(SKY_MAX_PIXEL_STATIC, Math.round(data.skymap.npixel))));
            }
          }
        }

        if (optionsRes.ok) {
          const options = (await optionsRes.json()) as OptionsResponse;
          const groupedFromApi =
            Array.isArray(options.instrument_groups) && options.instrument_groups.length > 0
              ? options.instrument_groups
                  .map((group) => ({
                    label: typeof group.label === "string" ? group.label : "",
                    items: Array.isArray(group.items) ? group.items.filter((name): name is string => typeof name === "string") : [],
                  }))
                  .filter((group) => group.label.length > 0 && group.items.length > 0)
              : [];
          const fallbackInstruments = Array.isArray(options.instruments) ? options.instruments : [];

          if (groupedFromApi.length > 0) {
            const groupedKnown = new Set(groupedFromApi.flatMap((group) => group.items));
            const ungrouped = fallbackInstruments.filter((name) => !groupedKnown.has(name));
            const merged = [...groupedFromApi];
            if (ungrouped.length > 0) {
              merged.push({ label: "Other", items: ungrouped });
            }
            setInstrumentGroups(merged);
          } else if (fallbackInstruments.length > 0) {
            setInstrumentGroups(groupInstrumentsByType(fallbackInstruments));
          } else {
            setInstrumentGroups([]);
          }

          if (Array.isArray(options.instruments)) {
            const allowed = new Set(options.instruments);
            setLcInstruments((prev) => prev.filter((name) => allowed.has(name)));
            setSedInstruments((prev) => prev.filter((name) => allowed.has(name)));
          }
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
    try {
      const storedLc = parseStoredObsGroups(window.localStorage.getItem(OBS_STORAGE_LC_KEY), true);
      if (storedLc !== null) {
        setLcObsGroups(storedLc);
      }
      const storedSed = parseStoredObsGroups(window.localStorage.getItem(OBS_STORAGE_SED_KEY), false);
      if (storedSed !== null) {
        setSedObsGroups(storedSed);
      }
    } catch {
      // Ignore localStorage errors (private mode / quota / unavailable).
    } finally {
      setObsStorageReady(true);
    }
  }, [bootReady]);

  useEffect(() => {
    if (!obsStorageReady) return;
    try {
      window.localStorage.setItem(OBS_STORAGE_LC_KEY, JSON.stringify(lcObsGroups));
    } catch {
      // Ignore localStorage write errors.
    }
  }, [lcObsGroups, obsStorageReady]);

  useEffect(() => {
    if (!obsStorageReady) return;
    try {
      window.localStorage.setItem(OBS_STORAGE_SED_KEY, JSON.stringify(sedObsGroups));
    } catch {
      // Ignore localStorage write errors.
    }
  }, [obsStorageReady, sedObsGroups]);

  useEffect(() => {
    try {
      const stored = parseStoredBookmarks(window.localStorage.getItem(BOOKMARKS_STORAGE_KEY));
      if (stored !== null) {
        setBookmarks(stored.slice(0, MAX_BOOKMARKS));
      }
    } catch {
      // Ignore localStorage read errors.
    }
  }, []);

  useEffect(() => {
    try {
      window.localStorage.setItem(BOOKMARKS_STORAGE_KEY, JSON.stringify(bookmarks));
    } catch {
      // Ignore localStorage write errors.
    }
  }, [bookmarks]);

  useEffect(() => {
    if (!bootReady || typeof window === "undefined") return;
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
  }, [applySnapshot, bootReady]);

  useEffect(() => {
    if (!urlStateReady || typeof window === "undefined") return;
    if (sliderInteracting) return;
    if (urlSyncTimerRef.current !== null) {
      window.clearTimeout(urlSyncTimerRef.current);
      urlSyncTimerRef.current = null;
    }
    const link = buildShareLink(currentSnapshot);
    if (!link) return;
    urlSyncTimerRef.current = window.setTimeout(() => {
      urlSyncTimerRef.current = null;
      const current = window.location.href;
      if (current === link) return;
      window.history.replaceState({}, "", link);
    }, 120);

    return () => {
      if (urlSyncTimerRef.current !== null) {
        window.clearTimeout(urlSyncTimerRef.current);
        urlSyncTimerRef.current = null;
      }
    };
  }, [buildShareLink, currentSnapshot, sliderInteracting, urlStateReady]);

  useEffect(() => {
    return () => {
      if (setupStatusTimerRef.current !== null) {
        window.clearTimeout(setupStatusTimerRef.current);
        setupStatusTimerRef.current = null;
      }
    };
  }, []);

  useEffect(() => {
    if (bookmarks.length === 0) {
      setCompareBookmarkId("");
      setCompareEnabled(false);
      return;
    }
    setCompareBookmarkId((prev) => (bookmarks.some((bookmark) => bookmark.id === prev) ? prev : bookmarks[0].id));
  }, [bookmarks]);

  useEffect(() => {
    setActiveLcObsTab((prev) => {
      if (lcObsGroups.length === 0) return 0;
      return Math.min(prev, lcObsGroups.length - 1);
    });
  }, [lcObsGroups.length]);

  useEffect(() => {
    setActiveSedObsTab((prev) => {
      if (sedObsGroups.length === 0) return 0;
      return Math.min(prev, sedObsGroups.length - 1);
    });
  }, [sedObsGroups.length]);

  const warnings = useMemo(() => (resultMode === mode ? result?.meta?.warnings ?? [] : []), [mode, result, resultMode]);

  const handlePlotRelayout = useCallback((event: Record<string, unknown>) => {
    const current = axisRangesRef.current[mode];
    const next: AxisRanges = { ...current };
    let changed = false;

    const axes: AxisName[] = ["xaxis", "yaxis", "yaxis2"];
    for (const axis of axes) {
      if (event[`${axis}.autorange`] === true && next[axis] !== null) {
        next[axis] = null;
        changed = true;
      }
      const range = parseRelayoutRange(event, axis);
      if (range) {
        const prev = next[axis];
        if (!rangesEqual(prev, range)) {
          next[axis] = range;
          changed = true;
        }
      }
    }

    if (changed) {
      axisRangesRef.current[mode] = next;
      setZoomRevision((v) => v + 1);
    }
  }, [mode]);

  useEffect(() => {
    if (resultMode !== mode) return;
    const layout = result?.figure?.layout;
    if (!layout) return;

    const current = axisRangesRef.current[mode];
    const next: AxisRanges = { ...current };
    let changed = false;

    const axes: AxisName[] = ["xaxis", "yaxis", "yaxis2"];
    for (const axis of axes) {
      if (next[axis] !== null) continue;
      const axisObj = layout[axis];
      if (axisObj && typeof axisObj === "object") {
        const range = parseAxisRange((axisObj as { range?: unknown }).range);
        if (range) {
          next[axis] = range;
          changed = true;
        }
      }
    }

    if (changed) {
      axisRangesRef.current[mode] = next;
      setZoomRevision((v) => v + 1);
    }
  }, [mode, result, resultMode]);

  useEffect(() => {
    if (resultMode !== mode) return;
    const layout = result?.figure?.layout;
    if (!layout) return;

    const nextYSig = axisSignature(layout, "yaxis");
    const nextY2Sig = axisSignature(layout, "yaxis2");
    const prevSig = axisSignaturesRef.current[mode];

    let clearY = false;
    let clearY2 = false;

    if (prevSig.yaxis !== null && prevSig.yaxis !== nextYSig) {
      clearY = true;
    }
    if (prevSig.yaxis2 !== null && prevSig.yaxis2 !== nextY2Sig) {
      clearY2 = true;
    }

    axisSignaturesRef.current[mode] = { yaxis: nextYSig, yaxis2: nextY2Sig };

    if (!clearY && !clearY2) return;

    const zoom = axisRangesRef.current[mode];
    const nextZoom: AxisRanges = { ...zoom };
    let changed = false;

    if (clearY && nextZoom.yaxis !== null) {
      nextZoom.yaxis = null;
      changed = true;
    }
    if (clearY2 && nextZoom.yaxis2 !== null) {
      nextZoom.yaxis2 = null;
      changed = true;
    }

    if (changed) {
      axisRangesRef.current[mode] = nextZoom;
      setZoomRevision((v) => v + 1);
    }
  }, [mode, result, resultMode]);

  useEffect(() => {
    if (typeof window === "undefined") return;
    const workspace = workspaceRef.current;
    if (!workspace) return;
    let observedPlotWrap: HTMLElement | null = null;

    const resolveMeasuredWidth = () => {
      const plotWrap = workspace.querySelector(".plot-wrap") as HTMLElement | null;
      return plotWrap?.clientWidth ?? workspace.clientWidth;
    };

    const updateWidth = () => {
      const next = resolveMeasuredWidth();
      if (next <= 0) return;
      setPlotWidthPx((prev) => (Math.abs(prev - next) > 1 ? next : prev));
    };

    const attachPlotObserver = (observer: ResizeObserver) => {
      const nextPlotWrap = workspace.querySelector(".plot-wrap") as HTMLElement | null;
      if (nextPlotWrap === observedPlotWrap) return;
      if (observedPlotWrap) {
        observer.unobserve(observedPlotWrap);
      }
      observedPlotWrap = nextPlotWrap;
      if (observedPlotWrap) {
        observer.observe(observedPlotWrap);
      }
    };

    updateWidth();
    const resizeObserver = new ResizeObserver(() => {
      attachPlotObserver(resizeObserver);
      updateWidth();
    });
    resizeObserver.observe(workspace);
    attachPlotObserver(resizeObserver);

    const mutationObserver = new MutationObserver(() => {
      attachPlotObserver(resizeObserver);
      updateWidth();
    });
    mutationObserver.observe(workspace, { childList: true, subtree: true });

    window.addEventListener("resize", updateWidth);
    window.addEventListener("orientationchange", updateWidth);

    return () => {
      if (observedPlotWrap) {
        resizeObserver.unobserve(observedPlotWrap);
      }
      resizeObserver.disconnect();
      mutationObserver.disconnect();
      window.removeEventListener("resize", updateWidth);
      window.removeEventListener("orientationchange", updateWidth);
    };
  }, []);

  const displayFigure = useMemo(() => {
    if (resultMode !== mode) return null;
    const figure = result?.figure;
    if (!figure?.data) return null;
    const rawData = Array.isArray(figure.data) ? [...figure.data] : [];
    const baseData =
      mode === "skymap"
        ? rawData
        : rawData.map((trace) => {
            if (!trace || typeof trace !== "object") return trace;
            return remapScientificHoverTemplate(trace as Record<string, unknown>);
          });
    let mergedData = baseData;

    if (mode !== "skymap" && compareEnabled && compareResultMode === mode && selectedCompareBookmark) {
      const compareData = compareResult?.figure?.data;
      if (Array.isArray(compareData) && compareData.length > 0) {
        const compareLabel = selectedCompareBookmark.name;
        const overlayData = compareData.map((trace, idx) => {
          if (!trace || typeof trace !== "object") return trace;
          const sourceTrace = remapScientificHoverTemplate(trace as Record<string, unknown>);
          const nextTrace = { ...sourceTrace };
          const name = typeof nextTrace.name === "string" ? nextTrace.name : `trace ${idx + 1}`;
          nextTrace.name = `${name} (${compareLabel})`;
          nextTrace.legendgroup =
            typeof nextTrace.legendgroup === "string" ? `cmp-${nextTrace.legendgroup}` : `cmp-${idx + 1}`;
          nextTrace.opacity =
            typeof nextTrace.opacity === "number" ? Math.min(0.9, nextTrace.opacity) : 0.85;
          if (typeof nextTrace.fill === "string" && nextTrace.fill !== "none") {
            nextTrace.fill = "none";
          }

          const traceType = typeof nextTrace.type === "string" ? nextTrace.type : "";
          if (traceType === "" || traceType === "scatter" || traceType === "scattergl") {
            const lineObj =
              nextTrace.line && typeof nextTrace.line === "object"
                ? { ...(nextTrace.line as Record<string, unknown>) }
                : {};
            lineObj.dash = "dot";
            nextTrace.line = lineObj;
          }
          return nextTrace;
        });
        mergedData = [...baseData, ...overlayData];
      }
    }

    const layout = { ...(figure.layout ?? {}) } as Record<string, unknown>;
    const zoom = axisRangesRef.current[mode];
    const axes: AxisName[] = ["xaxis", "yaxis", "yaxis2"];

    for (const axis of axes) {
      const saved = zoom[axis];
      if (!saved) continue;
      const axisObj = layout[axis];
      const nextAxis = axisObj && typeof axisObj === "object" ? { ...(axisObj as Record<string, unknown>) } : {};
      nextAxis.range = [saved[0], saved[1]];
      nextAxis.autorange = false;
      layout[axis] = nextAxis;
    }

    if (mode === "skymap") {
      // Remove backend animation buttons (Pause/Play) per UI requirement.
      delete layout.updatemenus;
      const xAxisObj =
        layout.xaxis && typeof layout.xaxis === "object" ? { ...(layout.xaxis as Record<string, unknown>) } : {};
      const yAxisObj =
        layout.yaxis && typeof layout.yaxis === "object" ? { ...(layout.yaxis as Record<string, unknown>) } : {};
      const xRange = parseAxisRange((xAxisObj as { range?: unknown }).range);
      const yRange = parseAxisRange((yAxisObj as { range?: unknown }).range);

      if (xRange && yRange) {
        const xCenter = 0.5 * (xRange[0] + xRange[1]);
        const yCenter = 0.5 * (yRange[0] + yRange[1]);
        const span = Math.max(Math.abs(xRange[1] - xRange[0]), Math.abs(yRange[1] - yRange[0]));
        const half = span * 0.5;
        xAxisObj.range = [xCenter - half, xCenter + half];
        yAxisObj.range = [yCenter - half, yCenter + half];
        xAxisObj.autorange = false;
        yAxisObj.autorange = false;
      }

      xAxisObj.constrain = "domain";
      yAxisObj.constrain = "domain";
      yAxisObj.scaleanchor = "x";
      yAxisObj.scaleratio = 1;
      layout.xaxis = xAxisObj;
      layout.yaxis = yAxisObj;
    }

    if (mode !== "skymap") {
      const legendFontSize = legendFontSizeForWidth(plotWidthPx);
      const legendObj =
        layout.legend && typeof layout.legend === "object" ? { ...(layout.legend as Record<string, unknown>) } : {};
      const legendFontObj =
        legendObj.font && typeof legendObj.font === "object" ? { ...(legendObj.font as Record<string, unknown>) } : {};
      legendFontObj.size = legendFontSize;
      legendObj.font = legendFontObj;
      layout.legend = legendObj;
    }

    // Keep UI interactions (zoom/pan) stable across realtime data updates.
    layout.uirevision = `${mode}-zoom-state`;
    return { ...figure, data: mergedData, layout };
  }, [compareEnabled, compareResult, compareResultMode, mode, plotWidthPx, result, resultMode, selectedCompareBookmark, zoomRevision]);
  const deferredFigure = useDeferredValue(displayFigure);

  const figureCaption = useMemo(() => {
    const title =
      mode === "lightcurve"
        ? "Multi-band GRB afterglow light curves"
        : mode === "spectrum"
          ? "GRB afterglow spectra"
          : "GRB afterglow sky images";
    const mediumLabel =
      shared.medium_type === "ISM" ? "uniform ISM" : shared.medium_type === "Wind" ? "stellar-wind medium" : "wind-bubble medium";
    const mediumParam =
      shared.medium_type === "ISM"
        ? (
          <>
            n<sub>ism</sub>={formatParamValueNode(shared.n_ism)} cm<sup>-3</sup>
          </>
        )
        : shared.medium_type === "Wind"
          ? (
            <>
              A<sub>*</sub>={formatParamValueNode(shared.A_star)}
            </>
          )
          : (
            <>
              A<sub>*</sub>={formatParamValueNode(shared.A_star)}, n<sub>floor</sub>={formatParamValueNode(shared.n_ism)} cm<sup>-3</sup>
            </>
          );
    const reverseShockTerms: ReactNode[] = [];
    if (shared.enable_rvs) {
      reverseShockTerms.push(
        <>
          jet duration={formatParamValueNode(shared.duration)} s
        </>,
      );
      reverseShockTerms.push(
        <>
          ε<sub>e,r</sub>={formatParamValueNode(shared.eps_e_r)}
        </>,
      );
      reverseShockTerms.push(
        <>
          ε<sub>B,r</sub>={formatParamValueNode(shared.eps_B_r)}
        </>,
      );
      reverseShockTerms.push(
        <>
          p<sub>r</sub>={formatParamValueNode(shared.p_r)}
        </>,
      );
      reverseShockTerms.push(
        <>
          ξ<sub>e,r</sub>={formatParamValueNode(shared.xi_e_r)}
        </>,
      );
      if (shared.rvs_ssc) reverseShockTerms.push("SSC_r enabled");
      if (shared.rvs_kn) reverseShockTerms.push("KN_r enabled");
    }

    return (
      <>
        {title} generated with VegasAfterglow for a {shared.jet_type.toLowerCase()} jet in a {mediumLabel}, with E
        <sub>iso</sub>={formatParamValueNode(shared.E_iso)} erg, Γ<sub>0</sub>={formatParamValueNode(shared.Gamma0)}, θ
        <sub>c</sub>={formatParamValueNode(shared.theta_c)} rad, {mediumParam}, ε<sub>e</sub>={formatParamValueNode(shared.eps_e)}, ε
        <sub>B</sub>={formatParamValueNode(shared.eps_B)}, p={formatParamValueNode(shared.p)}, ξ<sub>e</sub>=
        {formatParamValueNode(shared.xi_e)}, d<sub>L</sub>={formatParamValueNode(shared.d_L_mpc)} Mpc, and θ<sub>obs</sub>=
        {formatParamValueNode(shared.theta_obs)} rad
        {shared.enable_rvs ? (
          <>
            ; reverse shock with{" "}
            {reverseShockTerms.map((item, idx) => (
              <span key={`rvs-node-${idx}`}>
                {idx > 0 ? ", " : ""}
                {item}
              </span>
            ))}
          </>
        ) : null}
        .
      </>
    );
  }, [mode, shared]);

  const clearSavedXAxis = useCallback(
    (targetMode: Mode) => {
      const current = axisRangesRef.current[targetMode];
      if (current.xaxis === null) return;
      axisRangesRef.current[targetMode] = { ...current, xaxis: null };
      if (mode === targetMode) {
        setZoomRevision((v) => v + 1);
      }
    },
    [mode],
  );

  const clearSavedYAxes = useCallback(
    (targetMode: Mode) => {
      const current = axisRangesRef.current[targetMode];
      if (current.yaxis === null && current.yaxis2 === null) return;
      axisRangesRef.current[targetMode] = { ...current, yaxis: null, yaxis2: null };
      if (mode === targetMode) {
        setZoomRevision((v) => v + 1);
      }
    },
    [mode],
  );

  function updateObsGroup(isLc: boolean, index: number, patch: Partial<ObservationGroup>) {
    if (isLc) {
      setLcObsGroups((prev) => prev.map((g, i) => (i === index ? { ...g, ...patch } : g)));
    } else {
      setSedObsGroups((prev) => prev.map((g, i) => (i === index ? { ...g, ...patch } : g)));
    }
  }

  function removeObsGroup(isLc: boolean, index: number) {
    if (isLc) {
      setLcObsGroups((prev) => prev.filter((_, i) => i !== index));
      setActiveLcObsTab((prev) => {
        if (index < prev) return prev - 1;
        if (index === prev) return Math.max(0, prev - 1);
        return prev;
      });
    } else {
      setSedObsGroups((prev) => prev.filter((_, i) => i !== index));
      setActiveSedObsTab((prev) => {
        if (index < prev) return prev - 1;
        if (index === prev) return Math.max(0, prev - 1);
        return prev;
      });
    }
  }

  function addObsGroup(isLc: boolean) {
    if (isLc) {
      setLcObsGroups((prev) => [...prev, defaultObsGroup(true, prev.length + 1)]);
    } else {
      setSedObsGroups((prev) => [...prev, defaultObsGroup(false, prev.length + 1)]);
    }
  }

  function toggleInstrument(isLc: boolean, name: string, checked: boolean) {
    const setter = isLc ? setLcInstruments : setSedInstruments;
    setter((prev) => {
      if (checked) {
        return prev.includes(name) ? prev : [...prev, name];
      }
      return prev.filter((v) => v !== name);
    });
  }

  async function handleObservationUpload(isLc: boolean, index: number, file: File) {
    setError("");
    try {
      const parsed = await parseObservationUpload(file);
      updateObsGroup(isLc, index, { text: parsed });
    } catch (err) {
      const message = err instanceof Error ? err.message : "Failed to parse observation file.";
      setError(`Upload failed: ${message}`);
    }
  }

  function renderObservationEditor(isLc: boolean) {
    const groups = isLc ? lcObsGroups : sedObsGroups;
    const activeTab = isLc ? activeLcObsTab : activeSedObsTab;
    const setActiveTab = isLc ? setActiveLcObsTab : setActiveSedObsTab;
    const xOptions = isLc ? TIME_UNIT_OPTIONS : FREQ_UNIT_OPTIONS;
    const xName = isLc ? "t" : "nu";
    const xUnitLabel = isLc ? "t unit" : "ν unit";
    const rowLabel = isLc ? "Rows: t, value, error" : "Rows: ν, value, error";
    const obsHelpText = `One row per line. Columns: ${xName}, value, error (optional). Separated by space, tab, or comma.`;
    const obsPlaceholder = `${xName}  value  error\n1e4  0.5  0.1\n1e5  0.3  0.05`;
    const activeIndex = groups.length === 0 ? -1 : Math.max(0, Math.min(activeTab, groups.length - 1));
    const activeGroup = activeIndex >= 0 ? groups[activeIndex] : null;

    return (
      <details className="sb-expander">
        <summary>Observation Data</summary>
        <button
          className="sb-small-btn"
          type="button"
          onClick={() => {
            const nextIndex = groups.length;
            addObsGroup(isLc);
            setActiveTab(nextIndex);
          }}
        >
          + Add group
        </button>
        {groups.length === 0 ? <p className="sb-muted">No groups added.</p> : null}
        {groups.length > 0 ? (
          <div className="obs-tabs" role="tablist" aria-label={isLc ? "Light curve data groups" : "Spectrum data groups"}>
            {groups.map((group, idx) => (
              <button
                key={`${group.legend}-${idx}`}
                type="button"
                role="tab"
                aria-selected={idx === activeIndex}
                className={`obs-tab${idx === activeIndex ? " active" : ""}`}
                onClick={() => setActiveTab(idx)}
              >
                {group.legend.trim() || `data ${idx + 1}`}
              </button>
            ))}
          </div>
        ) : null}
        {activeGroup ? (
          <div className="obs-card" key={`${activeGroup.legend}-${activeIndex}`}>
            <div className="obs-card-row obs-card-row-meta">
              <label className="sb-field">
                <span className="sb-label">Legend</span>
                <input
                  value={activeGroup.legend}
                  onChange={(e) => updateObsGroup(isLc, activeIndex, { legend: e.target.value })}
                />
              </label>
              <label className="sb-checkbox-inline">
                <input
                  type="checkbox"
                  checked={activeGroup.visible}
                  onChange={(e) => updateObsGroup(isLc, activeIndex, { visible: e.target.checked })}
                />
                Visible
              </label>
            </div>
            <div className="obs-card-row obs-card-row-units">
              <label className="sb-field">
                <span className="sb-label">{xUnitLabel}</span>
                <select
                  value={activeGroup.x_unit}
                  onChange={(e) => updateObsGroup(isLc, activeIndex, { x_unit: e.target.value })}
                >
                  {xOptions.map((u) => (
                    <option key={u} value={u}>
                      {u}
                    </option>
                  ))}
                </select>
              </label>
              <label className="sb-field">
                <span className="sb-label">flux/flux den unit</span>
                <select
                  value={activeGroup.y_unit}
                  onChange={(e) => updateObsGroup(isLc, activeIndex, { y_unit: e.target.value })}
                >
                  {Y_UNIT_OPTIONS.map((u) => (
                    <option key={u} value={u}>
                      {u}
                    </option>
                  ))}
                </select>
              </label>
            </div>
            <label className="sb-field">
              <span className="sb-label sb-label-with-help">
                <span>{rowLabel}</span>
                <HelpHint text={obsHelpText} ariaLabel="Observation input help" />
              </span>
              <textarea
                rows={4}
                value={activeGroup.text}
                placeholder={obsPlaceholder}
                onChange={(e) => updateObsGroup(isLc, activeIndex, { text: e.target.value })}
              />
            </label>
            <div className="obs-card-actions">
              <label className="sb-upload-btn">
                Upload file
                <input
                  type="file"
                  accept=".csv,.txt,.dat,.xls,.xlsx,text/csv,text/plain,application/vnd.ms-excel,application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                  onChange={(e) => {
                    const file = e.target.files?.[0];
                    if (file) {
                      void handleObservationUpload(isLc, activeIndex, file);
                    }
                    e.currentTarget.value = "";
                  }}
                />
              </label>
              <button className="sb-danger-btn" type="button" onClick={() => removeObsGroup(isLc, activeIndex)}>
                Remove
              </button>
            </div>
            <p className="sb-muted sb-upload-hint">Supports .csv, .txt, .xls, .xlsx</p>
          </div>
        ) : null}
      </details>
    );
  }

  function renderDistanceObserverControls() {
    const formatZPlain = (value: number): string => {
      if (!Number.isFinite(value)) return "";
      if (value === 0) return "0";
      const abs = Math.abs(value);
      const exponent = Math.floor(Math.log10(abs));
      const decimals = Math.max(0, Math.min(12, 2 - exponent));
      return value.toFixed(decimals).replace(/\.?0+$/, "");
    };

    const convertedText =
      distanceDriver === "dL" ? (
        <>
          z<span className="sb-distance-converted-eq"> = </span>
          {formatZPlain(redshiftFromLuminosityDistanceMpc(shared.d_L_mpc))}
        </>
      ) : (
        <>
          d<sub>L</sub> (Mpc)
          <span className="sb-distance-converted-eq"> = </span>
          {normalizeDistanceMpc(luminosityDistanceMpcFromRedshift(shared.z))}
        </>
      );

    return (
      <div className="sb-stack sb-row-gap-top">
        <div className="sb-distance-panel">
          {distanceLinked ? (
            <div className="sb-distance-row sb-distance-row-linked">
              <div className="sb-distance-entry">
                <select
                  className="sb-distance-unit"
                  value={distanceDriver}
                  onChange={(e) => setDistanceDriver(e.target.value as DistanceDriver)}
                >
                  <option value="dL">dₗ (Mpc)</option>
                  <option value="z">z</option>
                </select>
                <input
                  className="sb-distance-number sb-distance-value"
                  type="number"
                  min={0}
                  max={distanceDriver === "dL" ? 1e7 : 20}
                  step={distanceDriver === "dL" ? 1 : 0.001}
                  value={distanceDriver === "dL" ? shared.d_L_mpc : formatZPlain(shared.z)}
                  onChange={(e) => {
                    const numeric = Number(e.target.value);
                    if (distanceDriver === "dL") {
                      setDistanceMpc(numeric);
                    } else {
                      setDistanceRedshift(numeric);
                    }
                  }}
                />
              </div>
              <span className="sb-muted sb-distance-converted">{convertedText}</span>
              <label className="sb-checkbox-inline sb-distance-unlock">
                <input type="checkbox" checked={!distanceLinked} onChange={(e) => setDistanceLinked(!e.target.checked)} />
                unlock
              </label>
            </div>
          ) : (
            <div className="sb-distance-row">
              <div className="sb-distance-unlocked-inline">
                <span className="sb-distance-prefix">
                  d<sub>L</sub> (Mpc)
                </span>
                <input
                  className="sb-distance-number sb-distance-value"
                  type="number"
                  min={0}
                  max={1e7}
                  step={1}
                  value={shared.d_L_mpc}
                  onChange={(e) => setDistanceMpc(Number(e.target.value))}
                />
                <span className="sb-distance-prefix">z</span>
                <input
                  className="sb-distance-number sb-distance-value"
                  type="number"
                  min={0}
                  max={20}
                  step={0.001}
                  value={formatZPlain(shared.z)}
                  onChange={(e) => setDistanceRedshift(Number(e.target.value))}
                />
              </div>
              <label className="sb-checkbox-inline sb-distance-unlock">
                <input type="checkbox" checked={!distanceLinked} onChange={(e) => setDistanceLinked(!e.target.checked)} />
                unlock
              </label>
            </div>
          )}
        </div>
        <SliderField
          label={
            <>
              θ<sub>obs</sub> (rad)
            </>
          }
          min={0}
          max={1.57}
          step={0.01}
          value={shared.theta_obs}
          onChange={(v) => setSharedField("theta_obs", v)}
        />
      </div>
    );
  }

  function renderFluxUnitControl() {
    return (
      <label className="sb-field">
        <span className="sb-label">Flux</span>
        <select
          value={shared.flux_unit}
          onChange={(e) => setSharedField("flux_unit", e.target.value as SharedParams["flux_unit"])}
        >
          <option value="mJy">mJy</option>
          <option value="Jy">Jy</option>
          <option value="uJy">uJy</option>
          <option value="cgs">cgs</option>
          <option value="AB mag">AB mag</option>
        </select>
      </label>
    );
  }

  function renderPlotUnitControls(timeDisabled: boolean, showTime = true) {
    if (!showTime) {
      return renderFluxUnitControl();
    }

    return (
      <div className="sb-row-2">
        {renderFluxUnitControl()}
        <label className="sb-field">
          <span className="sb-label">Time</span>
          <select
            value={shared.time_unit}
            disabled={timeDisabled}
            onChange={(e) => setSharedField("time_unit", e.target.value as SharedParams["time_unit"])}
          >
            {TIME_UNIT_OPTIONS.map((u) => (
              <option key={u} value={u}>
                {u}
              </option>
            ))}
          </select>
        </label>
      </div>
    );
  }

  function renderSharedControls(options?: { hideDistanceObserver?: boolean; hidePlotUnits?: boolean }) {
    const isSky = mode === "skymap";
    const isLc = mode === "lightcurve";

    return (
      <>
        {!options?.hideDistanceObserver ? <div className="sb-stack">{renderDistanceObserverControls()}</div> : null}

        {!isSky && !options?.hidePlotUnits ? renderPlotUnitControls(!isLc, isLc) : null}

        <details className="sb-expander sb-major" open>
          <summary>Jet</summary>
          <div className="sb-group-content">
            <div className="sb-row-2 sb-row-jet">
              <div className="sb-field">
                <span className="sb-label">Type</span>
                <div className="mode-radios sb-choice-radios">
                  {JET_TYPE_OPTIONS.map((jetType) => (
                    <label key={jetType} className="mode-radio">
                      <input
                        type="radio"
                        name="jet-type"
                        checked={shared.jet_type === jetType}
                        onChange={() => setSharedField("jet_type", jetType)}
                      />
                      <span>{jetType}</span>
                    </label>
                  ))}
                </div>
              </div>
              <SliderField
                label={
                  <>
                    θ<sub>c</sub>
                  </>
                }
                min={0.01}
                max={1.0}
                step={0.01}
                value={shared.theta_c}
                onChange={(v) => setSharedField("theta_c", v)}
              />
            </div>

            <div className="sb-row-2">
              <LogSliderField
                label={
                  <>
                    log10(E<sub>iso</sub> (erg))
                  </>
                }
                minExp={48}
                maxExp={57}
                step={0.1}
                value={shared.E_iso}
                defaultExp={52}
                onChange={(v) => setSharedField("E_iso", v)}
              />
              <LogSliderField
                label={
                  <>
                    log10(Γ<sub>0</sub>)
                  </>
                }
                minExp={1}
                maxExp={3.5}
                step={0.05}
                value={shared.Gamma0}
                defaultExp={Math.log10(300)}
                onChange={(v) => setSharedField("Gamma0", v)}
              />
            </div>

            <label className="sb-checkbox-inline">
              <input
                type="checkbox"
                checked={shared.spreading}
                onChange={(e) => setSharedField("spreading", e.target.checked)}
              />
              Spreading
            </label>

            {shared.jet_type === "Power-law" ? (
              <div className="sb-row-2">
                <SliderField
                  label={
                    <>
                      k<sub>e</sub>
                    </>
                  }
                  min={0.5}
                  max={10}
                  step={0.1}
                  value={shared.k_e}
                  onChange={(v) => setSharedField("k_e", v)}
                />
                <SliderField
                  label={
                    <>
                      k<sub>g</sub>
                    </>
                  }
                  min={0.5}
                  max={10}
                  step={0.1}
                  value={shared.k_g}
                  onChange={(v) => setSharedField("k_g", v)}
                />
              </div>
            ) : null}

            {shared.jet_type === "Two-component" ? (
              <>
                <div className="sb-row-2">
                  <SliderField
                    label={
                      <>
                        θ<sub>w</sub>
                      </>
                    }
                    min={0.05}
                    max={1.5}
                    step={0.01}
                    value={shared.theta_w}
                    onChange={(v) => setSharedField("theta_w", v)}
                  />
                  <LogSliderField
                    label={
                      <>
                        log10(E<sub>iso,w</sub>)
                      </>
                    }
                    minExp={48}
                    maxExp={55}
                    step={0.1}
                    value={shared.E_iso_w}
                    defaultExp={51}
                    onChange={(v) => setSharedField("E_iso_w", v)}
                  />
                </div>
                <LogSliderField
                  label={
                    <>
                      log10(Γ<sub>0,w</sub>)
                    </>
                  }
                  minExp={1}
                  maxExp={3}
                  step={0.05}
                  value={shared.Gamma0_w}
                  defaultExp={2}
                  onChange={(v) => setSharedField("Gamma0_w", v)}
                />
              </>
            ) : null}
          </div>
        </details>

        <details className="sb-expander sb-major">
          <summary>Medium</summary>
          <div className="sb-group-content">
            <div className="sb-row-2 sb-row-medium">
              <div className="sb-field">
                <span className="sb-label">Type</span>
                <div className="mode-radios sb-choice-radios">
                  {MEDIUM_TYPE_OPTIONS.map((mediumType) => (
                    <label key={mediumType.value} className="mode-radio">
                      <input
                        type="radio"
                        name="medium-type"
                        checked={shared.medium_type === mediumType.value}
                        onChange={() => setSharedField("medium_type", mediumType.value)}
                      />
                      <span>{mediumType.label}</span>
                    </label>
                  ))}
                </div>
              </div>
              {shared.medium_type === "ISM" ? (
                <LogSliderField
                  label={
                    <>
                      log10(n<sub>ism</sub> (cm⁻³))
                    </>
                  }
                  minExp={-5}
                  maxExp={5}
                  step={0.1}
                  value={shared.n_ism}
                  defaultExp={0}
                  onChange={(v) => setSharedField("n_ism", v)}
                />
              ) : (
                <LogSliderField
                  label={
                    <>
                      log10(A<sub>*</sub>)
                    </>
                  }
                  minExp={-3}
                  maxExp={2}
                  step={0.1}
                  value={shared.A_star}
                  defaultExp={-1}
                  onChange={(v) => setSharedField("A_star", v)}
                />
              )}
            </div>

            {shared.medium_type !== "ISM" ? (
              <SliderField
                label={
                  <>
                    k<sub>m</sub>
                  </>
                }
                min={0}
                max={4}
                step={0.1}
                value={shared.k_m}
                onChange={(v) => setSharedField("k_m", v)}
              />
            ) : null}
            {shared.medium_type === "Wind bubble" ? (
              <LogSliderField
                label={
                  <>
                    log10(n<sub>floor</sub> (cm⁻³))
                  </>
                }
                minExp={-8}
                maxExp={5}
                step={0.1}
                value={shared.n_ism}
                defaultExp={0}
                onChange={(v) => setSharedField("n_ism", v)}
              />
            ) : null}
          </div>
        </details>

        <details className="sb-expander sb-major">
          <summary>Radiation</summary>
          <div className="sb-group-content">
            <div className="sb-row-2">
              <LogSliderField
                label={
                  <>
                    log10(ε<sub>e</sub>)
                  </>
                }
                minExp={-5}
                maxExp={0}
                step={0.05}
                value={shared.eps_e}
                defaultExp={-1}
                onChange={(v) => setSharedField("eps_e", v)}
              />
              <LogSliderField
                label={
                  <>
                    log10(ε<sub>B</sub>)
                  </>
                }
                minExp={-6}
                maxExp={0}
                step={0.05}
                value={shared.eps_B}
                defaultExp={-3}
                onChange={(v) => setSharedField("eps_B", v)}
              />
            </div>
            <div className="sb-row-2">
              <SliderField label="p" min={2.01} max={3.0} step={0.01} value={shared.p} onChange={(v) => setSharedField("p", v)} />
              <LogSliderField
                label={
                  <>
                    log10(ξ<sub>e</sub>)
                  </>
                }
                minExp={-3}
                maxExp={0}
                step={0.05}
                value={shared.xi_e}
                defaultExp={0}
                onChange={(v) => setSharedField("xi_e", v)}
              />
            </div>
            <div className="sb-row-2 sb-check-row">
              <label className="sb-checkbox-inline">
                <input type="checkbox" checked={shared.ssc} onChange={(e) => setSharedField("ssc", e.target.checked)} />
                SSC
              </label>
              <label className="sb-checkbox-inline">
                <input type="checkbox" checked={shared.kn} onChange={(e) => setSharedField("kn", e.target.checked)} />
                KN
              </label>
            </div>
          </div>
        </details>

        <details className="sb-expander">
          <summary>Reverse Shock</summary>
          <label className="sb-checkbox-inline">
            <input
              type="checkbox"
              checked={shared.enable_rvs}
              onChange={(e) => setSharedField("enable_rvs", e.target.checked)}
            />
            Enable
          </label>
          {shared.enable_rvs ? (
            <>
              <LogSliderField
                label="log10(Jet Duration (s))"
                minExp={-1}
                maxExp={4}
                step={0.05}
                value={shared.duration}
                defaultExp={0}
                onChange={(v) => setSharedField("duration", v)}
              />
              <div className="sb-row-2">
                <LogSliderField
                  label={
                    <>
                      log10(ε<sub>e,r</sub>)
                    </>
                  }
                  minExp={-5}
                  maxExp={0}
                  step={0.05}
                  value={shared.eps_e_r}
                  defaultExp={-1}
                  onChange={(v) => setSharedField("eps_e_r", v)}
                />
                <LogSliderField
                  label={
                    <>
                      log10(ε<sub>B,r</sub>)
                    </>
                  }
                  minExp={-6}
                  maxExp={0}
                  step={0.05}
                  value={shared.eps_B_r}
                  defaultExp={-3}
                  onChange={(v) => setSharedField("eps_B_r", v)}
                />
              </div>
              <div className="sb-row-2">
                <SliderField
                  label={
                    <>
                      p<sub>r</sub>
                    </>
                  }
                  min={2.01}
                  max={3}
                  step={0.01}
                  value={shared.p_r}
                  onChange={(v) => setSharedField("p_r", v)}
                />
                <LogSliderField
                  label={
                    <>
                      log10(ξ<sub>e,r</sub>)
                    </>
                  }
                  minExp={-3}
                  maxExp={0}
                  step={0.05}
                  value={shared.xi_e_r}
                  defaultExp={0}
                  onChange={(v) => setSharedField("xi_e_r", v)}
                />
              </div>
              <div className="sb-row-2 sb-check-row">
                <label className="sb-checkbox-inline">
                  <input type="checkbox" checked={shared.rvs_ssc} onChange={(e) => setSharedField("rvs_ssc", e.target.checked)} />
                  SSC
                </label>
                <label className="sb-checkbox-inline">
                  <input type="checkbox" checked={shared.rvs_kn} onChange={(e) => setSharedField("rvs_kn", e.target.checked)} />
                  KN
                </label>
              </div>
            </>
          ) : null}
        </details>

        {!isSky ? (
          <details className="sb-expander">
            <summary>Instrument Sensitivities</summary>
            <div className="instrument-grid">
              {instrumentGroups.length === 0 ? <p className="sb-muted">No instrument list loaded.</p> : null}
              {instrumentGroups.map((group) => (
                <section key={group.label} className="instrument-group">
                  <p className="instrument-group-title">{group.label}</p>
                  <div className="instrument-group-grid">
                    {group.items.map((name) => {
                      const selected = isLc ? lcInstruments.includes(name) : sedInstruments.includes(name);
                      return (
                        <label key={name} className="sb-checkbox-inline instrument-option">
                          <input
                            type="checkbox"
                            checked={selected}
                            onChange={(e) => toggleInstrument(isLc, name, e.target.checked)}
                          />
                          {name}
                        </label>
                      );
                    })}
                  </div>
                </section>
              ))}
            </div>
          </details>
        ) : null}

        {!isSky ? renderObservationEditor(isLc) : null}

        <details className="sb-expander">
          <summary>Resolutions</summary>
          <div className="sb-stack sb-resolution-ppd-stack">
            <SliderField
              label={
                <>
                  φ ppd
                </>
              }
              min={0.05}
              max={0.5}
              step={0.05}
              decimals={2}
              value={shared.res_phi}
              onChange={(v) => setSharedField("res_phi", v)}
            />
            <SliderField
              label={
                <>
                  θ ppd
                </>
              }
              min={0.1}
              max={2}
              step={0.05}
              decimals={2}
              value={shared.res_theta}
              onChange={(v) => setSharedField("res_theta", v)}
            />
            <SliderField
              label={
                <>
                  t ppd
                </>
              }
              min={1}
              max={30}
              step={0.5}
              decimals={1}
              value={shared.res_t}
              onChange={(v) => setSharedField("res_t", v)}
            />
          </div>
        </details>
      </>
    );
  }

  const fullComputationSpec = useMemo<ComputationSpec>(() => {
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
  }, [
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
  ]);

  const computationSpec = useMemo<ComputationSpec>(
    () => buildInteractiveSpec(fullComputationSpec, ENABLE_INTERACTIVE_DOWNSAMPLE && sliderInteracting),
    [fullComputationSpec, sliderInteracting],
  );

  function renderBookmarkPanel() {
    const compareSupported = mode !== "skymap";
    return (
      <details className="sb-expander">
        <summary>Saved Setups</summary>
        <div className="bookmark-help-row">
          <HelpHint
            text="Save current parameters locally, then use each setup's Link button to share a URL so others can restore the same parameters."
            ariaLabel="Saved setup sharing help"
          />
        </div>
        <div className="bookmark-save-row">
          <input
            value={bookmarkNameDraft}
            onChange={(e) => setBookmarkNameDraft(e.target.value)}
            placeholder="Setup name (optional)"
          />
          <button className="sb-small-btn" type="button" onClick={saveCurrentBookmark}>
            Save Setup
          </button>
        </div>
        {setupStatusText ? <p className="sb-muted bookmark-status">{setupStatusText}</p> : null}

        {bookmarks.length === 0 ? <p className="sb-muted">No saved setups yet.</p> : null}
        {bookmarks.length > 0 ? (
          <div className="bookmark-list">
            {bookmarks.map((bookmark) => (
              <div className="bookmark-row" key={bookmark.id}>
                <span className="bookmark-name" title={bookmark.name}>
                  {bookmark.name}
                </span>
                <div className="bookmark-row-actions">
                  <button className="sb-small-btn" type="button" onClick={() => loadBookmarkById(bookmark.id)}>
                    Load
                  </button>
                  <button className="sb-small-btn" type="button" onClick={() => void copyShareLinkForSnapshot(bookmark.snapshot)}>
                    URL link
                  </button>
                  <button className="sb-danger-btn" type="button" onClick={() => removeBookmarkById(bookmark.id)}>
                    Remove
                  </button>
                </div>
              </div>
            ))}
          </div>
        ) : null}

        <label className="sb-checkbox-inline">
          <input
            type="checkbox"
            checked={compareEnabled}
            disabled={bookmarks.length === 0 || !compareSupported}
            onChange={(e) => setCompareEnabled(e.target.checked)}
          />
          Overlay saved setup
        </label>
        {!compareSupported ? <p className="sb-muted">Compare overlay supports Light Curve and Spectrum only.</p> : null}

        {compareEnabled && compareSupported ? (
          <>
            <label className="sb-field">
              <span className="sb-label">Setup to overlay</span>
              <select value={compareBookmarkId} onChange={(e) => setCompareBookmarkId(e.target.value)}>
                {bookmarks.map((bookmark) => (
                  <option key={bookmark.id} value={bookmark.id}>
                    {bookmark.name}
                  </option>
                ))}
              </select>
            </label>
            {compareRunning ? <p className="sb-muted">Updating overlay...</p> : null}
            {compareError ? <p className="sb-compare-error">{compareError}</p> : null}
          </>
        ) : null}
      </details>
    );
  }

  function renderModeSpecific() {
    const skyPixelMax = skyAnimate ? SKY_MAX_PIXEL_ANIMATE : SKY_MAX_PIXEL_STATIC;
    const skyPixelOptions = SKY_PIXEL_OPTIONS.filter((v) => v <= skyPixelMax);

    if (mode === "lightcurve") {
      return (
        <>
          <div className="sb-stack">
            <label className="sb-field">
              <span className="sb-label sb-label-with-help">
                <span>Frequencies</span>
                <HelpHint text={FREQ_HELP_TEXT} ariaLabel="Frequency input help" />
              </span>
              <input
                value={lcFreq}
                onChange={(e) => setLcFreq(e.target.value)}
                placeholder="e.g. 1e9, R, 1keV, XRT, [0.3keV,10keV]"
              />
            </label>
            {renderDistanceObserverControls()}
          </div>

          <div className="sb-stack">
            <div className="sb-row-2 sb-row-gap-top">
              <LogSliderField
                label={
                  <>
                    log10(t<sub>min</sub> (s))
                  </>
                }
                minExp={-1}
                maxExp={6}
                step={0.05}
                value={lcTMin}
                defaultExp={0}
                onChange={setLcTMin}
              />
              <LogSliderField
                label={
                  <>
                    log10(t<sub>max</sub> (s))
                  </>
                }
                minExp={3}
                maxExp={10}
                step={0.05}
                value={lcTMax}
                defaultExp={8}
                onChange={setLcTMax}
              />
            </div>
            <SliderField
              label={
                <>
                  t points
                </>
              }
              min={50}
              max={300}
              step={10}
              decimals={0}
              value={shared.num_t}
              onChange={(v) => setSharedField("num_t", Math.round(v))}
            />
            {renderPlotUnitControls(false)}
          </div>

          {renderSharedControls({ hideDistanceObserver: true, hidePlotUnits: true })}
        </>
      );
    }

    if (mode === "spectrum") {
      return (
        <>
          <div className="sb-stack">
            <label className="sb-field">
              <span className="sb-label">
                t<sub>obs</sub> (s)
              </span>
              <input
                value={sedTimesDraft}
                onChange={(e) => setSedTimesDraft(e.target.value)}
                onBlur={() => setSedTimes(sedTimesDraft)}
                placeholder="e.g. 1e3, 1e4, 1e5"
              />
            </label>
            {renderDistanceObserverControls()}
          </div>

          <div className="sb-stack">
            <div className="sb-row-2 sb-row-gap-top">
              <LogSliderField
                label={
                  <>
                    log10(ν<sub>min</sub> (Hz))
                  </>
                }
                minExp={6}
                maxExp={20}
                step={0.05}
                value={sedNuMin}
                defaultExp={8}
                onChange={setSedNuMin}
              />
              <LogSliderField
                label={
                  <>
                    log10(ν<sub>max</sub> (Hz))
                  </>
                }
                minExp={10}
                maxExp={35}
                step={0.05}
                value={sedNuMax}
                defaultExp={20}
                onChange={setSedNuMax}
              />
            </div>
            <div className="sb-row-2 sb-row-fit-right">
              <SliderField
                label={
                  <>
                    ν points
                  </>
                }
                min={50}
                max={500}
                step={10}
                decimals={0}
                value={sedNumNu}
                onChange={(v) => setSedNumNu(Math.round(v))}
              />
              <label className="sb-checkbox-inline sb-nufnu-inline">
                <input
                  type="checkbox"
                  checked={sedNuFNu}
                  disabled={shared.flux_unit === "AB mag"}
                  onChange={(e) => setSedNuFNu(e.target.checked)}
                />
                ν F<sub>ν</sub>
              </label>
            </div>
            <div className="sb-row-2">
              <label className="sb-field">
                <span className="sb-label">ν unit</span>
                <select value={sedFreqUnit} onChange={(e) => setSedFreqUnit(e.target.value)}>
                  {FREQ_UNIT_OPTIONS.map((u) => (
                    <option key={u} value={u}>
                      {u}
                    </option>
                  ))}
                </select>
              </label>
              {renderFluxUnitControl()}
            </div>
            {shared.flux_unit === "AB mag" ? <p className="sb-muted">ν F<sub>ν</sub> is disabled for AB mag.</p> : null}
          </div>

          {renderSharedControls({ hideDistanceObserver: true, hidePlotUnits: true })}
        </>
      );
    }

    return (
      <>
        <label className="sb-checkbox-inline">
          <input
            type="checkbox"
            checked={skyAnimate}
            onChange={(e) => {
              const nextAnimate = e.target.checked;
              setSkyAnimate(nextAnimate);
              if (nextAnimate && skyNpixel > SKY_MAX_PIXEL_ANIMATE) {
                setSkyNpixel(SKY_MAX_PIXEL_ANIMATE);
              }
            }}
          />
          Animate
        </label>

        {skyAnimate ? (
          <>
            <div className="sb-row-2">
              <LogSliderField
                label={
                  <>
                    log10(t<sub>min</sub> (s))
                  </>
                }
                minExp={2}
                maxExp={9}
                step={0.05}
                value={skyTMin}
                defaultExp={4}
                onChange={setSkyTMin}
              />
              <LogSliderField
                label={
                  <>
                    log10(t<sub>max</sub> (s))
                  </>
                }
                minExp={2}
                maxExp={9}
                step={0.05}
                value={skyTMax}
                defaultExp={7}
                onChange={setSkyTMax}
              />
            </div>
            <SliderField
              label="Frames"
              min={3}
              max={SKY_MAX_N_FRAMES}
              step={1}
              decimals={0}
              value={skyNFrames}
              onChange={(v) => setSkyNFrames(Math.max(3, Math.min(SKY_MAX_N_FRAMES, Math.round(v))))}
            />
          </>
        ) : (
          <LogSliderField
            label={
              <>
                log10(t<sub>obs</sub> (s))
              </>
            }
            minExp={2}
            maxExp={9}
            step={0.05}
            value={skyTObs}
            defaultExp={6}
            onChange={setSkyTObs}
          />
        )}

        <label className="sb-field">
          <span className="sb-label sb-label-with-help">
            <span>
              ν<sub>obs</sub>
            </span>
            <HelpHint text={FREQ_HELP_TEXT} ariaLabel="Sky frequency input help" />
          </span>
          <input value={skyNuInput} onChange={(e) => setSkyNuInput(e.target.value)} placeholder="e.g. 1e9, 1GHz, 1keV" />
        </label>

        <div className="sb-row-2 sb-row-gap-top">
          <LogSliderField
            label="log10(FOV (μas))"
            minExp={1}
            maxExp={5}
            step={0.05}
            value={skyFov}
            defaultExp={Math.log10(500)}
            onChange={setSkyFov}
          />
          <label className="sb-field">
            <span className="sb-label">Pixels</span>
            <select value={String(skyNpixel)} onChange={(e) => setSkyNpixel(Number(e.target.value))}>
              {skyPixelOptions.map((v) => (
                <option key={v} value={v}>
                  {v}
                </option>
              ))}
            </select>
          </label>
        </div>

        {renderSharedControls()}
      </>
    );
  }

  const runComputation = useCallback(async (spec: ComputationSpec) => {
    const requestMode = spec.endpoint;
    const controller = new AbortController();
    activeRequestRef.current = controller;
    const requestSeq = ++requestSeqRef.current;

    setError("");

    try {
      const response = await fetchFromApi(spec.endpoint, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          ...spec.payload,
          include_figure: true,
          include_exports: false,
        }),
        signal: controller.signal,
      });

      if (!response.ok) {
        const details = await response.text();
        throw new Error(details || `HTTP ${response.status}`);
      }

      const data = (await response.json()) as RunResponse;
      if (requestSeq !== requestSeqRef.current) return;

      startFigureTransition(() => {
        setResult(data);
        setResultMode(requestMode);
      });
    } catch (err) {
      if (controller.signal.aborted) return;
      if (err instanceof DOMException && err.name === "AbortError") return;
      if (requestSeq !== requestSeqRef.current) return;

      const message = err instanceof Error ? err.message : "Unknown error";
      const lower = message.toLowerCase();
      const isNetworkFailure = lower.includes("load failed") || lower.includes("failed to fetch");
      setError(
        isNetworkFailure ? `Cannot reach backend. Checked ${apiCandidates.join(" or ")}.` : message,
      );
    } finally {
      if (activeRequestRef.current === controller) {
        activeRequestRef.current = null;
      }
    }
  }, [apiCandidates, fetchFromApi, startFigureTransition]);

  const runCompareComputation = useCallback(async (spec: ComputationSpec) => {
    compareActiveRequestRef.current?.abort();
    const controller = new AbortController();
    compareActiveRequestRef.current = controller;
    const requestSeq = ++compareRequestSeqRef.current;
    setCompareRunning(true);
    setCompareError("");

    try {
      const response = await fetchFromApi(spec.endpoint, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          ...spec.payload,
          include_figure: true,
          include_exports: false,
        }),
        signal: controller.signal,
      });
      if (!response.ok) {
        const details = await response.text();
        throw new Error(details || `HTTP ${response.status}`);
      }

      const data = (await response.json()) as RunResponse;
      if (requestSeq !== compareRequestSeqRef.current) return;
      setCompareResult(data);
      setCompareResultMode(spec.endpoint);
    } catch (err) {
      if (controller.signal.aborted) return;
      if (err instanceof DOMException && err.name === "AbortError") return;
      if (requestSeq !== compareRequestSeqRef.current) return;
      const message = err instanceof Error ? err.message : "Compare request failed";
      setCompareError(message);
      setCompareResult(null);
      setCompareResultMode(null);
    } finally {
      if (compareActiveRequestRef.current === controller) {
        compareActiveRequestRef.current = null;
      }
      if (requestSeq === compareRequestSeqRef.current) {
        setCompareRunning(false);
      }
    }
  }, [fetchFromApi]);

  const flushQueuedComputation = useCallback(async () => {
    if (requestInFlightRef.current) return;
    if (!pendingSpecRef.current) return;

    requestInFlightRef.current = true;
    setRunning(true);
    try {
      let nextSpec: ComputationSpec | null = pendingSpecRef.current;
      while (nextSpec) {
        pendingSpecRef.current = null;
        await runComputation(nextSpec);
        nextSpec = pendingSpecRef.current;
      }
    } finally {
      requestInFlightRef.current = false;
      setRunning(false);
    }
  }, [runComputation]);

  useEffect(() => {
    if (!bootReady) return;
    if (isEmptyLightcurveFrequencySpec(computationSpec)) {
      pendingSpecRef.current = null;
      if (runTimerRef.current !== null) {
        window.clearTimeout(runTimerRef.current);
        runTimerRef.current = null;
      }
      activeRequestRef.current?.abort();
      setRunning(false);
      setError("");
      startFigureTransition(() => {
        setResult(null);
        setResultMode(null);
      });
      return;
    }
    pendingSpecRef.current = computationSpec;
    if (runTimerRef.current !== null) {
      window.clearTimeout(runTimerRef.current);
    }
    const debounceMs = sliderInteracting ? AUTO_RUN_DEBOUNCE_SLIDING_MS : AUTO_RUN_DEBOUNCE_IDLE_MS;
    runTimerRef.current = window.setTimeout(() => {
      runTimerRef.current = null;
      void flushQueuedComputation();
    }, debounceMs);
    return () => {
      if (runTimerRef.current !== null) {
        window.clearTimeout(runTimerRef.current);
        runTimerRef.current = null;
      }
    };
  }, [bootReady, computationSpec, flushQueuedComputation, sliderInteracting]);

  useEffect(() => {
    if (compareTimerRef.current !== null) {
      window.clearTimeout(compareTimerRef.current);
      compareTimerRef.current = null;
    }

    if (!bootReady || !compareSpec) {
      compareActiveRequestRef.current?.abort();
      setCompareRunning(false);
      if (!compareEnabled || mode === "skymap") {
        setCompareResult(null);
        setCompareResultMode(null);
        setCompareError("");
      }
      return;
    }

    if (isEmptyLightcurveFrequencySpec(compareSpec)) {
      compareActiveRequestRef.current?.abort();
      setCompareRunning(false);
      setCompareResult(null);
      setCompareResultMode(null);
      setCompareError("");
      return;
    }

    const debounceMs = sliderInteracting ? AUTO_RUN_DEBOUNCE_SLIDING_MS : AUTO_RUN_DEBOUNCE_IDLE_MS;
    compareTimerRef.current = window.setTimeout(() => {
      compareTimerRef.current = null;
      void runCompareComputation(compareSpec);
    }, debounceMs);

    return () => {
      if (compareTimerRef.current !== null) {
        window.clearTimeout(compareTimerRef.current);
        compareTimerRef.current = null;
      }
    };
  }, [bootReady, compareEnabled, compareSpec, mode, runCompareComputation, sliderInteracting]);

  useEffect(() => {
    return () => {
      if (runTimerRef.current !== null) {
        window.clearTimeout(runTimerRef.current);
        runTimerRef.current = null;
      }
      pendingSpecRef.current = null;
      activeRequestRef.current?.abort();
      if (compareTimerRef.current !== null) {
        window.clearTimeout(compareTimerRef.current);
        compareTimerRef.current = null;
      }
      compareActiveRequestRef.current?.abort();
    };
  }, []);

  useEffect(() => {
    if (!running) {
      setShowRunning(false);
      return;
    }
    const timer = window.setTimeout(() => {
      setShowRunning(true);
    }, 120);
    return () => window.clearTimeout(timer);
  }, [running]);

  useEffect(() => {
    const sliderSelector = ".sb-slider input[type='range']";
    const markInteracting = () => {
      if (sliderInteractingRef.current) return;
      sliderInteractingRef.current = true;
      setSliderInteracting(true);
    };
    const clearInteracting = () => {
      if (!sliderInteractingRef.current) return;
      sliderInteractingRef.current = false;
      setSliderInteracting(false);
    };

    const handlePointerDown = (event: PointerEvent) => {
      const target = event.target as Element | null;
      if (!target?.closest(sliderSelector)) return;
      markInteracting();
    };

    window.addEventListener("pointerdown", handlePointerDown, true);
    window.addEventListener("pointerup", clearInteracting, true);
    window.addEventListener("pointercancel", clearInteracting, true);
    window.addEventListener("blur", clearInteracting);

    return () => {
      window.removeEventListener("pointerdown", handlePointerDown, true);
      window.removeEventListener("pointerup", clearInteracting, true);
      window.removeEventListener("pointercancel", clearInteracting, true);
      window.removeEventListener("blur", clearInteracting);
    };
  }, []);

  useEffect(() => {
    if (typeof window === "undefined" || typeof document === "undefined") return;
    const root = document.documentElement;
    let rafId: number | null = null;
    let settleTimerShort: number | null = null;
    let settleTimerLong: number | null = null;

    const applyViewportHeight = () => {
      const vh = window.innerHeight * 0.01;
      root.style.setProperty("--app-vh", `${vh}px`);
    };

    const updateViewportHeight = () => {
      if (rafId !== null) return;
      rafId = window.requestAnimationFrame(() => {
        rafId = null;
        applyViewportHeight();
        if (settleTimerShort !== null) {
          window.clearTimeout(settleTimerShort);
        }
        if (settleTimerLong !== null) {
          window.clearTimeout(settleTimerLong);
        }
        // Mobile browsers can report transient viewport sizes during rotation.
        settleTimerShort = window.setTimeout(applyViewportHeight, 120);
        settleTimerLong = window.setTimeout(applyViewportHeight, 320);
      });
    };

    updateViewportHeight();
    window.addEventListener("resize", updateViewportHeight);
    window.addEventListener("orientationchange", updateViewportHeight);
    const viewport = window.visualViewport;
    viewport?.addEventListener("resize", updateViewportHeight);

    return () => {
      window.removeEventListener("resize", updateViewportHeight);
      window.removeEventListener("orientationchange", updateViewportHeight);
      viewport?.removeEventListener("resize", updateViewportHeight);
      if (rafId !== null) {
        window.cancelAnimationFrame(rafId);
      }
      if (settleTimerShort !== null) {
        window.clearTimeout(settleTimerShort);
      }
      if (settleTimerLong !== null) {
        window.clearTimeout(settleTimerLong);
      }
      root.style.removeProperty("--app-vh");
    };
  }, []);

  useEffect(() => {
    if (typeof window === "undefined" || typeof document === "undefined") return;
    const root = document.documentElement;
    const syncViewport = () => {
      const vh = window.innerHeight * 0.01;
      root.style.setProperty("--app-vh", `${vh}px`);
    };
    const triggerResize = () => {
      window.dispatchEvent(new Event("resize"));
    };

    syncViewport();
    triggerResize();

    // Sidebar open/close transition updates plot container size after animation.
    const settleTimer = window.setTimeout(() => {
      syncViewport();
      triggerResize();
    }, 220);

    return () => {
      window.clearTimeout(settleTimer);
    };
  }, [sidebarOpen]);

  useEffect(() => {
    if (typeof document === "undefined") return;
    const prev = document.body.style.overflow;
    if (sidebarOpen) {
      document.body.style.overflow = "hidden";
    } else {
      document.body.style.overflow = "";
    }
    return () => {
      document.body.style.overflow = prev;
    };
  }, [sidebarOpen]);

  useEffect(() => {
    if (sedTimesDraft === sedTimes) return;
    const timer = window.setTimeout(() => {
      setSedTimes(sedTimesDraft);
    }, SPECTRUM_TEXT_COMMIT_DEBOUNCE_MS);
    return () => window.clearTimeout(timer);
  }, [sedTimesDraft, sedTimes]);

  useEffect(() => {
    if (!skyAnimate) return;
    if (skyNpixel > SKY_MAX_PIXEL_ANIMATE) {
      setSkyNpixel(SKY_MAX_PIXEL_ANIMATE);
    }
  }, [skyAnimate, skyNpixel]);

  useEffect(() => {
    clearSavedXAxis("lightcurve");
  }, [clearSavedXAxis, lcTMin, lcTMax, shared.time_unit]);

  useEffect(() => {
    clearSavedXAxis("spectrum");
  }, [clearSavedXAxis, sedNuMin, sedNuMax, sedFreqUnit]);

  useEffect(() => {
    clearSavedYAxes("spectrum");
  }, [clearSavedYAxes, sedNuFNu, shared.flux_unit]);

  useEffect(() => {
    // Sky image axis extents must refresh when FOV/grid changes.
    clearSavedXAxis("skymap");
    clearSavedYAxes("skymap");
  }, [clearSavedXAxis, clearSavedYAxes, skyFov, skyNpixel]);

  useEffect(() => {
    // Mode switch should start from fresh axes/zoom for that mode.
    axisRangesRef.current[mode] = { xaxis: null, yaxis: null, yaxis2: null };
    setZoomRevision((v) => v + 1);
  }, [mode]);

  const downloadCurrentExport = useCallback(
    async (kind: DownloadKind) => {
      setError("");
      setDownloading(kind);
      try {
        let content = result?.exports?.[kind];
        if (!content) {
          const response = await fetchFromApi(fullComputationSpec.endpoint, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({
              ...fullComputationSpec.payload,
              include_figure: false,
              include_exports: true,
              export_kinds: [kind],
            }),
          });
          if (!response.ok) {
            const details = await response.text();
            throw new Error(details || `HTTP ${response.status}`);
          }
          const data = (await response.json()) as RunResponse;
          content = data.exports?.[kind];
          if (data.exports) {
            if (kind !== "gif") {
              setResult((prev) => (prev ? { ...prev, exports: { ...(prev.exports ?? {}), ...data.exports } } : prev));
            }
          }
        }

        if (!content) {
          throw new Error(`No ${kind.toUpperCase()} export returned by backend.`);
        }
        if (kind === "csv") {
          downloadText(content, `afterglow_${mode}.csv`, "text/csv");
        } else if (kind === "json") {
          downloadText(content, `afterglow_${mode}.json`, "application/json");
        } else {
          downloadBase64(content, "afterglow_skymap.gif", "image/gif");
        }
      } catch (err) {
        const message = err instanceof Error ? err.message : "Download failed";
        setError(message);
      } finally {
        setDownloading(null);
      }
    },
    [fetchFromApi, fullComputationSpec.endpoint, fullComputationSpec.payload, mode, result?.exports],
  );

  const copyBibtex = useCallback(async () => {
    const copied = await copyTextToClipboard(BIBTEX_TEXT);
    setCiteLinkText(copied ? "copied" : "copy failed");
    window.setTimeout(() => setCiteLinkText("cite"), 1200);
  }, [copyTextToClipboard]);

  const apiStatusRows = useMemo(
    () => {
      const activeStatus = activeApiKey ? apiStatusByKey[activeApiKey] : undefined;
      let activeDisplayKey = activeApiKey;
      if (!apiStatusDisplaySet.has(activeDisplayKey) && activeStatus?.region) {
        const mapped = apiEndpoints.find((endpoint) => {
          if (!apiStatusDisplaySet.has(endpoint.key)) return false;
          return apiStatusByKey[endpoint.key]?.region === activeStatus.region;
        });
        if (mapped) {
          activeDisplayKey = mapped.key;
        }
      }

      return apiEndpoints.filter((endpoint) => apiStatusDisplaySet.has(endpoint.key)).map((endpoint) => {
        const status = apiStatusByKey[endpoint.key];
        const isActive = activeDisplayKey === endpoint.key;
        let statusText = "checking...";
        let statusClass = "pending";
        if (status) {
          if (status.up) {
            statusText = `${status.latencyMs ?? "--"} ms`;
            statusClass = "up";
          } else if (status.statusCode !== null) {
            statusText = `down (${status.statusCode})`;
            statusClass = "down";
          } else {
            statusText = "down";
            statusClass = "down";
          }
        }
        return {
          ...endpoint,
          isActive,
          regionText: status?.locationLabel ?? status?.region ?? "region unknown",
          statusClass,
          statusText,
        };
      });
    },
    [activeApiKey, apiEndpoints, apiStatusByKey, apiStatusDisplaySet],
  );
  const activeApiStatusRow = useMemo(() => {
    if (apiStatusRows.length === 0) return null;
    return apiStatusRows.find((row) => row.isActive) ?? apiStatusRows[0];
  }, [apiStatusRows]);
  const otherApiStatusRows = useMemo(() => {
    if (!activeApiStatusRow) return apiStatusRows;
    return apiStatusRows.filter((row) => row.key !== activeApiStatusRow.key);
  }, [activeApiStatusRow, apiStatusRows]);
  const probeOtherServersOnce = useCallback(() => {
    const cancelledRef = { current: false };
    const targets = otherApiStatusRows
      .map((row) => apiEndpoints.find((endpoint) => endpoint.key === row.key))
      .filter((endpoint): endpoint is ApiEndpoint => Boolean(endpoint));
    if (targets.length === 0) return;
    void Promise.all(targets.map((endpoint) => probeEndpoint(endpoint, cancelledRef)));
  }, [apiEndpoints, otherApiStatusRows, probeEndpoint]);

  return (
    <main className={`app-shell ${sidebarOpen ? "sidebar-open" : ""}`}>
      {!sidebarOpen ? (
        <button className="sidebar-toggle" onClick={() => setSidebarOpen(true)} aria-label="Open controls">
          Controls
        </button>
      ) : null}
      {!sidebarOpen && activeApiStatusRow ? (
        <div className="sidebar-status-mini" aria-live="polite">
          <span className="sidebar-api-working">Working Server:</span>
          <span className="sidebar-status-mini-region">{activeApiStatusRow.regionText}</span>
          <span className={`sidebar-api-pill ${activeApiStatusRow.statusClass}`}>{activeApiStatusRow.statusText}</span>
        </div>
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
        <div className="sidebar-api-status" aria-live="polite">
          {activeApiStatusRow ? (
            <>
              <div key={activeApiStatusRow.key} className="sidebar-api-row sidebar-api-row-primary active">
                <span className="sidebar-api-name">
                  <span className="sidebar-api-working">Working Server:</span>
                  <span className="sidebar-api-location">{activeApiStatusRow.regionText}</span>
                  <span className={`sidebar-api-pill ${activeApiStatusRow.statusClass}`}>
                    {activeApiStatusRow.statusText}
                  </span>
                </span>
                {otherApiStatusRows.length > 0 ? (
                  <details
                    className="sidebar-api-dropdown-inline"
                    open={showAllServerStatus}
                    onToggle={(event) => {
                      const open = (event.currentTarget as HTMLDetailsElement).open;
                      setShowAllServerStatus(open);
                      if (open) {
                        probeOtherServersOnce();
                      }
                    }}
                  >
                    <summary>Other</summary>
                    <div className="sidebar-api-dropdown-menu">
                      {otherApiStatusRows.map((row) => (
                        <div key={row.key} className="sidebar-api-row">
                          <span className="sidebar-api-name">
                            <span className="sidebar-api-location">{row.regionText}</span>
                            <span className={`sidebar-api-pill ${row.statusClass}`}>{row.statusText}</span>
                          </span>
                        </div>
                      ))}
                    </div>
                  </details>
                ) : null}
              </div>
            </>
          ) : null}
        </div>
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

        {renderBookmarkPanel()}
        {renderModeSpecific()}

        <div className="sidebar-footer">
          <div className="sidebar-downloads">
            {mode !== "skymap" ? (
              <button
                disabled={downloading !== null || !displayFigure?.data}
                onClick={() => void downloadCurrentExport("csv")}
              >
                {downloading === "csv" ? "Preparing Data CSV..." : "Download Data (CSV)"}
              </button>
            ) : (
              <button
                disabled={downloading !== null || !displayFigure?.data}
                onClick={() => void downloadCurrentExport("gif")}
              >
                {downloading === "gif" ? "Preparing GIF..." : "Save GIF"}
              </button>
            )}
            <button
              disabled={downloading !== null || !displayFigure?.data}
              onClick={() => void downloadCurrentExport("json")}
            >
              {downloading === "json" ? "Preparing Data JSON..." : "Download Data (JSON)"}
            </button>
          </div>
          <p className="sb-footer-text">
            VegasAfterglow v{appVersion}
            <br />
            This interactive tool provides a subset of VegasAfterglow features. For MCMC fitting, custom jet/medium
            models, and the full API, see the{" "}
            <a
              className="sb-footer-link"
              href="https://github.com/YihanWangAstro/VegasAfterglow"
              target="_blank"
              rel="noreferrer"
            >
              Python package
            </a>
            .
          </p>
          <p className="sb-footer-text">
            VegasAfterglow Webtool is currently supported by personal funding from the developers. If this tool is
            helpful, any support is appreciated, including sharing it with others or{" "}
            <button className="sb-footer-link sb-footer-link-btn" type="button" onClick={() => void copyBibtex()}>
              {citeLinkText}
            </button>{" "}
            our work in research.
          </p>
        </div>

        {showRunning || isFigurePending ? <div className="status">{`Running ${mode} ...`}</div> : null}
        {error ? <div className="error">{error}</div> : null}
      </aside>

      {sidebarOpen ? (
        <button className="sidebar-backdrop" onClick={() => setSidebarOpen(false)} aria-label="Close sidebar" />
      ) : null}

      <section ref={workspaceRef} className="workspace">
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
            <p className="muted">Adjust parameters in the sidebar to render plot.</p>
          </div>
        )}

      </section>
    </main>
  );
}
