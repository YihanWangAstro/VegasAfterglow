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

const FALLBACK_SHARED = {
  d_L_mpc: 100,
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
const INTERACTIVE_LIGHTCURVE_NUM_T_MAX = 120;
const INTERACTIVE_SPECTRUM_NUM_NU_MAX = 120;
const INTERACTIVE_SKY_PIXEL_MAX_STATIC = 256;
const INTERACTIVE_SKY_PIXEL_MAX_ANIMATE = 256;
const INTERACTIVE_SKY_FRAMES_MAX = 8;
const AUTO_RUN_DEBOUNCE_IDLE_MS = 10;
const AUTO_RUN_DEBOUNCE_SLIDING_MS = 20;
const SLIDER_COMMIT_INTERVAL_MS = 30;
const SPECTRUM_TEXT_COMMIT_DEBOUNCE_MS = 220;
const AXIS_EPS = 1e-6;
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
const MEDIUM_TYPE_OPTIONS: SharedParams["medium_type"][] = ["ISM", "Wind", "Wind bubble"];
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

function safeLog10(value: number, fallback: number): number {
  return Number.isFinite(value) && value > 0 ? Math.log10(value) : fallback;
}

function sliderFillStyle(value: number, min: number, max: number): CSSProperties {
  const span = max - min;
  const raw = span > 0 ? ((value - min) / span) * 100 : 0;
  const bounded = Math.max(0, Math.min(100, raw));
  return { "--fill-percent": `${bounded}%` } as CSSProperties;
}

function legendFontSizeForWidth(plotWidthPx: number): number {
  if (plotWidthPx <= 440) return 9;
  if (plotWidthPx <= 620) return 10;
  if (plotWidthPx <= 860) return 11;
  return 12;
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

function normalizeShared(shared: SharedParams, mode: Mode): SharedParams {
  const next = { ...shared };

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

  const handleStart = () => {
    draggingRef.current = true;
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

  const handleStart = () => {
    draggingRef.current = true;
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

  const [mode, setMode] = useState<Mode>("lightcurve");
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const [shared, setShared] = useState<SharedParams>(FALLBACK_SHARED);

  const [error, setError] = useState<string>("");
  const [running, setRunning] = useState(false);
  const [showRunning, setShowRunning] = useState(false);
  const [result, setResult] = useState<RunResponse | null>(null);
  const [resultMode, setResultMode] = useState<Mode | null>(null);
  const [plotWidthPx, setPlotWidthPx] = useState(960);
  const [isFigurePending, startFigureTransition] = useTransition();
  const [downloading, setDownloading] = useState<DownloadKind | null>(null);
  const [appVersion, setAppVersion] = useState(DEFAULT_APP_VERSION);
  const [citeLinkText, setCiteLinkText] = useState("cite");
  const [bootReady, setBootReady] = useState(false);
  const [obsStorageReady, setObsStorageReady] = useState(false);
  const [zoomRevision, setZoomRevision] = useState(0);
  const requestSeqRef = useRef(0);
  const requestInFlightRef = useRef(false);
  const pendingSpecRef = useRef<ComputationSpec | null>(null);
  const runTimerRef = useRef<number | null>(null);
  const activeRequestRef = useRef<AbortController | null>(null);
  const [sliderInteracting, setSliderInteracting] = useState(false);
  const sliderInteractingRef = useRef(false);
  const workspaceRef = useRef<HTMLElement | null>(null);
  const axisRangesRef = useRef<Record<Mode, AxisRanges>>({
    lightcurve: { xaxis: null, yaxis: null, yaxis2: null },
    spectrum: { xaxis: null, yaxis: null, yaxis2: null },
    skymap: { xaxis: null, yaxis: null, yaxis2: null },
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

  const setSharedField = useCallback(<K extends keyof SharedParams>(key: K, value: SharedParams[K]) => {
    setShared((prev) => ({ ...prev, [key]: value }));
  }, []);

  const fetchFromApi = useCallback(async (path: string, init?: RequestInit): Promise<Response> => {
    const cleanPath = path.replace(/^\/+/, "");
    let lastErr: unknown = null;

    // Prefer direct browser->API calls to avoid an extra frontend proxy hop; keep proxy as fallback.
    const directTargets = apiCandidates.map((base) => `${base}/api/${cleanPath}`);
    const targets = [...directTargets, `/api-proxy/${cleanPath}`];
    for (const target of targets) {
      try {
        const response = await fetch(target, init);
        // If proxy backend is unavailable, continue to the next candidate.
        if (target.startsWith("/api-proxy/") && [502, 503, 504].includes(response.status)) {
          lastErr = new Error(`Proxy unavailable (${response.status})`);
          continue;
        }
        return response;
      } catch (err) {
        lastErr = err;
      }
    }
    throw lastErr instanceof Error ? lastErr : new Error("Load failed");
  }, [apiCandidates]);

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
    if (typeof window === "undefined") return;
    const workspace = workspaceRef.current;
    if (!workspace) return;

    const updateWidth = () => {
      const next = workspace.clientWidth;
      if (next <= 0) return;
      setPlotWidthPx((prev) => (Math.abs(prev - next) > 1 ? next : prev));
    };

    updateWidth();
    const observer = new ResizeObserver(updateWidth);
    observer.observe(workspace);
    window.addEventListener("resize", updateWidth);

    return () => {
      observer.disconnect();
      window.removeEventListener("resize", updateWidth);
    };
  }, []);

  const displayFigure = useMemo(() => {
    if (resultMode !== mode) return null;
    const figure = result?.figure;
    if (!figure?.data) return null;

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
    return { ...figure, layout };
  }, [mode, plotWidthPx, result, resultMode, zoomRevision]);
  const deferredFigure = useDeferredValue(displayFigure);

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

  function renderSharedControls() {
    const isSky = mode === "skymap";
    const isLc = mode === "lightcurve";

    return (
      <>
        <div className="sb-stack">
          <label className="sb-field">
            <span className="sb-label">
              d<sub>L</sub> (Mpc)
            </span>
            <input
              type="number"
              min={0.1}
              max={1e6}
              step={10}
              value={shared.d_L_mpc}
              onChange={(e) => setSharedField("d_L_mpc", Number(e.target.value))}
            />
          </label>
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

        {!isSky ? (
          <div className="sb-row-2">
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
              <label className="sb-field">
                <span className="sb-label">Time</span>
              <select
                value={shared.time_unit}
                disabled={!isLc}
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
        ) : null}

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
                    <label key={mediumType} className="mode-radio">
                      <input
                        type="radio"
                        name="medium-type"
                        checked={shared.medium_type === mediumType}
                        onChange={() => setSharedField("medium_type", mediumType)}
                      />
                      <span>{mediumType}</span>
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
                label="log10(Duration (s))"
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
          <div className="sb-row-3">
            <SliderField
              label={
                <>
                  φ ppd
                </>
              }
              min={0.05}
              max={1}
              step={0.05}
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
              step={0.1}
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
              max={20}
              step={0.5}
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
    () => buildInteractiveSpec(fullComputationSpec, sliderInteracting),
    [fullComputationSpec, sliderInteracting],
  );

  function renderModeSpecific() {
    const skyPixelMax = skyAnimate ? SKY_MAX_PIXEL_ANIMATE : SKY_MAX_PIXEL_STATIC;
    const skyPixelOptions = SKY_PIXEL_OPTIONS.filter((v) => v <= skyPixelMax);

    if (mode === "lightcurve") {
      return (
        <>
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
          {renderSharedControls()}
        </>
      );
    }

    if (mode === "spectrum") {
      return (
        <>
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

          <div className="sb-stack">
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
            <label className="sb-checkbox-inline sb-check-alone">
              <input
                type="checkbox"
                checked={sedNuFNu}
                disabled={shared.flux_unit === "AB mag"}
                onChange={(e) => setSedNuFNu(e.target.checked)}
              />
              ν F<sub>ν</sub>
            </label>
            {shared.flux_unit === "AB mag" ? <p className="sb-muted">ν F<sub>ν</sub> is disabled for AB mag.</p> : null}
          </div>

          {renderSharedControls()}
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
        isNetworkFailure
          ? `Cannot reach backend. Checked /api-proxy and ${apiCandidates.join(" or ")}.`
          : message,
      );
    } finally {
      if (activeRequestRef.current === controller) {
        activeRequestRef.current = null;
      }
    }
  }, [apiCandidates, fetchFromApi, startFigureTransition]);

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
    return () => {
      if (runTimerRef.current !== null) {
        window.clearTimeout(runTimerRef.current);
        runTimerRef.current = null;
      }
      pendingSpecRef.current = null;
      activeRequestRef.current?.abort();
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
    let copied = false;
    try {
      await navigator.clipboard.writeText(BIBTEX_TEXT);
      copied = true;
    } catch {
      try {
        const ta = document.createElement("textarea");
        ta.value = BIBTEX_TEXT;
        ta.style.position = "fixed";
        ta.style.left = "-9999px";
        document.body.appendChild(ta);
        ta.select();
        copied = document.execCommand("copy");
        document.body.removeChild(ta);
      } catch {
        copied = false;
      }
    }

    setCiteLinkText(copied ? "copied" : "copy failed");
    window.setTimeout(() => setCiteLinkText("cite"), 1200);
  }, []);

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

        {renderModeSpecific()}

        <div className="sidebar-footer">
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
            If you find VegasAfterglow useful in your research, we would be grateful if you could{" "}
            <button className="sb-footer-link sb-footer-link-btn" type="button" onClick={() => void copyBibtex()}>
              {citeLinkText}
            </button>{" "}
            our work.
          </p>
          <div className="sidebar-downloads">
            {mode !== "skymap" ? (
              <button
                disabled={downloading !== null || !displayFigure?.data}
                onClick={() => void downloadCurrentExport("csv")}
              >
                {downloading === "csv" ? "Preparing CSV..." : "CSV"}
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
              {downloading === "json" ? "Preparing JSON..." : "JSON"}
            </button>
          </div>
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
          <PlotFigure figure={deferredFigure} mode={mode} onRelayout={handlePlotRelayout} />
        ) : (
          <div className="workspace-empty">
            <p className="muted">Adjust parameters in the sidebar to render plot.</p>
          </div>
        )}

      </section>
    </main>
  );
}
