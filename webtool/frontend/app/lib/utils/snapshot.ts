import {
  FALLBACK_SHARED,
  INSTRUMENT_GROUP_ORDER,
  INSTRUMENT_TYPE_BY_NAME,
  INTERACTIVE_LIGHTCURVE_NUM_T_MAX,
  INTERACTIVE_SKY_FRAMES_MAX,
  INTERACTIVE_SKY_PIXEL_MAX_ANIMATE,
  INTERACTIVE_SKY_PIXEL_MAX_STATIC,
  INTERACTIVE_SPECTRUM_NUM_NU_MAX,
  SKY_DEFAULT_N_FRAMES,
  SKY_MAX_N_FRAMES,
  SKY_MAX_PIXEL_STATIC,
  URL_STATE_MAX_CHARS,
  URL_STATE_VERSION,
} from "../constants";
import type {
  BookmarkEntry,
  ComputationSpec,
  DistanceDriver,
  InstrumentGroup,
  Mode,
  ObservationGroup,
  OptionsResponse,
  SharedParams,
  UiSnapshot,
} from "../types";
import {
  clampInt,
  clampRangeValue,
  normalizeDistanceMpc,
  normalizeRedshift,
} from "./math";
import { defaultObsGroup } from "./obs";

export function cloneSnapshot(snapshot: UiSnapshot): UiSnapshot {
  return JSON.parse(JSON.stringify(snapshot)) as UiSnapshot;
}

function asFiniteNumber(value: unknown, fallback: number): number {
  const next = Number(value);
  return Number.isFinite(next) ? next : fallback;
}

function asString(value: unknown, fallback: string): string {
  return typeof value === "string" ? value : fallback;
}

function asBoolean(value: unknown, fallback: boolean): boolean {
  return typeof value === "boolean" ? value : fallback;
}

export function asStringArray(value: unknown): string[] {
  if (!Array.isArray(value)) return [];
  return value.filter((item): item is string => typeof item === "string");
}

export function setIfString(value: unknown, apply: (value: string) => void): void {
  if (typeof value === "string") {
    apply(value);
  }
}

export function setIfNumber(value: unknown, apply: (value: number) => void): void {
  if (typeof value === "number") {
    apply(value);
  }
}

export function setIfBoolean(value: unknown, apply: (value: boolean) => void): void {
  if (typeof value === "boolean") {
    apply(value);
  }
}

export function setIfArray<T>(value: unknown, apply: (value: T[]) => void): void {
  if (Array.isArray(value)) {
    apply(value as T[]);
  }
}

export function normalizeShared(shared: SharedParams, mode: Mode): SharedParams {
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

export function groupedInstrumentsFromOptions(options: OptionsResponse): InstrumentGroup[] {
  const groupedFromApi =
    Array.isArray(options.instrument_groups) && options.instrument_groups.length > 0
      ? options.instrument_groups
          .map((group) => ({
            label: typeof group.label === "string" ? group.label : "",
            items: asStringArray(group.items),
          }))
          .filter((group) => group.label.length > 0 && group.items.length > 0)
      : [];

  const fallbackInstruments = asStringArray(options.instruments);
  if (groupedFromApi.length === 0) {
    return fallbackInstruments.length > 0 ? groupInstrumentsByType(fallbackInstruments) : [];
  }

  const groupedKnown = new Set(groupedFromApi.flatMap((group) => group.items));
  const ungrouped = fallbackInstruments.filter((name) => !groupedKnown.has(name));
  if (ungrouped.length === 0) return groupedFromApi;
  return [...groupedFromApi, { label: "Other", items: ungrouped }];
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
        ...(typeof group.freq === "string" ? { freq: group.freq } : {}),
      } satisfies ObservationGroup;
    })
    .filter((group): group is ObservationGroup => group !== null);
}

export function defaultUiSnapshot(): UiSnapshot {
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

export function normalizeUiSnapshot(value: unknown): UiSnapshot | null {
  if (!value || typeof value !== "object") return null;
  const defaults = defaultUiSnapshot();
  const record = value as Record<string, unknown>;

  const modeCandidate = record.mode;
  const mode: Mode =
    modeCandidate === "lightcurve" || modeCandidate === "spectrum" || modeCandidate === "skymap"
      ? modeCandidate
      : defaults.mode;
  const distanceLinked = asBoolean(record.distance_linked, defaults.distance_linked);
  const distanceDriver: DistanceDriver =
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
      num_nu: clampRangeValue(asFiniteNumber(sedRaw.num_nu, defaults.spectrum.num_nu), 50, 2000),
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
      n_frames: clampRangeValue(asFiniteNumber(skyRaw.n_frames, defaults.skymap.n_frames), 3, SKY_MAX_N_FRAMES),
      nu_input: asString(skyRaw.nu_input, defaults.skymap.nu_input),
      fov: asFiniteNumber(skyRaw.fov, defaults.skymap.fov),
      npixel: clampRangeValue(asFiniteNumber(skyRaw.npixel, defaults.skymap.npixel), 64, SKY_MAX_PIXEL_STATIC),
    },
  };
}

export function makeShareableSnapshot(snapshot: UiSnapshot): UiSnapshot {
  const shared = cloneSnapshot(snapshot);
  shared.lightcurve.observation_groups = [];
  shared.spectrum.observation_groups = [];
  return shared;
}

export function encodeUrlState(snapshot: UiSnapshot): string | null {
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

export function decodeUrlState(raw: string): UiSnapshot | null {
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

export function parseStoredBookmarks(raw: string | null): BookmarkEntry[] | null {
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

export function buildComputationSpecFromSnapshot(snapshot: UiSnapshot, mode: Mode): ComputationSpec {
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

export function groupInstrumentsByType(instrumentNames: string[]): InstrumentGroup[] {
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

export function buildInteractiveSpec(base: ComputationSpec, interactive: boolean): ComputationSpec {
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

export function isEmptyPrimaryInputSpec(spec: ComputationSpec): boolean {
  if (spec.endpoint === "lightcurve") {
    const raw = spec.payload.frequencies_input;
    return typeof raw !== "string" || raw.trim().length === 0;
  }

  if (spec.endpoint === "spectrum") {
    const raw = spec.payload.t_snapshots_input;
    return typeof raw !== "string" || raw.trim().length === 0;
  }

  return false;
}

export function clampTabIndex(current: number, length: number): number {
  if (length === 0) return 0;
  return Math.min(current, length - 1);
}
