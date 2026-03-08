import type { ObsEntry, ObservationGroup } from "../types";
import { FLUX_SCALES_CGS, TIME_SCALES_S } from "./plot-builders";

export function defaultObsGroup(isLc: boolean, idx: number): ObservationGroup {
  return {
    legend: `data ${idx}`,
    x_unit: isLc ? "day" : "Hz",
    y_unit: "mJy",
    text: "",
    visible: true,
  };
}

export function parseStoredObsGroups(raw: string | null, isLc: boolean): ObservationGroup[] | null {
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
          ...(typeof group.freq === "string" && group.freq.trim() ? { freq: group.freq.trim() } : {}),
          ...(typeof group.shift === "number" && group.shift !== 1 ? { shift: group.shift } : {}),
        };
      })
      .filter((group): group is ObservationGroup => group !== null);

    return normalized;
  } catch {
    return null;
  }
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

const FREQ_SCALES: Record<string, number> = { Hz: 1, GHz: 1e9, keV: 2.417989242e17, MeV: 2.417989242e20 };

function magToCgs(mag: number, magErr: number): [number, number] {
  const fNu = 10 ** (-(mag + 48.6) / 2.5);
  const fErr = magErr ? fNu * (Math.LN10 / 2.5) * Math.abs(magErr) : 0;
  return [fNu, fErr];
}

/**
 * Parse observation groups into ObsEntry[] ready for the plot builder.
 * Replicates the backend's parse_observation_rows + _build_obs_raw in one pass.
 */
export function parseObsToEntries(groups: ObservationGroup[], isLc: boolean): ObsEntry[] {
  const xScales = isLc ? TIME_SCALES_S : FREQ_SCALES;
  const xDefault = isLc ? "day" : "Hz";
  const fnuGroups = new Map<string, [number, number, number][]>();
  const fbandGroups = new Map<string, [number, number, number][]>();

  for (const group of groups) {
    if (!group.visible || !group.text.trim()) continue;
    const label = group.legend.trim() || "data";
    const xUnit = group.x_unit || xDefault;
    const yUnit = group.y_unit;
    const xScale = xScales[xUnit] ?? 1;

    for (const line of group.text.trim().split(/\r?\n/)) {
      const parts = line.trim().split(/[,\s\t]+/);
      if (parts.length < 2) continue;
      const xVal = Number(parts[0]);
      const yVal = Number(parts[1]);
      if (!Number.isFinite(xVal) || !Number.isFinite(yVal) || xVal <= 0) continue;
      let errVal = parts.length > 2 ? Number(parts[2]) : 0;
      if (!Number.isFinite(errVal)) errVal = 0;

      const xPhys = xVal * xScale;
      if (yUnit === "erg/cm\u00b2/s") {
        const arr = fbandGroups.get(label) ?? [];
        arr.push([xPhys, yVal, Math.abs(errVal)]);
        fbandGroups.set(label, arr);
      } else if (yUnit === "AB mag") {
        const [fCgs, eCgs] = magToCgs(yVal, errVal);
        const arr = fnuGroups.get(label) ?? [];
        arr.push([xPhys, fCgs, eCgs]);
        fnuGroups.set(label, arr);
      } else {
        const scale = FLUX_SCALES_CGS[yUnit] ?? 1;
        const arr = fnuGroups.get(label) ?? [];
        arr.push([xPhys, yVal * scale, Math.abs(errVal) * scale]);
        fnuGroups.set(label, arr);
      }
    }
  }

  const labelSet = new Set<string>();
  fnuGroups.forEach((_, k) => labelSet.add(k));
  fbandGroups.forEach((_, k) => labelSet.add(k));
  const allLabels = Array.from(labelSet);
  return allLabels.map((label) => ({
    label,
    fnu: fnuGroups.get(label) ?? [],
    fband: fbandGroups.get(label) ?? [],
  }));
}

export async function parseObservationUpload(file: File): Promise<string> {
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
