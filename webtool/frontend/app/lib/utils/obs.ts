import type { ObservationGroup } from "../types";

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
        };
      })
      .filter((group): group is ObservationGroup => group !== null);

    return normalized;
  } catch {
    return null;
  }
}

export function compactObservationGroups(groups: ObservationGroup[]): ObservationGroup[] {
  return groups
    .filter((group) => group.visible && group.text.trim().length > 0)
    .map((group) => ({
      legend: group.legend.trim() || "data",
      x_unit: group.x_unit,
      y_unit: group.y_unit,
      text: group.text.trim(),
      visible: true,
      ...(group.freq?.trim() ? { freq: group.freq.trim() } : {}),
    }));
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
