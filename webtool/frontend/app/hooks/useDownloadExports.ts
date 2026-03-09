import { useCallback } from "react";
import { DOWNLOAD_TEXT_META } from "../lib/constants";
import type { DownloadKind, LcPlotData, Mode, RunResponse, SedPlotData, SkymapPlotData } from "../lib/types";
import {
  COMP_ORDER,
  FLUX_SCALES_CGS,
  FREQ_DISP_SCALES,
  TIME_SCALES_S,
  formatBandLabel,
  formatFreqLabel,
  formatTimeLabel,
} from "../lib/utils/plot-builders";
import { downloadText } from "../lib/utils/storage";

export type PresentationParams = {
  fluxUnit: string;
  timeUnit: string;
  sedFreqUnit: string;
};

// ---------------------------------------------------------------------------
// Local export builders — generate CSV/JSON from plot_data without backend.
// ---------------------------------------------------------------------------

function sci(v: number): string {
  return v.toExponential(6);
}

function exportLcCsv(pd: LcPlotData, fluxUnit: string, timeUnit: string): string {
  const tScale = TIME_SCALES_S[timeUnit] ?? 86400;
  const fScale = FLUX_SCALES_CGS[fluxUnit] ?? 1e-26;

  const cols = [`t(${timeUnit})`];
  const ptComps = pd.pt ? COMP_ORDER.filter((n) => n in pd.pt!.components) : [];
  if (pd.pt) {
    for (const nu of pd.pt.freq_hz) {
      const label = formatFreqLabel(nu);
      for (const comp of ptComps) cols.push(`F_${comp}(${label})[${fluxUnit}]`);
    }
  }
  const bandCompsList = pd.bands.map((band) => COMP_ORDER.filter((n) => n in band.components));
  for (let b = 0; b < pd.bands.length; b++) {
    const bl = formatBandLabel(pd.bands[b].nu_min, pd.bands[b].nu_max, pd.bands[b].name);
    for (const comp of bandCompsList[b]) cols.push(`F_${comp}(${bl})[erg/cm2/s]`);
  }

  const rows = [cols.join(",")];
  for (let j = 0; j < pd.times_s.length; j++) {
    const row = [sci(pd.times_s[j] / tScale)];
    if (pd.pt) {
      for (let i = 0; i < pd.pt.freq_hz.length; i++) {
        for (const comp of ptComps) row.push(sci(pd.pt.components[comp][i][j] / fScale));
      }
    }
    for (let b = 0; b < pd.bands.length; b++) {
      for (const comp of bandCompsList[b]) row.push(sci(pd.bands[b].components[comp][j]));
    }
    rows.push(row.join(","));
  }
  return rows.join("\n") + "\n";
}

function exportLcJson(pd: LcPlotData, fluxUnit: string, timeUnit: string): string {
  const tScale = TIME_SCALES_S[timeUnit] ?? 86400;
  const fScale = FLUX_SCALES_CGS[fluxUnit] ?? 1e-26;

  const obj: Record<string, unknown> = {
    units: { time: timeUnit, flux_density: fluxUnit },
    times: pd.times_s.map((t) => t / tScale),
    frequencies_Hz: pd.pt?.freq_hz ?? [],
    flux_density: {} as Record<string, Record<string, number[]>>,
  };

  if (pd.pt) {
    const fd = obj.flux_density as Record<string, Record<string, number[]>>;
    for (let i = 0; i < pd.pt.freq_hz.length; i++) {
      const label = formatFreqLabel(pd.pt.freq_hz[i]);
      fd[label] = {};
      for (const [comp, byFreq] of Object.entries(pd.pt.components)) {
        fd[label][comp] = byFreq[i].map((f) => f / fScale);
      }
    }
  }

  if (pd.bands.length > 0) {
    (obj.units as Record<string, string>).band_flux = "erg/cm2/s";
    const bands: Record<string, unknown> = {};
    for (const band of pd.bands) {
      const bl = formatBandLabel(band.nu_min, band.nu_max, band.name);
      bands[bl] = {
        nu_min_Hz: band.nu_min,
        nu_max_Hz: band.nu_max,
        flux: { ...band.components },
      };
    }
    obj.bands = bands;
  }

  return JSON.stringify(obj, null, 2);
}

function exportSedCsv(pd: SedPlotData, fluxUnit: string, freqUnit: string): string {
  const fScale = FLUX_SCALES_CGS[fluxUnit] ?? 1e-26;
  const nuScale = FREQ_DISP_SCALES[freqUnit] ?? 1;

  const comps = COMP_ORDER.filter((n) => n in pd.components);
  const cols = [`nu(${freqUnit})`];
  for (const t of pd.t_snapshots_s) {
    const tl = formatTimeLabel(t);
    for (const comp of comps) cols.push(`F_${comp}(t=${tl})[${fluxUnit}]`);
  }

  const rows = [cols.join(",")];
  for (let i = 0; i < pd.freq_hz.length; i++) {
    const row = [sci(pd.freq_hz[i] / nuScale)];
    for (let j = 0; j < pd.t_snapshots_s.length; j++) {
      for (const comp of comps) row.push(sci(pd.components[comp][j][i] / fScale));
    }
    rows.push(row.join(","));
  }
  return rows.join("\n") + "\n";
}

function exportSedJson(pd: SedPlotData, fluxUnit: string, freqUnit: string): string {
  const fScale = FLUX_SCALES_CGS[fluxUnit] ?? 1e-26;
  const nuScale = FREQ_DISP_SCALES[freqUnit] ?? 1;

  const obj: Record<string, unknown> = {
    units: { frequency: freqUnit, flux_density: fluxUnit },
    frequencies: pd.freq_hz.map((nu) => nu / nuScale),
    frequencies_Hz: pd.freq_hz,
    t_snapshots_s: pd.t_snapshots_s,
    flux_density: {} as Record<string, Record<string, number[]>>,
  };

  const fd = obj.flux_density as Record<string, Record<string, number[]>>;
  for (let j = 0; j < pd.t_snapshots_s.length; j++) {
    const tl = formatTimeLabel(pd.t_snapshots_s[j]);
    fd[tl] = {};
    for (const [comp, byT] of Object.entries(pd.components)) {
      fd[tl][comp] = byT[j].map((f) => f / fScale);
    }
  }

  return JSON.stringify(obj, null, 2);
}

function exportSkymapJson(pd: SkymapPlotData): string {
  const obj = {
    t_obs_s: pd.t_obs_s,
    nu_obs_Hz: pd.nu_obs_hz,
    fov_uas: pd.extent_uas[1] - pd.extent_uas[0],
    extent_uas: pd.extent_uas,
    units: { image: "log10(erg/cm2/s/Hz/sr)", extent: "uas" },
    nx: pd.nx,
    ny: pd.ny,
    frames_b64f32: pd.frames_b64f32,
  };
  return JSON.stringify(obj, null, 2);
}

// ---------------------------------------------------------------------------
// Hook
// ---------------------------------------------------------------------------

type Args = {
  mode: Mode;
  result: RunResponse | null;
  presentationParams: PresentationParams;
};

export function useDownloadExports({ mode, result, presentationParams }: Args) {
  const downloadCurrentExport = useCallback(
    (kind: DownloadKind) => {
      const plotData = result?.plot_data;
      if (!plotData) return;

      const { fluxUnit, timeUnit, sedFreqUnit } = presentationParams;
      let content: string;

      if (mode === "lightcurve") {
        const pd = plotData as LcPlotData;
        content = kind === "csv" ? exportLcCsv(pd, fluxUnit, timeUnit) : exportLcJson(pd, fluxUnit, timeUnit);
      } else if (mode === "spectrum") {
        const pd = plotData as SedPlotData;
        content = kind === "csv" ? exportSedCsv(pd, fluxUnit, sedFreqUnit) : exportSedJson(pd, fluxUnit, sedFreqUnit);
      } else {
        // skymap — only JSON
        if (kind === "csv") return;
        content = exportSkymapJson(plotData as SkymapPlotData);
      }

      const meta = DOWNLOAD_TEXT_META[kind];
      if (meta) downloadText(content, `afterglow_${mode}.${meta.ext}`, meta.mime);
    },
    [mode, presentationParams, result],
  );

  return { downloadCurrentExport };
}
