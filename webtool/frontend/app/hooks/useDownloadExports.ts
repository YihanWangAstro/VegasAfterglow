import { useCallback, useState, type Dispatch, type SetStateAction } from "react";
import { DOWNLOAD_TEXT_META } from "../lib/constants";
import type { ComputationSpec, DownloadKind, Mode, RunResponse } from "../lib/types";
import { downloadBase64, downloadText } from "../lib/utils/storage";

type Args = {
  mode: Mode;
  fullComputationSpec: ComputationSpec;
  postRunRequest: (
    spec: ComputationSpec,
    extraPayload: Record<string, unknown>,
    signal?: AbortSignal,
  ) => Promise<RunResponse>;
  exportsCache?: Record<string, string>;
  setResult: Dispatch<SetStateAction<RunResponse | null>>;
  setError: Dispatch<SetStateAction<string>>;
};

export function useDownloadExports({ mode, fullComputationSpec, postRunRequest, exportsCache, setResult, setError }: Args) {
  const [downloading, setDownloading] = useState<DownloadKind | null>(null);

  const downloadCurrentExport = useCallback(
    async (kind: DownloadKind) => {
      setError("");
      setDownloading(kind);
      try {
        let content = exportsCache?.[kind];
        if (!content) {
          const data = await postRunRequest(fullComputationSpec, {
            include_figure: false,
            include_exports: true,
            export_kinds: [kind],
          });
          content = data.exports?.[kind];
          if (data.exports && kind !== "gif") {
            setResult((prev) => (prev ? { ...prev, exports: { ...(prev.exports ?? {}), ...data.exports } } : prev));
          }
        }

        if (!content) {
          throw new Error(`No ${kind.toUpperCase()} export returned by backend.`);
        }
        if (kind === "gif") {
          downloadBase64(content, "afterglow_skymap.gif", "image/gif");
        } else {
          const meta = DOWNLOAD_TEXT_META[kind];
          downloadText(content, `afterglow_${mode}.${meta.ext}`, meta.mime);
        }
      } catch (err) {
        const message = err instanceof Error ? err.message : "Download failed";
        setError(message);
      } finally {
        setDownloading(null);
      }
    },
    [exportsCache, fullComputationSpec, mode, postRunRequest, setError, setResult],
  );

  return { downloading, downloadCurrentExport };
}
