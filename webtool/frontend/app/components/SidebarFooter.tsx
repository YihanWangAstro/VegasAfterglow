import type { DownloadKind, Mode } from "../lib/types";

type Props = {
  mode: Mode;
  downloading: DownloadKind | null;
  hasFigureData: boolean;
  appVersion: string;
  citeLinkText: string;
  onDownload: (kind: DownloadKind) => void;
  onCopyBibtex: () => void;
};

export function SidebarFooter({
  mode,
  downloading,
  hasFigureData,
  appVersion,
  citeLinkText,
  onDownload,
  onCopyBibtex,
}: Props) {
  return (
    <div className="sidebar-footer">
      <div className="sidebar-downloads">
        {mode !== "skymap" ? (
          <button disabled={downloading !== null || !hasFigureData} onClick={() => onDownload("csv")}>
            {downloading === "csv" ? "Preparing Data CSV..." : "Download Data (CSV)"}
          </button>
        ) : (
          <button disabled={downloading !== null || !hasFigureData} onClick={() => onDownload("gif")}>
            {downloading === "gif" ? "Preparing GIF..." : "Save GIF"}
          </button>
        )}
        <button disabled={downloading !== null || !hasFigureData} onClick={() => onDownload("json")}>
          {downloading === "json" ? "Preparing Data JSON..." : "Download Data (JSON)"}
        </button>
      </div>
      <p className="sb-footer-text">
        VegasAfterglow v{appVersion}
        <br />
        This interactive tool provides a subset of VegasAfterglow features. For MCMC fitting, custom jet/medium
        models, and the full API, see the{" "}
        <a className="sb-footer-link" href="https://github.com/YihanWangAstro/VegasAfterglow" target="_blank" rel="noreferrer">
          Python package
        </a>
        .
      </p>
      <p className="sb-footer-text">
        VegasAfterglow Webtool is currently supported by personal funding from the developers. If this tool is helpful,
        any support is appreciated, including sharing it with others or{" "}
        <button className="sb-footer-link sb-footer-link-btn" type="button" onClick={onCopyBibtex}>
          {citeLinkText}
        </button>{" "}
        our work in research.
      </p>
    </div>
  );
}
