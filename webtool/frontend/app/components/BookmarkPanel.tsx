import { HelpHint } from "./HelpHint";
import type { BookmarkEntry } from "../lib/types";

type Props = {
  mode: "lightcurve" | "spectrum" | "skymap";
  bookmarkNameDraft: string;
  setBookmarkNameDraft: (value: string) => void;
  saveCurrentBookmark: () => void;
  setupStatusText: string;
  bookmarks: BookmarkEntry[];
  loadBookmarkById: (id: string) => void;
  copyShareLinkForSnapshot: (snapshot: BookmarkEntry["snapshot"]) => Promise<void>;
  removeBookmarkById: (id: string) => void;
  compareEnabled: boolean;
  setCompareEnabled: (value: boolean) => void;
  compareBookmarkId: string;
  setCompareBookmarkId: (value: string) => void;
  compareRunning: boolean;
  compareError: string;
};

export function BookmarkPanel({
  mode,
  bookmarkNameDraft,
  setBookmarkNameDraft,
  saveCurrentBookmark,
  setupStatusText,
  bookmarks,
  loadBookmarkById,
  copyShareLinkForSnapshot,
  removeBookmarkById,
  compareEnabled,
  setCompareEnabled,
  compareBookmarkId,
  setCompareBookmarkId,
  compareRunning,
  compareError,
}: Props) {
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
