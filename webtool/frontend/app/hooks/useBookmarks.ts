import { useCallback, useState } from "react";
import { MAX_BOOKMARKS, URL_STATE_PARAM } from "../lib/constants";
import type { BookmarkEntry, UiSnapshot } from "../lib/types";
import { cloneSnapshot, encodeUrlState, makeShareableSnapshot } from "../lib/utils/snapshot";

export function useBookmarks() {
  const [bookmarks, setBookmarks] = useState<BookmarkEntry[]>([]);
  const [bookmarkNameDraft, setBookmarkNameDraft] = useState("");
  const [setupStatusText, setSetupStatusText] = useState("");

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

  const removeBookmarkById = useCallback((bookmarkId: string) => {
    setBookmarks((prev) => prev.filter((bookmark) => bookmark.id !== bookmarkId));
  }, []);

  const buildShareLink = useCallback((snapshot: UiSnapshot): string | null => {
    if (typeof window === "undefined") return null;
    const encoded = encodeUrlState(makeShareableSnapshot(snapshot));
    if (!encoded) return null;
    const url = new URL(window.location.href);
    url.searchParams.set(URL_STATE_PARAM, encoded);
    return url.toString();
  }, []);

  return {
    bookmarks,
    setBookmarks,
    bookmarkNameDraft,
    setBookmarkNameDraft,
    setupStatusText,
    setSetupStatusText,
    saveSnapshotAsBookmark,
    removeBookmarkById,
    buildShareLink,
  };
}
