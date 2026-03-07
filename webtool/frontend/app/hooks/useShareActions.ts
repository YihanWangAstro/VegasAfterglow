import { useCallback, useEffect, useRef, type Dispatch, type SetStateAction } from "react";
import type { UiSnapshot } from "../lib/types";
import { clearTimeoutRef } from "../lib/utils/async";

type Args = {
  buildShareLink: (snapshot: UiSnapshot) => string | null;
  setSetupStatusText: Dispatch<SetStateAction<string>>;
  setCiteLinkText: Dispatch<SetStateAction<string>>;
  bibtexText: string;
};

export function useShareActions({ buildShareLink, setSetupStatusText, setCiteLinkText, bibtexText }: Args) {
  const setupStatusTimerRef = useRef<number | null>(null);
  const citeTimerRef = useRef<number | null>(null);

  const copyTextToClipboard = useCallback(async (text: string): Promise<boolean> => {
    if (!navigator.clipboard?.writeText) return false;
    try {
      await navigator.clipboard.writeText(text);
      return true;
    } catch {
      return false;
    }
  }, []);

  const showSetupStatus = useCallback(
    (text: string, timeoutMs = 1500) => {
      setSetupStatusText(text);
      clearTimeoutRef(setupStatusTimerRef);
      setupStatusTimerRef.current = window.setTimeout(() => {
        setupStatusTimerRef.current = null;
        setSetupStatusText("");
      }, timeoutMs);
    },
    [setSetupStatusText],
  );

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

  const copyBibtex = useCallback(async () => {
    const copied = await copyTextToClipboard(bibtexText);
    setCiteLinkText(copied ? "copied" : "copy failed");
    clearTimeoutRef(citeTimerRef);
    citeTimerRef.current = window.setTimeout(() => {
      citeTimerRef.current = null;
      setCiteLinkText("cite");
    }, 1200);
  }, [bibtexText, copyTextToClipboard, setCiteLinkText]);

  useEffect(() => {
    return () => {
      clearTimeoutRef(setupStatusTimerRef);
      clearTimeoutRef(citeTimerRef);
    };
  }, []);

  return { copyShareLinkForSnapshot, copyBibtex };
}
