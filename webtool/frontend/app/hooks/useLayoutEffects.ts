import { useEffect, type Dispatch, type RefObject, type SetStateAction } from "react";

type PlotWidthObserverArgs = {
  workspaceRef: RefObject<HTMLElement | null>;
  setPlotWidthPx: Dispatch<SetStateAction<number>>;
};

export function usePlotWidthObserver({ workspaceRef, setPlotWidthPx }: PlotWidthObserverArgs) {
  useEffect(() => {
    if (typeof window === "undefined") return;
    const workspace = workspaceRef.current;
    if (!workspace) return;
    let observedPlotWrap: HTMLElement | null = null;

    const resolveMeasuredWidth = () => {
      const plotWrap = workspace.querySelector(".plot-wrap") as HTMLElement | null;
      return plotWrap?.clientWidth ?? workspace.clientWidth;
    };

    const updateWidth = () => {
      const next = resolveMeasuredWidth();
      if (next <= 0) return;
      setPlotWidthPx((prev) => (Math.abs(prev - next) > 1 ? next : prev));
    };

    const attachPlotObserver = (observer: ResizeObserver) => {
      const nextPlotWrap = workspace.querySelector(".plot-wrap") as HTMLElement | null;
      if (nextPlotWrap === observedPlotWrap) return;
      if (observedPlotWrap) {
        observer.unobserve(observedPlotWrap);
      }
      observedPlotWrap = nextPlotWrap;
      if (observedPlotWrap) {
        observer.observe(observedPlotWrap);
      }
    };

    updateWidth();
    const resizeObserver = new ResizeObserver(() => {
      attachPlotObserver(resizeObserver);
      updateWidth();
    });
    resizeObserver.observe(workspace);
    attachPlotObserver(resizeObserver);

    const mutationObserver = new MutationObserver(() => {
      attachPlotObserver(resizeObserver);
      updateWidth();
    });
    mutationObserver.observe(workspace, { childList: true, subtree: true });

    window.addEventListener("resize", updateWidth);
    window.addEventListener("orientationchange", updateWidth);

    return () => {
      if (observedPlotWrap) {
        resizeObserver.unobserve(observedPlotWrap);
      }
      resizeObserver.disconnect();
      mutationObserver.disconnect();
      window.removeEventListener("resize", updateWidth);
      window.removeEventListener("orientationchange", updateWidth);
    };
  }, [setPlotWidthPx, workspaceRef]);
}

export function useViewportLayout(sidebarOpen: boolean) {
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
}
