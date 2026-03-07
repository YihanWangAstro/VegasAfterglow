import { useCallback, useEffect, useMemo, useState } from "react";
import { API_STATUS_REFRESH_DIRECT_MS } from "../lib/constants";
import type { ApiEndpoint, ApiHealthStatus } from "../lib/types";
import { stripTrailingSlash, stringFromHealthPayload } from "../lib/utils/api";

export function useApiStatus() {
  const apiCandidates = useMemo(() => {
    const configured = process.env.NEXT_PUBLIC_API_URL;
    const host = typeof window !== "undefined" ? window.location.hostname : "";
    const isLocal = host === "localhost" || host === "127.0.0.1";
    const fallback = isLocal
      ? host === "127.0.0.1"
        ? "http://127.0.0.1:8000"
        : "http://localhost:8000"
      : "https://api.vegasafterglow.com";
    const primary = stripTrailingSlash((configured ?? fallback).trim());
    const candidates = [primary];
    if (primary.includes("localhost")) {
      candidates.push(primary.replace("localhost", "127.0.0.1"));
    }
    if (primary.includes("127.0.0.1")) {
      candidates.push(primary.replace("127.0.0.1", "localhost"));
    }
    return Array.from(new Set(candidates));
  }, []);

  const apiStatusExtraCandidates = useMemo(() => {
    const extraRaw = process.env.NEXT_PUBLIC_API_STATUS_URLS ?? "";
    return extraRaw
      .split(",")
      .map((item) => stripTrailingSlash(item.trim()))
      .filter((item) => item.length > 0);
  }, []);

  const apiStatusCandidates = useMemo(() => {
    return Array.from(new Set([...apiStatusExtraCandidates, ...apiCandidates]));
  }, [apiCandidates, apiStatusExtraCandidates]);

  const apiStatusDisplaySet = useMemo(() => {
    return apiStatusExtraCandidates.length > 0 ? new Set(apiStatusExtraCandidates) : new Set(apiStatusCandidates);
  }, [apiStatusCandidates, apiStatusExtraCandidates]);

  const [activeApiKey, setActiveApiKey] = useState<string>("");
  const [apiStatusByKey, setApiStatusByKey] = useState<Record<string, ApiHealthStatus>>({});

  const apiEndpoints = useMemo<ApiEndpoint[]>(
    () =>
      apiStatusCandidates.map((base) => ({
        key: base,
        healthTarget: `${base}/api/health`,
      })),
    [apiStatusCandidates],
  );

  const displayApiEndpoints = useMemo(
    () => apiEndpoints.filter((endpoint) => apiStatusDisplaySet.has(endpoint.key)),
    [apiEndpoints, apiStatusDisplaySet],
  );

  const resolveApiKeyFromTarget = useCallback(
    (target: string): string | null => {
      for (const base of apiCandidates) {
        if (target.startsWith(`${base}/api/`)) return base;
      }
      return null;
    },
    [apiCandidates],
  );

  const noteApiTarget = useCallback(
    (target: string) => {
      const key = resolveApiKeyFromTarget(target);
      if (!key) return;
      setActiveApiKey((prev) => (prev === key ? prev : key));
    },
    [resolveApiKeyFromTarget],
  );

  const updateApiStatus = useCallback((key: string, next: ApiHealthStatus) => {
    setApiStatusByKey((prev) => {
      const current = prev[key];
      if (
        current &&
        current.up === next.up &&
        current.latencyMs === next.latencyMs &&
        current.statusCode === next.statusCode &&
        current.region === next.region &&
        current.locationLabel === next.locationLabel
      ) {
        return prev;
      }
      return { ...prev, [key]: next };
    });
  }, []);

  const probeEndpoint = useCallback(
    async (endpoint: ApiEndpoint, cancelledRef?: { current: boolean }): Promise<void> => {
      const started = performance.now();
      try {
        const response = await fetch(endpoint.healthTarget, { cache: "no-store" });
        if (cancelledRef?.current) return;
        if (!response.ok) {
          updateApiStatus(endpoint.key, {
            up: false,
            latencyMs: null,
            statusCode: response.status,
            region: null,
            locationLabel: null,
          });
          return;
        }
        let region: string | null = null;
        let locationLabel: string | null = null;
        try {
          const payload = (await response.clone().json()) as unknown;
          region = stringFromHealthPayload(payload, "region");
          locationLabel = stringFromHealthPayload(payload, "location_label");
        } catch {
          region = null;
          locationLabel = null;
        }
        const elapsed = Math.max(1, Math.round(performance.now() - started));
        updateApiStatus(endpoint.key, {
          up: true,
          latencyMs: elapsed,
          statusCode: response.status,
          region,
          locationLabel,
        });
      } catch {
        if (cancelledRef?.current) return;
        updateApiStatus(endpoint.key, {
          up: false,
          latencyMs: null,
          statusCode: null,
          region: null,
          locationLabel: null,
        });
      }
    },
    [updateApiStatus],
  );

  useEffect(() => {
    setApiStatusByKey((prev) => {
      const allowed = new Set(apiEndpoints.map((endpoint) => endpoint.key));
      let changed = false;
      const next: Record<string, ApiHealthStatus> = {};
      for (const [key, value] of Object.entries(prev)) {
        if (!allowed.has(key)) {
          changed = true;
          continue;
        }
        next[key] = value;
      }
      return changed ? next : prev;
    });
    if (activeApiKey && !apiEndpoints.some((endpoint) => endpoint.key === activeApiKey)) {
      setActiveApiKey("");
    }
  }, [activeApiKey, apiEndpoints]);

  useEffect(() => {
    const cancelledRef = { current: false };
    const endpointMap = new Map(apiEndpoints.map((endpoint) => [endpoint.key, endpoint] as const));

    const resolveProbeTargets = (): ApiEndpoint[] => {
      const activeEndpoint = activeApiKey ? endpointMap.get(activeApiKey) ?? null : null;
      const fallbackDisplay = displayApiEndpoints[0] ?? apiEndpoints[0] ?? null;
      const targets = [activeEndpoint, fallbackDisplay].filter(
        (endpoint): endpoint is ApiEndpoint => endpoint !== null,
      );
      return targets.filter((endpoint, idx, arr) => arr.findIndex((item) => item.key === endpoint.key) === idx);
    };

    const probeActiveTargets = async () => {
      for (const endpoint of resolveProbeTargets()) {
        await probeEndpoint(endpoint, cancelledRef);
      }
    };

    void probeActiveTargets();
    const directTimer = window.setInterval(() => {
      void probeActiveTargets();
    }, API_STATUS_REFRESH_DIRECT_MS);

    return () => {
      cancelledRef.current = true;
      window.clearInterval(directTimer);
    };
  }, [activeApiKey, apiEndpoints, displayApiEndpoints, probeEndpoint]);

  useEffect(() => {
    const cancelledRef = { current: false };
    if (displayApiEndpoints.length === 0) return;
    void Promise.all(displayApiEndpoints.map((endpoint) => probeEndpoint(endpoint, cancelledRef)));
    return () => {
      cancelledRef.current = true;
    };
  }, [displayApiEndpoints, probeEndpoint]);

  const activeDisplayApiKey = useMemo(() => {
    if (displayApiEndpoints.length === 0) return "";
    if (activeApiKey && apiStatusDisplaySet.has(activeApiKey)) return activeApiKey;
    return displayApiEndpoints[0]?.key ?? "";
  }, [activeApiKey, apiStatusDisplaySet, displayApiEndpoints]);

  const apiStatusRows = useMemo(
    () =>
      displayApiEndpoints.map((endpoint) => {
        const status = apiStatusByKey[endpoint.key];
        const statusClass = !status ? "pending" : status.up ? "up" : "down";
        const statusText = !status
          ? "checking..."
          : status.up
            ? `${status.latencyMs ?? "--"} ms`
            : status.statusCode !== null
              ? `down (${status.statusCode})`
              : "down";
        return {
          ...endpoint,
          isActive: endpoint.key === activeDisplayApiKey,
          regionText: status?.locationLabel ?? status?.region ?? "region unknown",
          statusClass,
          statusText,
        };
      }),
    [activeDisplayApiKey, apiStatusByKey, displayApiEndpoints],
  );

  const activeApiStatusRow = apiStatusRows.find((row) => row.isActive) ?? apiStatusRows[0] ?? null;
  const otherApiStatusRows = activeApiStatusRow
    ? apiStatusRows.filter((row) => row.key !== activeApiStatusRow.key)
    : apiStatusRows;

  const probeOtherServersOnce = useCallback(() => {
    if (otherApiStatusRows.length === 0) return;
    const cancelledRef = { current: false };
    void Promise.all(otherApiStatusRows.map((row) => probeEndpoint(row, cancelledRef)));
  }, [otherApiStatusRows, probeEndpoint]);

  const fetchFromApi = useCallback(
    async (path: string, init?: RequestInit): Promise<Response> => {
      const cleanPath = path.replace(/^\/+/, "");
      let lastErr: unknown = null;

      const directTargets = apiCandidates.map((base) => `${base}/api/${cleanPath}`);
      for (const target of directTargets) {
        try {
          const response = await fetch(target, init);
          noteApiTarget(target);
          return response;
        } catch (err) {
          lastErr = err;
        }
      }
      throw lastErr instanceof Error ? lastErr : new Error("Load failed");
    },
    [apiCandidates, noteApiTarget],
  );

  return {
    apiCandidates,
    fetchFromApi,
    activeApiStatusRow,
    otherApiStatusRows,
    probeOtherServersOnce,
  };
}
