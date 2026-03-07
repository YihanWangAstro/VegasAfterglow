export function stripTrailingSlash(url: string): string {
  return url.replace(/\/+$/, "");
}

export function stringFromHealthPayload(payload: unknown, key: string): string | null {
  if (!payload || typeof payload !== "object") return null;
  const value = (payload as Record<string, unknown>)[key];
  return typeof value === "string" ? value.trim() || null : null;
}
