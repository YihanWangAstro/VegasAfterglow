import type { ReactNode } from "react";
import { C_KM_S, H0_KM_S_MPC } from "../constants";

export function safeLog10(value: number, fallback: number): number {
  return Number.isFinite(value) && value > 0 ? Math.log10(value) : fallback;
}

export function formatParamValueNode(value: number): ReactNode {
  if (!Number.isFinite(value)) return String(value);
  if (value === 0) return "0";
  const abs = Math.abs(value);
  const exponent = Math.floor(Math.log10(abs));
  const mantissa = value / Math.pow(10, exponent);
  const mantissaText = mantissa.toFixed(2).replace(/\.?0+$/, "");
  if (exponent === 0) return mantissaText;
  if (mantissaText === "1") {
    return (
      <>
        10<sup>{exponent}</sup>
      </>
    );
  }
  if (mantissaText === "-1") {
    return (
      <>
        -10<sup>{exponent}</sup>
      </>
    );
  }
  return (
    <>
      {mantissaText} x 10<sup>{exponent}</sup>
    </>
  );
}

export function clampRangeValue(value: number, min: number, max: number): number {
  return Math.max(min, Math.min(max, value));
}

export function luminosityDistanceMpcFromRedshift(z: number): number {
  const zSafe = Math.max(0, z);
  return (C_KM_S / H0_KM_S_MPC) * zSafe * (1 + zSafe);
}

export function redshiftFromLuminosityDistanceMpc(dLmpc: number): number {
  const dSafe = Math.max(0, dLmpc);
  const x = (dSafe * H0_KM_S_MPC) / C_KM_S;
  return (-1 + Math.sqrt(1 + 4 * x)) / 2;
}

export function normalizeDistanceMpc(value: number): number {
  if (!Number.isFinite(value)) return 0;
  return Math.max(0, Math.trunc(value));
}

export function roundToSignificant(value: number, digits: number): number {
  if (!Number.isFinite(value) || value === 0) return value;
  const abs = Math.abs(value);
  const exponent = Math.floor(Math.log10(abs));
  const scale = Math.pow(10, exponent - digits + 1);
  return Math.round(value / scale) * scale;
}

export function normalizeRedshift(value: number): number {
  if (!Number.isFinite(value)) return 0;
  return Math.max(0, roundToSignificant(value, 3));
}

export function formatZPlain(value: number): string {
  if (!Number.isFinite(value)) return "";
  if (value === 0) return "0";
  const abs = Math.abs(value);
  const exponent = Math.floor(Math.log10(abs));
  const decimals = Math.max(0, Math.min(12, 2 - exponent));
  return value.toFixed(decimals).replace(/\.?0+$/, "");
}

export function nearlyEqualRelative(a: number, b: number, relTol = 1e-6): boolean {
  return Math.abs(a - b) <= relTol * Math.max(1, Math.abs(a), Math.abs(b));
}

export function clampInt(value: unknown, min: number, max: number): number {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return min;
  return Math.max(min, Math.min(max, Math.round(numeric)));
}
