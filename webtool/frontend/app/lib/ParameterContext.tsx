"use client";

import { createContext, useContext } from "react";
import type { ParameterActions } from "../hooks/useParameterState";
import type { ParameterState } from "../hooks/useParameterState";
import type { SharedParams } from "./types";

export type ParameterContextValue = {
  state: ParameterState;
  actions: ParameterActions;
  setSharedField: <K extends keyof SharedParams>(key: K, value: SharedParams[K]) => void;
  setDistanceMpc: (value: number) => void;
  setDistanceRedshift: (value: number) => void;
};

export const ParameterContext = createContext<ParameterContextValue | null>(null);

export function useParameters(): ParameterContextValue {
  const ctx = useContext(ParameterContext);
  if (!ctx) throw new Error("useParameters must be used within a ParameterContext.Provider");
  return ctx;
}
