import { useCallback, useEffect } from "react";
import type { DistanceDriver, SharedParams } from "../lib/types";
import {
  luminosityDistanceMpcFromRedshift,
  nearlyEqualRelative,
  normalizeDistanceMpc,
  normalizeRedshift,
  redshiftFromLuminosityDistanceMpc,
} from "../lib/utils/math";

type Args = {
  distanceLinked: boolean;
  distanceDriver: DistanceDriver;
  shared: SharedParams;
  setShared: (value: SharedParams | ((prev: SharedParams) => SharedParams)) => void;
};

export function useDistanceLinking({ distanceLinked, distanceDriver, shared, setShared }: Args) {
  const setSharedField = useCallback(<K extends keyof SharedParams>(key: K, value: SharedParams[K]) => {
    setShared((prev) => ({ ...prev, [key]: value }));
  }, [setShared]);

  const setDistanceMpc = useCallback(
    (nextMpcRaw: number) => {
      if (!Number.isFinite(nextMpcRaw)) return;
      const nextMpc = normalizeDistanceMpc(nextMpcRaw);
      setShared((prev) => {
        const next: SharedParams = { ...prev, d_L_mpc: nextMpc };
        if (distanceLinked) {
          next.z = normalizeRedshift(redshiftFromLuminosityDistanceMpc(nextMpc));
        }
        return next;
      });
    },
    [distanceLinked, setShared],
  );

  const setDistanceRedshift = useCallback(
    (nextZRaw: number) => {
      if (!Number.isFinite(nextZRaw)) return;
      const nextZ = Math.max(0, nextZRaw);
      setShared((prev) => {
        const next: SharedParams = { ...prev, z: nextZ };
        if (distanceLinked) {
          next.d_L_mpc = normalizeDistanceMpc(luminosityDistanceMpcFromRedshift(nextZ));
        }
        return next;
      });
    },
    [distanceLinked, setShared],
  );

  useEffect(() => {
    if (!distanceLinked) return;
    setShared((prev) => {
      if (distanceDriver === "dL") {
        const nextMpc = normalizeDistanceMpc(prev.d_L_mpc);
        const nextZ = normalizeRedshift(redshiftFromLuminosityDistanceMpc(nextMpc));
        if (nextMpc !== prev.d_L_mpc) {
          return { ...prev, d_L_mpc: nextMpc, z: nextZ };
        }
        if (nearlyEqualRelative(prev.z, nextZ)) return prev;
        return { ...prev, z: nextZ };
      }
      const nextDL = normalizeDistanceMpc(luminosityDistanceMpcFromRedshift(prev.z));
      if (prev.d_L_mpc === nextDL) return prev;
      return { ...prev, d_L_mpc: nextDL };
    });
  }, [distanceDriver, distanceLinked, setShared, shared.d_L_mpc, shared.z]);

  return { setSharedField, setDistanceMpc, setDistanceRedshift };
}
