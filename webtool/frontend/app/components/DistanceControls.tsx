import { SliderField } from "./SliderField";
import { formatZPlain, luminosityDistanceMpcFromRedshift, normalizeDistanceMpc, redshiftFromLuminosityDistanceMpc } from "../lib/utils/math";
import type { DistanceDriver, SharedParams } from "../lib/types";

type Props = {
  distanceDriver: DistanceDriver;
  distanceLinked: boolean;
  shared: SharedParams;
  setDistanceLinked: (next: boolean) => void;
  setDistanceDriver: (next: DistanceDriver) => void;
  setDistanceMpc: (next: number) => void;
  setDistanceRedshift: (next: number) => void;
  setSharedField: <K extends keyof SharedParams>(key: K, value: SharedParams[K]) => void;
};

export function DistanceControls({
  distanceDriver,
  distanceLinked,
  shared,
  setDistanceLinked,
  setDistanceDriver,
  setDistanceMpc,
  setDistanceRedshift,
  setSharedField,
}: Props) {
  const convertedText =
    distanceDriver === "dL" ? (
      <>
        z<span className="sb-distance-converted-eq"> = </span>
        {formatZPlain(redshiftFromLuminosityDistanceMpc(shared.d_L_mpc))}
      </>
    ) : (
      <>
        d<sub>L</sub> (Mpc)
        <span className="sb-distance-converted-eq"> = </span>
        {normalizeDistanceMpc(luminosityDistanceMpcFromRedshift(shared.z))}
      </>
    );
  const unlockToggle = (
    <label className="sb-checkbox-inline sb-distance-unlock">
      <input type="checkbox" checked={!distanceLinked} onChange={(e) => setDistanceLinked(!e.target.checked)} />
      unlock
    </label>
  );

  return (
    <div className="sb-stack sb-row-gap-top">
      <div className="sb-distance-panel">
        {distanceLinked ? (
          <div className="sb-distance-row sb-distance-row-linked">
            <div className="sb-distance-entry">
              <select className="sb-distance-unit" value={distanceDriver} onChange={(e) => setDistanceDriver(e.target.value as DistanceDriver)}>
                <option value="dL">dₗ (Mpc)</option>
                <option value="z">z</option>
              </select>
              <input
                className="sb-distance-number sb-distance-value"
                type="number"
                min={0}
                max={distanceDriver === "dL" ? 1e7 : 20}
                step={distanceDriver === "dL" ? 1 : 0.001}
                value={distanceDriver === "dL" ? shared.d_L_mpc : formatZPlain(shared.z)}
                onChange={(e) => {
                  const numeric = Number(e.target.value);
                  if (distanceDriver === "dL") {
                    setDistanceMpc(numeric);
                  } else {
                    setDistanceRedshift(numeric);
                  }
                }}
              />
            </div>
            <span className="sb-muted sb-distance-converted">{convertedText}</span>
            {unlockToggle}
          </div>
        ) : (
          <div className="sb-distance-row">
            <div className="sb-distance-unlocked-inline">
              <span className="sb-distance-prefix">
                d<sub>L</sub> (Mpc)
              </span>
              <input
                className="sb-distance-number sb-distance-value"
                type="number"
                min={0}
                max={1e7}
                step={1}
                value={shared.d_L_mpc}
                onChange={(e) => setDistanceMpc(Number(e.target.value))}
              />
              <span className="sb-distance-prefix">z</span>
              <input
                className="sb-distance-number sb-distance-value"
                type="number"
                min={0}
                max={20}
                step={0.001}
                value={formatZPlain(shared.z)}
                onChange={(e) => setDistanceRedshift(Number(e.target.value))}
              />
            </div>
            {unlockToggle}
          </div>
        )}
      </div>
      <SliderField
        label={
          <>
            θ<sub>obs</sub> (rad)
          </>
        }
        min={0}
        max={1.57}
        step={0.01}
        value={shared.theta_obs}
        onChange={(v) => setSharedField("theta_obs", v)}
      />
    </div>
  );
}
