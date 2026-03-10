import type { ReactNode } from "react";
import {
  FREQ_UNIT_OPTIONS,
  JET_ENERGY_ROW,
  JET_POWERLAW_ROW,
  JET_THETA_C_SPEC,
  JET_TWOCOMP_GAMMA_SPEC,
  JET_TWOCOMP_ROW,
  JET_TYPE_OPTIONS,
  MEDIUM_ISM_SPEC,
  MEDIUM_KM_SPEC,
  MEDIUM_TYPE_OPTIONS,
  MEDIUM_WIND_FLOOR_SPEC,
  MEDIUM_WIND_SPEC,
  RADIATION_SLIDER_ROWS,
  RESOLUTION_SLIDERS,
  REVERSE_SHOCK_DURATION_SPEC,
  REVERSE_SHOCK_SLIDER_ROWS,
  TIME_UNIT_OPTIONS,
  Y_UNIT_OPTIONS,
} from "../lib/constants";
import type { InstrumentGroup, Mode, ObservationGroup, SelectOption, SharedParams, SharedSliderSpec } from "../lib/types";
import { defaultObsGroup, parseObservationUpload } from "../lib/utils/obs";
import type { UpdateValue } from "../hooks/useParameterState";
import { LogSliderField, SliderField } from "./SliderField";
import { ObservationEditor } from "./ObservationEditor";

type Props = {
  mode: Mode;
  shared: SharedParams;
  setSharedField: <K extends keyof SharedParams>(key: K, value: SharedParams[K]) => void;
  hideDistanceObserver?: boolean;
  hidePlotUnits?: boolean;
  distanceControls?: ReactNode;
  plotUnitControls?: ReactNode;
  instrumentGroups: InstrumentGroup[];
  lcInstruments: string[];
  sedInstruments: string[];
  setLcInstruments: (value: UpdateValue<string[]>) => void;
  setSedInstruments: (value: UpdateValue<string[]>) => void;
  lcObsGroups: ObservationGroup[];
  sedObsGroups: ObservationGroup[];
  activeLcObsTab: number;
  activeSedObsTab: number;
  setLcObsGroups: (value: UpdateValue<ObservationGroup[]>) => void;
  setSedObsGroups: (value: UpdateValue<ObservationGroup[]>) => void;
  setActiveLcObsTab: (value: UpdateValue<number>) => void;
  setActiveSedObsTab: (value: UpdateValue<number>) => void;
  setError: (value: string) => void;
  lcCurveOptions?: SelectOption[];
  sedCurveOptions?: SelectOption[];
};

export function SharedControls({
  mode,
  shared,
  setSharedField,
  hideDistanceObserver,
  hidePlotUnits,
  distanceControls,
  plotUnitControls,
  instrumentGroups,
  lcInstruments,
  sedInstruments,
  setLcInstruments,
  setSedInstruments,
  lcObsGroups,
  sedObsGroups,
  activeLcObsTab,
  activeSedObsTab,
  setLcObsGroups,
  setSedObsGroups,
  setActiveLcObsTab,
  setActiveSedObsTab,
  setError,
  lcCurveOptions,
  sedCurveOptions,
}: Props) {
  const isSky = mode === "skymap";
  const isLc = mode === "lightcurve";

  function obsSetters(isLcMode: boolean) {
    return {
      setGroups: isLcMode ? setLcObsGroups : setSedObsGroups,
      setActiveTab: isLcMode ? setActiveLcObsTab : setActiveSedObsTab,
    };
  }

  function updateObsGroup(isLcMode: boolean, index: number, patch: Partial<ObservationGroup>) {
    obsSetters(isLcMode).setGroups((prev) => prev.map((group, i) => (i === index ? { ...group, ...patch } : group)));
  }

  function removeObsGroup(isLcMode: boolean, index: number) {
    const { setGroups, setActiveTab } = obsSetters(isLcMode);
    setGroups((prev) => prev.filter((_, i) => i !== index));
    setActiveTab((prev) => {
      if (index < prev) return prev - 1;
      if (index === prev) return Math.max(0, prev - 1);
      return prev;
    });
  }

  function addObsGroup(isLcMode: boolean) {
    obsSetters(isLcMode).setGroups((prev) => [...prev, defaultObsGroup(isLcMode, prev.length + 1)]);
  }

  function toggleInstrument(isLcMode: boolean, name: string, checked: boolean) {
    (isLcMode ? setLcInstruments : setSedInstruments)((prev) => {
      if (checked) return prev.includes(name) ? prev : [...prev, name];
      return prev.filter((value) => value !== name);
    });
  }

  async function handleObservationUpload(isLcMode: boolean, index: number, file: File) {
    setError("");
    try {
      const parsed = await parseObservationUpload(file);
      updateObsGroup(isLcMode, index, { text: parsed });
    } catch (err) {
      const message = err instanceof Error ? err.message : "Failed to parse observation file.";
      setError(`Upload failed: ${message}`);
    }
  }

  function renderObservationEditor(isLcMode: boolean) {
    const groups = isLcMode ? lcObsGroups : sedObsGroups;
    const activeTab = isLcMode ? activeLcObsTab : activeSedObsTab;
    const setActiveTab = isLcMode ? setActiveLcObsTab : setActiveSedObsTab;
    const xOptions = isLcMode ? TIME_UNIT_OPTIONS : FREQ_UNIT_OPTIONS;
    const xName = isLcMode ? "t" : "nu";
    const xUnitLabel = isLcMode ? "t unit" : "ν unit";
    const rowLabel = isLcMode ? "Rows: t, value, error" : "Rows: ν, value, error";
    const obsHelpText = `One row per line. Columns: ${xName}, value, error (optional). Separated by space, tab, or comma.`;
    const obsPlaceholder = `${xName}  value  error\n1e4  0.5  0.1\n1e5  0.3  0.05`;

    return (
      <ObservationEditor
        isLc={isLcMode}
        groups={groups}
        activeTab={activeTab}
        setActiveTab={(value) => setActiveTab(value)}
        xOptions={xOptions}
        xUnitLabel={xUnitLabel}
        rowLabel={rowLabel}
        obsHelpText={obsHelpText}
        obsPlaceholder={obsPlaceholder}
        updateObsGroup={(index, patch) => updateObsGroup(isLcMode, index, patch)}
        removeObsGroup={(index) => removeObsGroup(isLcMode, index)}
        addObsGroup={() => addObsGroup(isLcMode)}
        handleObservationUpload={(index, file) => handleObservationUpload(isLcMode, index, file)}
        yUnitOptions={Y_UNIT_OPTIONS}
        curveOptions={isLcMode ? lcCurveOptions : sedCurveOptions}
      />
    );
  }

  function renderSharedSlider(spec: SharedSliderSpec) {
    if (spec.kind === "log") {
      return (
        <LogSliderField
          key={String(spec.key)}
          label={spec.label}
          minExp={spec.minExp}
          maxExp={spec.maxExp}
          step={spec.step}
          value={shared[spec.key] as number}
          defaultExp={spec.defaultExp}
          onChange={(value) => setSharedField(spec.key, value as SharedParams[typeof spec.key])}
        />
      );
    }
    return (
      <SliderField
        key={String(spec.key)}
        label={spec.label}
        min={spec.min}
        max={spec.max}
        step={spec.step}
        decimals={spec.decimals}
        value={shared[spec.key] as number}
        onChange={(value) => setSharedField(spec.key, value as SharedParams[typeof spec.key])}
      />
    );
  }

  function renderSharedSliderRow(specs: SharedSliderSpec[], className = "sb-row-2") {
    return <div className={className}>{specs.map((spec) => renderSharedSlider(spec))}</div>;
  }

  return (
    <>
      {!hideDistanceObserver ? <div className="sb-stack">{distanceControls}</div> : null}

      {!isSky && !hidePlotUnits ? plotUnitControls : null}

      <details className="sb-expander sb-major" open>
        <summary>Jet</summary>
        <div className="sb-group-content">
          <div className="sb-row-2 sb-row-jet">
            <div className="sb-field">
              <span className="sb-label">Type</span>
              <div className="mode-radios sb-choice-radios">
                {JET_TYPE_OPTIONS.map((jetType) => (
                  <label key={jetType} className="mode-radio">
                    <input
                      type="radio"
                      name="jet-type"
                      checked={shared.jet_type === jetType}
                      onChange={() => setSharedField("jet_type", jetType)}
                    />
                    <span>{jetType}</span>
                  </label>
                ))}
              </div>
            </div>
            {renderSharedSlider(JET_THETA_C_SPEC)}
          </div>

          {renderSharedSliderRow(JET_ENERGY_ROW)}

          <label className="sb-checkbox-inline">
            <input type="checkbox" checked={shared.spreading} onChange={(e) => setSharedField("spreading", e.target.checked)} />
            Spreading
          </label>

          {shared.jet_type === "Power-law" ? renderSharedSliderRow(JET_POWERLAW_ROW) : null}

          {shared.jet_type === "Two-component" ? (
            <>
              {renderSharedSliderRow(JET_TWOCOMP_ROW)}
              {renderSharedSlider(JET_TWOCOMP_GAMMA_SPEC)}
            </>
          ) : null}
        </div>
      </details>

      <details className="sb-expander sb-major">
        <summary>Medium</summary>
        <div className="sb-group-content">
          <div className="sb-row-2 sb-row-medium">
            <div className="sb-field">
              <span className="sb-label">Type</span>
              <div className="mode-radios sb-choice-radios">
                {MEDIUM_TYPE_OPTIONS.map((mediumType) => (
                  <label key={mediumType.value} className="mode-radio">
                    <input
                      type="radio"
                      name="medium-type"
                      checked={shared.medium_type === mediumType.value}
                      onChange={() => setSharedField("medium_type", mediumType.value)}
                    />
                    <span>{mediumType.label}</span>
                  </label>
                ))}
              </div>
            </div>
            {renderSharedSlider(shared.medium_type === "ISM" ? MEDIUM_ISM_SPEC : MEDIUM_WIND_SPEC)}
          </div>

          {shared.medium_type !== "ISM" ? renderSharedSlider(MEDIUM_KM_SPEC) : null}
          {shared.medium_type === "Wind bubble" ? renderSharedSlider(MEDIUM_WIND_FLOOR_SPEC) : null}
        </div>
      </details>

      <details className="sb-expander sb-major">
        <summary>Radiation</summary>
        <div className="sb-group-content">
          {RADIATION_SLIDER_ROWS.map((row, idx) => (
            <div className="sb-row-2" key={`rad-row-${idx}`}>
              {row.map((spec) => renderSharedSlider(spec))}
            </div>
          ))}
          <div className="sb-row-2 sb-check-row">
            <label className="sb-checkbox-inline">
              <input type="checkbox" checked={shared.ssc} onChange={(e) => setSharedField("ssc", e.target.checked)} />
              SSC
            </label>
            <label className="sb-checkbox-inline">
              <input type="checkbox" checked={shared.kn} onChange={(e) => setSharedField("kn", e.target.checked)} />
              KN
            </label>
          </div>
        </div>
      </details>

      <details className="sb-expander">
        <summary>Reverse Shock</summary>
        <label className="sb-checkbox-inline">
          <input type="checkbox" checked={shared.enable_rvs} onChange={(e) => setSharedField("enable_rvs", e.target.checked)} />
          Enable
        </label>
        {shared.enable_rvs ? (
          <>
            {renderSharedSlider(REVERSE_SHOCK_DURATION_SPEC)}
            {REVERSE_SHOCK_SLIDER_ROWS.map((row, idx) => (
              <div className="sb-row-2" key={`rvs-row-${idx}`}>
                {row.map((spec) => renderSharedSlider(spec))}
              </div>
            ))}
            <div className="sb-row-2 sb-check-row">
              <label className="sb-checkbox-inline">
                <input type="checkbox" checked={shared.rvs_ssc} onChange={(e) => setSharedField("rvs_ssc", e.target.checked)} />
                SSC
              </label>
              <label className="sb-checkbox-inline">
                <input type="checkbox" checked={shared.rvs_kn} onChange={(e) => setSharedField("rvs_kn", e.target.checked)} />
                KN
              </label>
            </div>
          </>
        ) : null}
      </details>

      {!isSky ? (
        <details className="sb-expander">
          <summary>Instrument Sensitivities</summary>
          <div className="instrument-grid">
            {instrumentGroups.length === 0 ? <p className="sb-muted">No instrument list loaded.</p> : null}
            {instrumentGroups.map((group) => (
              <section key={group.label} className="instrument-group">
                <p className="instrument-group-title">{group.label}</p>
                <div className="instrument-group-grid">
                  {group.items.map((name) => {
                    const selected = isLc ? lcInstruments.includes(name) : sedInstruments.includes(name);
                    return (
                      <label key={name} className="sb-checkbox-inline instrument-option">
                        <input
                          type="checkbox"
                          checked={selected}
                          onChange={(e) => toggleInstrument(isLc, name, e.target.checked)}
                        />
                        {name}
                      </label>
                    );
                  })}
                </div>
              </section>
            ))}
          </div>
        </details>
      ) : null}

      {!isSky ? renderObservationEditor(isLc) : null}

      <details className="sb-expander">
        <summary>Resolutions</summary>
        <div className="sb-stack sb-resolution-ppd-stack">{RESOLUTION_SLIDERS.map((spec) => renderSharedSlider(spec))}</div>
      </details>
    </>
  );
}
