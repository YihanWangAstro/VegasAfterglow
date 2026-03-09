import { useEffect, useMemo, useRef, useState } from "react";
import { HelpHint } from "./HelpHint";
import { sliderFillStyle } from "../lib/utils/slider";
import type { ObservationGroup, SelectOption } from "../lib/types";

const SHIFT_MIN = -10;
const SHIFT_MAX = 10;
const SHIFT_STEP = 0.001;

function superscript(n: number): string {
  const map: Record<string, string> = { "-": "\u207B", "0": "\u2070", "1": "\u00B9", "2": "\u00B2", "3": "\u00B3", "4": "\u2074", "5": "\u2075", "6": "\u2076", "7": "\u2077", "8": "\u2078", "9": "\u2079" };
  return String(n).split("").map((c) => map[c] ?? c).join("");
}

function formatShiftLabel(shift: number): string {
  if (shift === 1) return "1";
  const exp = Math.floor(Math.log10(shift));
  const leading = Math.round(shift / Math.pow(10, exp));
  if (exp === 0) return String(leading);
  if (leading === 1) return `10${superscript(exp)}`;
  return `${leading}\u00b710${superscript(exp)}`;
}

function shiftToSlider(s: number): number {
  return Math.log10(Math.max(1e-10, Math.min(1e10, s)));
}

function sliderToShift(v: number): number {
  const exp = Math.floor(v);
  const leading = Math.round(Math.pow(10, v - exp));
  const clamped = Math.max(1, Math.min(9, leading));
  return clamped * Math.pow(10, exp);
}

function renderStringOptions(options: readonly string[]) {
  return options.map((item) => (
    <option key={item} value={item}>
      {item}
    </option>
  ));
}

type Props = {
  isLc: boolean;
  groups: ObservationGroup[];
  activeTab: number;
  setActiveTab: (value: number) => void;
  xOptions: readonly string[];
  xUnitLabel: string;
  rowLabel: string;
  obsHelpText: string;
  obsPlaceholder: string;
  updateObsGroup: (index: number, patch: Partial<ObservationGroup>) => void;
  removeObsGroup: (index: number) => void;
  addObsGroup: () => void;
  handleObservationUpload: (index: number, file: File) => Promise<void>;
  yUnitOptions: readonly string[];
  curveOptions?: SelectOption[];
};

export function ObservationEditor({
  isLc,
  groups,
  activeTab,
  setActiveTab,
  xOptions,
  xUnitLabel,
  rowLabel,
  obsHelpText,
  obsPlaceholder,
  updateObsGroup,
  removeObsGroup,
  addObsGroup,
  handleObservationUpload,
  yUnitOptions,
  curveOptions,
}: Props) {
  // Clear stale freq bindings when available curves change.
  // Uses refs for groups/updateObsGroup to avoid re-firing on every parent render.
  const groupsRef = useRef(groups);
  groupsRef.current = groups;
  const updateRef = useRef(updateObsGroup);
  updateRef.current = updateObsGroup;

  useEffect(() => {
    if (!curveOptions) return;
    groupsRef.current.forEach((group, idx) => {
      if (group.freq && !curveOptions.some((o) => o.value === group.freq)) {
        updateRef.current(idx, { freq: undefined });
      }
    });
  }, [curveOptions]);

  const activeIndex = groups.length === 0 ? -1 : Math.max(0, Math.min(activeTab, groups.length - 1));
  const activeGroup = activeIndex >= 0 ? groups[activeIndex] : null;

  const [shiftDraft, setShiftDraft] = useState(() => shiftToSlider(activeGroup?.shift ?? 1));
  useEffect(() => { setShiftDraft(shiftToSlider(activeGroup?.shift ?? 1)); }, [activeIndex, activeGroup?.shift]);
  const shiftFillStyle = useMemo(() => sliderFillStyle(shiftDraft, SHIFT_MIN, SHIFT_MAX), [shiftDraft]);

  const hasCurves = curveOptions && curveOptions.length > 0;
  const freqMatchesCurve = activeGroup?.freq && hasCurves && curveOptions.some((o) => o.value === activeGroup.freq);
  const isCustom = !activeGroup?.freq || !freqMatchesCurve;

  // Collect freqs already bound by other groups so we can exclude them from the dropdown.
  const takenFreqs = new Set(
    groups
      .filter((_, i) => i !== activeIndex)
      .map((g) => g.freq)
      .filter((f): f is string => !!f),
  );
  const availableCurveOptions = hasCurves
    ? curveOptions.filter((o) => !takenFreqs.has(o.value))
    : undefined;

  const legendField = activeGroup ? (
    <label className="sb-field">
      <span className="sb-label">Legend</span>
      <input value={activeGroup.legend} onChange={(e) => updateObsGroup(activeIndex, { legend: e.target.value })} />
    </label>
  ) : null;

  return (
    <details className="sb-expander">
      <summary>Observation Data</summary>
      <button
        className="sb-small-btn"
        type="button"
        onClick={() => {
          const nextIndex = groups.length;
          addObsGroup();
          setActiveTab(nextIndex);
        }}
      >
        + Add group
      </button>
      {groups.length === 0 ? <p className="sb-muted">No groups added.</p> : null}
      {groups.length > 0 ? (
        <div className="obs-tabs" role="tablist" aria-label={isLc ? "Light curve data groups" : "Spectrum data groups"}>
          {groups.map((group, idx) => (
            <button
              key={`${group.legend}-${idx}`}
              type="button"
              role="tab"
              aria-selected={idx === activeIndex}
              className={`obs-tab${idx === activeIndex ? " active" : ""}`}
              onClick={() => setActiveTab(idx)}
            >
              {group.legend.trim() || `data ${idx + 1}`}
            </button>
          ))}
        </div>
      ) : null}
      {activeGroup ? (
        <div className="obs-card" key={activeIndex}>
          <div className="obs-card-row obs-card-row-meta">
            {hasCurves ? (
              <>
                <label className="sb-field">
                  <span className="sb-label">Group</span>
                  <select
                    value={freqMatchesCurve ? activeGroup.freq! : ""}
                    onChange={(e) => {
                      const val = e.target.value;
                      if (val) {
                        const opt = curveOptions.find((o) => o.value === val);
                        updateObsGroup(activeIndex, { freq: val, legend: opt?.label ?? val });
                      } else {
                        updateObsGroup(activeIndex, { freq: undefined });
                      }
                    }}
                  >
                    <option value="">— custom —</option>
                    {(availableCurveOptions ?? []).map((opt) => (
                      <option key={opt.value} value={opt.value}>{opt.label}</option>
                    ))}
                  </select>
                </label>
                {isCustom ? legendField : null}
              </>
            ) : (
              legendField
            )}
            <label className="sb-checkbox-inline">
              <input
                type="checkbox"
                checked={activeGroup.visible}
                onChange={(e) => updateObsGroup(activeIndex, { visible: e.target.checked })}
              />
              Visible
            </label>
          </div>
          <div className="obs-card-row obs-card-row-units">
            <label className="sb-field">
              <span className="sb-label">{xUnitLabel}</span>
              <select value={activeGroup.x_unit} onChange={(e) => updateObsGroup(activeIndex, { x_unit: e.target.value })}>
                {renderStringOptions(xOptions)}
              </select>
            </label>
            <label className="sb-field">
              <span className="sb-label">flux/flux den unit</span>
              <select value={activeGroup.y_unit} onChange={(e) => updateObsGroup(activeIndex, { y_unit: e.target.value })}>
                {renderStringOptions(yUnitOptions)}
              </select>
            </label>
            <label className="sb-field sb-slider" style={{ gridColumn: "1 / -1" }}>
              <span className="sb-label">y shift: {formatShiftLabel(activeGroup.shift ?? 1)}</span>
              <div className="sb-slider-track">
                <input
                  type="range"
                  min={SHIFT_MIN}
                  max={SHIFT_MAX}
                  step={SHIFT_STEP}
                  value={shiftDraft}
                  style={shiftFillStyle}
                  onInput={(e) => {
                    const raw = parseFloat((e.target as HTMLInputElement).value);
                    setShiftDraft(raw);
                    const snapped = sliderToShift(raw);
                    updateObsGroup(activeIndex, { shift: snapped });
                  }}
                  onChange={(e) => {
                    const raw = parseFloat((e.target as HTMLInputElement).value);
                    setShiftDraft(raw);
                    const snapped = sliderToShift(raw);
                    updateObsGroup(activeIndex, { shift: snapped });
                  }}
                />
                <span className="sb-value">{formatShiftLabel(activeGroup.shift ?? 1)}</span>
              </div>
            </label>
          </div>
          <label className="sb-field">
            <span className="sb-label sb-label-with-help">
              <span>{rowLabel}</span>
              <HelpHint text={obsHelpText} ariaLabel="Observation input help" />
            </span>
            <textarea
              rows={4}
              value={activeGroup.text}
              placeholder={obsPlaceholder}
              onChange={(e) => updateObsGroup(activeIndex, { text: e.target.value })}
            />
          </label>
          <div className="obs-card-actions">
            <label className="sb-upload-btn">
              Upload file
              <input
                type="file"
                accept=".csv,.txt,.dat,.xls,.xlsx,text/csv,text/plain,application/vnd.ms-excel,application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                onChange={(e) => {
                  const file = e.target.files?.[0];
                  if (file) {
                    void handleObservationUpload(activeIndex, file);
                  }
                  e.currentTarget.value = "";
                }}
              />
            </label>
            <button className="sb-danger-btn" type="button" onClick={() => removeObsGroup(activeIndex)}>
              Remove
            </button>
          </div>
          <p className="sb-muted sb-upload-hint">Supports .csv, .txt, .xls, .xlsx</p>
        </div>
      ) : null}
    </details>
  );
}
