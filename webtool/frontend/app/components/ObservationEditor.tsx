import { HelpHint } from "./HelpHint";
import type { ObservationGroup } from "../lib/types";

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
  xName: string;
  xUnitLabel: string;
  rowLabel: string;
  obsHelpText: string;
  obsPlaceholder: string;
  updateObsGroup: (index: number, patch: Partial<ObservationGroup>) => void;
  removeObsGroup: (index: number) => void;
  addObsGroup: () => void;
  handleObservationUpload: (index: number, file: File) => Promise<void>;
  yUnitOptions: readonly string[];
};

export function ObservationEditor({
  isLc,
  groups,
  activeTab,
  setActiveTab,
  xOptions,
  xName,
  xUnitLabel,
  rowLabel,
  obsHelpText,
  obsPlaceholder,
  updateObsGroup,
  removeObsGroup,
  addObsGroup,
  handleObservationUpload,
  yUnitOptions,
}: Props) {
  const activeIndex = groups.length === 0 ? -1 : Math.max(0, Math.min(activeTab, groups.length - 1));
  const activeGroup = activeIndex >= 0 ? groups[activeIndex] : null;

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
        <div className="obs-card" key={`${activeGroup.legend}-${activeIndex}`}>
          <div className="obs-card-row obs-card-row-meta">
            <label className="sb-field">
              <span className="sb-label">Legend</span>
              <input value={activeGroup.legend} onChange={(e) => updateObsGroup(activeIndex, { legend: e.target.value })} />
            </label>
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
