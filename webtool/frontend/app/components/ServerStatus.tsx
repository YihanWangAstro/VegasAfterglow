type ServerRow = {
  key: string;
  regionText: string;
  statusClass: string;
  statusText: string;
};

type ActiveServerRow = ServerRow & {
  isActive: boolean;
};

type Props = {
  activeApiStatusRow: ActiveServerRow | null;
  otherApiStatusRows: ServerRow[];
  probeOtherServersOnce: () => void;
};

export function ServerStatus({ activeApiStatusRow, otherApiStatusRows, probeOtherServersOnce }: Props) {
  return (
    <div className="sidebar-api-status" aria-live="polite">
      {activeApiStatusRow ? (
        <div className="sidebar-api-row sidebar-api-row-primary">
          <div className="sidebar-api-primary">
            <span className="sidebar-api-name">
              <span className="sidebar-api-working">Working Server:</span>
              <span className="sidebar-api-location">{activeApiStatusRow.regionText}</span>
              <span className={`sidebar-api-pill ${activeApiStatusRow.statusClass}`}>{activeApiStatusRow.statusText}</span>
            </span>
          </div>
          {otherApiStatusRows.length > 0 ? (
            <details
              className="sidebar-api-dropdown-inline"
              onToggle={(event) => {
                if ((event.currentTarget as HTMLDetailsElement).open) {
                  probeOtherServersOnce();
                }
              }}
            >
              <summary>Other</summary>
              <div className="sidebar-api-dropdown-menu">
                {otherApiStatusRows.map((row) => (
                  <div key={row.key} className="sidebar-api-row">
                    <span className="sidebar-api-name">
                      <span className="sidebar-api-location">{row.regionText}</span>
                      <span className={`sidebar-api-pill ${row.statusClass}`}>{row.statusText}</span>
                    </span>
                  </div>
                ))}
              </div>
            </details>
          ) : null}
        </div>
      ) : null}
      <p className="sidebar-api-note">
        Mainland China users can temporarily use{" "}
        <a href="http://8.136.116.255" target="_blank" rel="noreferrer">
          8.136.116.255
        </a>
        {" "}for lower latency.
      </p>
    </div>
  );
}
