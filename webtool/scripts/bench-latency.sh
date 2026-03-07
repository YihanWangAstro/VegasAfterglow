#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Benchmark VegasAfterglow webtool API latency.

Usage:
  bench-latency.sh [--direct-base URL] [--proxy-base URL] [--runs N] [--warmup N]

Examples:
  webtool/scripts/bench-latency.sh --direct-base https://api.vegasafterglow.com
  webtool/scripts/bench-latency.sh --direct-base https://api.vegasafterglow.com --proxy-base https://www.vegasafterglow.com
EOF
}

DIRECT_BASE=""
PROXY_BASE=""
RUNS=24
WARMUP=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    --direct-base)
      DIRECT_BASE="${2:-}"
      shift 2
      ;;
    --proxy-base)
      PROXY_BASE="${2:-}"
      shift 2
      ;;
    --runs)
      RUNS="${2:-}"
      shift 2
      ;;
    --warmup)
      WARMUP="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -z "$DIRECT_BASE" && -z "$PROXY_BASE" ]]; then
  echo "Set at least one of --direct-base or --proxy-base." >&2
  usage
  exit 2
fi

if ! [[ "$RUNS" =~ ^[0-9]+$ ]] || ! [[ "$WARMUP" =~ ^[0-9]+$ ]]; then
  echo "--runs and --warmup must be non-negative integers." >&2
  exit 2
fi

trim_slash() {
  local v="$1"
  echo "${v%/}"
}

DIRECT_BASE="$(trim_slash "$DIRECT_BASE")"
PROXY_BASE="$(trim_slash "$PROXY_BASE")"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
BENCH_SEED="$(( $(date +%s) % 10000 ))"

make_payload() {
  local endpoint="$1"
  case "$endpoint" in
    lightcurve)
      cat <<'JSON'
{"shared":{"d_L_mpc":100,"theta_obs":0.25,"flux_unit":"mJy","time_unit":"s","jet_type":"Top-hat","theta_c":0.1,"E_iso":1e52,"Gamma0":300,"spreading":false,"duration":1,"k_e":2,"k_g":2,"theta_w":0.3,"E_iso_w":1e51,"Gamma0_w":100,"medium_type":"ISM","n_ism":1,"A_star":0.1,"k_m":2,"eps_e":0.1,"eps_B":0.001,"p":2.3,"xi_e":1,"ssc":false,"kn":false,"enable_rvs":false,"eps_e_r":0.1,"eps_B_r":0.001,"p_r":2.3,"xi_e_r":1,"rvs_ssc":false,"rvs_kn":false,"num_t":100,"res_phi":0.1,"res_theta":0.25,"res_t":10},"frequencies_input":"1e9, R, 1keV","t_min":1,"t_max":1e8,"selected_instruments":[],"observation_groups":[],"include_figure":true,"include_exports":false}
JSON
      ;;
    spectrum)
      cat <<'JSON'
{"shared":{"d_L_mpc":100,"theta_obs":0.25,"flux_unit":"mJy","time_unit":"s","jet_type":"Top-hat","theta_c":0.1,"E_iso":1e52,"Gamma0":300,"spreading":false,"duration":1,"k_e":2,"k_g":2,"theta_w":0.3,"E_iso_w":1e51,"Gamma0_w":100,"medium_type":"ISM","n_ism":1,"A_star":0.1,"k_m":2,"eps_e":0.1,"eps_B":0.001,"p":2.3,"xi_e":1,"ssc":false,"kn":false,"enable_rvs":false,"eps_e_r":0.1,"eps_B_r":0.001,"p_r":2.3,"xi_e_r":1,"rvs_ssc":false,"rvs_kn":false,"num_t":100,"res_phi":0.1,"res_theta":0.25,"res_t":10},"t_snapshots_input":"1e3,1e4,1e5,1e6","nu_min":1e8,"nu_max":1e20,"num_nu":200,"freq_unit":"Hz","show_nufnu":false,"selected_instruments":[],"observation_groups":[],"include_figure":true,"include_exports":false}
JSON
      ;;
    skymap)
      cat <<'JSON'
{"shared":{"d_L_mpc":100,"theta_obs":0.25,"flux_unit":"cgs","time_unit":"s","jet_type":"Top-hat","theta_c":0.1,"E_iso":1e52,"Gamma0":300,"spreading":false,"duration":1,"k_e":2,"k_g":2,"theta_w":0.3,"E_iso_w":1e51,"Gamma0_w":100,"medium_type":"ISM","n_ism":1,"A_star":0.1,"k_m":2,"eps_e":0.1,"eps_B":0.001,"p":2.3,"xi_e":1,"ssc":false,"kn":false,"enable_rvs":false,"eps_e_r":0.1,"eps_B_r":0.001,"p_r":2.3,"xi_e_r":1,"rvs_ssc":false,"rvs_kn":false,"num_t":100,"res_phi":0.1,"res_theta":0.25,"res_t":10},"animate":false,"t_obs":1e6,"nu_input":"1e9","fov":500,"npixel":256,"include_figure":true,"include_exports":false}
JSON
      ;;
    *)
      echo "Unsupported endpoint: $endpoint" >&2
      return 2
      ;;
  esac
}

run_one() {
  local label="$1"
  local url="$2"
  local endpoint="$3"
  local out_csv="$4"
  local base_payload
  base_payload="$(make_payload "$endpoint")"

  local i
  for ((i = 0; i < WARMUP; i += 1)); do
    local payload
    payload="$(python3 - <<'PY' "$base_payload" "$endpoint" "$i" "$label" "$BENCH_SEED"
import json
import sys

raw, endpoint, idx, label, seed = sys.argv[1:6]
i = int(idx)
offset = (1000 if label == "proxy" else 0) + int(seed)
obj = json.loads(raw)
shared = obj.get("shared", {})
shared["theta_obs"] = 0.2 + 0.005 * (i + offset)
if endpoint == "skymap":
    obj["t_obs"] = 1e6 * (1.0 + 0.01 * (i + offset))
elif endpoint == "spectrum":
    obj["nu_min"] = 1e8 * (1.0 + 0.01 * (i + offset))
print(json.dumps(obj, separators=(",", ":")))
PY
)"
    curl -sS -o /dev/null -X POST "$url" \
      -H "Content-Type: application/json" \
      --data "$payload" >/dev/null
  done

  for ((i = 0; i < RUNS; i += 1)); do
    local payload
    payload="$(python3 - <<'PY' "$base_payload" "$endpoint" "$i" "$label" "$BENCH_SEED"
import json
import sys

raw, endpoint, idx, label, seed = sys.argv[1:6]
i = int(idx)
offset = (1000 if label == "proxy" else 0) + int(seed)
obj = json.loads(raw)
shared = obj.get("shared", {})
shared["theta_obs"] = 0.2 + 0.003 * (i + 1 + offset)
if endpoint == "skymap":
    obj["t_obs"] = 1e6 * (1.0 + 0.007 * (i + 1 + offset))
elif endpoint == "spectrum":
    obj["nu_min"] = 1e8 * (1.0 + 0.008 * (i + 1 + offset))
print(json.dumps(obj, separators=(",", ":")))
PY
)"
    local body="$TMP_DIR/body-${label}-${endpoint}-${i}.json"
    local metrics="$TMP_DIR/metrics-${label}-${endpoint}-${i}.txt"
    curl -sS -o "$body" -w '%{http_code} %{time_total} %{size_download}\n' \
      -X POST "$url" \
      -H "Content-Type: application/json" \
      --data "$payload" >"$metrics"

    python3 - <<'PY' "$metrics" "$body" "$out_csv" "$i" "$label" "$endpoint"
import json
import math
import sys
from pathlib import Path

metrics_path, body_path, out_csv, run_idx, label, endpoint = sys.argv[1:7]
http_code, total_s, size_download = Path(metrics_path).read_text().strip().split()
body = Path(body_path).read_text()
compute_s = float("nan")
try:
    payload = json.loads(body)
    compute_s = float(payload.get("meta", {}).get("compute_seconds", float("nan")))
except Exception:
    pass
rtt_ms = float(total_s) * 1000.0
compute_ms = compute_s * 1000.0 if math.isfinite(compute_s) else float("nan")
overhead_ms = rtt_ms - compute_ms if math.isfinite(compute_ms) else float("nan")
with open(out_csv, "a", encoding="utf-8") as f:
    f.write(f"{label},{endpoint},{run_idx},{http_code},{rtt_ms:.3f},{compute_ms:.3f},{overhead_ms:.3f},{size_download}\n")
PY
  done
}

CSV_OUT="$TMP_DIR/bench.csv"
echo "path,endpoint,run,http_code,rtt_ms,compute_ms,overhead_ms,size_download" >"$CSV_OUT"

if [[ -n "$DIRECT_BASE" ]]; then
  run_one "direct" "${DIRECT_BASE}/api/lightcurve" "lightcurve" "$CSV_OUT"
  run_one "direct" "${DIRECT_BASE}/api/spectrum" "spectrum" "$CSV_OUT"
  run_one "direct" "${DIRECT_BASE}/api/skymap" "skymap" "$CSV_OUT"
fi

if [[ -n "$PROXY_BASE" ]]; then
  run_one "proxy" "${PROXY_BASE}/api-proxy/lightcurve" "lightcurve" "$CSV_OUT"
  run_one "proxy" "${PROXY_BASE}/api-proxy/spectrum" "spectrum" "$CSV_OUT"
  run_one "proxy" "${PROXY_BASE}/api-proxy/skymap" "skymap" "$CSV_OUT"
fi

python3 - <<'PY' "$CSV_OUT"
import csv
import math
import statistics
import sys
from collections import defaultdict

csv_path = sys.argv[1]
groups = defaultdict(list)

with open(csv_path, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        key = (row["path"], row["endpoint"])
        groups[key].append(row)

def pct(values, q):
    if not values:
        return float("nan")
    idx = int(round((len(values) - 1) * q))
    idx = max(0, min(len(values) - 1, idx))
    return sorted(values)[idx]

print("Latency benchmark summary")
print("=========================")
for (path, endpoint), rows in sorted(groups.items()):
    ok_rows = [r for r in rows if r["http_code"] == "200"]
    rtt = [float(r["rtt_ms"]) for r in ok_rows]
    comp = [float(r["compute_ms"]) for r in ok_rows if math.isfinite(float(r["compute_ms"]))]
    overhead = [float(r["overhead_ms"]) for r in ok_rows if math.isfinite(float(r["overhead_ms"]))]
    sizes = [int(r["size_download"]) for r in ok_rows]
    print(f"- {path}/{endpoint}: ok={len(ok_rows)}/{len(rows)}")
    if not ok_rows:
        continue
    print(f"  rtt(ms):      p50={pct(rtt, 0.50):.1f} p90={pct(rtt, 0.90):.1f} p95={pct(rtt, 0.95):.1f}")
    if comp:
        print(f"  compute(ms):  p50={pct(comp, 0.50):.1f} p90={pct(comp, 0.90):.1f} p95={pct(comp, 0.95):.1f}")
    if overhead:
        print(f"  overhead(ms): p50={pct(overhead, 0.50):.1f} p90={pct(overhead, 0.90):.1f} p95={pct(overhead, 0.95):.1f}")
    if sizes:
        print(f"  payload:      p50={pct(sizes, 0.50)} bytes")
PY

echo
echo "Raw CSV: $CSV_OUT"
