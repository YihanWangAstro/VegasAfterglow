# webtool Deployment Runbook (Vercel + Cloud Run + Global API)

This runbook matches your current project:

- `PROJECT_ID=project-819bd021-ada6-4176-b14`
- frontend: Vercel (`webtool/frontend`)
- backend: FastAPI on Cloud Run (`webtool-api`)
- recommended global API domain: `api.vegasafterglow.com`

This setup is optimized for realtime slider updates:

- frontend calls backend API directly first (no extra proxy hop)
- Next.js `/api-proxy/*` is kept as fallback only
- multi-region backend + global LB reduces user RTT for US/Asia users

## 0) One-time local tools

```bash
brew install --cask google-cloud-sdk
npm install -g vercel
```

## 1) First-time auth and project config

```bash
gcloud auth login
gcloud config set project project-819bd021-ada6-4176-b14

vercel login
```

Enable required Google APIs:

```bash
gcloud services enable \
  run.googleapis.com \
  cloudbuild.googleapis.com \
  artifactregistry.googleapis.com \
  compute.googleapis.com \
  certificatemanager.googleapis.com
```

## 2) Vercel project setup (one-time)

In Vercel dashboard:

1. import this repository
2. set `Root Directory` = `webtool/frontend`
3. build command = `npm run build`
4. install command = `npm install --legacy-peer-deps`

(`webtool/frontend/vercel.json` already sets install/build commands.)

## 3) Deploy backend to multi regions

Run from repo root:

```bash
cd /Users/yihanwang/Repositories/afterglow
```

Set common origins (must include every frontend origin you use):

```bash
export ALLOWED_ORIGINS="https://www.vegasafterglow.com,https://vegasafterglow.com,https://vegasafterglow.vercel.app"
```

Deploy US West:

```bash
PROJECT_ID=project-819bd021-ada6-4176-b14 \
REGION=us-west2 \
SERVICE_NAME=webtool-api \
ARTIFACT_REPO=containers \
MIN_INSTANCES=0 \
ALLOWED_ORIGINS="$ALLOWED_ORIGINS" \
bash webtool/scripts/deploy-cloudrun.sh
```

Deploy US East:

```bash
PROJECT_ID=project-819bd021-ada6-4176-b14 \
REGION=us-east4 \
SERVICE_NAME=webtool-api \
ARTIFACT_REPO=containers \
MIN_INSTANCES=0 \
ALLOWED_ORIGINS="$ALLOWED_ORIGINS" \
bash webtool/scripts/deploy-cloudrun.sh
```

Deploy Europe:

```bash
PROJECT_ID=project-819bd021-ada6-4176-b14 \
REGION=europe-west4 \
SERVICE_NAME=webtool-api \
ARTIFACT_REPO=containers \
MIN_INSTANCES=0 \
ALLOWED_ORIGINS="$ALLOWED_ORIGINS" \
bash webtool/scripts/deploy-cloudrun.sh
```

Deploy Asia (Hong Kong):

```bash
PROJECT_ID=project-819bd021-ada6-4176-b14 \
REGION=asia-east2 \
SERVICE_NAME=webtool-api \
ARTIFACT_REPO=containers \
MIN_INSTANCES=0 \
ALLOWED_ORIGINS="$ALLOWED_ORIGINS" \
bash webtool/scripts/deploy-cloudrun.sh
```

Deploy Asia (Tokyo):

```bash
PROJECT_ID=project-819bd021-ada6-4176-b14 \
REGION=asia-northeast1 \
SERVICE_NAME=webtool-api \
ARTIFACT_REPO=containers \
MIN_INSTANCES=0 \
ALLOWED_ORIGINS="$ALLOWED_ORIGINS" \
bash webtool/scripts/deploy-cloudrun.sh
```

Check all service URLs:

```bash
gcloud run services describe webtool-api --region us-west2 --format='value(status.url)'
gcloud run services describe webtool-api --region us-east4 --format='value(status.url)'
gcloud run services describe webtool-api --region europe-west4 --format='value(status.url)'
gcloud run services describe webtool-api --region asia-east2 --format='value(status.url)'
gcloud run services describe webtool-api --region asia-northeast1 --format='value(status.url)'
```

## 4) Create global HTTPS load balancer for API

All commands below run once unless noted.

Create serverless NEGs:

Note: for serverless NEG, Cloud Run region is derived from `--region`.
If your `gcloud` does not recognize `--cloud-run-region`, do not pass it.

```bash
gcloud compute network-endpoint-groups create webtool-api-usw2-neg \
  --region=us-west2 \
  --network-endpoint-type=serverless \
  --cloud-run-service=webtool-api

gcloud compute network-endpoint-groups create webtool-api-asiae2-neg \
  --region=asia-east2 \
  --network-endpoint-type=serverless \
  --cloud-run-service=webtool-api

gcloud compute network-endpoint-groups create webtool-api-use4-neg \
  --region=us-east4 \
  --network-endpoint-type=serverless \
  --cloud-run-service=webtool-api

gcloud compute network-endpoint-groups create webtool-api-euw4-neg \
  --region=europe-west4 \
  --network-endpoint-type=serverless \
  --cloud-run-service=webtool-api

gcloud compute network-endpoint-groups create webtool-api-asiane1-neg \
  --region=asia-northeast1 \
  --network-endpoint-type=serverless \
  --cloud-run-service=webtool-api
```

Create global backend service:

```bash
gcloud compute backend-services create webtool-api-global-bs \
  --global \
  --load-balancing-scheme=EXTERNAL_MANAGED \
  --protocol=HTTP \
  --port-name=http \
  --timeout=30s
```

Attach all NEGs:

```bash
gcloud compute backend-services add-backend webtool-api-global-bs \
  --global \
  --network-endpoint-group=webtool-api-usw2-neg \
  --network-endpoint-group-region=us-west2

gcloud compute backend-services add-backend webtool-api-global-bs \
  --global \
  --network-endpoint-group=webtool-api-asiae2-neg \
  --network-endpoint-group-region=asia-east2

gcloud compute backend-services add-backend webtool-api-global-bs \
  --global \
  --network-endpoint-group=webtool-api-use4-neg \
  --network-endpoint-group-region=us-east4

gcloud compute backend-services add-backend webtool-api-global-bs \
  --global \
  --network-endpoint-group=webtool-api-euw4-neg \
  --network-endpoint-group-region=europe-west4

gcloud compute backend-services add-backend webtool-api-global-bs \
  --global \
  --network-endpoint-group=webtool-api-asiane1-neg \
  --network-endpoint-group-region=asia-northeast1
```

Create URL map + certificate + HTTPS proxy + global IP + forwarding rule:

```bash
gcloud compute url-maps create webtool-api-urlmap \
  --default-service=webtool-api-global-bs

gcloud compute ssl-certificates create webtool-api-cert \
  --domains=api.vegasafterglow.com \
  --global

gcloud compute target-https-proxies create webtool-api-https-proxy \
  --url-map=webtool-api-urlmap \
  --ssl-certificates=webtool-api-cert

gcloud compute addresses create webtool-api-lb-ip --global

gcloud compute forwarding-rules create webtool-api-https-fr \
  --global \
  --target-https-proxy=webtool-api-https-proxy \
  --ports=443 \
  --address=webtool-api-lb-ip
```

Get LB IP:

```bash
gcloud compute addresses describe webtool-api-lb-ip --global --format='value(address)'
```

Set it as a shell var (recommended):

```bash
export GLOBAL_LB_IP="$(gcloud compute addresses describe webtool-api-lb-ip --global --format='value(address)')"
echo "$GLOBAL_LB_IP"
```

### 4.1 DNS setup (critical)

At your DNS provider for `vegasafterglow.com`, create this record:

- Type: `A`
- Name/Host: `api`
- Value/Content: `<GLOBAL_LB_IP>`
- TTL: `Auto` (or `300`)

Notes:

- If using Cloudflare, set `Proxy status = DNS only` (gray cloud) while cert is provisioning.
- Do not create a conflicting `CNAME` for `api`.
- If you have an old `A` record for `api`, replace it with `GLOBAL_LB_IP`.

### 4.2 Verify DNS propagation

```bash
dig +short api.vegasafterglow.com A
dig @8.8.8.8 +short api.vegasafterglow.com A
dig +short api.vegasafterglow.com AAAA
```

Expected:

- `A` output contains exactly `GLOBAL_LB_IP`
- `AAAA` is empty unless you also configured IPv6 forwarding

### 4.3 Verify managed cert status

```bash
gcloud compute ssl-certificates describe webtool-api-cert --global \
  --format='yaml(managed.status,managed.domainStatus,managed.domains)'
```

Expected final state:

- `managed.status: ACTIVE`
- `managed.domainStatus.api.vegasafterglow.com: ACTIVE`

If it stays `PROVISIONING` for more than 1-2 hours:

1. Re-check DNS A record points to `GLOBAL_LB_IP`
2. Ensure Cloudflare proxy is off (DNS only)
3. Wait for DNS propagation, then check cert status again

## 5) Configure frontend API endpoint (Vercel env)

Use the global API domain (not a single-region Cloud Run URL):

```bash
cd /Users/yihanwang/Repositories/afterglow/webtool/frontend

vercel env rm NEXT_PUBLIC_API_URL production --yes || true
vercel env rm NEXT_PUBLIC_API_URL preview --yes || true

printf 'https://api.vegasafterglow.com' | vercel env add NEXT_PUBLIC_API_URL production
printf 'https://api.vegasafterglow.com' | vercel env add NEXT_PUBLIC_API_URL preview

vercel env rm NEXT_PUBLIC_API_STATUS_URLS production --yes || true
vercel env rm NEXT_PUBLIC_API_STATUS_URLS preview --yes || true

printf 'https://webtool-api-qlcoh6gldq-wl.a.run.app,https://webtool-api-qlcoh6gldq-uk.a.run.app,https://webtool-api-qlcoh6gldq-ez.a.run.app,https://webtool-api-qlcoh6gldq-df.a.run.app,https://webtool-api-qlcoh6gldq-an.a.run.app' | vercel env add NEXT_PUBLIC_API_STATUS_URLS production
printf 'https://webtool-api-qlcoh6gldq-wl.a.run.app,https://webtool-api-qlcoh6gldq-uk.a.run.app,https://webtool-api-qlcoh6gldq-ez.a.run.app,https://webtool-api-qlcoh6gldq-df.a.run.app,https://webtool-api-qlcoh6gldq-an.a.run.app' | vercel env add NEXT_PUBLIC_API_STATUS_URLS preview
```

Deploy frontend:

```bash
vercel --prod
```

## 6) Verify end-to-end

Health check via global API:

```bash
curl https://api.vegasafterglow.com/api/health
```

Expected:

```json
{"status":"ok","service":"fastapi"}
```

Then open your frontend and verify:

1. Light Curve realtime slider updates
2. Spectrum realtime slider updates
3. Sky Image render + GIF download

## 7) Future updates

### 7.1 Backend code update

The build uses Docker layer caching: if only app code changed (not `requirements*.txt`), the
VegasAfterglow C++ compilation layer is reused and the build finishes in ~1 min instead of ~5 min.

Redeploy all active regions:

```bash
cd /Users/yihanwang/Repositories/afterglow

PROJECT_ID=project-819bd021-ada6-4176-b14 REGION=us-west2 SERVICE_NAME=webtool-api ARTIFACT_REPO=containers MIN_INSTANCES=0 ALLOWED_ORIGINS="$ALLOWED_ORIGINS" bash webtool/scripts/deploy-cloudrun.sh
PROJECT_ID=project-819bd021-ada6-4176-b14 REGION=us-east4 SERVICE_NAME=webtool-api ARTIFACT_REPO=containers MIN_INSTANCES=0 ALLOWED_ORIGINS="$ALLOWED_ORIGINS" bash webtool/scripts/deploy-cloudrun.sh
PROJECT_ID=project-819bd021-ada6-4176-b14 REGION=europe-west4 SERVICE_NAME=webtool-api ARTIFACT_REPO=containers MIN_INSTANCES=0 ALLOWED_ORIGINS="$ALLOWED_ORIGINS" bash webtool/scripts/deploy-cloudrun.sh
PROJECT_ID=project-819bd021-ada6-4176-b14 REGION=asia-east2 SERVICE_NAME=webtool-api ARTIFACT_REPO=containers MIN_INSTANCES=0 ALLOWED_ORIGINS="$ALLOWED_ORIGINS" bash webtool/scripts/deploy-cloudrun.sh
PROJECT_ID=project-819bd021-ada6-4176-b14 REGION=asia-northeast1 SERVICE_NAME=webtool-api ARTIFACT_REPO=containers MIN_INSTANCES=0 ALLOWED_ORIGINS="$ALLOWED_ORIGINS" bash webtool/scripts/deploy-cloudrun.sh
```

No Vercel env change is needed if API domain is unchanged.

### 7.2 Frontend-only update

```bash
cd /Users/yihanwang/Repositories/afterglow/webtool/frontend
vercel --prod
```

### 7.3 Change API domain or endpoint

Update Vercel `NEXT_PUBLIC_API_URL` (production + preview) and redeploy frontend.

## 8) Cost / cold start knobs

- `MIN_INSTANCES=0`: cheapest, but cold start can happen after idle period.
- `MIN_INSTANCES=1`: lower latency, but not free (charged idle baseline).

Current script defaults:

- `cpu=1`, `memory=2Gi`, `concurrency=8`, `max-instances=10`

For light traffic and cost control, keep `MIN_INSTANCES=0`.

## 9) Build caching details

The Dockerfile uses a multi-stage build to keep the runtime image small:

- **builder stage** — installs build tools (cmake, ninja, gcc) and compiles VegasAfterglow from source
- **runtime stage** — copies only the Python venv; no build tools in the final image

Dependency layers:

- local VegasAfterglow source tree — compiled in the builder stage
- `requirements.txt` — pure-Python deps (fastapi, uvicorn, plotly, etc.)

`cloudbuild.yaml` passes `--cache-from :latest` to Docker so Cloud Build can reuse unchanged
layers across builds. On each successful build, both a commit-tagged image and `:latest` are
pushed to Artifact Registry.

**Cache invalidation rules:**

| What changed | Layers rebuilt |
| --- | --- |
| `webtool/backend/app/**` only | none (copy layer only) |
| `requirements.txt` | pure-Python deps layer |
| `pyproject.toml`, `src/**`, `include/**`, `pybind/**`, `external/**`, `VegasAfterglow/**` | VegasAfterglow C++ build + all layers above |
| Base image (`python:3.12-slim`) | everything |
