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

## 3) Deploy backend to two regions

Run from repo root:

```bash
cd /Users/yihanwang/Repositories/afterglow
```

Set common origins (replace custom domain if different):

```bash
export ALLOWED_ORIGINS="https://vegasafterglow.vercel.app,https://vegasafterglow.com"
```

Deploy US region:

```bash
PROJECT_ID=project-819bd021-ada6-4176-b14 \
REGION=us-west2 \
SERVICE_NAME=webtool-api \
ARTIFACT_REPO=containers \
MIN_INSTANCES=0 \
ALLOWED_ORIGINS="$ALLOWED_ORIGINS" \
bash webtool/scripts/deploy-cloudrun.sh
```

Deploy Asia region (prefer `asia-east2`; if unavailable, use `asia-northeast1`):

```bash
PROJECT_ID=project-819bd021-ada6-4176-b14 \
REGION=asia-east2 \
SERVICE_NAME=webtool-api \
ARTIFACT_REPO=containers \
MIN_INSTANCES=0 \
ALLOWED_ORIGINS="$ALLOWED_ORIGINS" \
bash webtool/scripts/deploy-cloudrun.sh
```

Check both service URLs:

```bash
gcloud run services describe webtool-api --region us-west2 --format='value(status.url)'
gcloud run services describe webtool-api --region asia-east2 --format='value(status.url)'
```

## 4) Create global HTTPS load balancer for API

All commands below run once unless noted.

Create serverless NEGs:

```bash
gcloud compute network-endpoint-groups create webtool-api-usw2-neg \
  --region=us-west2 \
  --network-endpoint-type=serverless \
  --cloud-run-service=webtool-api \
  --cloud-run-region=us-west2

gcloud compute network-endpoint-groups create webtool-api-asiae2-neg \
  --region=asia-east2 \
  --network-endpoint-type=serverless \
  --cloud-run-service=webtool-api \
  --cloud-run-region=asia-east2
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

Attach both NEGs:

```bash
gcloud compute backend-services add-backend webtool-api-global-bs \
  --global \
  --network-endpoint-group=webtool-api-usw2-neg \
  --network-endpoint-group-region=us-west2

gcloud compute backend-services add-backend webtool-api-global-bs \
  --global \
  --network-endpoint-group=webtool-api-asiae2-neg \
  --network-endpoint-group-region=asia-east2
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

At your DNS provider, create:

- type `A`
- host `api`
- value `<GLOBAL_LB_IP>`

Wait for SSL to become active:

```bash
gcloud compute ssl-certificates describe webtool-api-cert --global --format='value(managed.status)'
```

Expected status: `ACTIVE`.

## 5) Configure frontend API endpoint (Vercel env)

Use the global API domain (not a single-region Cloud Run URL):

```bash
cd /Users/yihanwang/Repositories/afterglow/webtool/frontend

vercel env rm NEXT_PUBLIC_API_URL production --yes || true
vercel env rm NEXT_PUBLIC_API_URL preview --yes || true

printf '%s\n' "https://api.vegasafterglow.com" | vercel env add NEXT_PUBLIC_API_URL production
printf '%s\n' "https://api.vegasafterglow.com" | vercel env add NEXT_PUBLIC_API_URL preview
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

Redeploy both regions:

```bash
cd /Users/yihanwang/Repositories/afterglow

PROJECT_ID=project-819bd021-ada6-4176-b14 REGION=us-west2 SERVICE_NAME=webtool-api ARTIFACT_REPO=containers MIN_INSTANCES=0 ALLOWED_ORIGINS="$ALLOWED_ORIGINS" bash webtool/scripts/deploy-cloudrun.sh
PROJECT_ID=project-819bd021-ada6-4176-b14 REGION=asia-east2 SERVICE_NAME=webtool-api ARTIFACT_REPO=containers MIN_INSTANCES=0 ALLOWED_ORIGINS="$ALLOWED_ORIGINS" bash webtool/scripts/deploy-cloudrun.sh
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

- `cpu=1`, `memory=2Gi`, `concurrency=20`, `max-instances=10`

For light traffic and cost control, keep `MIN_INSTANCES=0`.
