# webtool Deployment (Step-by-Step)

This is the exact deployment runbook for this repository.

Constants used in this project:
- `PROJECT_ID=project-819bd021-ada6-4176-b14`
- `REGION=us-west2`
- `SERVICE_NAME=webtool-api`
- `ARTIFACT_REPO=containers`
- `REPO_ROOT=/Users/yihanwang/Repositories/afterglow`

Preferred backend URL (Cloud Run service URL):
- `https://webtool-api-476587253129.us-west2.run.app`

## 0) Prerequisites

Run these once on your machine:

```bash
brew install --cask google-cloud-sdk
npm install -g vercel
```

## 1) Login and project setup

```bash
gcloud auth login
gcloud config set project project-819bd021-ada6-4176-b14
gcloud services enable run.googleapis.com cloudbuild.googleapis.com artifactregistry.googleapis.com

vercel login
```

## 2) Vercel project settings (one-time)

In Vercel dashboard for this project:
1. Set `Root Directory` to `webtool/frontend`
2. Set `Install Command` to `npm install --legacy-peer-deps`

`webtool/frontend/vercel.json` already matches this.

## 3) Deploy frontend first (get public frontend URL)

```bash
cd /Users/yihanwang/Repositories/afterglow/webtool/frontend
vercel --prod
```

Copy the resulting frontend URL, for example:
- `https://vegasafterglow-xxxxx.vercel.app`

Set it in shell for next step:

```bash
export FRONTEND_URL="https://vegasafterglow-xxxxx.vercel.app"
```

## 4) Deploy backend (Cloud Run)

```bash
cd /Users/yihanwang/Repositories/afterglow

PROJECT_ID=project-819bd021-ada6-4176-b14 \
REGION=us-west2 \
SERVICE_NAME=webtool-api \
ARTIFACT_REPO=containers \
MIN_INSTANCES=0 \
ALLOWED_ORIGINS="$FRONTEND_URL" \
bash webtool/scripts/deploy-cloudrun.sh
```

## 5) Confirm backend service URL

Use the Cloud Run **Service URL** (stable), not random aliases.

```bash
gcloud run services describe webtool-api \
  --project project-819bd021-ada6-4176-b14 \
  --region us-west2 \
  --format='value(status.url)'
```

Set it in shell:

```bash
export BACKEND_URL="https://webtool-api-476587253129.us-west2.run.app"
```

## 6) Set Vercel env vars exactly

Run from frontend directory:

```bash
cd /Users/yihanwang/Repositories/afterglow/webtool/frontend
```

If vars already exist, remove then re-add:

```bash
vercel env rm NEXT_PUBLIC_API_URL production --yes || true
vercel env rm BACKEND_URL production --yes || true
```

Add values:

```bash
printf '%s\n' "$BACKEND_URL" | vercel env add NEXT_PUBLIC_API_URL production
printf '%s\n' "$BACKEND_URL" | vercel env add BACKEND_URL production
```

## 7) Redeploy frontend (apply env changes)

```bash
cd /Users/yihanwang/Repositories/afterglow/webtool/frontend
vercel --prod
```

## 8) Verify

Backend health:

```bash
curl "$BACKEND_URL/api/health"
```

Expected:

```json
{"status":"ok","service":"fastapi"}
```

Open frontend URL and verify:
1. Light Curve computes successfully
2. Spectrum computes successfully
3. Sky Image computes successfully

## Future updates

### A) Update backend only

```bash
cd /Users/yihanwang/Repositories/afterglow

PROJECT_ID=project-819bd021-ada6-4176-b14 \
REGION=us-west2 \
SERVICE_NAME=webtool-api \
ARTIFACT_REPO=containers \
MIN_INSTANCES=0 \
ALLOWED_ORIGINS="$FRONTEND_URL" \
bash webtool/scripts/deploy-cloudrun.sh
```

If backend URL did not change, no Vercel env update needed.

### B) Update frontend only

```bash
cd /Users/yihanwang/Repositories/afterglow/webtool/frontend
vercel --prod
```

### C) Update frontend env to a new backend URL

```bash
cd /Users/yihanwang/Repositories/afterglow/webtool/frontend

export BACKEND_URL="https://<new-cloud-run-service-url>"
vercel env rm NEXT_PUBLIC_API_URL production --yes || true
vercel env rm BACKEND_URL production --yes || true
printf '%s\n' "$BACKEND_URL" | vercel env add NEXT_PUBLIC_API_URL production
printf '%s\n' "$BACKEND_URL" | vercel env add BACKEND_URL production
vercel --prod
```

## Troubleshooting used in this repo

1. `COPY webtool/backend/requirements.txt ... file does not exist` in Cloud Build:
`/.gcloudignore` must include:
- `!webtool/backend/requirements.txt`

2. Cloud Run quota error about maxScale / total CPU:
Use the repo defaults in `webtool/scripts/deploy-cloudrun.sh`:
- `--cpu 1`
- `--max-instances 10`

3. Vercel install fails on eslint peer dependency:
- `webtool/frontend/vercel.json` uses `npm install --legacy-peer-deps`
