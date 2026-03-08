#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

PROJECT_ID="${PROJECT_ID:?Set PROJECT_ID}"
REGION="${REGION:-us-central1}"
SERVICE_NAME="${SERVICE_NAME:-webtool-api}"
ARTIFACT_REPO="${ARTIFACT_REPO:-containers}"
ALLOWED_ORIGINS="${ALLOWED_ORIGINS:?Set ALLOWED_ORIGINS (comma-separated, e.g. https://app.vercel.app)}"
IMAGE_TAG="${IMAGE_TAG:-$(git rev-parse --short HEAD 2>/dev/null || date +%Y%m%d%H%M%S)}"
VA_VERSION="$(python -m setuptools_scm 2>/dev/null || git describe --tags --always 2>/dev/null || echo 0.0.0)"
VA_VERSION="${VA_VERSION%.dev*}"
IMAGE_URI="${REGION}-docker.pkg.dev/${PROJECT_ID}/${ARTIFACT_REPO}/${SERVICE_NAME}:${IMAGE_TAG}"
IMAGE_LATEST="${REGION}-docker.pkg.dev/${PROJECT_ID}/${ARTIFACT_REPO}/${SERVICE_NAME}:latest"
MIN_INSTANCES="${MIN_INSTANCES:-0}"
CONCURRENCY="${CONCURRENCY:-8}"

if ! gcloud artifacts repositories describe "$ARTIFACT_REPO" \
  --location "$REGION" \
  --project "$PROJECT_ID" >/dev/null 2>&1; then
  gcloud artifacts repositories create "$ARTIFACT_REPO" \
    --repository-format=docker \
    --location="$REGION" \
    --project "$PROJECT_ID"
fi

gcloud builds submit \
  --project "$PROJECT_ID" \
  --config webtool/backend/cloudrun/cloudbuild.yaml \
  --substitutions "_IMAGE_URI=${IMAGE_URI},_IMAGE_LATEST=${IMAGE_LATEST},_VA_VERSION=${VA_VERSION}" \
  .

gcloud run deploy "$SERVICE_NAME" \
  --project "$PROJECT_ID" \
  --region "$REGION" \
  --image "$IMAGE_URI" \
  --allow-unauthenticated \
  --port 8080 \
  --execution-environment gen2 \
  --cpu 1 \
  --memory 2Gi \
  --concurrency "$CONCURRENCY" \
  --timeout 300 \
  --min-instances "$MIN_INSTANCES" \
  --max-instances 10 \
  --set-env-vars "^|^ALLOWED_ORIGINS=${ALLOWED_ORIGINS}|SERVER_REGION=${REGION}"

SERVICE_URL="$(gcloud run services describe "$SERVICE_NAME" --project "$PROJECT_ID" --region "$REGION" --format='value(status.url)')"
echo "Cloud Run deployed: ${SERVICE_URL}"
