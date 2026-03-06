#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

APP_NAME="${APP_NAME:?Set APP_NAME (e.g. webtool-api-prod)}"
PRIMARY_REGION="${PRIMARY_REGION:-ord}"
ALLOWED_ORIGINS="${ALLOWED_ORIGINS:?Set ALLOWED_ORIGINS (comma-separated)}"

sed \
  -e "s/__APP_NAME__/${APP_NAME}/g" \
  -e "s/__PRIMARY_REGION__/${PRIMARY_REGION}/g" \
  webtool/backend/fly.toml.example > webtool/backend/fly.toml

flyctl apps create "$APP_NAME" --machines >/dev/null 2>&1 || true
flyctl secrets set ALLOWED_ORIGINS="$ALLOWED_ORIGINS" --app "$APP_NAME"
flyctl deploy --config webtool/backend/fly.toml --app "$APP_NAME"

APP_URL="https://${APP_NAME}.fly.dev"
echo "Fly deployed: ${APP_URL}"
