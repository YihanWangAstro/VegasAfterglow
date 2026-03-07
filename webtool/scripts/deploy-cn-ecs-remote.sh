#!/usr/bin/env bash
set -euo pipefail

APP_ROOT="${APP_ROOT:-/opt/afterglow}"
RELEASE_ID="${RELEASE_ID:?Set RELEASE_ID}"
RELEASE_DIR="${RELEASE_DIR:-${APP_ROOT}/releases/${RELEASE_ID}}"
CURRENT_DIR="${CURRENT_DIR:-${APP_ROOT}/current}"

BACKEND_IMAGE_NAME="${BACKEND_IMAGE_NAME:-webtool-api}"
CONTAINER_NAME="${CONTAINER_NAME:-webtool-api}"
BACKEND_HOST_PORT="${BACKEND_HOST_PORT:-8080}"
BACKEND_CONTAINER_PORT="${BACKEND_CONTAINER_PORT:-8080}"
SERVER_REGION="${SERVER_REGION:-cn-hangzhou}"

PUBLIC_BASE_URL="${PUBLIC_BASE_URL:?Set PUBLIC_BASE_URL (for example: http://8.136.116.255)}"
NEXT_PUBLIC_API_URL="${NEXT_PUBLIC_API_URL:-$PUBLIC_BASE_URL}"
NEXT_PUBLIC_API_STATUS_URLS="${NEXT_PUBLIC_API_STATUS_URLS:-$NEXT_PUBLIC_API_URL}"
SERVER_NAME="${SERVER_NAME:-_}"
GLOBAL_FRONTEND_ORIGINS="${GLOBAL_FRONTEND_ORIGINS:-https://www.vegasafterglow.com,https://vegasafterglow.com}"
ALLOWED_ORIGINS="${ALLOWED_ORIGINS:-${PUBLIC_BASE_URL},${GLOBAL_FRONTEND_ORIGINS}}"

FRONTEND_SERVICE="${FRONTEND_SERVICE:-vegasafterglow-frontend}"
FRONTEND_PORT="${FRONTEND_PORT:-3000}"
NGINX_SITE_NAME="${NGINX_SITE_NAME:-vegasafterglow-cn}"
KEEP_RELEASES="${KEEP_RELEASES:-3}"

render_server_name="$SERVER_NAME"
if [[ "$render_server_name" == "_" ]]; then
  render_server_name="${PUBLIC_BASE_URL#*://}"
  render_server_name="${render_server_name%%/*}"
fi

if [[ $(id -u) -ne 0 ]]; then
  echo "Run deploy as root or use SSH_TARGET=root@host." >&2
  exit 1
fi

for required in docker npm systemctl nginx curl sed; do
  if ! command -v "$required" >/dev/null 2>&1; then
    echo "Missing required command on server: $required" >&2
    exit 1
  fi
done

mkdir -p "${APP_ROOT}/releases"
ln -sfnT "$RELEASE_DIR" "$CURRENT_DIR"

echo "Using uploaded backend image: ${BACKEND_IMAGE_NAME}:${RELEASE_ID}"
if ! docker image inspect "${BACKEND_IMAGE_NAME}:${RELEASE_ID}" >/dev/null 2>&1; then
  echo "Expected uploaded backend image not found: ${BACKEND_IMAGE_NAME}:${RELEASE_ID}" >&2
  exit 1
fi
docker tag "${BACKEND_IMAGE_NAME}:${RELEASE_ID}" "${BACKEND_IMAGE_NAME}:current"

if docker container inspect "$CONTAINER_NAME" >/dev/null 2>&1; then
  echo "Stopping existing backend container: ${CONTAINER_NAME}"
  docker rm -f "$CONTAINER_NAME"
fi

echo "Starting backend container: ${CONTAINER_NAME}"
docker run -d \
  --name "$CONTAINER_NAME" \
  --restart unless-stopped \
  -p "127.0.0.1:${BACKEND_HOST_PORT}:${BACKEND_CONTAINER_PORT}" \
  -e ALLOWED_ORIGINS="$ALLOWED_ORIGINS" \
  -e SERVER_REGION="$SERVER_REGION" \
  "${BACKEND_IMAGE_NAME}:current"

echo "Waiting for backend health check"
backend_health_url="http://127.0.0.1:${BACKEND_HOST_PORT}/api/health"
for attempt in $(seq 1 20); do
  if curl -fsS "$backend_health_url" >/dev/null; then
    break
  fi

  if (( attempt == 20 )); then
    echo "Backend health check failed after ${attempt} attempts." >&2
    docker logs --tail 200 "$CONTAINER_NAME" >&2 || true
    exit 1
  fi

  sleep 1
done

cd "${CURRENT_DIR}/webtool/frontend"
echo "Installing frontend dependencies"
npm ci --legacy-peer-deps

echo "Building frontend bundle"
NEXT_PUBLIC_API_URL="$NEXT_PUBLIC_API_URL" \
NEXT_PUBLIC_API_STATUS_URLS="$NEXT_PUBLIC_API_STATUS_URLS" \
npm run build

echo "Restarting frontend service: ${FRONTEND_SERVICE}"
sed \
  -e "s#__FRONTEND_DIR__#${CURRENT_DIR}/webtool/frontend#g" \
  "${CURRENT_DIR}/webtool/deploy/alicloud/systemd/vegasafterglow-frontend.service" \
  >"/etc/systemd/system/${FRONTEND_SERVICE}.service"

systemctl daemon-reload
systemctl enable "$FRONTEND_SERVICE" >/dev/null 2>&1 || true
systemctl restart "$FRONTEND_SERVICE"

echo "Waiting for frontend service"
frontend_local_url="http://127.0.0.1:${FRONTEND_PORT}/"
for attempt in $(seq 1 20); do
  if curl -fsSI "$frontend_local_url" >/dev/null; then
    break
  fi

  if (( attempt == 20 )); then
    echo "Frontend service failed to become ready after ${attempt} attempts." >&2
    systemctl status "$FRONTEND_SERVICE" --no-pager >&2 || true
    journalctl -u "$FRONTEND_SERVICE" -n 120 --no-pager >&2 || true
    exit 1
  fi

  sleep 1
done

echo "Reloading Nginx"
sed \
  -e "s/__PUBLIC_SERVER_NAME__/${render_server_name}/g" \
  -e "s/__FRONTEND_PORT__/${FRONTEND_PORT}/g" \
  -e "s/__BACKEND_PORT__/${BACKEND_HOST_PORT}/g" \
  "${CURRENT_DIR}/webtool/deploy/alicloud/nginx/vegasafterglow-cn-single-host.conf.template" \
  >"/etc/nginx/conf.d/${NGINX_SITE_NAME}.conf"

nginx -t
systemctl reload nginx

echo "Waiting for Nginx frontend route"
nginx_check_url="http://127.0.0.1/"
if [[ "$SERVER_NAME" == "_" ]]; then
  nginx_check_header=()
else
  nginx_check_header=(-H "Host: ${SERVER_NAME}")
fi

for attempt in $(seq 1 20); do
  if curl -fsSI "${nginx_check_header[@]}" "$nginx_check_url" >/dev/null; then
    break
  fi

  if (( attempt == 20 )); then
    echo "Nginx frontend route failed after ${attempt} attempts." >&2
    systemctl status nginx --no-pager >&2 || true
    exit 1
  fi

  sleep 1
done

if [[ "$KEEP_RELEASES" =~ ^[0-9]+$ ]] && (( KEEP_RELEASES > 0 )); then
  mapfile -t release_dirs < <(find "${APP_ROOT}/releases" -mindepth 1 -maxdepth 1 -type d | sort)
  if ((${#release_dirs[@]} > KEEP_RELEASES)); then
    delete_count=$((${#release_dirs[@]} - KEEP_RELEASES))
    for old_dir in "${release_dirs[@]:0:${delete_count}}"; do
      if [[ "$old_dir" != "$RELEASE_DIR" ]]; then
        rm -rf "$old_dir"
      fi
    done
  fi
fi

echo "ECS deploy complete."
echo "Release: ${RELEASE_ID}"
echo "App URL: ${PUBLIC_BASE_URL}"
echo "Backend health: http://127.0.0.1:${BACKEND_HOST_PORT}/api/health"
