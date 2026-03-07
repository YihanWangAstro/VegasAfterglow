#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"
. "${ROOT_DIR}/webtool/scripts/_ssh_common.sh"

SSH_TARGET="${SSH_TARGET:?Set SSH_TARGET (for example: root@8.136.116.255)}"
APP_ROOT="${APP_ROOT:-/opt/afterglow}"
SSH_OPTS="${SSH_OPTS:-}"
BACKEND_IMAGE_NAME="${BACKEND_IMAGE_NAME:-webtool-api}"
FRONTEND_SERVICE="${FRONTEND_SERVICE:-vegasafterglow-frontend}"
NGINX_SITE_NAME="${NGINX_SITE_NAME:-vegasafterglow-cn}"
SERVER_REGION="${SERVER_REGION:-cn-hangzhou}"
KEEP_RELEASES="${KEEP_RELEASES:-3}"
SERVER_NAME="${SERVER_NAME:-_}"
BACKEND_IMAGE_PLATFORM="${BACKEND_IMAGE_PLATFORM:-linux/amd64}"

target_host="${SSH_TARGET##*@}"
PUBLIC_BASE_URL="${PUBLIC_BASE_URL:-http://${target_host}}"
NEXT_PUBLIC_API_URL="${NEXT_PUBLIC_API_URL:-$PUBLIC_BASE_URL}"
NEXT_PUBLIC_API_STATUS_URLS="${NEXT_PUBLIC_API_STATUS_URLS:-$NEXT_PUBLIC_API_URL}"
GLOBAL_FRONTEND_ORIGINS="${GLOBAL_FRONTEND_ORIGINS:-https://www.vegasafterglow.com,https://vegasafterglow.com}"
ALLOWED_ORIGINS="${ALLOWED_ORIGINS:-${PUBLIC_BASE_URL},${GLOBAL_FRONTEND_ORIGINS}}"

git_ref="$(git rev-parse --short HEAD 2>/dev/null || echo local)"
git_dirty=""
if [[ -n "$(git status --porcelain 2>/dev/null || true)" ]]; then
  git_dirty="-dirty"
fi
RELEASE_ID="${RELEASE_ID:-$(date +%Y%m%d%H%M%S)-${git_ref}${git_dirty}}"
RELEASE_DIR="${APP_ROOT}/releases/${RELEASE_ID}"

if ! command -v docker >/dev/null 2>&1; then
  echo "This deploy flow requires a local Docker installation on your Mac." >&2
  exit 1
fi

echo "Preparing remote release: ${RELEASE_DIR}"
ssh_run "$SSH_TARGET" "mkdir -p $(printf '%q' "$RELEASE_DIR")"

tmp_dir="$(mktemp -d)"
archive_path="${tmp_dir}/webtool-${RELEASE_ID}.tar.gz"
remote_archive="/tmp/webtool-${RELEASE_ID}.tar.gz"
trap 'rm -rf "$tmp_dir"' EXIT

backend_image_archive="${tmp_dir}/backend-image-${RELEASE_ID}.tar.gz"
remote_backend_image_archive="/tmp/backend-image-${RELEASE_ID}.tar.gz"

echo "Creating local deploy archive ..."
COPYFILE_DISABLE=1 COPY_EXTENDED_ATTRIBUTES_DISABLE=1 tar \
  --no-mac-metadata \
  --no-xattrs \
  --no-acls \
  --no-fflags \
  --exclude-vcs \
  --exclude='.DS_Store' \
  --exclude='webtool/frontend/node_modules' \
  --exclude='webtool/frontend/.next' \
  --exclude='webtool/backend/.venv' \
  -czf "$archive_path" \
  pyproject.toml \
  README.md \
  LICENSE \
  CMakeLists.txt \
  include \
  src \
  pybind \
  external \
  VegasAfterglow \
  webtool \
  2> >(grep -v 'LIBARCHIVE.xattr\.' >&2)

echo "Building backend image locally for ${BACKEND_IMAGE_PLATFORM} ..."
docker build \
  --platform "$BACKEND_IMAGE_PLATFORM" \
  -f webtool/backend/Dockerfile \
  -t "${BACKEND_IMAGE_NAME}:${RELEASE_ID}" \
  .

echo "Saving backend image archive ..."
docker save "${BACKEND_IMAGE_NAME}:${RELEASE_ID}" | gzip > "$backend_image_archive"

echo "Syncing deploy bundle to ${SSH_TARGET} ..."
scp_run "$archive_path" "${SSH_TARGET}:$(printf '%q' "$remote_archive")"
ssh_run "$SSH_TARGET" \
  "tar -xzf $(printf '%q' "$remote_archive") -C $(printf '%q' "$RELEASE_DIR") && rm -f $(printf '%q' "$remote_archive")"

echo "Syncing backend image to ${SSH_TARGET} ..."
scp_run "$backend_image_archive" "${SSH_TARGET}:$(printf '%q' "$remote_backend_image_archive")"
ssh_run "$SSH_TARGET" \
  "gzip -dc $(printf '%q' "$remote_backend_image_archive") | docker load && rm -f $(printf '%q' "$remote_backend_image_archive")"

remote_script="${RELEASE_DIR}/webtool/scripts/deploy-cn-ecs-remote.sh"

echo "Running remote deploy steps ..."
ssh_run "$SSH_TARGET" \
  "APP_ROOT=$(printf '%q' "$APP_ROOT") \
RELEASE_ID=$(printf '%q' "$RELEASE_ID") \
RELEASE_DIR=$(printf '%q' "$RELEASE_DIR") \
PUBLIC_BASE_URL=$(printf '%q' "$PUBLIC_BASE_URL") \
NEXT_PUBLIC_API_URL=$(printf '%q' "$NEXT_PUBLIC_API_URL") \
NEXT_PUBLIC_API_STATUS_URLS=$(printf '%q' "$NEXT_PUBLIC_API_STATUS_URLS") \
SERVER_NAME=$(printf '%q' "$SERVER_NAME") \
ALLOWED_ORIGINS=$(printf '%q' "$ALLOWED_ORIGINS") \
GLOBAL_FRONTEND_ORIGINS=$(printf '%q' "$GLOBAL_FRONTEND_ORIGINS") \
BACKEND_IMAGE_NAME=$(printf '%q' "$BACKEND_IMAGE_NAME") \
FRONTEND_SERVICE=$(printf '%q' "$FRONTEND_SERVICE") \
NGINX_SITE_NAME=$(printf '%q' "$NGINX_SITE_NAME") \
SERVER_REGION=$(printf '%q' "$SERVER_REGION") \
KEEP_RELEASES=$(printf '%q' "$KEEP_RELEASES") \
bash $(printf '%q' "$remote_script")"

echo "Local sync complete."
echo "Remote release: ${RELEASE_DIR}"
