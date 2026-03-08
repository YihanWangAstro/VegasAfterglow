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
API_SERVER_NAME="${API_SERVER_NAME:-}"
BACKEND_IMAGE_PLATFORM="${BACKEND_IMAGE_PLATFORM:-linux/amd64}"
REGISTRY_MIRROR="${REGISTRY_MIRROR:-}"
PIP_INDEX_URL="${PIP_INDEX_URL:-https://pypi.tuna.tsinghua.edu.cn/simple}"
APT_MIRROR="${APT_MIRROR:-mirrors.aliyun.com}"
CMAKE_BUILD_PARALLEL_LEVEL="${CMAKE_BUILD_PARALLEL_LEVEL:-1}"

target_host="${SSH_TARGET##*@}"
PUBLIC_BASE_URL="${PUBLIC_BASE_URL:-http://${target_host}}"
public_scheme="${PUBLIC_BASE_URL%%://*}"
if [[ "$public_scheme" == "$PUBLIC_BASE_URL" ]]; then
  public_scheme="http"
fi
if [[ -z "${NEXT_PUBLIC_API_URL:-}" && -n "$API_SERVER_NAME" ]]; then
  NEXT_PUBLIC_API_URL="${public_scheme}://${API_SERVER_NAME}"
else
  NEXT_PUBLIC_API_URL="${NEXT_PUBLIC_API_URL:-$PUBLIC_BASE_URL}"
fi
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

# Sync external/ (boost, xtensor headers — ~36k files, rarely changes) to a shared
# cache so repeat deploys skip re-uploading unchanged files.
EXTERNAL_CACHE="${APP_ROOT}/cache/external"
echo "Syncing external/ to shared cache ..."
ssh_run "$SSH_TARGET" "mkdir -p $(printf '%q' "$EXTERNAL_CACHE")"
rsync_run \
  --delete \
  "$ROOT_DIR/external/" \
  "${SSH_TARGET}:${EXTERNAL_CACHE}/"

echo "Syncing source to ${SSH_TARGET}:${RELEASE_DIR} ..."
# Only upload paths needed by the Dockerfile and frontend build.
rsync_run \
  --delete \
  --exclude='node_modules/' \
  --exclude='.next/' \
  --exclude='.venv/' \
  --exclude='__pycache__/' \
  --exclude='*.pyc' \
  --exclude='.DS_Store' \
  --include='.dockerignore' \
  --include='pyproject.toml' \
  --include='README.md' \
  --include='LICENSE' \
  --include='CMakeLists.txt' \
  --include='include/***' \
  --include='src/***' \
  --include='pybind/***' \
  --include='VegasAfterglow/***' \
  --include='webtool/' \
  --include='webtool/backend/***' \
  --include='webtool/frontend/***' \
  --include='webtool/deploy/***' \
  --include='webtool/scripts/***' \
  --exclude='*' \
  "$ROOT_DIR/" \
  "${SSH_TARGET}:${RELEASE_DIR}/"

# Hard-link external/ from cache into this release (instant, no extra disk space).
ssh_run "$SSH_TARGET" "cp -al $(printf '%q' "$EXTERNAL_CACHE") $(printf '%q' "$RELEASE_DIR/external")"

VA_VERSION="$(python -m setuptools_scm 2>/dev/null || git describe --tags --always 2>/dev/null || echo 0.0.0)"
VA_VERSION="${VA_VERSION%.dev*}"

echo "Building backend image on ${SSH_TARGET} (VegasAfterglow ${VA_VERSION}) ..."
# REGISTRY_MIRROR prefixes Docker Hub image refs with a CN mirror (e.g. docker.m.daocloud.io/).
# PIP_INDEX_URL points pip at a CN PyPI mirror. Docker layer cache means only changed layers rebuild.
build_cmd="docker build"
build_cmd+=" --build-arg VEGASAFTERGLOW_VERSION='${VA_VERSION}'"
[[ -n "$REGISTRY_MIRROR" ]] && build_cmd+=" --build-arg REGISTRY_MIRROR='${REGISTRY_MIRROR}'"
[[ -n "$PIP_INDEX_URL" ]]   && build_cmd+=" --build-arg PIP_INDEX_URL='${PIP_INDEX_URL}'"
[[ -n "$APT_MIRROR" ]]      && build_cmd+=" --build-arg APT_MIRROR='${APT_MIRROR}'"
[[ -n "$CMAKE_BUILD_PARALLEL_LEVEL" ]] && build_cmd+=" --build-arg CMAKE_BUILD_PARALLEL_LEVEL='${CMAKE_BUILD_PARALLEL_LEVEL}'"
build_cmd+=" --platform '${BACKEND_IMAGE_PLATFORM}'"
build_cmd+=" -f '${RELEASE_DIR}/webtool/backend/Dockerfile'"
build_cmd+=" -t '${BACKEND_IMAGE_NAME}:${RELEASE_ID}'"
build_cmd+=" '${RELEASE_DIR}'"

ssh_run "$SSH_TARGET" "$build_cmd"
echo "Docker build succeeded."

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
API_SERVER_NAME=$(printf '%q' "$API_SERVER_NAME") \
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
