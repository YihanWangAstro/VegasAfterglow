#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
. "${ROOT_DIR}/webtool/scripts/_ssh_common.sh"

SSH_TARGET="${SSH_TARGET:?Set SSH_TARGET (for example: root@8.136.116.255)}"
APP_ROOT="${APP_ROOT:-/opt/afterglow}"
SWAP_SIZE_GB="${SWAP_SIZE_GB:-4}"
SSH_OPTS="${SSH_OPTS:-}"

ssh_run "$SSH_TARGET" \
  "APP_ROOT=$(printf '%q' "$APP_ROOT") SWAP_SIZE_GB=$(printf '%q' "$SWAP_SIZE_GB") bash -s" <<'EOF'
set -euo pipefail

if [[ $(id -u) -ne 0 ]]; then
  echo "Run bootstrap as root or use SSH_TARGET=root@host." >&2
  exit 1
fi

install_node_20() {
  if command -v node >/dev/null 2>&1; then
    local current_major
    current_major="$(node -v | sed -E 's/^v([0-9]+).*/\1/')"
    if [[ "$current_major" == "20" ]]; then
      return
    fi
  fi

  if command -v dnf >/dev/null 2>&1; then
    curl -fsSL https://rpm.nodesource.com/setup_20.x | bash -
    dnf install -y nodejs
    return
  fi

  if command -v apt-get >/dev/null 2>&1; then
    apt-get update
    apt-get install -y ca-certificates curl gnupg
    install -m 0755 -d /etc/apt/keyrings
    curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg
    chmod a+r /etc/apt/keyrings/nodesource.gpg
    . /etc/os-release
    echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_20.x nodistro main" \
      >/etc/apt/sources.list.d/nodesource.list
    apt-get update
    apt-get install -y nodejs
    return
  fi

  echo "Unsupported package manager for Node.js install." >&2
  exit 1
}

install_docker_dnf() {
  if command -v docker >/dev/null 2>&1; then
    return
  fi
  rm -f /etc/yum.repos.d/docker*.repo
  dnf -y remove docker-ce containerd.io docker-ce-rootless-extras docker-buildx-plugin docker-ce-cli docker-compose-plugin || true
  wget -O /etc/yum.repos.d/docker-ce.repo http://mirrors.cloud.aliyuncs.com/docker-ce/linux/centos/docker-ce.repo
  sed -i 's|https://mirrors.aliyun.com|http://mirrors.cloud.aliyuncs.com|g' /etc/yum.repos.d/docker-ce.repo
  dnf -y install dnf-plugin-releasever-adapter --repo alinux3-plus || true
  dnf -y install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
}

install_docker_apt() {
  if command -v docker >/dev/null 2>&1; then
    return
  fi
  apt-get update
  apt-get install -y curl wget nginx tar docker.io
}

if command -v dnf >/dev/null 2>&1; then
  dnf -y install curl wget nginx tar
  install_docker_dnf
elif command -v apt-get >/dev/null 2>&1; then
  install_docker_apt
else
  echo "Unsupported package manager. Expected dnf or apt-get." >&2
  exit 1
fi

install_node_20

systemctl enable docker
systemctl restart docker
systemctl enable --now nginx

mkdir -p "$APP_ROOT/releases"

if [[ "${SWAP_SIZE_GB}" != "0" ]] && ! swapon --show | grep -q '^/swapfile '; then
  fallocate -l "${SWAP_SIZE_GB}G" /swapfile || dd if=/dev/zero of=/swapfile bs=1M count=$((SWAP_SIZE_GB * 1024))
  chmod 600 /swapfile
  mkswap /swapfile
  swapon /swapfile
  grep -q '^/swapfile ' /etc/fstab || echo '/swapfile none swap sw 0 0' >> /etc/fstab
fi

echo "Bootstrap complete."
docker --version
node -v
npm -v
free -h
swapon --show || true
EOF
