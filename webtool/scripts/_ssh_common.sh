#!/usr/bin/env bash

# Keep SSH alive during long operations (C++ builds can take 5-10 min with no output).
_SSH_KEEPALIVE="-o ServerAliveInterval=30 -o ServerAliveCountMax=10"

ssh_run() {
  # shellcheck disable=SC2086
  ssh $_SSH_KEEPALIVE ${SSH_OPTS:-} "$@"
}

scp_run() {
  # shellcheck disable=SC2086
  scp $_SSH_KEEPALIVE ${SSH_OPTS:-} "$@"
}

rsync_run() {
  # shellcheck disable=SC2086
  rsync -az --progress -e "ssh $_SSH_KEEPALIVE ${SSH_OPTS:-}" "$@"
}
