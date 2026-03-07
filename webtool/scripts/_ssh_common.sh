#!/usr/bin/env bash

ssh_run() {
  if [[ -n "${SSH_OPTS:-}" && -n "${SSH_OPTS//[[:space:]]/}" ]]; then
    # shellcheck disable=SC2086
    ssh $SSH_OPTS "$@"
  else
    ssh "$@"
  fi
}

scp_run() {
  if [[ -n "${SSH_OPTS:-}" && -n "${SSH_OPTS//[[:space:]]/}" ]]; then
    # shellcheck disable=SC2086
    scp $SSH_OPTS "$@"
  else
    scp "$@"
  fi
}
