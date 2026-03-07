# China Mainland Deployment

The mainland deployment path is now standardized on the SSH workflow in
[DEPLOY_ECS_SSH.md](./DEPLOY_ECS_SSH.md).

Use that document for:

- one-time ECS bootstrap
- repeat deploys from your Mac over SSH
- local backend image builds on your Mac
- frontend, backend, systemd, and Nginx updates in one command

## Current recommended topology

- one ECS host in mainland China
- backend container on `127.0.0.1:8080`
- frontend service on `127.0.0.1:3000`
- Nginx on `80/443`

## Compliance prerequisites

- mainland ECS with public bandwidth enabled
- ICP filing before public mainland domain cutover
- domain ownership and DNS control

If you are still in internal testing, deploy with the ECS public IP first and do
not switch public DNS until ICP and TLS are ready.

## Domain strategy

Suggested rollout:

1. Deploy and validate on the ECS public IP first.
2. Reserve the production mainland hostnames:
   - `www.vegasafterglow.cn`
   - `api.vegasafterglow.cn`
3. After ICP approval, issue certificates with `certbot --nginx`.
4. Redeploy with:
   - `PUBLIC_BASE_URL=https://www.vegasafterglow.cn`
   - `SERVER_NAME=www.vegasafterglow.cn`
   - `API_SERVER_NAME=api.vegasafterglow.cn`

## Pre-filing prep already done

The SSH deploy flow now supports both:

- current IP-based single-host mode
- later split-host mode for `www.vegasafterglow.cn` + `api.vegasafterglow.cn`

So you do not need a second deployment path after ICP approval. Re-run the same
deploy command with the domain variables when DNS and TLS are ready.

## Notes

- The old split deployment path was removed to keep a single mainland deploy workflow.
