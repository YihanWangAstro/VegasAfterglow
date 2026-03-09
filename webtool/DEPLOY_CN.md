# ECS Deployment via Mac SSH

This flow is the closest ECS equivalent to the existing Cloud Run deploy:

- local machine sends the current repo state to the server
- backend image is built locally on your Mac and uploaded over SSH
- ECS builds the Next.js frontend locally
- ECS restarts Docker, systemd, and Nginx in place

This is the recommended mainland deployment path.

The deploy script uploads the **current local deploy bundle**, so local edits in
the shipped paths are included, similar to `gcloud builds submit .`.

## 1) One-time bootstrap

Prerequisites:

- You can SSH from your Mac to the ECS host as `root`
- The security group allows `22/TCP` from your current public IP
- Docker is installed and running on your Mac

Run from the repo root on your Mac:

```bash
cd /Users/yihanwang/Repositories/afterglow

SSH_TARGET=root@8.136.116.255 \
bash webtool/scripts/bootstrap-cn-ecs-via-ssh.sh
```

This installs:

- Docker
- Node.js 20
- Nginx
- curl / tar
- optional swap (default: 4 GiB)

If you already configured the server manually, you can skip this step.

## 2) Deploy to ECS over SSH

### 2.1 Deploy with the server IP first

This is the simplest first deploy and works before DNS / ICP / TLS:

```bash
cd /Users/yihanwang/Repositories/afterglow

SSH_TARGET=root@8.136.116.255 \
PUBLIC_BASE_URL=http://8.136.116.255 \
SERVER_NAME=_ \
bash webtool/scripts/deploy-cn-ecs-via-ssh.sh
```

If Docker Hub is blocked (common on CN mainland), add a registry mirror:

```bash
REGISTRY_MIRROR=docker.m.daocloud.io/ \
SSH_TARGET=root@8.136.116.255 \
PUBLIC_BASE_URL=http://8.136.116.255 \
SERVER_NAME=_ \
bash webtool/scripts/deploy-cn-ecs-via-ssh.sh
```

What it does:

1. copies the current local deploy bundle to `/opt/afterglow/releases/<release-id>`
2. updates `/opt/afterglow/current`
3. builds the backend image locally on your Mac
4. runs the backend container on `127.0.0.1:8080`
5. runs `npm ci` + `npm run build` for `webtool/frontend`
6. installs / updates the frontend systemd service
7. installs / updates an Nginx config that serves:
   - `/` -> Next.js on `127.0.0.1:3000`
   - `/api/*` -> FastAPI on `127.0.0.1:8080`

After the script succeeds, open:

```text
http://8.136.116.255
```

## 3) Repeat deploys

After you change code locally, rerun the same command:

```bash
cd /Users/yihanwang/Repositories/afterglow

REGISTRY_MIRROR=docker.m.daocloud.io/ \
SSH_TARGET=root@8.136.116.255 \
PUBLIC_BASE_URL=http://8.136.116.255 \
SERVER_NAME=_ \
bash webtool/scripts/deploy-cn-ecs-via-ssh.sh
```

Each deploy creates a new release directory and keeps the latest 3 by default.

## 4) Important environment knobs

Defaults are chosen so the IP-based deploy works immediately.

- `SSH_TARGET`
  - required
  - example: `root@8.136.116.255`
- `PUBLIC_BASE_URL`
  - default: `http://<ssh-host>`
  - example: `http://8.136.116.255`
- `SERVER_NAME`
  - default: `_`
  - set a real domain later, for example: `cn.vegasafterglow.com`
- `API_SERVER_NAME`
  - default: empty
  - optional separate API hostname, for example: `api.vegasafterglow.cn`
  - if set, the deploy writes a split Nginx config for frontend and API hosts
- `ALLOWED_ORIGINS`
  - default: `<PUBLIC_BASE_URL>,https://www.vegasafterglow.com,https://vegasafterglow.com`
- `REGISTRY_MIRROR`
  - default: empty (Docker Hub)
  - set to a CN Docker mirror if Docker Hub is blocked, for example: `docker.m.daocloud.io/`
  - **must include the trailing slash**
- `PIP_INDEX_URL`
  - default: `https://pypi.tuna.tsinghua.edu.cn/simple` (Tsinghua mirror)
  - set to empty string to use upstream PyPI
- `KEEP_RELEASES`
  - default: `3`
- `RELEASE_ID`
  - default: timestamp + local git short SHA

## 5) Later: switch from IP to domain

After DNS and ICP are ready, redeploy with your frontend domain:

```bash
cd /Users/yihanwang/Repositories/afterglow

REGISTRY_MIRROR=docker.m.daocloud.io/ \
SSH_TARGET=root@8.136.116.255 \
PUBLIC_BASE_URL=https://www.vegasafterglow.cn \
SERVER_NAME=www.vegasafterglow.cn \
API_SERVER_NAME=api.vegasafterglow.cn \
ALLOWED_ORIGINS=https://www.vegasafterglow.cn,https://www.vegasafterglow.com,https://vegasafterglow.com \
bash webtool/scripts/deploy-cn-ecs-via-ssh.sh
```

Then issue TLS certificates on the server:

```bash
ssh root@8.136.116.255
certbot --nginx -d www.vegasafterglow.cn -d api.vegasafterglow.cn
nginx -t && systemctl reload nginx
```

## 6) Verification

On the server:

```bash
docker ps
curl -fsS http://127.0.0.1:8080/api/health
systemctl status vegasafterglow-frontend --no-pager
nginx -t
```

From your Mac:

```bash
curl -I http://8.136.116.255
curl -fsS http://8.136.116.255/api/health
```

After the domain cutover:

```bash
curl -I https://www.vegasafterglow.cn
curl -fsS https://api.vegasafterglow.cn/api/health
```
