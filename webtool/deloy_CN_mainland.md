# China Mainland Deployment Plan (Alibaba Cloud)

This plan adds a China mainland mirror without breaking the current global stack.

- Current global stack stays as-is:
  - frontend: `www.vegasafterglow.com` (Vercel)
  - API: `api.vegasafterglow.com` (GCP Global LB + Cloud Run)
- New mainland entry:
  - frontend: `cn.vegasafterglow.com`
  - API: `api-cn.vegasafterglow.com`

## 1) Scope and rollout strategy

1. Keep current global production unchanged.
2. Deploy mainland mirror on Alibaba Cloud first.
3. Validate with internal tests.
4. Publicly announce `cn.vegasafterglow.com` as the mainland-optimized entry.
5. Optional later: geo-routing at DNS/CDN layer.

## 2) Compliance prerequisites (mandatory for mainland hosting)

1. ICP filing for the domain used in mainland hosting.
2. Public security filing (if required by your provider/policy).
3. Domain control at DNS provider (currently Namecheap).

If compliance is not completed, use Hong Kong/overseas only and do not cut traffic to mainland infra.

## 3) Recommended architecture (Phase 1)

Single ECS host in mainland (simple and fast to launch):

- Region: `cn-hangzhou` (primary)
- Nginx reverse proxy + TLS
- Next.js frontend process on `127.0.0.1:3000`
- FastAPI backend container on `127.0.0.1:8080`

Optional Phase 2 (higher availability):

- Add second region ECS (for example `cn-beijing`)
- Add Alibaba SLB/ALB in front of both nodes
- Add active health checks and failover

## 4) Step-by-step deployment

## 4.1 Prepare Alibaba resources

1. Create ECS (Ubuntu 22.04+), public IP, and security group.
2. Open ports: `22`, `80`, `443`.
3. Create ACR namespace and registry endpoint (for example in `cn-hangzhou`).

## 4.2 Build and push backend image (from your existing v2.0.1 tag pin)

Run on local machine:

```bash
cd /Users/yihanwang/Repositories/afterglow

export CN_REGISTRY="registry.cn-hangzhou.aliyuncs.com/<namespace>"
export BACKEND_IMAGE="$CN_REGISTRY/webtool-api:v2.0.1"

docker login "$CN_REGISTRY"
docker build -f webtool/backend/Dockerfile -t "$BACKEND_IMAGE" .
docker push "$BACKEND_IMAGE"
```

Note: `webtool/backend/requirements.txt` is already pinned to:

- `VegasAfterglow @ git+https://github.com/YihanWangAstro/VegasAfterglow.git@v2.0.1`

## 4.3 Configure ECS host base software

```bash
ssh root@<ECS_PUBLIC_IP>
apt update
apt install -y docker.io nginx certbot python3-certbot-nginx git curl
systemctl enable --now docker
systemctl enable --now nginx
```

## 4.4 Run backend container

```bash
docker pull "$BACKEND_IMAGE"
docker rm -f webtool-api || true

docker run -d \
  --name webtool-api \
  --restart unless-stopped \
  -p 127.0.0.1:8080:8080 \
  -e ALLOWED_ORIGINS="https://cn.vegasafterglow.com,https://www.vegasafterglow.com,https://vegasafterglow.com" \
  -e SERVER_REGION="cn-mainland" \
  "$BACKEND_IMAGE"

curl -s http://127.0.0.1:8080/api/health
```

## 4.5 Deploy frontend process (Next.js)

```bash
mkdir -p /opt/afterglow
cd /opt/afterglow
git clone https://github.com/YihanWangAstro/VegasAfterglow.git repo || true
cd repo
git fetch --tags
git checkout main
git pull

cd webtool/frontend
npm ci --legacy-peer-deps

NEXT_PUBLIC_API_URL="https://api-cn.vegasafterglow.com" \
NEXT_PUBLIC_API_STATUS_URLS="https://api-cn.vegasafterglow.com" \
npm run build
```

Create systemd service:

```bash
cat >/etc/systemd/system/vegasafterglow-frontend.service <<'EOF'
[Unit]
Description=VegasAfterglow CN Frontend (Next.js)
After=network.target

[Service]
Type=simple
WorkingDirectory=/opt/afterglow/repo/webtool/frontend
Environment=NODE_ENV=production
ExecStart=/usr/bin/npm run start -- --hostname 127.0.0.1 --port 3000
Restart=always
RestartSec=3
User=root

[Install]
WantedBy=multi-user.target
EOF

systemctl daemon-reload
systemctl enable --now vegasafterglow-frontend
systemctl status vegasafterglow-frontend --no-pager
```

## 4.6 Configure Nginx reverse proxy

Create `/etc/nginx/conf.d/vegasafterglow-cn.conf`:

```nginx
map $http_upgrade $connection_upgrade {
  default upgrade;
  ''      close;
}

server {
  listen 80;
  server_name cn.vegasafterglow.com;

  location / {
    proxy_pass http://127.0.0.1:3000;
    proxy_http_version 1.1;
    proxy_set_header Host $host;
    proxy_set_header X-Real-IP $remote_addr;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
    proxy_set_header X-Forwarded-Proto $scheme;
    proxy_set_header Upgrade $http_upgrade;
    proxy_set_header Connection $connection_upgrade;
  }
}

server {
  listen 80;
  server_name api-cn.vegasafterglow.com;

  location / {
    proxy_pass http://127.0.0.1:8080;
    proxy_http_version 1.1;
    proxy_set_header Host $host;
    proxy_set_header X-Real-IP $remote_addr;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
    proxy_set_header X-Forwarded-Proto $scheme;
  }
}
```

Then:

```bash
nginx -t
systemctl reload nginx
```

## 4.7 DNS records

At Namecheap, add:

1. `A` record: `cn` -> `<ECS_PUBLIC_IP>`
2. `A` record: `api-cn` -> `<ECS_PUBLIC_IP>`

Verify:

```bash
dig +short cn.vegasafterglow.com A @8.8.8.8
dig +short api-cn.vegasafterglow.com A @8.8.8.8
```

## 4.8 TLS certificates

After DNS resolves:

```bash
certbot --nginx -d cn.vegasafterglow.com -d api-cn.vegasafterglow.com
nginx -t && systemctl reload nginx
```

## 5) Validation checklist

1. Backend health:

```bash
curl -i https://api-cn.vegasafterglow.com/api/health
```

2. Frontend:

```bash
curl -I https://cn.vegasafterglow.com
```

3. Browser test:
- open `https://cn.vegasafterglow.com`
- light curve slider should update
- no CORS errors in console
- server status should show mainland API endpoint

## 6) Cutover and rollback

Cutover:

1. Keep global entry unchanged.
2. Add a visible mainland entry link on main site:
   - `Mainland optimized entry: https://cn.vegasafterglow.com`

Rollback:

1. If CN mirror fails, keep global site running.
2. Set `cn` and `api-cn` DNS to maintenance page or disable access.
3. Revert backend image tag:

```bash
docker rm -f webtool-api
docker run ... <previous_image_tag>
```

## 7) Future updates

## 7.1 Backend update

1. Build and push new image tag to ACR.
2. Pull and restart container on ECS.
3. Verify `/api/health` and one compute endpoint.

## 7.2 Frontend update

```bash
cd /opt/afterglow/repo
git pull
cd webtool/frontend
npm ci --legacy-peer-deps
NEXT_PUBLIC_API_URL="https://api-cn.vegasafterglow.com" \
NEXT_PUBLIC_API_STATUS_URLS="https://api-cn.vegasafterglow.com" \
npm run build
systemctl restart vegasafterglow-frontend
```

## 8) Recommended next hardening tasks

1. Add daily backup of Nginx config + systemd files.
2. Add basic monitoring:
   - process up/down
   - `/api/health` probe
   - disk usage and memory alerts
3. Add second ECS node + SLB/ALB for failover.
