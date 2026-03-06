# webtool (FastAPI + Next.js)

`webtool/` is the FastAPI + Next.js replacement for the Streamlit `webapp/`, with the same three modes:
- Light Curve
- Spectrum
- Sky Image

It is designed for real-time parameter tuning, mobile-phone compatible interaction, and observational data input/upload workflows.

The backend exposes compute + plotting + export APIs, and the frontend renders Plotly figures and downloads CSV/JSON.

## 1) Run backend

```bash
cd webtool/backend
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

`pip install -r requirements.txt` will also build/install local `VegasAfterglow` (`-e ../..`), which can take a few minutes on first run.

Health check:

```bash
curl http://localhost:8000/api/health
```

Useful API endpoints:
- `GET /api/options`
- `GET /api/defaults`
- `POST /api/lightcurve`
- `POST /api/spectrum`
- `POST /api/skymap`

## 2) Run frontend

Open another terminal:

```bash
cd webtool/frontend
cp .env.local.example .env.local
npm install
npm run dev
```

Open `http://localhost:3000`.

## Environment variables

- Backend: `webtool/backend/.env.example`
  - `ALLOWED_ORIGINS`: comma-separated CORS origins (default: `http://localhost:3000,http://127.0.0.1:3000`)
- Frontend: `webtool/frontend/.env.local.example`
  - `NEXT_PUBLIC_API_URL`: browser fallback backend URL (default: `http://127.0.0.1:8000`)
  - `BACKEND_URL`: Next.js server-side rewrite target for `/api-proxy/*` (default: `http://127.0.0.1:8000`)

## Notes

- The backend reuses the same physical model logic as the Streamlit app.
- Sky Image output is interactive Plotly heatmap/animation in the Next.js UI.
- Production deploy configs are in [DEPLOY.md](./DEPLOY.md):
  - Vercel frontend config: `webtool/frontend/vercel.json`
  - Cloud Run script: `webtool/scripts/deploy-cloudrun.sh`
  - Fly.io script/template: `webtool/scripts/deploy-fly.sh`, `webtool/backend/fly.toml.example`
