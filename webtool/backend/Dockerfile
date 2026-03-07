# ── builder: compile VegasAfterglow C++ extension ──────────────────────────────
FROM python:3.12-slim AS builder

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ninja-build \
    git \
 && rm -rf /var/lib/apt/lists/*

RUN python -m venv /venv
ENV PATH="/venv/bin:$PATH"

# VegasAfterglow (slow C++ build — own layer so app changes skip recompilation).
COPY webtool/backend/requirements-base.txt /tmp/requirements-base.txt
RUN pip install --upgrade pip setuptools wheel \
 && pip install -r /tmp/requirements-base.txt

# Remaining pure-Python deps (fast).
COPY webtool/backend/requirements.txt /tmp/requirements.txt
RUN pip install -r /tmp/requirements.txt

# ── runtime: lean image without build tools ──────────────────────────────────
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

COPY --from=builder /venv /venv
ENV PATH="/venv/bin:$PATH"

WORKDIR /app
COPY webtool/backend/app /app/app

ENV PORT=8080
EXPOSE 8080

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8080"]
