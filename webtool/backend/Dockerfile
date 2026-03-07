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

# Local VegasAfterglow source so CN deployments do not depend on GitHub access.
COPY pyproject.toml README.md LICENSE CMakeLists.txt /tmp/vegasafterglow/
COPY include /tmp/vegasafterglow/include
COPY src /tmp/vegasafterglow/src
COPY pybind /tmp/vegasafterglow/pybind
COPY external /tmp/vegasafterglow/external
COPY VegasAfterglow /tmp/vegasafterglow/VegasAfterglow

RUN pip install --upgrade pip setuptools wheel scikit-build-core pybind11 numpy setuptools-scm \
 && pip install --no-build-isolation /tmp/vegasafterglow

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
