# ── builder: compile VegasAfterglow C++ extension ──────────────────────────────
ARG REGISTRY_MIRROR=
FROM ${REGISTRY_MIRROR}python:3.12-slim AS builder

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# Optional CN Debian mirror (e.g. mirrors.aliyun.com) for faster apt in mainland China.
ARG APT_MIRROR=
RUN if [ -n "$APT_MIRROR" ]; then \
      sed -i "s|deb.debian.org|${APT_MIRROR}|g" /etc/apt/sources.list.d/debian.sources; \
    fi \
 && apt-get update && apt-get install -y --no-install-recommends \
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

# setuptools-scm needs .git to resolve the version; pass it via build arg instead.
ARG VEGASAFTERGLOW_VERSION=0.0.0
ENV SETUPTOOLS_SCM_PRETEND_VERSION=${VEGASAFTERGLOW_VERSION}

# Limit parallel compilation to avoid OOM on small instances (e.g. 1-2 GB RAM ECS).
# Empty default = unlimited (uses all cores). CN deploy sets this to 1.
ARG CMAKE_BUILD_PARALLEL_LEVEL=
ENV CMAKE_BUILD_PARALLEL_LEVEL=${CMAKE_BUILD_PARALLEL_LEVEL}

# Optional PyPI mirror for CN mainland deploys (e.g. https://pypi.tuna.tsinghua.edu.cn/simple).
ARG PIP_INDEX_URL=

RUN pip install --upgrade pip setuptools wheel scikit-build-core pybind11 numpy setuptools-scm \
 && pip install --no-build-isolation /tmp/vegasafterglow

# Remaining pure-Python deps (fast).
COPY webtool/backend/requirements.txt /tmp/requirements.txt
RUN pip install -r /tmp/requirements.txt

# ── runtime: lean image without build tools ──────────────────────────────────
FROM ${REGISTRY_MIRROR}python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

COPY --from=builder /venv /venv
ENV PATH="/venv/bin:$PATH"

WORKDIR /app
COPY webtool/backend/app /app/app

ENV PORT=8080
EXPOSE 8080

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8080"]
