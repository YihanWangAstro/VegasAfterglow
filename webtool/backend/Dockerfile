FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ninja-build \
    git \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# API + VegasAfterglow dependencies (VegasAfterglow pinned to GitHub tag in requirements).
COPY webtool/backend/requirements.txt /tmp/requirements.txt
RUN python -m pip install --upgrade pip setuptools wheel \
 && python -m pip install -r /tmp/requirements.txt

# FastAPI service.
COPY webtool/backend/app /app/app

ENV PORT=8080
EXPOSE 8080

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8080"]
