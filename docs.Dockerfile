FROM python:3.10-slim AS build

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV POETRY_VERSION=1.8.3

SHELL [ "/bin/bash", "-euxo", "pipefail", "-c" ]

RUN apt-get update && \
    apt-get install --no-install-recommends -y gcc libc-dev && \
    pip install --no-cache-dir "poetry==$POETRY_VERSION" && \
    poetry config virtualenvs.create false && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

COPY pyproject.toml poetry.lock ./

RUN poetry export --without-hashes --with dev --format requirements.txt --output requirements.txt

WORKDIR /build

RUN python -m venv .venv && \
    ./.venv/bin/pip install --no-cache-dir -r /tmp/requirements.txt

COPY mkdocs.yml .
COPY ./docs/ ./docs/
COPY ./primerdriver/ ./primerdriver/

RUN ./.venv/bin/mkdocs build

FROM bitnami/nginx:1.27.2-debian-12-r2

WORKDIR /app

COPY --from=build /build/site/ ./
