FROM python:3.14-slim AS build

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV UV_LINK_MODE=copy
ENV UV_PYTHON_DOWNLOADS=never

SHELL [ "/bin/bash", "-euxo", "pipefail", "-c" ]

COPY --from=ghcr.io/astral-sh/uv:0.13.6 /uv /usr/local/bin/uv

WORKDIR /build

COPY pyproject.toml uv.lock ./

RUN uv sync --frozen --no-install-project --no-cache

COPY mkdocs.yml .
COPY ./docs/ ./docs/
COPY ./primerdriver/ ./primerdriver/

RUN ./.venv/bin/mkdocs build

FROM nginx/nginx-unprivileged:1.31-alpine

COPY --from=build /build/site/ /usr/share/nginx/html/
