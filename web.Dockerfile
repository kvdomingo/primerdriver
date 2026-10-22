FROM oven/bun:1.1-alpine AS build

WORKDIR /web

COPY ./ui/ ./

RUN bun install && bun run build

FROM python:3.14-slim AS prod

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV UV_LINK_MODE=copy
ENV UV_PYTHON_DOWNLOADS=never

ARG SHORT_SHA
ENV SHORT_SHA=$SHORT_SHA

SHELL [ "/bin/bash", "-euxo", "pipefail", "-c" ]

COPY --from=ghcr.io/astral-sh/uv:0.13.4 /uv /usr/local/bin/uv

WORKDIR /primerdriver

COPY pyproject.toml uv.lock ./

RUN uv sync --frozen --no-dev --no-install-project --no-cache

COPY ./primerx/ ./primerx/
COPY ./primerdriver/ ./primerdriver/
COPY ./*.py ./
COPY ./*.sh ./
COPY --from=build /web/dist ./ui/

ENTRYPOINT [ "/bin/bash", "-euxo", "pipefail", "-c" ]
CMD [ "/primerdriver/docker-entrypoint.sh" ]
