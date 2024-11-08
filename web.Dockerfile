FROM oven/bun:1.1-alpine AS build

WORKDIR /web

COPY ./ui/ ./

RUN bun install && bun run build

FROM python:3.10-slim AS base

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV POETRY_VERSION=1.8.3

ARG SHORT_SHA=$SHORT_SHA
ENV SHORT_SHA=$SHORT_SHA

SHELL [ "/bin/bash", "-euxo", "pipefail", "-c" ]

# hadolint ignore=DL3009
RUN apt-get update && \
    apt-get install --no-install-recommends -y gcc libc-dev

FROM base AS export

RUN pip install --no-cache-dir "poetry==$POETRY_VERSION" && \
    poetry config virtualenvs.create false && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

COPY pyproject.toml poetry.lock ./

RUN poetry export --without-hashes --without dev --format requirements.txt --output requirements.txt

FROM base AS prod

WORKDIR /tmp

COPY --from=export /tmp/requirements.txt .

WORKDIR /primerdriver

RUN python -m venv .venv && \
    ./.venv/bin/pip install --no-cache-dir -r /tmp/requirements.txt

COPY ./primerx/ ./primerx/
COPY ./primerdriver/ ./primerdriver/
COPY ./*.py ./
COPY ./*.sh ./
COPY --from=build /web/dist ./ui/

ENTRYPOINT [ "/bin/bash", "-euxo", "pipefail", "-c" ]
CMD [ "/primerdriver/docker-entrypoint.sh" ]
