FROM python:3.14

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV UV_LINK_MODE=copy
ENV UV_PYTHON_DOWNLOADS=never
ENV UV_PROJECT_ENVIRONMENT=/tmp/.venv

COPY --from=ghcr.io/astral-sh/uv:0.13.4 /uv /usr/local/bin/uv

WORKDIR /tmp

COPY pyproject.toml uv.lock ./

RUN apt-get update && \
    apt-get install --no-install-recommends -y upx-ucl && \
    uv sync --frozen --no-install-project --no-cache && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build

ENTRYPOINT [ "/bin/bash", "-euxo", "pipefail", "-c" ]
CMD [ "/tmp/.venv/bin/pyinstaller -F --clean --name primerdriver --add-data primerdriver/settings.json:primerdriver --add-data primerdriver/tables:primerdriver/tables --add-data primerdriver/expression_systems:primerdriver/expression_systems setup.py" ]
