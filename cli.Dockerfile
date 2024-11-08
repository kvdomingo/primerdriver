FROM python:3.10-bullseye

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV POETRY_VERSION=1.8.3
ENV POETRY_HOME=/opt/poetry
ENV PATH=${PATH}:${POETRY_HOME}/bin

WORKDIR /tmp

COPY pyproject.toml poetry.lock ./

RUN apt-get update && \
    apt-get install --no-install-recommends -y upx-ucl && \
    pip install --no-cache-dir "poetry==$POETRY_VERSION" && \
    poetry export --without-hashes --with dev --format requirements.txt --output requirements.txt && \
    python -m venv .venv && \
    /tmp/.venv/bin/pip install --no-cache-dir -r requirements.txt && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build

ENTRYPOINT [ "/bin/bash", "-euxo", "pipefail", "-c" ]
CMD [ "/tmp/.venv/bin/pyinstaller -F --clean --name primerdriver setup.py" ]
