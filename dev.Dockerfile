FROM python:3.10-bullseye

ENV DEBIAN_FRONTEND noninteractive
ENV PYTHONUNBUFFERED 1
ENV PYTHONDONTWRITEBYTECODE 1
ENV POETRY_VERSION 1.8.3

WORKDIR /primerdriver

SHELL [ "/bin/bash", "-euxo", "pipefail", "-c" ]

RUN pip install --no-cache-dir poetry==$POETRY_VERSION && \
    poetry config virtualenvs.create true && \
    poetry config virtualenvs.in-project true

ENTRYPOINT [ "/bin/bash", "-euxo", "pipefail", "-c" ]
