FROM python:3.10-bullseye

ENV PYTHONUNBUFFERED 1
ENV PYTHONDONTWRITEBYTECODE 1
ENV POETRY_VERSION 1.3.2

RUN apt-get update && apt-get install upx-ucl libgfortran-10-dev libquadmath0 -y

RUN pip install "poetry==$POETRY_VERSION"

RUN sh -c "$(curl --location https://taskfile.dev/install.sh)" -- -d -b /bin

WORKDIR /primerdriver

COPY pyproject.toml poetry.lock ./

RUN poetry config virtualenvs.create false && poetry install --no-interaction --no-ansi


ENTRYPOINT [ "/bin/task", "build" ]
