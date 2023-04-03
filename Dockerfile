FROM python:3.10-bullseye as base

ENV PYTHONUNBUFFERED 1
ENV PYTHONDONTWRITEBYTECODE 1
ENV POETRY_VERSION 1.3.2

RUN pip install poetry==$POETRY_VERSION && poetry config virtualenvs.create false

FROM base as dev

WORKDIR /tmp

COPY pyproject.toml poetry.lock ./

RUN poetry install

WORKDIR /primerdriver

ENTRYPOINT ["poetry", "run", "flask", "run", "-h", "0.0.0.0", "-p", "5000", "--reloSomead"]

FROM node:16-alpine as web_build

WORKDIR /web

COPY ./web/app/public ./public
COPY ./web/app/src ./src
COPY ./web/app/package.json ./web/app/yarn.lock ./

RUN yarn install && yarn build

FROM base as prod

ARG SHORT_SHA=$SHORT_SHA
ENV SHORT_SHA $SHORT_SHA

WORKDIR /tmp

COPY pyproject.toml poetry.lock ./

RUN poetry export -f requirements.txt | pip install --no-cache-dir -r /dev/stdin

WORKDIR /primerdriver

COPY ./primerdriver/ ./primerdriver/
COPY ./primerx/ ./primerx/
COPY ./*.py ./
COPY ./*.sh ./
COPY --from=web_build /web/build ./web/app/

RUN chmod +x ./docker-entrypoint.sh

EXPOSE $PORT

ENTRYPOINT [ "./docker-entrypoint.sh" ]
