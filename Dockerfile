FROM python:3.9.7-bullseye as base

ENV PYTHONUNBUFFERED 1

COPY requirements.txt /tmp/requirements.txt

RUN python -m pip install -U pip setuptools

RUN pip install --no-cache-dir -r /tmp/requirements.txt

FROM base as dev

COPY requirements.dev.txt /tmp/requirements.dev.txt
COPY requirements.txt /tmp/requirements.txt

RUN pip install --no-cache-dir -r /tmp/requirements.dev.txt

WORKDIR /primerdriver

ENTRYPOINT python manage.py migrate && \
            SHORT_SHA=$(git show --format="%h" --no-patch) gunicorn primerx.wsgi \
            -b 0.0.0.0:5000 \
            --workers 2 \
            --threads 4 \
            --log-file - \
            --capture-output \
            --reload

FROM node:16-alpine as build

WORKDIR /web

COPY ./web/app/ ./

RUN yarn install --prod

RUN yarn build

FROM base as prod

WORKDIR /primerdriver

COPY ./primerdriver/ ./primerdriver/
COPY ./primerx/ ./primerx/
COPY ./sdm/ ./sdm/
COPY ./manage.py ./manage.py
COPY ./runserver.sh ./runserver.sh
COPY --from=build /web/build ./web/app/

ARG SHORT_SHA=$SHORT_SHA

ENV SHORT_SHA $SHORT_SHA

ENTRYPOINT python manage.py collectstatic --noinput \
            python manage.py migrate \
            gunicorn primerx.wsgi -b 0.0.0.0:$PORT --workers 1 --threads 2 --log-file -
