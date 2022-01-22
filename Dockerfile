FROM python:3.9.7-bullseye as base

ENV PYTHONUNBUFFERED 1

COPY requirements.txt /tmp/requirements.txt

RUN python -m pip install -U pip setuptools

RUN pip install --no-cache-dir -r /tmp/requirements.txt

FROM base as dev

COPY requirements.dev.txt /tmp/requirements.dev.txt
COPY requirements.txt /tmp/requirements.txt

RUN pip install --no-cache-dir -r /tmp/requirements.dev.txt

RUN sed -i "s/'_headers'/'headers'/" /usr/local/lib/python3.9/site-packages/revproxy/utils.py
RUN sed -i "s/'_headers'/'headers'/" /usr/local/lib/python3.9/site-packages/revproxy/response.py

WORKDIR /primerdriver

ENTRYPOINT gunicorn primerx.wsgi -b 0.0.0.0:$PORT --log-file - --reload

FROM node:16-alpine as build

WORKDIR /web

COPY ./web/app/ ./

RUN yarn install --prod

RUN yarn build

FROM base as prod

WORKDIR /primerdriver

COPY ./pdcli/ ./pdcli/
COPY ./primerx/ ./primerx/
COPY ./sdm/ ./sdm/
COPY --from=build /web/build ./web/app/
COPY ./*.py ./

ENTRYPOINT gunicorn primerx.wsgi -b 0.0.0.0:$PORT --log-file -
