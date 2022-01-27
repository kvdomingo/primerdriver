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

ENTRYPOINT python manage.py migrate && gunicorn primerx.wsgi -b 0.0.0.0:$PORT --log-file - --reload

FROM base as make-linux

RUN apt-get update
RUN apt-get install upx-ucl libgfortran-10-dev libquadmath0 -y

ENV PYTHONDONTWRITEBYTECODE 1

COPY requirements.dev.txt /tmp/requirements.dev.txt
COPY requirements.txt /tmp/requirements.txt

RUN pip install --no-cache-dir -r /tmp/requirements.dev.txt

WORKDIR /primerdriver

ENTRYPOINT [ "sh", "build.sh" ]

FROM node:16-alpine as build

WORKDIR /web

COPY ./web/app/ ./

RUN yarn install --prod

RUN yarn build

FROM base as prod

WORKDIR /primerdriver

COPY primerdriver/ ./pdcli/
COPY ./primerx/ ./primerx/
COPY ./sdm/ ./sdm/
COPY --from=build /web/build ./web/app/
COPY ./*.py ./

RUN python manage.py collectstatic --noinput

ENTRYPOINT gunicorn primerx.wsgi -b 0.0.0.0:$PORT --log-file -
