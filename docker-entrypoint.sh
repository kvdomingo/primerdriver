#!/bin/bash

exec /primerdriver/.venv/bin/gunicorn -b 0.0.0.0:$PORT -c ./gunicorn.conf.py
