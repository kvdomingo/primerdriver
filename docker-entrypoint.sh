#!/bin/sh

gunicorn -b 0.0.0.0:$PORT -c ./gunicorn.conf.py