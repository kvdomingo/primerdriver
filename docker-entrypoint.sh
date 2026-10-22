#!/bin/bash

exec /primerdriver/.venv/bin/fastapi run primerx/app.py --port "$PORT"
