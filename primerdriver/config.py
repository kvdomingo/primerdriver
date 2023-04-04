import json
import os
from functools import lru_cache
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()

BASE_DIR = Path(__file__).parent.parent

PYTHON_ENV = (
    "development" if bool(os.environ.get("FLASK_DEBUG", False)) else "production"
)

SHORT_SHA = os.environ.get("SHORT_SHA", "")


@lru_cache
def get_settings():
    with open(BASE_DIR / "primerdriver" / "settings.json", "r") as f:
        return json.load(f)
