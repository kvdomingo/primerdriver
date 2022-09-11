import os
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()

BASE_DIR = Path(__file__).parent.parent

PYTHON_ENV = "development" if bool(os.environ.get("FLASK_DEBUG", False)) else "production"

SHORT_SHA = os.environ.get("SHORT_SHA", "")
