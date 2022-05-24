import os
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()

BASE_DIR = Path(__file__).parent.parent

PYTHON_ENV = os.environ.get("FLASK_ENV", "production")

SHORT_SHA = os.environ.get("SHORT_SHA", "")
