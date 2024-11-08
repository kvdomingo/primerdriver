import json
import os
from functools import lru_cache
from pathlib import Path
from typing import Any, Literal

from dotenv import load_dotenv
from pydantic import BaseSettings

load_dotenv()

BASE_DIR = Path(__file__).parent.parent

PYTHON_ENV = (
    "development" if bool(os.environ.get("FLASK_DEBUG", False)) else "production"
)

SHORT_SHA = os.environ.get("SHORT_SHA", "")


class Settings(BaseSettings):
    Tm_range_min: float
    Tm_range_max: float
    gc_range_min: float
    gc_range_max: float
    length_min: int
    length_max: int
    flank5_range_min: int
    flank5_range_max: int
    flank3_range_min: int
    flank3_range_max: int
    forward_overlap5: int
    forward_overlap3: int
    terminate_gc: bool
    center_mutation: bool
    primer_mode: Literal["overlapping", "complementary"]
    expression_system: str

    class Config:
        @staticmethod
        def json_config_settings_source(_):
            dev_loc = BASE_DIR / "primerdriver" / "settings.json"
            prod_loc = Path.home() / ".primerdriverrc.json"

            if dev_loc.resolve().exists():
                with open(dev_loc) as f:
                    return json.load(f)

            if prod_loc.resolve().exists():
                with open(prod_loc) as f:
                    return json.load(f)

            raise FileNotFoundError("Config JSON not found.")

        @classmethod
        def customise_sources(cls, init_settings, env_settings, file_secret_settings):
            return (
                init_settings,
                cls.json_config_settings_source,
                env_settings,
                file_secret_settings,
            )


@lru_cache
def get_settings():
    return Settings()


@lru_cache
def get_lookup_tables() -> dict[str, Any]:
    with open(BASE_DIR / "primerdriver" / "tables" / "lut.json") as f:
        return json.load(f)
