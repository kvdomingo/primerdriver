import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Literal

from dotenv import load_dotenv
from pydantic_settings import BaseSettings, JsonConfigSettingsSource

load_dotenv()

BASE_DIR = Path(__file__).parent.parent


def _settings_file() -> Path:
    for loc in (
        BASE_DIR / "primerdriver" / "settings.json",
        Path.home() / ".primerdriverrc.json",
    ):
        if loc.exists():
            return loc
    raise FileNotFoundError("Config JSON not found.")


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

    @classmethod
    def settings_customise_sources(
        cls,
        settings_cls,
        init_settings,
        env_settings,
        dotenv_settings,
        file_secret_settings,
    ):
        return (
            init_settings,
            JsonConfigSettingsSource(settings_cls, json_file=_settings_file()),
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
