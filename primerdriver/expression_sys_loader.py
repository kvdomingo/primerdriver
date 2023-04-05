import json

from primerdriver.config import BASE_DIR

EXPRESSION_SYSTEM_DIR = BASE_DIR / "primerdriver" / "expression_systems"


def eager_load():
    lookup = {}
    for filename in EXPRESSION_SYSTEM_DIR.glob("*.json"):
        name = filename.stem
        with open(EXPRESSION_SYSTEM_DIR / f"{name}.json", "r") as f:
            lookup[name] = json.load(f)
    return lookup


def lazy_load(name: str):
    with open(EXPRESSION_SYSTEM_DIR / f"{name}.json", "r") as f:
        return json.load(f)
