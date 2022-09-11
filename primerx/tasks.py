import json

from .cache import cache
from .config import BASE_DIR
from .log import logger


def on_ready():
    expression_systems = sorted(
        ["".join(f.name.split(".")[:-1]).strip() for f in (BASE_DIR / "primerdriver" / "expression system").glob("*")]
    )
    try:
        cache.set("expression_systems", json.dumps(expression_systems))
        logger.info("Successfully cached expression systems")
    except Exception as e:
        logger.exception(e)
