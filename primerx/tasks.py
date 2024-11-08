import json

from primerdriver import expression_sys_loader
from primerx.cache import cache
from primerx.log import logger


def on_ready():
    expression_systems = expression_sys_loader.eager_load()
    try:
        cache.set("expression_systems", json.dumps(expression_systems))
        logger.info("Successfully cached expression systems")
    except Exception as e:
        logger.exception(e)
