from flask_caching import Cache

from primerx.config import PYTHON_ENV

cache = Cache(
    config={
        "DEBUG": PYTHON_ENV == "development",
        "CACHE_TYPE": "SimpleCache",
        "CACHE_DEFAULT_TIMEOUT": 0,
    }
)
