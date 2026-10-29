from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse
from pandas import concat
from scalar_fastapi import get_scalar_api_reference

from primerdriver import __version__, expression_sys_loader
from primerdriver.checks import PrimerChecks
from primerdriver.config import get_settings
from primerdriver.exceptions import PrimerCheckError
from primerdriver.primer_design import OperationMode, PrimerDesign
from primerx.config import BASE_DIR, settings
from primerx.log import logger
from primerx.models import (
    DesignRequest,
    DesignResult,
    ExpressionSystemsResponse,
    MessageResponse,
    VersionResponse,
)

EXPRESSION_SYSTEMS = sorted(expression_sys_loader.eager_load())

app = FastAPI(
    title="PrimerDriver",
    version=__version__,
    docs_url=None,
    redoc_url=None,
    openapi_url="/api/openapi.json",
)


@app.exception_handler(PrimerCheckError)
def primer_check_error(_: Request, e: PrimerCheckError):
    logger.error(str(e))
    return JSONResponse({"detail": str(e)}, status_code=400)


@app.get("/api/docs", include_in_schema=False)
def docs():
    return get_scalar_api_reference(openapi_url=app.openapi_url, title=app.title)


@app.get("/api/health")
def health():
    raise RuntimeError("deliberately broken to test fast preview failure")


@app.get("/api/version")
def api_version() -> VersionResponse:
    return VersionResponse(
        program_version=__version__,
        web_version=f"web {settings.SHORT_SHA or __version__}",
    )


@app.get("/api/expressionsys")
def api_expression_systems() -> ExpressionSystemsResponse:
    return ExpressionSystemsResponse(data=EXPRESSION_SYSTEMS)


@app.post("/api")
def api(body: DesignRequest) -> DesignResult | MessageResponse:
    checks = PrimerChecks(body.sequence)
    sequence = (
        checks.is_valid_protein()
        if body.mode == OperationMode.PROTEIN
        else checks.is_valid_dna()
    )
    res = PrimerDesign(
        **body.model_dump(mode="json", exclude={"sequence", "settings"}),
        sequence=sequence,
        settings=body.settings or get_settings(),
    )
    res.main()
    if body.mode == OperationMode.CHARACTERIZATION:
        return res.df.to_dict()
    try:
        return concat([*res.df]).T.to_dict()
    except (TypeError, ValueError):
        return MessageResponse(data="No valid primers found")


if settings.PYTHON_ENV != "development":
    app.frontend("/", directory=BASE_DIR / "ui")
