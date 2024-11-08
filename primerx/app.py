import json
import os
from http import HTTPStatus

from flask import Flask, Response, jsonify, request, send_from_directory
from pandas import concat

from primerdriver import __version__
from primerdriver.checks import PrimerChecks
from primerdriver.config import Settings
from primerdriver.exceptions import PrimerCheckError
from primerdriver.primer_design import OperationMode, PrimerDesign
from primerx.cache import cache
from primerx.config import BASE_DIR, PYTHON_ENV, SHORT_SHA
from primerx.log import logger
from primerx.models import AdvancedSettings
from primerx.tasks import on_ready

app = Flask(__name__, static_url_path="", static_folder=BASE_DIR / "ui")
cache.init_app(app)

on_ready()


@app.route("/api/health")
def health():
    return {"status": "ok"}


@app.route("/api/version")
def api_version():
    return jsonify(
        {
            "program_version": __version__,
            "web_version": f"web {SHORT_SHA if SHORT_SHA else __version__}",
        }
    )


@app.route("/api/expressionsys")
def api_expression_systems():
    try:
        expression = json.loads(cache.get("expression_systems"))
        return jsonify({"data": sorted(expression.keys())})
    except json.JSONDecodeError as e:
        logger.error(str(e))
        return Response(status=HTTPStatus.INTERNAL_SERVER_ERROR)


@app.route("/api", methods=["POST"])
def api():
    data = request.json
    data["settings"] = AdvancedSettings(
        **(data.get("settings", None) or Settings().dict())
    )
    checks = PrimerChecks(data["sequence"])
    if data["mode"] == OperationMode.PROTEIN.value:
        try:
            data["sequence"] = checks.is_valid_protein()
        except PrimerCheckError as e:
            logger.error(str(e))
            return Response(str(e), status=HTTPStatus.BAD_REQUEST)
    else:
        try:
            data["sequence"] = checks.is_valid_dna()
        except PrimerCheckError as e:
            logger.error(str(e))
            return Response(str(e), status=HTTPStatus.BAD_REQUEST)
    res = PrimerDesign(**data)
    res.main()
    if data["mode"] != OperationMode.CHARACTERIZATION.value:
        try:
            df = concat([*res.df]).T
            out = df.to_dict()
        except (TypeError, ValueError):
            out = {"data": "No valid primers found"}
    else:
        df = res.df
        out = df.to_dict()
    return jsonify(out)


if PYTHON_ENV != "development":

    @app.route("/", defaults={"path": ""})
    @app.route("/<path:path>")
    def serve(path: str):
        if path != "" and os.path.exists(f"{app.static_folder}/{path}"):
            return send_from_directory(app.static_folder, path)
        return send_from_directory(app.static_folder, "index.html")

    @app.errorhandler(404)
    def not_found(_):
        return send_from_directory(app.static_folder, "index.html")
