import os
import json
from flask import Flask, request, Response, send_from_directory, jsonify
from http import HTTPStatus as status
from pandas import concat
from primerdriver import __version__
from primerdriver.checks import PrimerChecks
from primerdriver.exceptions import PrimerCheckError
from primerdriver.primer_design import PrimerDesign
from .cache import cache
from .config import BASE_DIR, SHORT_SHA, PYTHON_ENV
from .log import logger
from .tasks import on_ready

app = Flask(__name__, static_url_path="", static_folder=BASE_DIR / "web" / "app")
cache.init_app(app)

on_ready()


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
        return jsonify({"data": expression})
    except json.JSONDecodeError as e:
        logger.error(str(e))
        return Response(status=status.INTERNAL_SERVER_ERROR)


@app.route("/api", methods=["POST"])
def api():
    data = request.json
    checks = PrimerChecks(data["sequence"])
    if data["mode"] == "PRO":
        try:
            data["sequence"] = checks.check_valid_protein()
        except PrimerCheckError as e:
            logger.error(str(e))
            return Response(str(e), status=status.BAD_REQUEST)
    else:
        try:
            data["sequence"] = checks.check_valid_base()
        except PrimerCheckError as e:
            logger.error(str(e))
            return Response(str(e), status=status.BAD_REQUEST)
    res = PrimerDesign(**data)
    res.main()
    if data["mode"] != "CHAR":
        try:
            df = concat([*res.df]).T
            out = df.to_dict()
        except TypeError:
            out = {"data": "No valid primers found"}
    else:
        df = res.df
        out = df.to_dict()
    return jsonify(out)


if PYTHON_ENV != "development":

    @app.route("/", defaults={"path": ""})
    @app.route("/<path:path")
    def serve(path: str):
        if path != "" and os.path.exists(f"{app.static_folder}/{path}"):
            return send_from_directory(app.static_folder, path)
        return send_from_directory(app.static_folder, "index.html")

    @app.errorhandler(404)
    def not_found(_):
        return send_from_directory(app.static_folder, "index.html")
