import os
from pandas import concat
from rest_framework import status
from rest_framework.views import APIView
from rest_framework.response import Response
from django.conf import settings
from primerdriver.primerclass import PrimerDesign
from primerdriver.checks import PrimerChecks
from primerdriver.version import __version__

BASE_DIR = settings.BASE_DIR


class PrimerDriverAPIView(APIView):
    def post(self, request):
        data = request.data
        checks = PrimerChecks(data["sequence"])
        if data["mode"] == "PRO":
            data["sequence"] = checks.check_valid_protein()
        else:
            data["sequence"] = checks.check_valid_base()
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
        return Response(data=out, status=status.HTTP_200_OK)


class VersionView(APIView):
    def get(self, request):
        return Response(
            data={
                "program_version": str(__version__),
                "web_version": f'(web {os.environ.get("HEROKU_RELEASE_VERSION")})' if settings.ON_HEROKU else "",
            }
        )


class ExpressionSystemsView(APIView):
    def get(self, request):
        return Response(
            {
                "data": sorted(
                    [
                        "".join(f.split(".")[:-1]).strip()
                        for f in os.listdir(BASE_DIR / "primerdriver" / "expression system")
                    ]
                )
            }
        )
