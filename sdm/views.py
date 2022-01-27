import os
from pandas import concat
from rest_framework import status
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.permissions import AllowAny
from django.conf import settings
from primerdriver.primer_design import PrimerDesign
from primerdriver.checks import PrimerChecks
from primerdriver.version import __version__

BASE_DIR = settings.BASE_DIR


class PrimerDriverAPIView(APIView):
    permission_classes = [AllowAny]

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
    permission_classes = [AllowAny]

    def get(self, request):
        return Response(
            data={
                "program_version": str(__version__),
                "web_version": str(__version__),
            }
        )


class ExpressionSystemsView(APIView):
    permission_classes = [AllowAny]

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
