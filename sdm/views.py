from os import listdir
from json import loads, dumps
from pandas import concat
from django.http import HttpResponse, HttpResponseRedirect, JsonResponse, HttpResponseForbidden
from django.core.files import File
from django.shortcuts import render
from django.templatetags.static import static
from django.urls import reverse
from django.conf import settings
from django.utils.html import escapejs
from pdcli.primerclass import PrimerDesign
from pdcli.checks import PrimerChecks


def index(request):
    exp_systems = {
        "data": sorted([''.join(f.split('.')[:-1]).strip() for f in listdir(f'{settings.BASE_DIR}/pdcli/expression system')])
    }
    context = {
        "html_title": "Home",
        "active_page": "index",
        "expression": escapejs(dumps(exp_systems)),
        "settings": settings,
    }
    return render(request, "sdm/index.html.j2", context)


def api(request):
    if request.method == 'POST':
        data = dict()
        for k, v in request.POST.items():
            data[k] = v
        checks = PrimerChecks(data['sequence'])
        if data['mode'] == 'PRO':
            data['sequence'] = checks.check_valid_protein()
        else:
            data['sequence'] = checks.check_valid_base()
        res = PrimerDesign(**data)
        res.main()
        if data['mode'] != 'CHAR':
            df = concat([*res.df]).T
        else:
            df = res.df
        out = df.to_dict()
        return JsonResponse(out)
    else:
        return HttpResponseForbidden(f'{request.method} not allowed on /api')
