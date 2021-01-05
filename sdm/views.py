from os import listdir, environ
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
from pdcli.version import __version__


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
    return render(request, "sdm/index.html", context)


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
            try:
                df = concat([*res.df]).T
                out = df.to_dict()
            except TypeError:
                out = {'data': 'No valid primers found'}
        else:
            df = res.df
            out = df.to_dict()
        return JsonResponse(out)
    else:
        return HttpResponseForbidden(f'{request.method} not allowed on /api')


def api_version(request):
    if request.method == 'GET':
        return JsonResponse({
            'program_version': str(__version__),
            'web_version': f'(web {environ["HEROKU_RELEASE_VERSION"]})' if settings.ON_HEROKU else '',
        })
    else:
        return HttpResponseForbidden(f'{request.method} not allowed on /version')

def api_expressions(request):
    if request.method == 'GET':
        return JsonResponse({
            "data": sorted([''.join(f.split('.')[:-1]).strip() for f in listdir(f'{settings.BASE_DIR}/pdcli/expression system')])
        })
    else:
        return HttpResponseForbidden(f'{request.method} not allowed on /expressionsys')
