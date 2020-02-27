from os import listdir
from json import loads, dumps
from pandas import concat
from django.http import HttpResponseRedirect, JsonResponse
from django.core.files import File
from django.shortcuts import render, HttpResponse
from django.templatetags.static import static
from django.urls import reverse
from django.conf import settings
from pdcli.primerclass import PrimerDesign
from pdcli.checks import PrimerChecks


def index(request):
    context = {
        "html_title": "Home",
        "active_page": "index",
        "settings": settings,
    }
    return render(request, "sdm/index.html.jinja2", context)


def train(request):
    context = {
        "html_title": "Design",
        "active_page": "design",
        "settings": settings
    }
    return render(request, "sdm/train.html.jinja2", context)


def characterize(request):
    context = {
        "html_title": "Primer Characterization",
        "active_page": "design",
        "settings": settings
    }
    return render(request, "sdm/characterize.html.jinja2", context)


def dna_based(request):
    context = {
        "html_title": "DNA-based",
        "active_page": "design",
        "settings": settings
    }
    return render(request, "sdm/dnabased.html.jinja2", context)

def protein_based(request):
    exp_systems = sorted([''.join(f.split('.')[:-1]).strip() for f in listdir(f'{settings.BASE_DIR}/pdcli/expression system')])
    context = {
        "html_title": "Protein-based",
        "active_page": "design",
        "settings": settings,
        "expression_system": exp_systems
    }
    return render(request, "sdm/proteinbased.html.jinja2", context)


def result(request):
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
    out = df.to_json()
    return HttpResponse(out)
