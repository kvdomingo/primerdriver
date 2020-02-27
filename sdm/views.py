from json import loads, dumps
from django.http import HttpResponseRedirect, JsonResponse
from django.core.files import File
from django.shortcuts import render, HttpResponse
from django.templatetags.static import static
from django.urls import reverse
from django.conf import settings
from pdcli.primerclass import *


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


def result(request):
    data = dict()
    for k, v in request.POST.items():
        data[k] = v
    res = PrimerDesign(**data)
    res.main()
    out = res.df.to_json()
    return HttpResponse(out)
