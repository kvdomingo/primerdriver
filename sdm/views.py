from json import loads, dumps
from django.http import HttpResponseRedirect, JsonResponse
from django.core.files import File
from django.shortcuts import render, HttpResponse
from django.templatetags.static import static
from django.urls import reverse
from django.conf import settings
from django.views.decorators.csrf import csrf_exempt
from mimetypes import guess_type
from .forms import Characterize
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
        "html_title": "Driver",
        "active_page": "design",
        "settings": settings
    }
    return render(request, "sdm/train.html.jinja2", context)


def characterize(request):
    if request.method == "POST":
        form = Characterize(request.POST)
        if form.is_valid():
            return HttpResponseRedirect("sdm/characterize.html.jinja2")
    else:
        form = Characterize()
    context = {
        "html_title": "Driver",
        "active_page": "design",
        "settings": settings,
        "form": form
    }
    return render(request, "sdm/characterize.html.jinja2", context)


def result(request):
    data = dict()
    if request.POST.get('mode') == 'CHAR':
        for k, v in request.POST.items():
            data[k] = v
        res = PrimerDesign(**data)
        res.main()
        out = res.df.to_json()
    return HttpResponse(out)
