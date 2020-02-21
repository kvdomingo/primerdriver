from django.core.files import File
from django.shortcuts import render, HttpResponse
from django.templatetags.static import static
from django.conf import settings
from mimetypes import guess_type


def index(request):
    context = {
        "html_title": "Home",
        "active_page": "index",
        "settings": settings,
    }
    return render(request, "sdm/index.html.j2", context)