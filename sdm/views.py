from django.core.files import File
from django.shortcuts import render, HttpResponse
from django.templatetags.static import static
from mimetypes import guess_type


def index(request):
    context = {
        "html_title": "Home",
        "active_page": "index"
    }
    return render(request, "sdm/index.html.j2", context)


def download(request):
    filename = "primerdriver.py"
    path = static(f"sdm/media/public/{filename}")
    fl = File(open(path, "r"))
    mime_type = guess_type(fl)[0]
    response = HttpResponse(fl, content_type=mime_type)
    response["Content-Disposition"] = f"attachment; filename={filename}/"
    return response
