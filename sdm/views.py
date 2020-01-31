from django.shortcuts import render


def index(request):
    context = {
        "html_title": "Home",
        "active_page": "index"
    }
    return render(request, "sdm/index.html.j2", context)
