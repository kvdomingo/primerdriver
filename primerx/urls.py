from django.contrib import admin
from django.conf import settings
from django.urls import path, re_path, include
from django.shortcuts import render

urlpatterns = [
    path("admin/", admin.site.urls),
    path("api/", include("sdm.urls")),
]

if settings.PYTHON_ENV == "development":
    from revproxy.views import ProxyView

    urlpatterns.append(re_path(r"^(?P<path>.*)$", ProxyView.as_view(upstream="http://frontend:3000")))
else:
    urlpatterns.append(re_path(r"^.*/?$", lambda req: render(req, "index.html")))
