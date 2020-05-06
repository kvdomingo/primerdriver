from django.urls import path
from . import views


urlpatterns = [
    path("", views.index, name='index'),
    path("api", views.api, name='api'),
    path("api/version", views.api_version, name='version'),
    path("api/expressionsys", views.api_expressions, name='expressions'),
]
