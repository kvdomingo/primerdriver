from django.urls import path
from . import views


urlpatterns = [
    path("", views.api, name='api'),
    path("version", views.api_version, name='version'),
    path("expressionsys", views.api_expressions, name='expressions'),
]
