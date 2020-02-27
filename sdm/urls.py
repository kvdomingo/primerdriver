from django.urls import path
from . import views


urlpatterns = [
    path("", views.index, name='index'),
    path("design/", views.train, name="design"),
    path("characterize/", views.characterize, name="characterize"),
    path("result", views.result, name="result")
]
