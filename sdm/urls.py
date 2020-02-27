from django.urls import path
from . import views


urlpatterns = [
    path("", views.index, name='index'),
    path("design", views.train, name="design"),
    path("characterize", views.characterize, name="characterize"),
    path("dna-based", views.dna_based, name="dna"),
    path("protein-based", views.protein_based, name="protein"),
    path("result", views.result, name="result")
]
