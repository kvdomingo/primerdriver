from django.urls import path
from rest_framework.urlpatterns import format_suffix_patterns
from . import views


urlpatterns = [
    path("", views.index, name='index'),
    path("design/", views.CharacterizeList.as_view(), name="design"),
    path("design/<int:pk>", views.CharacterizeDetail.as_view())
]
