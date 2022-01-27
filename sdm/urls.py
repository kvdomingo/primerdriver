from django.urls import path
from . import views


urlpatterns = [
    path("", views.PrimerDriverAPIView.as_view(), name="api"),
    path("version", views.VersionView.as_view(), name="version"),
    path("expressionsys", views.ExpressionSystemsView.as_view(), name="expressions"),
]
