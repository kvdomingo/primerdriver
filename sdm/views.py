from django.core.files import File
from django.shortcuts import render, HttpResponse
from django.templatetags.static import static
from django.conf import settings
from mimetypes import guess_type
from rest_framework import generics
from .models import Characterize
from .serializers import CharacterizeSerializer


class CharacterizeList(generics.ListCreateAPIView):
    queryset = Characterize.objects.all()
    serializer_class = CharacterizeSerializer


class CharacterizeDetail(generics.RetrieveUpdateDestroyAPIView):
    queryset = Characterize.objects.all()
    serializer_class = CharacterizeSerializer
    

def index(request):
    context = {
        "html_title": "Home",
        "active_page": "index",
        "settings": settings,
    }
    return render(request, "sdm/index.html.j2", context)


def train(request):
    context = {
        "html_title": "Driver",
        "active_page": "design",
        "settings": settings
    }
    return render(request, "sdm/train.html.j2", context)