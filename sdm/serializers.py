from rest_framework import serializers
from .models import Characterize


class CharacterizeSerializer(serializers.ModelSerializer):
    class Meta:
        model = Characterize
        fields = (
            'sequence', 
            'mut_type', 
            'mismatch'
        )