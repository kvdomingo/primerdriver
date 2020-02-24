from django import forms


class Characterize(forms.Form):
    sequence = forms.CharField(label="Enter sequence", max_length=20000, widget=forms.Textarea)
    mut_type = forms.CharField(label="Mutation type", max_length=3)
    mismatch = forms.IntegerField(label="Mismatched bases", min_value=0)