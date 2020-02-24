from django.db import models


class Characterize(models.Model):
    sequence = models.TextField()
    mut_type = models.CharField(max_length=3)
    mismatch = models.IntegerField()

    def __str__(self):
        return f'{self.mismatch}-{self.mut_type} for {self.sequence}'