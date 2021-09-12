import io
from contextlib import redirect_stdout
from django.test import TestCase
from .primerclass import *


class CharacterizeTestCase(TestCase):
    databases = []
    pd = PrimerDesign('char', 'GATTACA', 'sub', 1, 'T', 'A', 4)

    def test_gc_content(self):
        """Calculate GC content"""
        self.assertAlmostEqual(
            self.pd.calculate_gc_content(self.pd.sequence),
            0.286,
            3
        )

    def test_calculate_mismatch(self):
        """Calculate percent mismatch"""
        self.assertAlmostEqual(
            self.pd.calculate_mismatch(self.pd.sequence, self.pd.mismatched_bases),
            1/7
        )

    def test_reverse_complement(self):
        """Calculate reverse complement"""
        self.assertEqual(
            self.pd.get_reverse_complement(self.pd.sequence),
            list('TGTAATC')
        )

    def test_gc_end(self):
        """Check if sequence ends in G/C"""
        self.assertFalse(self.pd.is_gc_end(self.pd.sequence))

    def test_calculate_Tm(self):
        """Calculate melting temperature"""
        self.assertAlmostEqual(
            self.pd.calculate_Tm(
                self.pd.sequence,
                self.pd.mutation_type,
                self.pd.replacement,
                0.286,
                self.pd.mismatched_bases
            ),
            -103.449,
            3
        )

    def test_characterization_is_dataframe(self):
        """Check if output is a DataFrame"""
        out_catch = io.StringIO()
        with redirect_stdout(out_catch):
            out = self.pd.characterize_primer(
                self.pd.sequence,
                self.pd.mutation_type,
                self.pd.replacement,
                self.pd.mismatched_bases
            )
        self.assertTrue(isinstance(out, DataFrame))
