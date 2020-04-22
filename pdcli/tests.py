import io, subprocess
from pandas import DataFrame
from contextlib import redirect_stdout
from django.test import TestCase
from .primerclass import *


class CommandLineTestCase(TestCase):
    def test_return_help_on_no_arguments(self):
        proc = subprocess.run(
            ['python', 'pdcli.py'],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        self.assertEqual(proc.returncode, 0)

    def test_no_exception_on_help(self):
        proc = subprocess.run(
            ['python', 'pdcli.py', '-h'],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        self.assertEqual(proc.returncode, 0)


class CharacterizeTestCase(TestCase):
    pd = PrimerDesign('char', 'GATTACA', 'sub', 1, 'T', 'A', 4)

    def test_gc_content(self):
        self.assertAlmostEqual(
            self.pd.calculate_gc_content(self.pd.sequence),
            0.286,
            3
        )

    def test_calculate_mismatch(self):
        self.assertAlmostEqual(
            self.pd.calculate_mismatch(self.pd.sequence, self.pd.mismatched_bases),
            1/7
        )

    def test_reverse_complement(self):
        self.assertEqual(
            self.pd.get_reverse_complement(self.pd.sequence),
            list('TGTAATC')
        )

    def test_gc_end(self):
        self.assertEqual(
            self.pd.is_gc_end(self.pd.sequence),
            False
        )

    def test_calculate_Tm(self):
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
        out_catch = io.StringIO()
        with redirect_stdout(out_catch):
            out = self.pd.characterize_primer(
                self.pd.sequence,
                self.pd.mutation_type,
                self.pd.replacement,
                self.pd.mismatched_bases
            )
        self.assertEqual(isinstance(out, DataFrame), True)
