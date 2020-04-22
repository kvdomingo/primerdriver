import subprocess
from pdcli.primerclass import PrimerDesign
from django.test import TestCase


class CommandLineTests(TestCase):
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


class CharacterizeTests(TestCase):
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
