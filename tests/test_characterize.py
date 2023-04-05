import os
from contextlib import redirect_stdout

from pandas import DataFrame
from pytest import approx

from primerdriver.primer_design import MutationType, PrimerDesign

pd = PrimerDesign("char", "GATTACA", "SUB", 1, "T", "A", 4)

decimal_tolerance = 1e-3


def test_gc_content():
    """Calculate GC content"""
    assert pd.calculate_gc_content(pd.sequence) == approx(0.286, abs=decimal_tolerance)


def test_calculate_mismatch():
    """Calculate percent mismatch"""
    assert pd.calculate_mismatch(pd.sequence, pd.mismatched_bases) == approx(
        1 / 7, abs=decimal_tolerance
    )


def test_reverse_complement():
    """Calculate reverse complement"""
    assert pd.get_reverse_complement(pd.sequence) == list("TGTAATC")


def test_gc_end():
    """Check if sequence ends in G/C"""
    assert not pd.is_gc_end(pd.sequence)


def test_calculate_Tm():
    """Calculate melting temperature"""
    assert pd.calculate_melting_temperature(
        pd.sequence,
        pd.mutation_type,
        pd.replacement,
        0.286,
        pd.mismatched_bases,
    ) == approx(-103.449, rel=decimal_tolerance)


def test_characterization_returns_dataframe():
    """Check if output is a DataFrame"""
    with open(os.devnull, "w") as f:
        with redirect_stdout(f):
            out = pd.characterize_primer(
                sequence=pd.sequence,
                mutation_type=MutationType(pd.mutation_type),
                replacement=pd.replacement,
                mismatched_bases=pd.mismatched_bases,
            )
    assert isinstance(out, DataFrame)
