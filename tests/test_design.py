import os
from contextlib import redirect_stdout
from random import Random

from pytest import raises

from primerdriver.config import get_settings
from primerdriver.exceptions import PrimerCheckError
from primerdriver.output_handler import single_command_saver
from primerdriver.primer_design import PrimerDesign

TEMPLATE = "".join(Random(0).choices("ACGT", k=200))

RELAXED = get_settings().model_copy(
    update={
        "Tm_range_min": 50,
        "Tm_range_max": 100,
        "gc_range_min": 20,
        "gc_range_max": 80,
        "terminate_gc": False,
    }
)


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def design(primer_mode="complementary", settings=RELAXED, **kwargs):
    pd = PrimerDesign(
        settings=settings.model_copy(update={"primer_mode": primer_mode}), **kwargs
    )
    with open(os.devnull, "w") as f, redirect_stdout(f):
        pd.main()
    return pd.df


def pairs(dfs):
    return [(df["Forward"].iloc[0], df["Reverse"].iloc[0]) for df in dfs]


def test_substitution_targets_given_position():
    target = TEMPLATE[99]
    replacement = "A" if target != "A" else "C"
    mutant = TEMPLATE[:99] + replacement + TEMPLATE[100:]
    dfs = design(
        mode="dna",
        sequence=TEMPLATE,
        mutation_type="S",
        target=target,
        replacement=replacement,
        position=100,
    )
    assert dfs
    for fwd, rev in pairs(dfs):
        assert fwd in mutant and fwd not in TEMPLATE
        assert rev == revcomp(fwd)


def test_overlapping_reverse_is_reverse_complement_of_mutant():
    mutant = TEMPLATE[:100] + "GGG" + TEMPLATE[100:]
    dfs = design(
        "overlapping",
        mode="dna",
        sequence=TEMPLATE,
        mutation_type="I",
        replacement="GGG",
        position=101,
    )
    assert dfs
    for df in dfs:
        fwd, rev = df["Forward"].iloc[0], df["Reverse"].iloc[0]
        assert fwd in mutant and revcomp(rev) in mutant
        assert rev != revcomp(fwd)
        fwd_Tm = float(df["Fwd melting temp"].iloc[0].split()[0])
        rev_Tm = float(df["Rev melting temp"].iloc[0].split()[0])
        assert abs(fwd_Tm - rev_Tm) <= 2


def test_deletion_complementary():
    target = TEMPLATE[99:102]
    mutant = TEMPLATE[:99] + TEMPLATE[102:]
    dfs = design(
        mode="dna", sequence=TEMPLATE, mutation_type="D", target=target, position=100
    )
    assert dfs
    assert all(fwd in mutant for fwd, _ in pairs(dfs))


def test_mutation_near_start_never_wraps():
    target = TEMPLATE[2]
    dfs = design(
        mode="dna",
        sequence=TEMPLATE,
        mutation_type="S",
        target=target,
        replacement="A" if target != "A" else "C",
        position=3,
    )
    assert dfs is None


def test_wrong_target_rejected():
    wrong = "A" if TEMPLATE[99] != "A" else "C"
    with raises(PrimerCheckError):
        design(
            mode="dna",
            sequence=TEMPLATE,
            mutation_type="S",
            target=wrong,
            replacement="G",
            position=100,
        )


def test_custom_settings_are_honored():
    strict = RELAXED.model_copy(update={"length_min": 30, "length_max": 32})
    dfs = design(
        settings=strict,
        mode="dna",
        sequence=TEMPLATE,
        mutation_type="I",
        replacement="T",
        position=100,
    )
    assert dfs
    assert all(30 <= len(fwd) <= 32 for fwd, _ in pairs(dfs))


def test_protein_multi_residue_substitution():
    protein = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"
    dfs = design(
        mode="pro",
        sequence=protein,
        mutation_type="S",
        target="QR",
        replacement="WW",
        position=9,
    )
    assert dfs
    assert all("TGGTGG" in fwd for fwd, _ in pairs(dfs))


def test_missing_mutation_arguments_rejected():
    for mutation_type, kwargs in (
        ("I", {}),
        ("D", {}),
        ("S", {"target": TEMPLATE[99]}),
    ):
        with raises(PrimerCheckError):
            design(
                mode="dna",
                sequence=TEMPLATE,
                mutation_type=mutation_type,
                position=100,
                **kwargs,
            )


def test_characterization_without_mismatch_count():
    pd = PrimerDesign("char", "GATTACAGATTACA", "S")
    with open(os.devnull, "w") as f, redirect_stdout(f):
        pd.main()
    assert pd.df["Fwd mismatch"].iloc[0] == "0.00%"


def test_save_characterization_csv(tmp_path):
    pd = PrimerDesign("char", "GATTACAGATTACA", "S", mismatched_bases=1)
    with open(os.devnull, "w") as f, redirect_stdout(f):
        pd.main()
    single_command_saver(pd.df, str(tmp_path / "out.csv"))
    assert (tmp_path / "out.csv").exists()
