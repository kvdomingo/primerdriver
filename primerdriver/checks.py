from __future__ import annotations

from json import load

from primerdriver.config import BASE_DIR, Settings, get_settings
from primerdriver.exceptions import PrimerCheckError

settings = get_settings()

TABLES_DIR = BASE_DIR / "primerdriver" / "tables"


class PrimerChecks:
    def __init__(self, sequence: str):
        """
        A set of validation checks to perform on the DNA/protein input before processing.

        Args:
            sequence: The input DNA/protein sequence.
        """
        self.sequence = sequence

    def is_valid_dna(self) -> str:
        """
        Check if the `self.sequence` contains valid bases [ATCG].

        Raises:
            primerdriver.checks.PrimerCheckError: DNA/protein sequence contains bases/amino acids.
        Returns:
             The input DNA sequence.
        """
        unique_bases = set(self.sequence.upper())
        true_bases = {"A", "C", "T", "G"}
        invalid_bases = unique_bases.difference(true_bases)
        if len(invalid_bases) > 0:
            raise PrimerCheckError(
                f"Sequence contains invalid bases: {', '.join(list(invalid_bases))}"
            )
        return self.sequence

    def is_valid_protein(self) -> str:
        """
        Check if the `self.sequence` contains valid amino acids.

        Raises:
            primerdriver.checks.PrimerCheckError: Protein sequence contains invalid amino acids.
        Returns:
             The input protein sequence.
        """
        unique_proteins = set(self.sequence.upper())
        with open(TABLES_DIR / "AAcompressed.json") as f:
            true_proteins = load(f)
        invalid_proteins = unique_proteins.difference(true_proteins.keys())
        if len(invalid_proteins) != 0:
            raise PrimerCheckError(
                f"Sequence contains invalid proteins: {', '.join(list(invalid_proteins))}"
            )
        return self.sequence

    def is_valid_sequence_length(self) -> None:
        """
        Check if the `self.sequence` is within the allowed processing length (40 <= sequence <= 8000).

        Raises:
             primerdriver.checks.PrimerCheckError: DNA/protein sequence is too short/long.
        """
        if len(self.sequence) < 40:
            raise PrimerCheckError("DNA sequence is too short")
        if len(self.sequence) > 8000:
            raise PrimerCheckError("DNA sequence is too long")

    def is_valid_gc_content(self) -> None:
        """
        Check if the `self.sequence` has valid %GC content (determined by settings.json).

        Raises:
             primerdriver.checks.PrimerCheckError: DNA/protein sequence has too little/too much GC content.
        """
        seq = list(self.sequence)
        gc = (seq.count("C") + seq.count("G")) / len(seq)
        if gc < 0.40:
            error_message = "GC content is less than 40%"
            raise PrimerCheckError(error_message)
        elif gc > 0.60:
            error_message = "GC content is greater than 60%"
            raise PrimerCheckError(error_message)


class SequenceChecks:
    def __init__(self, sequence: str | list[str], settings: Settings = settings):
        """
        Set of checks to perform on a generated primer.

        Args:
             sequence: The generated primer's DNA sequence.
             settings: The design settings to check against.
        """
        self.sequence = "".join(sequence).upper()
        self.settings = settings

    def is_valid_length(self) -> bool:
        return self.settings.length_min <= len(self.sequence) <= self.settings.length_max

    def is_valid_gc_content(self) -> bool:
        seq = list(self.sequence)
        gc = (seq.count("C") + seq.count("G")) / len(seq) * 100
        return self.settings.gc_range_min < gc < self.settings.gc_range_max

    def is_valid_melting_temp(self, melting_temp) -> bool:
        return self.settings.Tm_range_min < melting_temp < self.settings.Tm_range_max

    @staticmethod
    def are_melting_temps_close(fwd_Tm, rev_Tm) -> bool:
        return abs(fwd_Tm - rev_Tm) <= 2

    def is_gc_clamped(self) -> bool:
        return not self.settings.terminate_gc or (
            self.sequence[0] in ["C", "G"] and self.sequence[-1] in ["C", "G"]
        )
