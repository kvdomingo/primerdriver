from __future__ import annotations

from enum import Enum

from numpy import array
from pandas import DataFrame
from tabulate import tabulate

from primerdriver import expression_sys_loader
from primerdriver.checks import SequenceChecks
from primerdriver.config import BASE_DIR, get_lookup_tables, get_settings
from primerdriver.exceptions import PrimerCheckError
from primerdriver.log import logger

settings = get_settings()

TABLES_DIR = BASE_DIR / "primerdriver" / "tables"

EXPRESSION_SYS_DIR = BASE_DIR / "primerdriver" / "expression_systems"

LOOKUP_TABLES = get_lookup_tables()


class OperationMode(Enum):
    CHARACTERIZATION = "CHAR"
    DNA = "DNA"
    PROTEIN = "PRO"


class MutationType(Enum):
    SUBSTITUTION = "S"
    INSERTION = "I"
    DELETION = "D"


class PrimerMode(Enum):
    COMPLEMENTARY = "complementary"
    OVERLAPPING = "overlapping"


class PrimerDesign:
    df: DataFrame
    lut = LOOKUP_TABLES

    def __init__(
        self,
        mode: str,
        sequence,
        mutation_type,
        mismatched_bases=None,
        target=None,
        replacement=None,
        position=None,
        savename=None,
        print_buffer=20,
        settings=settings,
    ):
        self.settings = settings
        self.mode = OperationMode(mode.upper())
        self.sequence: str = sequence.upper()
        self.mutation_type = MutationType(mutation_type[0].upper())
        self.mismatched_bases: int | None = (
            int(mismatched_bases) if mismatched_bases is not None else None
        )
        self.replacement: str | None = (
            replacement.upper() if replacement is not None else None
        )
        self.position: int | None = int(position) - 1 if position is not None else None
        self.target: str | None = target.upper() if target is not None else None
        self.Tm_range: tuple[float, float] = (
            self.settings.Tm_range_min,
            self.settings.Tm_range_max,
        )
        self.length_range: tuple[int, int] = (
            self.settings.length_min,
            self.settings.length_max,
        )
        self.gc_range: tuple[float, float] = (
            self.settings.gc_range_min,
            self.settings.gc_range_max,
        )
        self.flank5_range: tuple[int, int] = (
            self.settings.flank5_range_min,
            self.settings.flank5_range_max,
        )
        self.flank3_range: tuple[int, int] = (
            self.settings.flank3_range_min,
            self.settings.flank3_range_max,
        )
        self.primer_mode = PrimerMode(self.settings.primer_mode)
        self.print_buffer: int = print_buffer
        self.savename = savename
        self.expression_system: dict[str, str] = expression_sys_loader.lazy_load(
            self.settings.expression_system
        )

    @staticmethod
    def calculate_gc_content(seq: str | list[str]) -> float:
        """
        Calculate the GC content of a sequence.

        Parameters:
          seq: The sequence to calculate the GC content of.

        Returns:
          The %GC content of the sequence.
        """
        return (seq.count("G") + seq.count("C")) / len(seq)

    @staticmethod
    def calculate_mismatch(seq: str | list[str], mismatched_bases: int) -> float:
        """
        Calculate the base mismatch percentage of a sequence.

        Parameters:
          seq: The sequence to calculate the mismatch percentage of.
          mismatched_bases: The number of mismatched bases in the sequence.

        Returns:
          The base mismatch percentage of the sequence.
        """
        return mismatched_bases / len(seq)

    def get_reverse_complement(self, seq: str | list[str]) -> list[str]:
        """
        Get the reverse complement of a sequence.

        Parameters:
          seq: The sequence to get the reverse complement of.

        Returns:
          The reverse complement of the sequence.
        """
        seq = list(seq)
        return [self.lut["complement"][b] for b in seq][::-1]

    @staticmethod
    def is_gc_end(sequence: str | list[str]) -> bool:
        """
        Check if a sequence starts/ends in G or C.

        Parameters:
          sequence: The sequence to check.

        Returns:
          True if the sequence starts/ends in G or C, False otherwise.
        """
        sequence = "".join(sequence)
        return (sequence.startswith("G") or sequence.startswith("C")) and (
            sequence.endswith("G") or sequence.endswith("C")
        )

    @staticmethod
    def calculate_melting_temperature(
        seq: str | list[str],
        mutation_type: MutationType,
        replacement: str | list[str] | None,
        gc_content: float,
        mismatch: float,
    ) -> float:
        """
        Calculate the melting temperature of a sequence.

        Parameters:
          seq: The sequence to calculate the melting temperature of.
          mutation_type: The type of mutation (substitution, insertion, or deletion).
          replacement: The replacement sequence for the mutation.
          gc_content: The %GC content of the sequence.
          mismatch: The base mismatch percentage of the sequence.

        Returns:
          The melting temperature of the sequence in degrees Celsius.
        """
        Tm = 81.5 + 0.41 * gc_content * 100
        if mutation_type == MutationType.SUBSTITUTION:
            return Tm - 675 / len(seq) - mismatch * 100
        # Indel N excludes inserted bases; deleted bases are already absent from the primer.
        if mutation_type == MutationType.INSERTION:
            return Tm - 675 / (len(seq) - len(replacement or ""))
        return Tm - 675 / len(seq)

    def characterize_primer(
        self,
        sequence: str | list[str],
        mutation_type: MutationType,
        replacement: str | list[str] | None,
        mismatched_bases: int,
        index: int | None = None,
        reverse: str | None = None,
    ) -> DataFrame:
        """
        Characterize a primer.

        Parameters:
          sequence: The sequence of the primer.
          mutation_type: The type of mutation (substitution, insertion, or deletion).
          replacement: The replacement sequence for the mutation.
          mismatched_bases: The number of mismatched bases in the sequence.
          index: The index of the primer in the list of primers.
          reverse: The reverse primer; defaults to the reverse complement of the sequence.

        Returns:
          A DataFrame containing the characteristics of the primer.
        """
        mol_weight = self.lut["mol_weight"]
        seq = list(sequence)
        rev = reverse or "".join(self.get_reverse_complement(sequence))

        forwards = {
            "sequence": sequence,
            "length": len(seq),
            "gc_content": self.calculate_gc_content(seq),
            "mismatch": self.calculate_mismatch(seq, mismatched_bases),
            "gc_end": self.is_gc_end(seq),
            "mol_weight": sum(float(mol_weight[b]) for b in seq),
        }
        forwards["Tm"] = self.calculate_melting_temperature(
            seq,
            mutation_type,
            replacement,
            forwards["gc_content"],
            forwards["mismatch"],
        )

        reverses = {
            "sequence": rev,
            "length": len(rev),
            "gc_content": self.calculate_gc_content(rev),
            "mismatch": self.calculate_mismatch(rev, mismatched_bases),
            "gc_end": self.is_gc_end(rev),
            "mol_weight": sum(float(mol_weight[b]) for b in rev),
        }
        reverses["Tm"] = self.calculate_melting_temperature(
            rev,
            mutation_type,
            replacement,
            reverses["gc_content"],
            reverses["mismatch"],
        )

        col = [
            "Forward",
            "Reverse",
            "Fwd length",
            "Rev length",
            "Fwd GC content",
            "Rev GC content",
            "Fwd melting temp",
            "Rev melting temp",
            "Fwd mol. weight",
            "Rev mol. weight",
            "Fwd mismatch",
            "Rev mismatch",
            "Fwd ends in G/C",
            "Rev ends in G/C",
        ]
        dat = [
            forwards["sequence"],
            reverses["sequence"],
            f"{forwards['length']} bp",
            f"{reverses['length']} bp",
            f"{forwards['gc_content'] * 100:.2f}%",
            f"{reverses['gc_content'] * 100:.2f}%",
            f"{forwards['Tm']:.2f} C",
            f"{reverses['Tm']:.2f} C",
            f"{forwards['mol_weight']:.2f} g/mol",
            f"{reverses['mol_weight']:.2f} g/mol",
            f"{forwards['mismatch'] * 100:.2f}%",
            f"{reverses['mismatch'] * 100:.2f}%",
            forwards["gc_end"],
            reverses["gc_end"],
        ]
        if index is None:
            index = 1
        if index - 1 < self.print_buffer:
            print(
                "\n",
                tabulate(
                    array([col, dat]).T, headers=[f"Primer {index}"], tablefmt="orgtbl"
                ),
                sep="",
            )
        elif index - 1 == self.print_buffer:
            logger.info("Too many results; truncating output...")
        dat = array([dat])
        df = DataFrame(data=dat, columns=col, index=[index])
        self.df = df
        return df

    def _flank_pairs(self):
        for f5 in range(self.flank5_range[0], self.flank5_range[1] + 1):
            for f3 in range(self.flank3_range[0], self.flank3_range[1] + 1):
                if self.settings.center_mutation and abs(f5 - f3) > 1:
                    continue
                yield f5, f3

    def _passing_melting_temp(
        self,
        candidate: str,
        mutation_type: MutationType,
        replacement: str,
        mismatched_bases: int,
    ) -> float | None:
        """Return the candidate's melting temperature if it passes all primer checks, else None."""
        melting_temp = self.calculate_melting_temperature(
            candidate,
            mutation_type,
            replacement,
            self.calculate_gc_content(candidate),
            self.calculate_mismatch(candidate, mismatched_bases),
        )
        sc = SequenceChecks(candidate, self.settings)
        if (
            sc.is_valid_gc_content()
            and sc.is_valid_melting_temp(melting_temp)
            and sc.is_gc_clamped()
            and sc.is_valid_length()
        ):
            return melting_temp
        return None

    def design_primers(  # noqa: C901
        self,
        template: str,
        mutation_type: MutationType,
        position: int,
        deleted_length: int,
        inserted: str,
        mismatched_bases: int,
    ) -> list[DataFrame] | None:
        """
        Generate mutagenic primer pairs.

        Parameters:
          template: The template DNA sequence.
          mutation_type: The type of mutation (substitution, insertion, or deletion).
          position: 0-based index in `template` where the mutation starts.
          deleted_length: The number of template bases removed at `position`.
          inserted: The bases placed at `position` in the mutant.
          mismatched_bases: The number of mismatched bases in the primer.

        Returns:
          A list of DataFrames, one per primer pair, or None if no valid primers were found.
        """
        mutant = template[:position] + inserted + template[position + deleted_length :]
        mutation_end = position + len(inserted)

        def region(f5: int, f3: int) -> str | None:
            start, end = position - f5, mutation_end + f3
            if start < 0 or end > len(mutant):
                return None
            return mutant[start:end]

        pairs = []
        for f5, f3 in self._flank_pairs():
            forward = region(f5, f3)
            if forward is None:
                continue
            fwd_Tm = self._passing_melting_temp(
                forward, mutation_type, inserted, mismatched_bases
            )
            if fwd_Tm is None:
                continue
            if self.primer_mode == PrimerMode.COMPLEMENTARY:
                pairs.append((forward, "".join(self.get_reverse_complement(forward))))
                continue
            if f5 < self.settings.forward_overlap5:
                continue
            # The reverse primer reaches further 5' and less far 3' than the forward, so their 3' ends don't anneal to each other.
            for r5, r3 in self._flank_pairs():
                if r5 <= f5 or r3 >= f3 or r3 < self.settings.forward_overlap3:
                    continue
                reverse_region = region(r5, r3)
                if reverse_region is None:
                    continue
                reverse = "".join(self.get_reverse_complement(reverse_region))
                rev_Tm = self._passing_melting_temp(
                    reverse, mutation_type, inserted, mismatched_bases
                )
                if rev_Tm is not None and SequenceChecks.are_melting_temps_close(
                    fwd_Tm, rev_Tm
                ):
                    pairs.append((forward, reverse))

        if not pairs:
            print("No valid primers found")
            return None
        print(f"\nGenerated primer pairs: {len(pairs)}")
        if self.mode == OperationMode.PROTEIN:
            print(f"Using expression system: {self.settings.expression_system}")
        return [
            self.characterize_primer(
                forward, mutation_type, inserted, mismatched_bases, i, reverse
            )
            for i, (forward, reverse) in enumerate(pairs, start=1)
        ]

    def _back_translate(self, protein: str) -> str:
        """
        Reconstruct a DNA template from a protein using one preferred codon per residue.

        The reconstruction matches a real gene only where that gene uses the expression system's
        preferred codons, so primers designed from it mismatch a non-optimized template at silent
        positions. Documented as a known limitation of protein mode; DNA mode avoids it.
        """
        return "".join(self.expression_system[aa] for aa in protein).replace("U", "T")

    def _resolve_mutation(self) -> tuple[str, str]:
        """Return the bases removed from, and inserted into, the template at `self.position`."""
        if self.mutation_type == MutationType.INSERTION:
            if not self.replacement:
                raise PrimerCheckError("Insertion requires a replacement sequence")
            deleted, inserted = "", self.replacement
        elif self.mutation_type == MutationType.DELETION:
            if not self.target:
                raise PrimerCheckError("Deletion requires a target sequence")
            deleted, inserted = self.target, ""
        else:
            if not self.replacement:
                raise PrimerCheckError("Substitution requires a replacement sequence")
            inserted = self.replacement
            deleted = (
                self.target
                or self.sequence[self.position : self.position + len(inserted)]
            )
        return deleted, inserted

    def mutagenesis(self) -> list[DataFrame] | None:
        """DNA or protein based primer design."""
        if self.position is None:
            raise PrimerCheckError("Mutation position is required")
        deleted, inserted = self._resolve_mutation()

        if not 0 <= self.position <= len(self.sequence) - len(deleted):
            raise PrimerCheckError(
                f"Position {self.position + 1} is outside the sequence"
            )
        found = self.sequence[self.position : self.position + len(deleted)]
        if found != deleted:
            raise PrimerCheckError(
                f"Target '{deleted}' does not match '{found}' at position {self.position + 1}"
            )

        template, position = self.sequence, self.position
        if self.mode == OperationMode.PROTEIN:
            template, deleted, inserted = (
                self._back_translate(s) for s in (template, deleted, inserted)
            )
            position *= 3
        if self.mismatched_bases is None:
            self.mismatched_bases = len(
                deleted if self.mutation_type == MutationType.DELETION else inserted
            )
        return self.design_primers(
            template,
            self.mutation_type,
            position,
            len(deleted),
            inserted,
            self.mismatched_bases,
        )

    def main(self):
        if self.mode == OperationMode.CHARACTERIZATION:
            df = self.characterize_primer(
                self.sequence,
                self.mutation_type,
                self.replacement,
                self.mismatched_bases or 0,
            )
        else:
            df = self.mutagenesis()
        self.df = df
