from json import load

from primerdriver.exceptions import PrimerCheckError
from primerdriver.log import logger


@logger.catch
class PrimerChecks:
    def __init__(self, sequence: str, no_interaction: bool = False):
        self.sequence = sequence
        self.no_interaction = no_interaction

    def check_valid_base(self) -> str:
        unique_bases = set(list(self.sequence.upper()))
        true_bases = {"A", "C", "T", "G"}
        invalid_bases = unique_bases.difference(true_bases)
        if len(invalid_bases) != 0:
            if self.no_interaction:
                logger.warning(
                    "Sequence contains invalid bases. Automatically removing..."
                )
                for b in invalid_bases:
                    self.sequence = self.sequence.upper().replace(b, "")
            else:
                raise PrimerCheckError(
                    f"Sequence contains invalid bases: {', '.join(list(invalid_bases))}"
                )
        return self.sequence

    def check_valid_protein(self) -> str:
        unique_proteins = set(list(self.sequence.upper()))
        with open("primerdriver/AAcompressed.json", "r") as f:
            true_proteins = load(f)
        invalid_proteins = unique_proteins.difference(true_proteins.keys())
        if len(invalid_proteins) != 0:
            if self.no_interaction:
                logger.warning(
                    "Sequence contains invalid proteins. Automatically removing..."
                )
                for b in invalid_proteins:
                    self.sequence = self.sequence.upper().replace(b, "")
            else:
                raise PrimerCheckError(
                    f"Sequence contains invalid proteins: {', '.join(list(invalid_proteins))}"
                )
        return self.sequence

    def check_sequence_length(self) -> None:
        if len(self.sequence) < 40:
            error_message = "DNA sequence is too short"
        elif len(self.sequence) > 8000:
            error_message = "DNA sequence is too long"
        else:
            return

        if self.no_interaction:
            logger.warning(error_message)
        else:
            raise PrimerCheckError(error_message)

    def check_gc_content(self) -> None:
        seq = list(self.sequence)
        gc = (seq.count("C") + seq.count("G")) / len(seq)
        if gc < 0.40:
            error_message = "GC content is less than 40%"
            if self.no_interaction:
                logger.warning(error_message)
            else:
                raise PrimerCheckError(error_message)
        elif gc > 0.60:
            error_message = "GC content is greater than 60%"
            if self.no_interaction:
                logger.warning(error_message)
            else:
                raise PrimerCheckError(error_message)


class SequenceChecks:
    def __init__(self, sequence: str | list[str]):
        self.sequence = sequence.upper()

    def check_sequence_length(self, length_range: tuple[int, int]) -> bool:
        return length_range[0] <= len(self.sequence) <= length_range[1]

    def check_gc_content(self, gc_range: [int, int]) -> bool:
        seq = list(self.sequence)
        gc = (seq.count("C") + seq.count("G")) / len(seq) * 100
        return not (gc < gc_range[0] or gc > gc_range[1])

    def check_Tm(self, Tm, Tm_range) -> bool:
        return not (Tm < Tm_range[0] or Tm > Tm_range[1])

    def check_close_Tm(self, fwd_Tm, rev_Tm) -> bool:
        return abs(fwd_Tm - rev_Tm) <= 2

    def check_ends_gc(self, terminate_gc) -> bool:
        return not terminate_gc or (
            self.sequence[0] in ["C", "G"] and self.sequence[-1] in ["C", "G"]
        )
