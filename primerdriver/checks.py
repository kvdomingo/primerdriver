from json import load
from primerx.log import logger
from .exceptions import PrimerCheckError


@logger.catch
class PrimerChecks:
    def __init__(self, sequence, no_interaction=False):
        self.sequence = sequence
        self.no_interaction = no_interaction

    def check_valid_base(self):
        unique_bases = set(list(self.sequence.upper()))
        true_bases = {"A", "C", "T", "G"}
        invalid_bases = unique_bases.difference(true_bases)
        if len(invalid_bases) != 0:
            if self.no_interaction:
                logger.warning("Sequence contains invalid bases. Automatically removing...")
                for b in invalid_bases:
                    self.sequence = self.sequence.upper().replace(b, "")
            else:
                raise PrimerCheckError(f"Sequence contains invalid bases: {', '.join(list(invalid_bases))}")
        return self.sequence

    def check_valid_protein(self):
        unique_prots = set(list(self.sequence.upper()))
        with open("primerdriver/AAcompressed.json", "r") as f:
            true_prots = load(f)
        invalid_prots = unique_prots.difference(true_prots.keys())
        if len(invalid_prots) != 0:
            if self.no_interaction:
                logger.warning("Sequence contains invalid proteins. Automatically removing...")
                for b in invalid_prots:
                    self.sequence = self.sequence.upper().replace(b, "")
            else:
                raise PrimerCheckError(f"Sequence contains invalid proteins: {', '.join(list(invalid_prots))}")
        return self.sequence

    def check_sequence_length(self):
        if len(self.sequence) < 40:
            error_message = "DNA sequence is too short"
            if self.no_interaction:
                logger.warning(error_message)
            else:
                raise PrimerCheckError(error_message)
        elif len(self.sequence) > 8000:
            error_message = "DNA sequence is too long"
            if self.no_interaction:
                logger.warning(error_message)
            else:
                raise PrimerCheckError(error_message)
        else:
            pass

    def check_gc_content(self):
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
    def __init__(self, sequence):
        self.sequence = sequence.upper()

    def check_sequence_length(self, length_range):
        if length_range[0] <= len(self.sequence) <= length_range[1]:
            return True
        else:
            return False

    def check_gc_content(self, gc_range):
        seq = list(self.sequence)
        gc = (seq.count("C") + seq.count("G")) / len(seq) * 100
        if gc < gc_range[0] or gc > gc_range[1]:
            return False
        else:
            return True

    def check_Tm(self, Tm, Tm_range):
        if Tm < Tm_range[0] or Tm > Tm_range[1]:
            return False
        else:
            return True

    def check_close_Tm(self, fwd_Tm, rev_Tm):
        if abs(fwd_Tm - rev_Tm) <= 2:
            return True
        else:
            return False

    def check_ends_gc(self, terminate_gc):
        if not terminate_gc or (self.sequence[0] in ["C", "G"] and self.sequence[-1] in ["C", "G"]):
            return True
        else:
            return False
