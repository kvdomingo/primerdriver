from json import load
from warnings import warn


class PrimerChecks:
    def __init__(self, sequence):
        self.sequence = sequence

    def check_valid_base(self):
        unique_bases = set(list(self.sequence.upper()))
        true_bases = {"A", "C", "T", "G"}
        invalid_bases = unique_bases.difference(true_bases)
        if len(invalid_bases) != 0:
            warn("Sequence contains invalid bases. Automatically removing...", Warning)
            for b in invalid_bases:
                self.sequence = self.sequence.upper().replace(b, "")
        return self.sequence

    def check_valid_protein(self):
        unique_prots = set(list(self.sequence.upper()))
        with open("pdcli/AAcompressed.json", "r") as f:
            true_prots = load(f)
        invalid_prots = unique_prots.difference(true_prots.keys())
        if len(invalid_prots) != 0:
            warn("Sequence contains invalid proteins. Automatically removing...", Warning)
            for b in invalid_prots:
                self.sequence = self.sequence.upper().replace(b, "")
        return self.sequence

    def check_sequence_length(self):
        if len(self.sequence) < 40:
            warn("DNA sequence is too short", Warning)
        elif len(self.sequence) > 8000:
            warn("DNA sequence is too long", Warning)

    def check_gc_content(self):
        seq = list(self.sequence)
        gc = (seq.count("C") + seq.count("G")) / len(seq)
        if gc < 0.40:
            warn("GC content is less than 40%", Warning)
        elif gc > 0.60:
            warn("GC content is greater than 60%", Warning)


class SequenceChecks:
    def __init__(self, sequence):
        self.sequence = sequence.upper()

    def check_sequence_length(self, length_range):
        if len(self.sequence) >= length_range[0] and len(self.sequence) <= length_range[1]:
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
