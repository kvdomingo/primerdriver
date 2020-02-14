from warnings import warn
from json import load
from pandas import DataFrame
from tabulate import tabulate
from numpy import array


class PrimerDesign:
    def __init__(
        self,
        mode,
        sequence,
        mutation_type,
        mismatched_bases=None,
        target=None,
        replacement=None,
        position=None,
        Tm_range=(75, 85),
        length_range=(25, 45),
        gc_range=(40, 60),
        flank5_range=(11, 21),
        flank3_range=(11, 21),
        terminate_gc=True,
        center_mutation=True,
        primer_mode=0,
        savename=None
    ):
        self.mode = mode.upper()
        self.sequence = sequence.upper()
        self.mutation_type = mutation_type[0].upper()
        self.mismatched_bases = int(mismatched_bases) if mismatched_bases is not None else None
        self.replacement = replacement.upper() if replacement is not None else None
        self.position = int(position) if position is not None else None
        self.target = target.upper() if target is not None else None
        self.Tm_range = Tm_range
        self.length_range = length_range
        self.gc_range = gc_range
        self.flank5_range = flank5_range
        self.flank3_range = flank3_range
        self.terminate_gc = terminate_gc
        self.center_mutation = center_mutation
        self.primer_mode = primer_mode
        self.savename = savename
        with open("pdcli/lut.json", "r", encoding="utf-8") as f:
            lut = load(f)
        self.lut = lut

    def calculate_gc_content(self, seq):
        return (seq.count('G') + seq.count('C'))/len(seq)

    def calculate_mismatch(self, seq):
        return self.mismatched_bases/len(seq)

    def get_reverse_complement(self, seq):
        seq = list(seq)
        return [self.lut["complement"][b] for b in seq][::-1]

    def is_gc_end(self, seq):
        return (self.sequence.startswith('G') or self.sequence.startswith('C')) and (self.sequence.endswith('G') or self.sequence.endswith('C'))

    def calculate_Tm(self, seq):
        gc_content = int(self.calculate_gc_content(seq)*100)
        mismatch = int(self.calculate_mismatch(seq)*100)
        if self.mutation_type == 'S':
            N = len(self.sequence)
            return 81.5 + 0.41*gc_content - 675/N - mismatch
        else:
            if self.mutation_type == 'I':
                N = len(self.sequence)
            else:
                N = len(self.sequence) - len(self.replacement)
            return 81.5 + 0.41*gc_content - 675/N

    def characterize_primer(self):
        mol_weight = self.lut["mol_weight"]
        complement_dict = self.lut["complement"]
        seq = list(self.sequence)
        rev = self.get_reverse_complement(self.sequence)
        primer_length = len(seq)
        gc_content = self.calculate_gc_content(seq)
        mismatch = self.calculate_mismatch(seq)
        melt_temp = self.calculate_Tm(seq)
        gc_end = self.is_gc_end(seq)
        molweight_fwd = sum(float(mol_weight[b])*2 for b in seq)
        molweight_rev = sum(float(mol_weight[b])*2 for b in rev)
        col = [
            'Length',
            'GC content',
            'Melting temp',
            'Mol. weight (fwd)',
            'Mol. weight (rev)',
            'Mismatch',
            'Ends in G/C'
        ]
        dat = [
            f'{primer_length} bp',
            f'{gc_content*100:.2f}%',
            f'{melt_temp:.2f} C',
            f'{molweight_fwd:.2f} g/mol',
            f'{molweight_rev:.2f} g/mol',
            f'{mismatch*100:.2f}%',
            gc_end
        ]
        print('\n', tabulate(
            array([col, dat]).T,
            headers=['Primer 1'],
            tablefmt='orgtbl'
        ), sep="")
        self.df = DataFrame(
            data=dat,
            columns=['Primer 1'],
            index=col
        )

    def dna_based(self):
        pass

    def protein_based(self):
        pass

    def cut(self):
        #sequence
        for f in range(*self.flank5_range):
            self.cut = self.sequence[self.position-1-f: self.position-1+f]
            SequenceChecks(self.cut)
            calculate_Tm()
        #find primer sequences with mutation w/ min max primer length
        #forward sequence with mutation at the middle
        #return value

    def main(self):
        if self.mode == 'CHAR':
            self.characterize_primer()
        elif self.mode == 'DNA':
            self.dna_based()
        elif self.mode == 'PRO':
            self.protein_based()


class PrimerChecks:
    def __init__(self, sequence):
        self.sequence = sequence

    def check_valid_base(self):
        unique_bases = set(list(self.sequence.upper()))
        true_bases = {'A', 'C', 'T', 'G'}
        invalid_bases = unique_bases.symmetric_difference(true_bases)
        print(invalid_bases)
        if len(invalid_bases) != 0:
            warn("Sequence contains invalid bases. Automatically removing...", Warning)
            for b in invalid_bases:
                self.sequence = self.sequence.upper().replace(b, "")
        else:
            return 0

    def check_sequence_length(self):
        if len(self.sequence) < 40:
            warn('DNA sequence is too short', Warning)
        elif len(self.sequence) > 8000:
            warn('DNA sequence is too long', Warning)

    def check_gc_content(self):
        seq = list(self.sequence)
        gc = (seq.count('C') + seq.count('G'))/len(seq)
        if gc < 0.40:
            warn("GC content is less than 40%", Warning)
        elif gc > 0.60:
            warn('GC content is greater than 60%', Warning)

class SequenceChecks:
    def __init__(self, sequence):
        self.sequence = sequence

    def check_sequence_length(self,length_range):
        if len(self.sequence) < length_range[0] and len(self.sequence) > length_range[1]:
            return False
        else:
            return True

    def check_gc_content(self,gc_range):
        seq = list(self.sequence)
        gc = (seq.count('C') + seq.count('G'))/len(seq)
        if gc < gc_range[0] and gc > gc_range[1]:
            return False
        else:
            return True

    def check_Tm(self,Tm_range):
        seq = list(self.sequence)
        Tm = (seq.count('C') + seq.count('G'))/len(seq)
        if Tm < Tm_range[0] and Tm > Tm_range[1]:
            return False
        else:
            return True
