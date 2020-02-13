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
        self.mol_weight = lut["mol_weight"]
        self.complement_dict = lut["complement"]

    def __len__(self):
        return len(self.forward)

    def calculate_Tm(self):
        gc_content = int(self.gc_content*100)
        mismatch = int(self.mismatch*100)
        if self.mutation_type == 'S':
            N = len(self.sequence)
            self.melt_temp = 81.5 + 0.41*gc_content - 675/N - mismatch
        else:
            if self.mutation_type == 'I':
                N = len(self.sequence)
            else:
                N = len(self.sequence) - len(self.replacement)
            self.melt_temp = 81.5 + 0.41*gc_content - 675/N

    def substitution(self):
        seq = list(self.sequence)
        if self.target != seq[self.position-1]:
            raise ValueError('Sequence position does not match target base')
        seq[self.position-1] = self.replacement
        self.forward = ''.join(seq)

    def insertion(self):
        if not self.replacement:
            raise RuntimeError('Insertion mutation requires replacement')
        seq = list(self.sequence)
        seq[self.position-1:self.position-1] = self.replacement
        self.forward = ''.join(seq)

    def deletion(self):
        seq = list(self.sequence)
        del seq[self.position-1]
        self.forward = ''.join(seq)

    def characterize_primer(self):
        mol_weight = self.lut["mol_weight"]
        complement_dict = self.lut["complement"]
        seq = list(self.sequence)
        rev = [complement_dict[b] for b in seq][::-1]
        primer_length = len(seq)
        self.gc_content = (seq.count('G') + seq.count('C'))/len(seq)
        self.mismatch = self.mismatched_bases/len(seq)
        self.calculate_Tm()
        gc_end = (self.sequence.startswith('G') or self.sequence.startswith('C')) and (self.sequence.endswith('G') or self.sequence.endswith('C'))
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
            f'{self.gc_content*100:.2f}%',
            f'{self.melt_temp:.2f} C',
            f'{molweight_fwd:.2f} g/mol',
            f'{molweight_rev:.2f} g/mol',
            f'{self.mismatch*100:.2f}%',
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

    def main(self):
        if self.mutation_type == 'S':
            self.substitution()
        elif self.mutation_type == 'I':
            self.insertion()
        elif self.mutation_type == 'D':
            self.deletion()
        else:
            raise NotImplementedError('Invalid mutation type')

        with open("pdcli/lut.json", "r", encoding="utf-8") as f:
            complement_dict = load(f)["complement"]
        self.rev_compl = ''.join([complement_dict[b] for b in list(self.forward[::-1])])
        self.gc_content = (self.forward.count('G') + self.forward.count('C'))/len(self.forward)
        self.mismatch = len(self.replacement)/len(self.forward)

        if self.mutation_type in ['S', 'SUB']:
            self.calculate_Tm()
        else:
            self.Tm_insdel()

        molweight_fwd = sum(float(self.mol_weight[b])*2 for b in self.forward)
        molweight_rev = sum(float(self.mol_weight[b])*2 for b in self.rev_compl)
        if self.mode not in ['C', 'CHAR']:
            col = [
                'Forward',
                'Reverse',
                'GC content',
                'Mol. weight (fwd)',
                'Mol. weight (rev)',
                'Melting temp.',
                'Length'
            ]
            dat = [
                f'{self.forward}',
                f'{self.rev_compl}',
                f'{self.gc_content*100:.2f}%',
                f'{molweight_fwd:.2f} g/mol',
                f'{molweight_rev:.2f} g/mol',
                f'{self.melt_temp:.2f} C',
                f'{len(self.forward)} bp'
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
