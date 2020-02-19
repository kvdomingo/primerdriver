from warnings import warn
from json import load
from pandas import DataFrame
from tabulate import tabulate
from numpy import array
from Bio.Seq import Seq


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
        savename=None,
        settings="pdcli/settings.json"
    ):
        with open(settings, "r") as f:
            settings = load(f)
        self.mode = mode.upper()
        self.sequence = sequence.upper()
        self.mutation_type = mutation_type[0].upper()
        self.mismatched_bases = int(mismatched_bases) if mismatched_bases is not None else None
        self.replacement = replacement.upper() if replacement is not None else None
        self.position = int(position) if position is not None else None
        self.target = target.upper() if target is not None else None
        self.Tm_range = (settings["Tm_range_min"], settings["Tm_range_max"])
        self.length_range = (settings["length_min"], settings["length_max"])
        self.gc_range = (settings["gc_range_min"], settings["gc_range_max"])
        self.flank5_range = (settings["flank5_range_min"], settings["flank5_range_max"])
        self.flank3_range = (settings["flank3_range_min"], settings["flank3_range_max"])
        self.terminate_gc = bool(settings["terminate_gc"])
        self.center_mutation = bool(settings["center_mutation"])
        self.primer_mode = settings["primer_mode"]
        self.expression_system = settings["expression_system"]
        self.savename = savename
        self.settings = settings
        with open("pdcli/lut.json", "r", encoding="utf-8") as f:
            lut = load(f)
        with open("pdcli/settings.json", "r") as f:
            settings = load(f)
        self.lut = lut

    def calculate_gc_content(self, seq):
        return (seq.count('G') + seq.count('C'))/len(seq)

    def calculate_mismatch(self, seq, mismatched_bases):
        return mismatched_bases/len(seq)

    def get_reverse_complement(self, seq):
        seq = list(seq)
        return [self.lut["complement"][b] for b in seq][::-1]

    def is_gc_end(self, sequence):
        sequence = ''.join(sequence)
        return (sequence.startswith('G') or sequence.startswith('C')) and (sequence.endswith('G') or sequence.endswith('C'))

    def calculate_Tm(self, seq, mutation_type, replacement, gc_content, mismatch):
        gc_content = int(gc_content*100)
        mismatch = int(mismatch*100)
        if mutation_type == 'S':
            N = len(seq)
            return 81.5 + 0.41*gc_content - 675/N - mismatch
        else:
            if mutation_type == 'I':
                N = len(seq)
            elif mutation_type == 'D':
                N = len(seq)
            else:
                N = len(seq) - len(replacement)
            return 81.5 + 0.41*gc_content - 675/N

    def characterize_primer(self, sequence, mutation_type, replacement, mismatched_bases, index=None, reverse=None):
        mol_weight = self.lut["mol_weight"]
        complement_dict = self.lut["complement"]
        seq = list(sequence)
        if self.primer_mode == "complementary":
            rev = ''.join(self.get_reverse_complement(sequence))
        else:
            rev = reverse
        primer_length = len(seq)
        gc_content = self.calculate_gc_content(seq)
        mismatch = self.calculate_mismatch(seq, mismatched_bases)
        melt_temp = self.calculate_Tm(seq, mutation_type, replacement, gc_content, mismatch)
        gc_end = self.is_gc_end(seq)
        molweight_fwd = sum(float(mol_weight[b])*2 for b in seq)
        molweight_rev = sum(float(mol_weight[b])*2 for b in rev)
        col = [
            'Forward',
            'Reverse',
            'Length',
            'GC content',
            'Melting temp',
            'Mol. weight (fwd)',
            'Mol. weight (rev)',
            'Mismatch',
            'Ends in G/C'
        ]
        dat = [
            sequence,
            rev,
            f'{primer_length} bp',
            f'{gc_content*100:.2f}%',
            f'{melt_temp:.2f} C',
            f'{molweight_fwd:.2f} g/mol',
            f'{molweight_rev:.2f} g/mol',
            f'{mismatch*100:.2f}%',
            gc_end
        ]
        if index == None:
            index = 1
        print('\n', tabulate(
            array([col, dat]).T,
            headers=[f'Primer {index}'],
            tablefmt='orgtbl'
        ), sep="")
        dat = array([dat])
        df = DataFrame(
            data=dat,
            columns=col,
            index=[f"Primer {index}"]
        )
        self.df = df
        return df

    def substitution(self, sequence, mutation_type, target, replacement, start_position, mismatched_bases):
        if self.primer_mode == 'complementary':
            valid_primers = []
            seq = list(sequence)
            seqlen = len(replacement)
            seq[start_position-1] = replacement
            for f5 in range(*self.flank5_range):
                for f3 in range(*self.flank3_range):
                    if abs(f5 - f3) > 1 and self.center_mutation:
                        continue
                    candidate1 = seq[start_position-1-f5 : start_position-1]
                    candidate3 = list(replacement)
                    candidate2 = seq[start_position+seqlen-1 : start_position+seqlen+f3]
                    candidate = candidate1 + candidate3 + candidate2
                    candidate = ''.join(candidate)
                    if len(candidate) == 0:
                        continue
                    gc_content = self.calculate_gc_content(candidate)
                    mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                    Tm = self.calculate_Tm(candidate, mutation_type, replacement, gc_content, mismatch)
                    sc = SequenceChecks(candidate)
                    valid_gc = sc.check_gc_content(self.gc_range)
                    valid_temp = sc.check_Tm(Tm, self.Tm_range)
                    valid_ends = sc.check_ends_gc(self.terminate_gc)
                    valid_length = sc.check_sequence_length(self.length_range)
                    if valid_gc and valid_temp and valid_ends and valid_length:
                        valid_primers.append(candidate)
        else:
            valid_primers = []
            forseq = list(sequence)
            revseq = list(sequence)[::-1]
            seqlen = len(replacement)
            forseq[start_position-1] = replacement
            for f5 in range(*self.flank5_range):
                for f3 in range(*self.flank3_range):
                    if abs(f5 - f3) > 1 and self.center_mutation:
                        continue
                    candidate1 = forseq[start_position-1-f5 : start_position-1]
                    candidate3 = list(replacement)
                    candidate2 = forseq[start_position+seqlen-1 : start_position+seqlen+f3]
                    candidate = candidate1 + candidate3 + candidate2
                    candidate = ''.join(candidate)
                    if len(candidate) == 0:
                        continue
                    gc_content = self.calculate_gc_content(candidate)
                    mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                    Tm = self.calculate_Tm(candidate, mutation_type, replacement, gc_content, mismatch)
                    sc = SequenceChecks(candidate)
                    valid_gc = sc.check_gc_content(self.gc_range)
                    valid_temp = sc.check_Tm(Tm,self.Tm_range)
                    valid_ends = sc.check_ends_gc(self.terminate_gc)
                    valid_length = sc.check_sequence_length(self.length_range)
                    if valid_gc and valid_temp and valid_ends and valid_length:
                        valid_primers.append(candidate)
            valid_reverse = []
            for primers in valid_primers:
                sequence = sequence[:start_position-1] + replacement + sequence[start_position+seqlen:]
                start = sequence.find(primers)
                end = start + len(primers)-1
                prilen = len(primers)
                while start < self.position-7 and end > self.position+seqlen+7:
                    start = start-1
                    end = end-1

                while start < self.position-7:
                    candidate = sequence[start:prilen]
                    candidate = [self.lut["complement"][b] for b in candidate]
                    if len(candidate) == 0:
                        continue
                    gc_content = self.calculate_gc_content(candidate)
                    mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                    Tm = self.calculate_Tm(candidate, mutation_type, replacement, gc_content, mismatch)
                    sc = SequenceChecks(primers)
                    valid_gc = sc.check_gc_content(self.gc_range)
                    valid_temp = sc.check_Tm(Tm,self.Tm_range)
                    valid_ends = sc.check_ends_gc(self.terminate_gc)
                    valid_length = sc.check_sequence_length(self.length_range)
                    if valid_gc and valid_temp and valid_ends and valid_length:
                        valid_reverse.append(primers)
                    start = start + 1

        if not len(valid_primers) > 0:
            print("No valid primers found")
            return
        else:
            df = []
            print(f"\nGenerated forward primers: {len(valid_primers)}")
            if self.primer_mode == "complementary":
                for i, p in enumerate(valid_primers):
                    df.append(self.characterize_primer(p, mutation_type, replacement, mismatched_bases, i+1))
            else:
                if len(valid_reverse) == 0:
                    print('No valid reverse primers found')
                    return
                for i, p in enumerate(valid_primers):
                    for j, rev in enumerate(valid_reverse):
                        for r in rev:
                            df.append(self.characterize_primer(p, mutation_type, replacement, mismatched_bases, i+1, r))
        return df

    def deletion(self, sequence, mutation_type, target, replacement, start_position, mismatched_bases):
        valid_primers = []
        seqlen = len(self.target)
        seq = list(sequence)
        for f5 in range(*self.flank5_range):
            for f3 in range(*self.flank3_range):
                if abs(f5 - f3) > 1 and self.center_mutation:
                    continue
                candidate1 = seq[start_position-1-f5 : start_position-1]
                candidate2 = seq[start_position+seqlen-1 : start_position+seqlen+f3]
                candidate = candidate1 + candidate2
                candidate = ''.join(candidate)
                if len(candidate) == 0:
                    continue
                sc = SequenceChecks(candidate)
                valid_gc = sc.check_gc_content(self.gc_range)
                valid_temp = sc.check_Tm(self.Tm_range)
                valid_ends = sc.check_ends_gc(self.terminate_gc)
                valid_length = sc.check_sequence_length(self.length_range)
                if valid_gc and valid_temp and valid_ends and valid_length:
                    valid_primers.append(candidate)
        if len(valid_primers) > 0:
            df = []
            print(f"\nGenerated primers: {len(valid_primers)}")
            for i, p in enumerate(valid_primers):
                df.append(self.characterize_primer(p, mutation_type, replacement, mismatched_bases, i+1))
        else:
            print("No valid primers found")
        return df

    def insertion(self, sequence, mutation_type, target, replacement, start_position, mismatched_bases):
        valid_primers = []
        seq = list(sequence)
        seqlen = len(replacement)
        seq[start_position-1] = replacement
        for f5 in range(*self.flank5_range):
            for f3 in range(*self.flank3_range):
                if abs(f5 - f3) > 1 and self.center_mutation:
                    continue
                candidate1 = seq[start_position-1-f5 : start_position-1]
                candidate3 = list(replacement)
                candidate2 = seq[start_position-1 : start_position+seqlen+f3]
                candidate = candidate1 + candidate3 + candidate2
                candidate = ''.join(candidate)
                if len(candidate) == 0:
                    continue
                sc = SequenceChecks(candidate)
                valid_gc = sc.check_gc_content(self.gc_range)
                valid_temp = sc.check_Tm(self.Tm_range)
                valid_ends = sc.check_ends_gc(self.terminate_gc)
                valid_length = sc.check_sequence_length(self.length_range)
                if valid_gc and valid_temp and valid_ends and valid_length:
                    valid_primers.append(candidate)
        if len(valid_primers) > 0:
            df = []
            print(f"\nGenerated primers: {len(valid_primers)}")
            for i, p in enumerate(valid_primers):
                df.append(self.characterize_primer(p, mutation_type, replacement, mismatched_bases, i+1))
        else:
            print("No valid primers found")
        return df

    def dna_based(self):
        if self.mutation_type in ['S', 'SUB']:
            if self.mismatched_bases is None:
                self.mismatched_bases = len(self.replacement)
            result = self.substitution(self.sequence, self.mutation_type, self.target, self.replacement, self.position, self.mismatched_bases)
        elif self.mutation_type in ['D', 'DEL']:
            if self.mismatched_bases is None:
                self.mismatched_bases = len(self.target)
            result = self.deletion(self.sequence, self.mutation_type, self.target, self.replacement, self.position, self.mismatched_bases)
        elif self.mutation_type in ['I', 'INS']:
            if self.mismatched_bases is None:
                self.mismatched_bases = len(self.replacement)
            result = self.insertion(self.sequence, self.mutation_type, self.target, self.replacement, self.position, self.mismatched_bases)
        return result

    def protein_based(self):
        if self.mutation_type in ['S', 'SUB']:
            if self.mismatched_bases is None:
                self.mismatched_bases = len(self.replacement)
            with open("pdcli/AAcompressed.json", "r") as f:
                cod1 = load(f)
            with open("pdcli/symnucbases.json", "r") as f:
                sym = load(f)
            rna = ''.join([cod1[b][0] for b in self.sequence])
            rna = ''.join([sym[b][0] for b in rna])
            dna = rna.replace('U', 'T')
            target = cod1[self.target][0]
            target = ''.join([sym[b][0] for b in target]).replace('U', 'T')
            replacement = cod1[self.replacement][0]
            replacement = ''.join([sym[b][0] for b in replacement]).replace('U', 'T')
            result = self.substitution(dna, self.mutation_type, target, replacement, self.position*3, self.mismatched_bases)
            return result

    def main(self):
        if self.mode == 'CHAR':
            df = self.characterize_primer(self.sequence, self.mutation_type, self.replacement, self.mismatched_bases)
        elif self.mode == 'DNA':
            df = self.dna_based()
        elif self.mode == 'PRO':
            df = self.protein_based()
        self.df = df


class PrimerChecks:
    def __init__(self, sequence):
        self.sequence = sequence

    def check_valid_base(self):
        unique_bases = set(list(self.sequence.upper()))
        true_bases = {'A', 'C', 'T', 'G'}
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
        self.sequence = sequence.upper()

    def check_sequence_length(self, length_range):
        if len(self.sequence) >= length_range[0] and len(self.sequence) <= length_range[1]:
            return True
        else:
            return False

    def check_gc_content(self, gc_range):
        seq = list(self.sequence)
        gc = (seq.count('C') + seq.count('G'))/len(seq) * 100
        if gc < gc_range[0] or gc > gc_range[1]:
            return False
        else:
            return True

    def check_Tm(self, Tm, Tm_range):
        if Tm < Tm_range[0] or Tm > Tm_range[1]:
            return False
        else:
            return True

    def check_ends_gc(self, terminate_gc):
        if not terminate_gc or (self.sequence[0] in ['C', 'G'] and self.sequence[-1] in ['C', 'G']):
            return True
        else:
            return False
