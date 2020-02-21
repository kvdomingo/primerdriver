from warnings import warn
from json import load
from pandas import DataFrame
from tabulate import tabulate
from numpy import array
from Bio.Seq import Seq
from .checks import *


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
        self.forward_overlap5 = settings["forward_overlap5"]
        self.forward_overlap3 = settings["forward_overlap3"]
        self.terminate_gc = bool(settings["terminate_gc"])
        self.center_mutation = bool(settings["center_mutation"])
        self.primer_mode = settings["primer_mode"]
        self.expression_system = settings["expression_system"]
        self.savename = savename
        self.settings = settings
        with open("pdcli/lut.json", "r", encoding="utf-8") as f:
            lut = load(f)
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
        if self.primer_mode == "complementary" or self.mode == "CHAR":
            rev = ''.join(self.get_reverse_complement(sequence))
        else:
            rev = reverse

        forwards = {}
        forwards["sequence"] = sequence
        forwards["length"] = len(seq)
        forwards["gc_content"] = self.calculate_gc_content(seq)
        forwards["mismatch"] = self.calculate_mismatch(seq, mismatched_bases)
        forwards["Tm"] = self.calculate_Tm(seq, mutation_type, replacement, forwards["gc_content"], forwards["mismatch"])
        forwards["gc_end"] = self.is_gc_end(seq)
        forwards["mol_weight"] = sum(float(mol_weight[b]) for b in seq)
        
        reverses = {}
        reverses["sequence"] = rev
        reverses["length"] = len(rev)
        reverses["gc_content"] = self.calculate_gc_content(rev)
        reverses["mismatch"] = self.calculate_mismatch(rev, mismatched_bases)
        reverses["Tm"] = self.calculate_Tm(seq, mutation_type, replacement, reverses["gc_content"], reverses["mismatch"])
        reverses["gc_end"] = self.is_gc_end(rev)
        reverses["mol_weight"] = sum(float(mol_weight[b]) for b in rev)

        col = [
            'Forward',
            'Reverse',
            'Fwd length',
            'Rev length',
            'Fwd GC content',
            'Rev GC content',
            'Fwd melting temp',
            'Rev melting temp',
            'Fwd mol. weight',
            'Rev mol. weight',
            'Fwd mismatch',
            'Rev mismatch',
            'Fwd ends in G/C',
            'Rev ends in G/C'
        ]
        dat = [
            forwards["sequence"],
            reverses["sequence"],
            f'{forwards["length"]} bp',
            f'{reverses["length"]} bp',
            f'{forwards["gc_content"]*100:.2f}%',
            f'{reverses["gc_content"]*100:.2f}%',
            f'{forwards["Tm"]:.2f} C',
            f'{reverses["Tm"]:.2f} C',
            f'{forwards["mol_weight"]:.2f} g/mol',
            f'{reverses["mol_weight"]:.2f} g/mol',
            f'{forwards["mismatch"]*100:.2f}%',
            f'{reverses["mismatch"]*100:.2f}%',
            forwards["gc_end"],
            reverses["gc_end"],
        ]
        if index == None:
            index = 1
        if index <= 10:
            print('\n', tabulate(
                array([col, dat]).T,
                headers=[f'Primer {index}'],
                tablefmt='orgtbl'
            ), sep="")
        elif index == 11:
            print("Too many results; truncating output...")
        dat = array([dat])
        df = DataFrame(
            data=dat,
            columns=col,
            index=[index]
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
            print(self.sequence, len(self.sequence))
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
                    valid_temp = sc.check_Tm(Tm, self.Tm_range)
                    valid_ends = sc.check_ends_gc(self.terminate_gc)
                    valid_length = sc.check_sequence_length(self.length_range)
                    if valid_gc and valid_temp and valid_ends and valid_length:
                        valid_primers.append(candidate)
            valid_reverse = []
            for primers in valid_primers:
                sequence = sequence[:start_position-1] + replacement + sequence[start_position-1+seqlen:]
                start = sequence.find(primers)
                end = start + len(primers)-1
                prilen = len(primers)
                while start < self.position-self.forward_overlap5 and end > self.position+seqlen+self.forward_overlap3:
                    start -= 1
                    end -= 1
                for i in range(self.position-self.forward_overlap5, end):
                    for j in range(self.flank3_range[1]):
                        candidate = sequence[start-j:end]
                        if len(candidate) == 0:
                            continue
                        gc_content = self.calculate_gc_content(candidate)
                        mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                        fwd_Tm = Tm
                        Tm = self.calculate_Tm(candidate, mutation_type, replacement, gc_content, mismatch)
                        sc = SequenceChecks(candidate)
                        valid_gc = sc.check_gc_content(self.gc_range)
                        valid_temp = sc.check_Tm(Tm, self.Tm_range)
                        valid_trange = sc.check_close_Tm(fwd_Tm, Tm)
                        valid_ends = sc.check_ends_gc(self.terminate_gc)
                        valid_length = sc.check_sequence_length(self.length_range)
                        if valid_gc and valid_temp and valid_ends and valid_length and valid_trange:
                            valid_reverse.append(candidate)
                    end += 1

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
                count = 1
                for i, p in enumerate(valid_primers):
                    for j, r in enumerate(valid_reverse):
                        df.append(self.characterize_primer(p, mutation_type, replacement, mismatched_bases, count, r))
                        count += 1
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