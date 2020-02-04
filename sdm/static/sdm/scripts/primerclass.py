from warnings import warn
from json import load


class PrimerDesign:
    def __init__(self, mode, sequence, mutation_type, destination, position, target=None):
        self.mode = mode.upper()
        self.sequence = sequence.upper()
        self.mutation_type = mutation_type.upper()
        self.destination = destination.upper()
        self.position = int(position)
        if target is not None:
            self.target = target.upper()

    def __len__(self):
        return len(self.forward)

    def substitution(self):
        seq = list(self.sequence)
        if self.target != seq[self.position-1]:
            raise ValueError('Sequence position does not match target base')
        seq[self.position-1] = self.destination
        self.forward = ''.join(seq)

    def insertion(self):
        if not self.destination:
            raise RuntimeError('Insertion mutation requires destination')
        seq = list(self.sequence)
        seq[self.position-1:self.position-1] = self.destination
        self.forward = ''.join(seq)

    def deletion(self):
        seq = list(self.sequence)
        del seq[self.position-1]
        self.forward = ''.join(seq)

    def Tm_subsitution(self):
        Tm = 81.5 + 0.41*int(self.gc_content*100) - 675/len(self.forward) - int(self.mismatch*100)
        self.melt_temp = Tm

    def Tm_insdel(self):
        Tm = 81.5 + 0.41*int(self.gc_content*100) - 675/len(self.forward)
        self.melt_temp = Tm

    def main(self):
        if self.mutation_type in ['S', 'SUB']:
            self.substitution()
        elif self.mutation_type in ['I', 'INS']:
            self.insertion()
        elif self.mutation_type in ['D', 'DEL']:
            self.deletion()
        else:
            raise NotImplementedError

        with open("lut.json", "r", encoding="utf-8") as f:
            complement_dict = load(f)["complement"]
        self.rev_compl = ''.join([complement_dict[b] for b in list(self.forward[::-1])])
        self.gc_content = (self.forward.count('G') + self.forward.count('C'))/len(self.forward)
        self.mismatch = len(self.destination)/len(self.forward)

        if self.mutation_type in ['S', 'SUB']:
            self.Tm_subsitution()
        else:
            self.Tm_insdel()


class PrimerChecks:
    def __init__(self, sequence):
        self.sequence = sequence

    def check_valid_base(self):
        unique_bases = set(list(self.sequence.upper()))
        true_bases = {'A', 'C', 'T', 'G'}
        if len(unique_bases.union(true_bases)) != 4:
            raise ValueError('DNA sequence contains invalid bases')
        else:
            return 0

    def check_sequence_length(self):
        if len(self.sequence) < 40:
            raise ValueError('DNA sequence is too short')
        elif len(self.sequence) > 8000:
            raise ValueError('DNA sequence is too long')
        else:
            return 0

    def check_gc_content(self):
        seq = list(self.sequence)
        gc = (seq.count('C') + seq.count('G'))/len(seq)
        if gc < 0.40:
            warn("GC content is less than 40%", Warning)
        elif gc > 0.60:
            warn('GC content is greater than 60%', Warning)
        else:
            return 0
