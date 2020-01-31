from argparse import ArgumentParser
from pandas import DataFrame


bases_dict = {
    'T': 'thymine',
    'A': 'adenine',
    'G': 'guanine',
    'C': 'cytosine'
}

complement_dict = {
    'T': 'A',
    'A': 'T',
    'G': 'C',
    'C': 'G'
}

protein_dict = {
    'A': 'alanine',
    'R': 'arginine',
    'N': 'asparagine',
    'D': 'aspartic acid',
    'C': 'cysteine',
    'Q': 'glutamine',
    'E': 'glutamic acid',
    'G': 'glycine',
    'H': 'histidine',
    'I': 'isoleucine',
    'L': 'leucine',
    'K': 'lysine',
    'M': 'methionine',
    'F': 'phenylalanine',
    'P': 'proline',
    'S': 'serine',
    'T': 'threonine',
    'W': 'tryptophan',
    'Y': 'tyrosine',
    'V': 'valine'
}


class PrimerDesign:
    def __init__(self, mode, sequence, mutation_type, destination, position, target=None):
        self.mode = mode
        self.sequence = sequence
        self.mutation_type = mutation_type.lower()
        self.target = target
        self.destination = destination
        self.position = position

        self.checks()

        if self.mutation_type in ['s', 'sub']:
            self.substitution()
        elif self.mutation_type in ['i', 'ins']:
            self.insertion()
        elif self.mutation_type in ['d', 'deletion']:
            self.deletion()
        else:
            raise NotImplementedError

    def checks(self):
        if self.mode.lower() == 'dna':
            unique_bases = set(list(self.sequence))
            true_bases = {'A', 'C', 'T', 'G'}
            if len(unique_bases.union(true_bases)) != 4:
                raise ValueError('DNA sequence contains invalid bases.')

    def substitution(self):
        seq = list(self.sequence)
        if self.target != seq[self.position-1]:
            raise ValueError('Sequence position does not match target base')
        seq[self.position-1] = self.destination
        self.out = ''.join(seq)

    def insertion(self):
        if not self.destination:
            raise RuntimeError('Insertion mutation requires destination')
        seq = list(self.sequence)
        seq[self.position-1:self.position-1] = self.destination
        self.out = ''.join(seq)

    def deletion(self):
        seq = list(self.sequence)
        del seq[self.position-1]
        self.out = ''.join(seq)


def main():
    parser = ArgumentParser()
    parser.add_argument('-s', '--sequence', help='Template DNA sequence', type=str)
    parser.add_argument('-m', '--mutation-type', help='Mutation type', type=str)
    parser.add_argument('-t', '--target', help='Target base', type=str)
    parser.add_argument('-d', '--destination', help='Replacement for target base', type=str)
    parser.add_argument('-p', '--position', help='Target base position', type=int)
    parser.add_argument('--mode', help="Choose between `dna` or `pro` mode", type=str)
    parser.add_argument('-i', '--interactive', help='Interactive mode', action='store_true')
    args = parser.parse_args()
    args_dict = dict()

    if args.interactive:
        print('')
        print('===================================')
        print('====                           ====')
        print('====      PrimerX v0.1.1       ====')
        print('====                           ====')
        print('===================================\n')
        args_dict['mode'] = input('Enter primer mode [dna/protein]: ')
        args_dict['sequence'] = input('Enter DNA sequence: ')
        args_dict['mutation_type'] = input('Enter mutation type [s/i/d]: ')
        if args_dict['mutation_type'].lower() in ['s', 'sub']:
            args_dict['target'] = input('Enter target base: ')
        else:
            args_dict['target'] = None
        args_dict['destination'] = input('Enter replacement for target base: ')
        args_dict['position'] = int(input('Enter position of target: '))
    else:
        args_dict['mode'] = args.mode
        args_dict['sequence'] = args.sequence
        args_dict['mutation_type'] = args.mutation_type
        if args_dict['mutation_type'].lower() in ['s', 'sub']:
            args_dict['target'] = args.target
        else:
            args_dict['target'] = None
        args_dict['destination'] = args.destination
        args_dict['position'] = args.position

    res = PrimerDesign(**args_dict)
    gccont = (res.out.count('G') + res.out.count('C'))/len(res.out)
    idx = ['Forward', 'Reverse', 'GC content', 'Length']
    revcomp = ''.join([complement_dict[b] for b in list(res.out[::-1])])
    dat = [res.out, revcomp, f'{gccont*100:.2f}%', f'{len(res.out)} bp']
    df = DataFrame(
        data=dat,
        index=idx,
        dtype=str,
    )
    print(df)
    return 0


if __name__ == '__main__':
    main()
