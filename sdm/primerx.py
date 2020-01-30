from argparse import ArgumentParser


def main():
    parser = ArgumentParser()
    parser.add_argument('-s', '--sequence', help='Template DNA sequence', type=str)
    parser.add_argument('-m', '--mutation-type', help='Mutation type', type=str)
    parser.add_argument('-t', '--target', help='Target base', type=str)
    parser.add_argument('-d', '--destination', help='Replacement for target base', type=str)
    parser.add_argument('-p', '--position', help='Target base position', type=int)
    parser.add_argument('-i', '--interactive', help='Interactive mode', action='store_true')
    args = parser.parse_args()

    if args.interactive:
        print('')
        print('===================================')
        print('====                           ====')
        print('====      PrimerX v0.1.0       ====')
        print('====                           ====')
        print('===================================\n')
        sequence = input('Enter DNA sequence: ')
        mutation_type = input('Enter mutation type [s/i/d]: ')
        if mutation_type.lower() == 's':
            target = input('Enter target base: ')
        destination = input('Enter replacement for target base: ')
        position = int(input('Enter position of target: '))

        if mutation_type.lower() == 's':
            out = list(sequence)
            if target != out[position-1]:
                raise ValueError('Sequence position does not match target base')
            out[position-1] = destination
            print(f"\nResult: {''.join(out)}")
            return 0

        return 0

    if args.mutation_type.lower() in ['sub', 's']:
        if not (args.target and args.destination):
            raise RuntimeError('Substitution mutation requires target and destination')
        out = list(args.sequence)
        if args.target != args.sequence[args.position-1]:
            raise ValueError('Sequence position does not match target base')
        out[args.position-1] = args.destination
        print(''.join(out))
        return 0
    elif args.mutation_type.lower() in ['ins', 'i']:
        if not args.destination:
            raise RuntimeError('Insertion mutation requires destination')
        out = list(args.sequence)
        out[args.position-1:args.position-1] = args.destination
        print(''.join(out))
        return 0
    elif args.mutation_type.lower() in ['del', 'd']:
        out = list(args.sequence)
        del out[args.position-1]
        print(''.join(out))
        return 0
    else:
        raise NotImplementedError


if __name__ == '__main__':
    main()
