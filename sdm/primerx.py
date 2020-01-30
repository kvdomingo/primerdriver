from argparse import ArgumentParser

def main():
    parser = ArgumentParser()
    parser.add_argument('-s', '--sequence', help='Template DNA sequence', required=True, type=str)
    parser.add_argument('-m', '--mutation-type', help='Mutation type', required=True, type=str)
    parser.add_argument('-t', '--target', help='Target base', required=True, type=str)
    parser.add_argument('-d', '--destination', help='Replacement for target base', required=True, type=str)
    parser.add_argument('-p', '--position', help='Target base position', required=True, type=int)
    args = parser.parse_args()

    if args.mutation_type.lower() == 'sub':
        out = list(args.sequence)
        if args.target != args.sequence[args.position-1]:
            raise ValueError('Sequence position does not math target base')
        out[args.position-1] = args.destination
        print(''.join(out))
        return 0


if __name__ == '__main__':
    main()
