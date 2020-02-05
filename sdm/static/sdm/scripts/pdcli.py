import os
from datetime import datetime
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from primerclass import *


def main():
    parser = ArgumentParser()
    parser.add_argument('--mode', help="Choose between 'dna' (DNA), 'pro' (protein), or 'char' (primer characterization) mode", type=str)
    parser.add_argument('-s', '--sequence', help='Template DNA sequence', type=str)
    parser.add_argument('-m', '--mutation-type', help='Mutation type', type=str)
    parser.add_argument('-t', '--target', help='Target base', type=str)
    parser.add_argument('-d', '--destination', help='Replacement for target base', type=str)
    parser.add_argument('-p', '--position', help='Target base position', type=int)
    parser.add_argument('-i', '--interactive', help='Interactive mode', action='store_true')
    args = parser.parse_args()
    args_dict = dict()

    if args.interactive:
        print('')
        print('====================================')
        print('======                        ======')
        print('===     PrimerDriver v0.1.2      ===')
        print('======                        ======')
        print('====================================\n')
        print('(c) 2020 Kenneth Domingo & Nomer Gutierrez\n')
        args_dict['mode'] = input('Enter primer mode [dna/pro/char]: ')
        if args_dict['mode'].upper() == 'DNA':
            args_dict['sequence'] = input('Enter DNA sequence: ')
            PrimerChecks(args_dict['sequence']).check_sequence_length()
            PrimerChecks(args_dict['sequence']).check_valid_base()
            args_dict['mutation_type'] = input('Enter mutation type [s/i/d]: ')
            if args_dict['mutation_type'].upper() in ['S', 'SUB']:
                args_dict['target'] = input('Enter target base: ')
                args_dict['destination'] = input('Enter replacement for target base: ')
                args_dict['position'] = int(input('Enter position of target: '))
            elif args_dict['mutation_type'].upper() in ['I', 'INS']:
                args_dict['target'] = None
                args_dict['destination'] = input('Enter insertion sequence: ')
                args_dict['position'] = int(input('Enter insertion position: '))
            elif args_dict['mutation_type'].upper() in ['I', 'INS']:
                args_dict['target'] = None
                args_dict['destination'] = input('Enter starting position to delete: ')
            else:
                raise ValueError("Invalid argument passed to 'MUTATION_TYPE'")
        elif args_dict['mode'].upper() == 'CHAR':
            args_dict['sequence'] = input('Enter primer sequence: ')
            args_dict['mutation_type'] = input('Enter mutation type [s/i/d]: ')
            args_dict['mismatched_bases'] = input('Enter number of mismatched bases: ')
        else:
            raise NotImplementedError(f"{args_dict['mode']} mode not implemented (yet).")

    else:
        args_dict['mode'] = args.mode
        if args.mode.upper() =='DNA':
            args_dict['sequence'] = args.sequence
            PrimerChecks(args.sequence).check_sequence_length()
            PrimerChecks(args.sequence).check_valid_base()
            args_dict['mutation_type'] = args.mutation_type
            args_dict['position'] = args.position
            args_dict['destination'] = args.destination
            if args.mutation_type.upper() in ['S', 'SUB']:
                args_dict['target'] = args.target
            else:
                raise ValueError("Invalid argument passed to 'MUTATION_TYPE'")
        else:
            raise NotImplementedError(f"{args_dict['mode']} mode not implemented (yet).")

    res = PrimerDesign(**args_dict)
    if args_dict['mode'].upper() == 'CHAR':
        res.characterize_primer()
        df = res.df
    else:
        res.main()
        idx = [
            'Forward',
            'Reverse',
            'GC content',
            'Melting temp',
            'Length'
        ]
        dat = [
            res.forward, res.rev_compl,
            f'{res.gc_content*100:.2f}%',
            f'{res.melt_temp:.2f} C',
            f'{len(res)} bp'
        ]
        df = DataFrame(
            data=dat,
            index=idx,
            dtype=str,
        )
    save = input("Save? [y/n] ")
    if save.upper() == "Y":
        while True:
            savename = input("Enter filename: ")
            if savename.endswith(".csv"):
                df.to_csv(savename)
                break
            elif savename.endswith(".html"):
                df.to_html(savename)
                break
            elif savename.endswith(".fasta"):
                seq_forward = SeqRecord(
                    Seq(res.forward, IUPAC.unambiguous_dna),
                    id=f"{datetime.now()}",
                    description=f"{''.join(savename.split('.')[:-1])} forward strand, {res.gc_content*100:.2f}% G/C, Tm = {res.melt_temp:.2f} C, {len(res)} bp"
                )
                seq_reverse = SeqRecord(
                    Seq(res.rev_compl, IUPAC.unambiguous_dna),
                    id=f"{datetime.now()}",
                    description=f"{''.join(savename.split('.')[:-1])} reverse strand, {res.gc_content*100:.2f}% G/C, Tm = {res.melt_temp:.2f} C, {len(res)} bp"
                )
                seq_primer = [seq_forward, seq_reverse]
                with open(savename, "w") as f:
                    SeqIO.write(seq_primer, f, "fasta")
                break
            else:
                print("Unsupported filetype. Supported filetypes are: .csv, .html, .fasta")
    return 0


if __name__ == '__main__':
    main()
