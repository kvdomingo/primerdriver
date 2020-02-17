import os
from datetime import datetime
from argparse import ArgumentParser
from pdcli.input_handler import *
from pdcli.output_handler import *


def main():
    parser = ArgumentParser()
    parser.add_argument('--mode', help="Choose between 'dna' (DNA), 'pro' (protein), or 'char' (primer characterization) mode", type=str)
    parser.add_argument('-s', '--sequence', help='Template DNA sequence', type=str)
    parser.add_argument('-m', '--mutation-type', help='Mutation type', type=str)
    parser.add_argument('-t', '--target', help='Target base', type=str)
    parser.add_argument('-r', '--replacement', help='Replacement for target base', type=str)
    parser.add_argument('-p', '--position', help='Target base position', type=int)
    parser.add_argument('-i', '--interactive', help='Interactive mode', action='store_true')
    parser.add_argument('--save', help='Filename to save to', type=str)
    args = parser.parse_args()

    if args.interactive:
        args_dict = interactive_handler()
    else:
        args_dict = singleCommand_handler(args)

    res = PrimerDesign(**args_dict)
    if args_dict['mode'].upper() == 'CHAR':
        res.characterize_primer(
            args_dict['sequence'],
            args_dict['mutation_type'],
            args_dict['mismatched_bases']
        )
        df = res.df
    elif args_dict['mode'].upper() == 'DNA':
        res.main()
        # idx = [
        #     'Forward',
        #     'Reverse',
        #     'GC content',
        #     'Melting temp',
        #     'Length'
        # ]
        # dat = [
        #     res.forward, res.rev_compl,
        #     f'{res.gc_content*100:.2f}%',
        #     f'{res.melt_temp:.2f} C',
        #     f'{len(res)} bp'
        # ]
        # df = DataFrame(
        #     data=dat,
        #     index=idx,
        #     dtype=str,
        # )
    else:
        raise NotImplementedError

    # if 'savename' in args_dict.keys():
    #     savename = args_dict['savename']
    # else:
    #     savename = None
    #
    # if args.interactive:
    #     interactive_saver(res, df)
    # else:
    #     if savename is not None:
    #         singleCommand_saver(res, df, savename)

    return 0


if __name__ == '__main__':
    main()
