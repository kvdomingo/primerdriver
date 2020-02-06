import os
from datetime import datetime
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from pdcli.primerclass import *
from pdcli.input_handler import *


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

    if args.interactive:
        args_dict = interactive_handler()
    else:
        args_dict = singleCommand_handler(args)

    res = PrimerDesign(**args_dict)
    if args_dict['mode'].upper() == 'CHAR':
        res.characterize_primer()
        df = res.df
    elif args_dict['mode'].upper() == 'DNA':
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
    else:
        raise NotImplementedError

    save = input("\nSave? [y/n] ")
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
