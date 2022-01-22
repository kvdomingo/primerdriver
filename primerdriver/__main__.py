import sys
from argparse import ArgumentParser
from .input_handler import *
from .output_handler import *


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "-M",
        "--mode",
        help="Choose between 'dna' (DNA), 'pro' (protein), or 'char' (primer characterization) mode",
        type=str,
    )
    parser.add_argument("-s", "--sequence", help="Template DNA sequence", type=str)
    parser.add_argument("-m", "--mutation-type", help="Mutation type", type=str)
    parser.add_argument("-t", "--target", help="Target base", type=str)
    parser.add_argument("-r", "--replacement", help="Replacement for target base", type=str)
    parser.add_argument("-p", "--position", help="Target base position", type=int)
    parser.add_argument("-i", "--interactive", help="Interactive mode", action="store_true")
    parser.add_argument("--save", help="Filename to save to", type=str)

    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        return 0

    args = parser.parse_args()

    if args.interactive:
        args_dict = interactive_handler()
    else:
        args_dict = singleCommand_handler(args)

    res = PrimerDesign(**args_dict)
    if args_dict["mode"].upper() == "CHAR":
        res.characterize_primer(args_dict["sequence"], args_dict["mutation_type"], [1], args_dict["mismatched_bases"])
    else:
        res.main()

    if args.interactive:
        interactive_saver(res.df)
    else:
        if res.savename is not None:
            singleCommand_saver(res, res.df, res.savename)

    return 0


if __name__ == "__main__":
    main()
