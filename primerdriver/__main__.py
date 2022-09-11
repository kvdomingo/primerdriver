import sys
from argparse import ArgumentParser

from primerx.log import logger

from .input_handler import interactive_handler, single_command_handler
from .output_handler import interactive_saver, single_command_saver
from .primer_design import MutationType, OperationMode, PrimerDesign
from .version import __version__


@logger.catch
def main():
    print(
        f"""
        ---.   .------------.
        ||||\\ /||||||||||||||\\    
       Primer · Driver
\\|||||||||||/ \\|||||||              
 `---------`   `------
     
PrimerDriver v{__version__}
(c) 2020 Kenneth V. Domingo & Numeriano Amer E. Gutierrez
    """
    )

    parser = ArgumentParser(prog="primerdriver")
    parser.add_argument(
        "-M",
        "--mode",
        help="Operation mode. Choose between 'dna' (DNA), 'pro' (protein), or 'char' (primer characterization) mode",
        type=str,
    )
    parser.add_argument("-s", "--sequence", help="Template DNA sequence", type=str)
    parser.add_argument(
        "-m",
        "--mutation-type",
        help="Mutation type. Choose between 'S' (substitution), 'I' (insertion), or 'D' (deletion)",
        type=str,
    )
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
        args_dict = single_command_handler(args)

    res = PrimerDesign(**args_dict)
    if OperationMode(args_dict["mode"].upper()) == OperationMode.CHARACTERIZATION:
        res.characterize_primer(
            sequence=args_dict["sequence"],
            mutation_type=MutationType(args_dict["mutation_type"].upper()),
            mismatched_bases=args_dict["mismatched_bases"],
            replacement=None,
        )
    else:
        res.main()

    if args.interactive:
        interactive_saver(res.df)
    else:
        if res.savename is not None:
            single_command_saver(res.df, res.savename)


if __name__ == "__main__":
    main()
