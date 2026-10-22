from __future__ import annotations

from datetime import datetime
from os import PathLike

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger
from pandas import DataFrame, concat


def interactive_saver(df: DataFrame | list[DataFrame]) -> None:
    save = input("\nSave? [y/n] ")
    if save.upper() != "Y":
        return
    while True:
        try:
            single_command_saver(df, input("Enter filename: "))
            break
        except ValueError as e:
            logger.error(str(e))


def single_command_saver(df: DataFrame | list[DataFrame], savename: str | PathLike):
    dfs = [df] if isinstance(df, DataFrame) else df
    savename = str(savename)
    if savename.endswith(".csv"):
        concat(dfs).to_csv(savename)
    elif savename.endswith(".html"):
        concat(dfs).to_html(savename)
    elif savename.endswith(".fasta"):
        seq_list = []
        for i, primer in enumerate(dfs, start=1):
            row = primer.iloc[0]
            for strand, prefix in (("Forward", "Fwd"), ("Reverse", "Rev")):
                seq_list.append(
                    SeqRecord(
                        Seq(row[strand], IUPAC.unambiguous_dna),
                        id=f"{datetime.now()}",
                        description=f"Primer {i} {strand.lower()} strand, {row[f'{prefix} GC content']} G/C, Tm = {row[f'{prefix} melting temp']}, {row[f'{prefix} length']}",
                    )
                )
        with open(savename, "w") as f:
            SeqIO.write(seq_list, f, "fasta")
    elif savename.endswith(".json"):
        concat(dfs).to_json(savename, indent=4)
    else:
        raise ValueError(
            "Unsupported filetype. Supported filetypes are: .csv, .html, .fasta, .json"
        )
