from datetime import datetime
from os import PathLike

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger
from pandas import DataFrame, concat


def interactive_saver(df: DataFrame) -> None:
    save = input("\nSave? [y/n] ")
    if save.upper() == "Y":
        while True:
            savename = input("Enter filename: ")
            if savename.endswith(".csv"):
                df = concat([*df])
                df.to_csv(savename)
                break
            elif savename.endswith(".html"):
                df = concat([*df])
                df.to_html(savename)
                break
            elif savename.endswith(".fasta"):
                seq_list = []
                for i in range(len(df)):
                    seq_list.append(
                        SeqRecord(
                            Seq(df[i]["Forward"].values[0], IUPAC.unambiguous_dna),
                            id=f"{datetime.now()}",
                            description=f"Primer {i+1} forward strand, {df[i]['Fwd GC content'].values[0]*100:.2s}% G/C, Tm = {df[i]['Fwd melting temp'].values[0]:.2s} C, {df[i]['Fwd length'].values[0]}",
                        )
                    )
                    seq_list.append(
                        SeqRecord(
                            Seq(df[i]["Reverse"].values[0], IUPAC.unambiguous_dna),
                            id=f"{datetime.now()}",
                            description=f"Primer {i+1} reverse strand, {df[i]['Rev GC content'].values[0]*100:.2s}% G/C, Tm = {df[i]['Rev melting temp'].values[0]:.2s} C, {df[i]['Rev length'].values[0]}",
                        )
                    )
                with open(savename, "w") as f:
                    SeqIO.write(seq_list, f, "fasta")
                break
            elif savename.endswith(".json"):
                df = concat([*df])
                df.to_json(savename, indent=4)
            else:
                logger.error(
                    "Unsupported filetype. Supported filetypes are: .csv, .html, .fasta, .json"
                )


def single_command_saver(df: DataFrame, savename: str | PathLike):
    if savename.endswith(".csv"):
        df = concat([*df])
        df.to_csv(savename)
    elif savename.endswith(".html"):
        df = concat([*df])
        df.to_html(savename)
    elif savename.endswith(".fasta"):
        seq_list = []
        for i in range(len(df)):
            seq_list.append(
                SeqRecord(
                    Seq(df[i]["Forward"].values[0], IUPAC.unambiguous_dna),
                    id=f"{datetime.now()}",
                    description=f"Primer {i+1} forward strand, {df[i]['Fwd GC content'].values[0]*100:.2s}% G/C, Tm = {df[i]['Fwd melting temp'].values[0]:.2s} C, {df[i]['Fwd length'].values[0]}",
                )
            )
            seq_list.append(
                SeqRecord(
                    Seq(df[i]["Reverse"].values[0], IUPAC.unambiguous_dna),
                    id=f"{datetime.now()}",
                    description=f"Primer {i+1} reverse strand, {df[i]['Rev GC content'].values[0]*100:.2s}% G/C, Tm = {df[i]['Rev melting temp'].values[0]:.2s} C, {df[i]['Rev length'].values[0]}",
                )
            )
        with open(savename, "w") as f:
            SeqIO.write(seq_list, f, "fasta")
    elif savename.endswith(".json"):
        df = concat([*df])
        df.to_json(savename, indent=4)
    else:
        raise ValueError(
            "Unsupported filetype. Supported filetypes are: .csv, .html, .fasta, .json"
        )
