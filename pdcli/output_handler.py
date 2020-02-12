from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def interactive_saver(res, df):
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

def singleCommand_saver(res, df, savename):
    if savename.endswith(".csv"):
        df.to_csv(savename)
        return 0
    elif savename.endswith(".html"):
        df.to_html(savename)
        return 0
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
    else:
        print("Unsupported filetype. Supported filetypes are: .csv, .html, .fasta")
        return 1
