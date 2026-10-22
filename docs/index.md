# ![PrimerDriver](assets/PrimerDriver_logo.png)

*PrimerDriver* is a user-friendly bioinformatics platform specifically designed to generate primers for site-directed
mutagenesis experiments into workflows with an easy-to-use command-line interface. "Traversing" *PrimerDriver* lanes,
users can generate possible mutagenic primers upon input of a DNA sequence (`DNA`), design primers for an array of
species codon expression systems through direct mutation of the amino acid (`PRO`), characterize and report
user-designed primers (`CHAR`).

!!! warning

    `PRO` mode back-translates the input protein into the selected organism's preferred codons and designs primers
    against that reconstruction, so the primers match a real template only if the gene was synthesized with exactly
    those codons. See [Software Features](software.md#choosing-primer-design-modes).
