"""This converts an AA alignment into a DNA alignment and adds the species names"""

import pandas as pd
from skbio import read, TabularMSA, Protein, DNA


def convert_seq(aa_align: str, dna_unalign):
    dna_align = ""
    for aa in aa_align:
        if aa == "-":
            dna_align += "---"
        else:
            dna_align += dna_unalign[:3]
            dna_unalign = dna_unalign[3:]
    if len(dna_unalign) != 3:
        raise ValueError
    return dna_align


dna_seqs = {s.metadata['id']: s for s in read("../final_seqs_dna.fasta", "fasta")}
df = pd.read_csv("../cluster_seqs.csv", sep=';', index_col=0)

for trim in ("trim6", "trim34", "trim5", "trim22"):
    aa_align = TabularMSA.read(f"Sequences/{trim}_align.fasta", "fasta", constructor=Protein)
    dna_aligns = []
    for seq in aa_align:
        dna_name = df.loc[seq.metadata['id'], "RNAId"]
        dna_unalign = dna_seqs[dna_name]
        seq_str = convert_seq(str(seq), str(dna_unalign))
        dna_aligns.append(DNA(seq_str, metadata={"id": dna_name}))
    TabularMSA(dna_aligns).write(format="fasta", file=f"Sequences/{trim}_align_dna.fasta")