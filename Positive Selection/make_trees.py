#!/usr/bin/env python
# coding: utf-8

# In[58]:


import re
import pandas as pd
from ete3 import Tree
from skbio import read, DNA


# In[59]:


t = Tree("../Final Tree/final.tree")
df = pd.read_csv("../cluster_seqs.csv", sep=";", index_col=0)
pat = re.compile(r"[XN]M_\d+\.\d+")


# In[62]:


for trim in [6, 34, 5, 22]:
    seqs = read(f"Sequences/trim{trim}_align_dna.fasta", format="fasta")
    seq_names = {re.search(pat, seq.metadata["id"]).group(): seq.metadata["id"] for seq in seqs}
    ids = df.loc[df["RNAId"].isin(seq_names), "RNAId"]
    sub_t = t.copy()
    sub_t.prune(ids.index, preserve_branch_length=True)
    for node in sub_t.iter_leaves():
        node.name = seq_names[ids[node.name]]
    sub_t.write(outfile=f"./Trees/trim{trim}.nwk")

