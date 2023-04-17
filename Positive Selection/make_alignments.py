#!/usr/bin/env python
# coding: utf-8

# In[27]:


import os
import pandas as pd
from skbio import read, write


# In[28]:


aa_seqs = list(read("../final_seqs.fasta", format="fasta"))
df = pd.read_csv("../cluster_seqs.csv", sep=';', index_col=0)
df = df[df["TreeStatus"].str.contains("INCLUDED")]


# In[29]:


trims = {"TRIM6": [],
         "TRIM34": [],
         "TRIM5": [],
         "TRIM22": []}
for seq in aa_seqs:
    seq_id = seq.metadata["id"]
    trim = df.loc[seq_id, "TreeAnnotation"]
    trims[trim].append(seq)


# In[31]:


for trim, seqs in trims.items():
    raw_seq_path = f"Sequences/{trim.lower()}_seqs.fasta"
    write((seq for seq in seqs), format="fasta", into=raw_seq_path)
    os.system(f"clustalo --force -i {raw_seq_path} -o Sequences/{trim.lower()}_align.fasta")

