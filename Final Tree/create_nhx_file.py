#! py -3.7
"""Add sequence information to the newick file and encode it as NHX

REMEMBER TO USE THE NEWICK REPRESENTATION FROM THE .iqtree FILE, NOT FROM THE .treefile
"""

import pandas as pd
from ete3 import Tree
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from add_info_to_file import guess_trim

NEWICK_PATH = "tree.tree"
CSV_SEQ_INFO_PATH = "../cluster_seqs.csv"
DEST_PATH = "treeNHX.tree"


def main():
	t = Tree(NEWICK_PATH, format=0)
	df = pd.read_csv(CSV_SEQ_INFO_PATH, sep=';', index_col=0)

	for node in t.iter_descendants():
		if node.is_leaf():
			features = {
				"species": df.loc[node.name, "Species"],
				"gene_id": df.loc[node.name, "GeneId"],
				"trim": guess_trim(df.loc[node.name, "ProtDef"])
			}
			node.add_features(**features)

	t.write(outfile=DEST_PATH, format=0, features=["species", "gene_id", "trim"])


if __name__ == "__main__":
	main()