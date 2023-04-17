import pandas as pd
import re
from optparse import OptionParser
from functools import partial


def guess_trim(seq_title):
	regex = re.compile(r"TRIM(\d+)|tripartite[- ]motif[^\d]*(\d+)", re.IGNORECASE)
	matches = re.findall(regex, seq_title)
	if not matches:
		match = re.search("uncharacterized", seq_title, re.IGNORECASE)
		if match is None:
			raise Exception(f"could not guess from \"{seq_title}\"")
		return "UNCHAR"
	return "-".join(set(f"TRIM{trim[0] if trim[0] else trim[1]}" for trim in matches))


def _sub_func(match, encl_chars, sep_char):
	seq = match.group()
	title = _df.loc[seq, "ProtDef"]
	locus = _df.loc[seq, "GeneId"]
	sp = _df.loc[seq, "Species"]
	sp = f"{sp[0]}. {sp[sp.index(' ') + 1:]}"
	sp = sp.replace(" ", "_")
	return f"{seq}{encl_chars[0]}{guess_trim(title)}{sep_char}{locus}{sep_char}{sp}{encl_chars[1]}"


def main():
	parser = OptionParser()
	parser.add_option("-b", action="store_true", help="create backup files beforehand")
	parser.add_option("-e", default="[]", help="the two characters in which the metadata will be enclosed [default: \"%default\"]")
	parser.add_option("-s", default="|", help="separation character for metadata [default: \"%default\"]")
	options, args = parser.parse_args()

	if not args:
		args = [input("Input name of the file in which the identifiers will be added: ")]

	with open(args[0]) as f:
		text = f.read()
	if options.b:
		with open(args[0] + ".old", 'w') as f:
			f.write(text)
	with open(args[0], 'w') as f:
		f.write(re.sub(r"[XN]P_[\d\.]+", partial(_sub_func, encl_chars=options.e, sep_char=options.s), text))


if __name__ == "__main__":
	_df = pd.read_csv("cluster_seqs.csv", sep=';', index_col=0)
	main()