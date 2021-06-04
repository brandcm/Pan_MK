# Tim Webster, University of Utah, 2021

from __future__ import print_function
import argparse
import collections
from itertools import groupby
import numpy as np
import sys
import textwrap


def parse_args():
	"""
	Parse command line arguments
	"""
	parser = argparse.ArgumentParser(description="")

	parser.add_argument(
		"--fasta", type=str, required=True,
		help="Full path to input fasta file (output of sfs_code)")

	parser.add_argument(
		"--outfile", type=str, required=True, help="Full path to and name of "
		"desired output file. Will overwrite if it exists.")

	args = parser.parse_args()

	return args


def trim_incomplete_codons(seq):
	"""
	Checks to see if final codon is three bases long, and if not, trims it.
	"""
	if len(seq[-1]) != 3:
		return seq[:-1]
	else:
		return seq


def main():
	args = parse_args()

	amino_acids = {
		"TGC": "C",
		"TGT": "C",
		"AGC": "S",
		"AGT": "S",
		"TCA": "S",
		"TCC": "S",
		"TCG": "S",
		"TCT": "S",
		"ACA": "T",
		"ACC": "T",
		"ACG": "T",
		"ACT": "T",
		"CCA": "P",
		"CCC": "P",
		"CCG": "P",
		"CCT": "P",
		"GCA": "A",
		"GCC": "A",
		"GCG": "A",
		"GCT": "A",
		"GGA": "G",
		"GGC": "G",
		"GGG": "G",
		"GGT": "G",
		"AAC": "N",
		"AAT": "N",
		"GAC": "D",
		"GAT": "D",
		"GAA": "E",
		"GAG": "E",
		"CAA": "Q",
		"CAG": "Q",
		"CAC": "H",
		"CAT": "H",
		"AGA": "R",
		"AGG": "R",
		"CGA": "R",
		"CGC": "R",
		"CGG": "R",
		"CGT": "R",
		"AAA": "K",
		"AAG": "K",
		"ATG": "M",
		"ATA": "I",
		"ATC": "I",
		"ATT": "I",
		"CTA": "L",
		"CTC": "L",
		"CTG": "L",
		"CTT": "L",
		"TTA": "L",
		"TTG": "L",
		"GTA": "V",
		"GTC": "V",
		"GTG": "V",
		"GTT": "V",
		"TTC": "F",
		"TTT": "F",
		"TAC": "Y",
		"TAT": "Y",
		"TGG": "W",
		"TAA": ".",
		"TAG": ".",
		"TGA": "."}

	t_repl_sil = {
		"TGC": [2.666, .334],
		"TGT": [2.666, .334],
		"AGC": [2.666, .334],
		"AGT": [2.666, .334],
		"TCA": [2.0, 1.0],
		"TCC": [2.0, 1.0],
		"TCG": [2.0, 1.0],
		"TCT": [2.0, 1.0],
		"ACA": [2.0, 1.0],
		"ACC": [2.0, 1.0],
		"ACG": [2.0, 1.0],
		"ACT": [2.0, 1.0],
		"CCA": [2.0, 1.0],
		"CCC": [2.0, 1.0],
		"CCG": [2.0, 1.0],
		"CCT": [2.0, 1.0],
		"GCA": [2.0, 1.0],
		"GCC": [2.0, 1.0],
		"GCG": [2.0, 1.0],
		"GCT": [2.0, 1.0],
		"GGA": [2.0, 1.0],
		"GGC": [2.0, 1.0],
		"GGG": [2.0, 1.0],
		"GGT": [2.0, 1.0],
		"AAC": [2.666, .334],
		"AAT": [2.666, .334],
		"GAC": [2.666, .334],
		"GAT": [2.666, .334],
		"GAA": [2.666, .334],
		"GAG": [2.666, .334],
		"CAA": [2.666, .334],
		"CAG": [2.666, .334],
		"CAC": [2.666, .334],
		"CAT": [2.666, .334],
		"AGA": [2.333, .667],
		"AGG": [2.333, .667],
		"CGA": [1.666, 1.334],
		"CGC": [2.0, 1.0],
		"CGG": [2.0, 1.0],
		"CGT": [2.0, 1.0],
		"AAA": [2.666, .334],
		"AAG": [2.666, .334],
		"ATG": [3.0, 0.0],
		"ATA": [2.333, .667],
		"ATC": [2.333, .667],
		"ATT": [2.333, .667],
		"CTA": [1.666, 1.334],
		"CTC": [2.0, 1.0],
		"CTG": [1.666, 1.334],
		"CTT": [2.0, 1.0],
		"TTA": [2.333, .667],
		"TTG": [2.333, .667],
		"GTA": [2.0, 1.0],
		"GTC": [2.0, 1.0],
		"GTG": [2.0, 1.0],
		"GTT": [2.0, 1.0],
		"TTC": [2.666, .334],
		"TTT": [2.666, .334],
		"TAC": [2.666, .334],
		"TAT": [2.666, .334],
		"TGG": [3.0, 0.0],
		"TAA": [2.333, .667],
		"TAG": [2.666, .334],
		"TGA": [2.666, .334]}

	if sys.version_info[0] < 3:
		sys.exit("Error. This script requires Python 3 or greater")

	else:
		# python 3 code
		with open(args.fasta, "r") as f:
			# Create dictionary to hold sequences
			sequence_dict = collections.OrderedDict()

			# Read and process Fasta
			# Code for parsing fasta based on code provided by Bret Pedersen: https://www.biostars.org/p/710/
			faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))

			counter = 0
			for header in faiter:
				# drop the ">" and grab the header
				header = next(header)[1:].strip()

				# join all sequence lines to one.
				seq = "".join(s.strip() for s in next(faiter))
				
				# ensure that no sites are soft-masked to ensure the dictionary works
				seq = seq.upper()

				# concatenate coding sequences for genes
				if header in sequence_dict:
					sequence_dict[header] += seq
				else:
					sequence_dict[header] = seq

				counter += 1
				if counter % 100 == 0:
					print("{} fasta records processed".format(counter))

		with open(args.outfile, "w") as o:
			o.write(
				"locus\traw_length\tclipped_length\tt_sil\tt_repl\n")
			for i in sequence_dict:
				print(i)
				seq_len_raw = len(sequence_dict[i])
				seq_codons = trim_incomplete_codons(
					textwrap.wrap(sequence_dict[i], 3))
				seq_clipped = "".join(seq_codons)
				tsil = np.sum([t_repl_sil[codon][1] for codon in seq_codons])
				trepl = np.sum([t_repl_sil[codon][0] for codon in seq_codons])

				o.write(
					"{}\t{}\t{}\t{}\t{}\n".format(
						i, seq_len_raw, len(seq_clipped), tsil, trepl))

if __name__ == "__main__":
	main()
