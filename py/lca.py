"""
Script to rank alignments based on LCA and proportion of reads mapped over 80% of the reference.

How it works:
1. Query JGI taxonomy server to get full taxonomies of the organisms in the reference fasta.
2. Do 1. again, but with the reported organism in each dataset.
3. Calculate the LCA between 1. and 2., and use this as the comparison point.
4. For each reference aligned to, get the LCA between it and 3. Multiply by # reads aligned > 80% 
5. Repeat for all parameter sets, for all datasets.
6. Pick the param set with the most min. number of scores.
"""

import sys
import os
import requests
from collections import OrderedDict
from argparse import ArgumentParser
from Bio import SeqIO
from sam_pairwise_parser import build_ref_seq_map



def get_args():
	parser = ArgumentParser(description = "Script to rank alignments based on LCA and proportion of reads mapped over 80% of the reference.")
	parser.add_argument("-r", "--reference-fasta", required = True, help = "Path to the fasta file containing all reference sequences.")
	parser.add_argument("-s", "--sample-dir", required = True, help = "Directory containing folders for each sample.")
	parser.add_argument("-t", "--dataset-taxonomies", required = True, help = "Tab delimited file matching dataset names to their reported taxonomies.")
	parser.add_argument("-o", "--output-dir", default = ".", help = "Path to the output directory. [current working directory]")
	args = parser.parse_args()
	return args

def query_jgi_taxa():
	pass

def main():
	args = get_args()
	dataset_to_reported_taxa = {}
	with open(args.dataset_taxonomies, "r") as f:
		for line in f:
			dataset_id, taxonomy = line.strip().split("\t")
			full_lineage = query_jgi_taxa(taxonomy)
			dataset_to_reported_taxa[dataset_id] = full_lineage