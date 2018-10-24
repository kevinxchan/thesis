"""
Parser for some info I need from sam files.

I need:

1. reference fasta file parsing. get species names from IDs (build dict)
2. dataset name (from dataset_names.txt)
3. parameters used (@PG header)
4. % reads mapped (samtools flagstat output / total # of reads)

example usage:
python sam_pairwise_parser.py -r /home/kchan/thesis/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta -d /home/kchan/thesis/processed/minimap2
"""

from Bio import SeqIO
import argparse
import os
import re

class ReferenceRecord:
	def __init__(self, id, taxonomy, sequence):
		self.id = id
		self.taxonomy = taxonomy
		self.sequence = sequence

def list_dir_abs(dir):
	for f in os.listdir(dir):
		yield os.path.abspath(os.path.join(dir, f))

def build_ref_seq_map(fasta_file):
	ret = {}
	with open(fasta_file, "rU") as f:
		for record in SeqIO.parse(f, "fasta"):
			_id = record.id
			taxonomy = record.description.split("organism=")[1].split(",")[0]
			seq = record.seq
			ret[_id] = ReferenceRecord(_id, taxonomy, seq)
	return ret

def get_total_percent_mapped(flagstats):
	percentage = ""
	f = open(flagstats, "r").readlines()
	match = re.match(r"(.)*(\(.*\))(.)*", f[4])
	if match:
		percentage = match.group(2)
		percentage = percentage.split(":")[0][1:].rstrip()
	return percentage if percentage != "N/A" else "0.0%"

def process_sample_folders(folders):
	for sample in folders:
		for f in list_dir_abs(sample):
			if f.endswith(".sam"):
				txt = os.path.basename(f).replace(".sam", ".txt")
				flagstats_file = os.path.join(sample, "flagstats", txt)
				total_percent_mapped = get_total_percent_mapped(flagstats_file)
				print total_percent_mapped
	return None

def get_args():
	parser = argparse.ArgumentParser(description = "Script for parsing FASTA reference files and SAM files for some stuff I need. " +
		"Outputs a summary table.")
	parser.add_argument("-r", "--reference-fasta", required = True, help = "Path to the fasta file containing all reference sequences.")
	parser.add_argument("-d", "--sample-dir", required = True, help = "Directory containing folders for each sample.")
	parser.add_argument("-o", "--output-dir", default = ".", help = "Output directory. Default: current working directory.")
	args = parser.parse_args()
	return args

def main():
	args = get_args()
	fasta_reference = args.reference_fasta
	ref_seqs = build_ref_seq_map(fasta_reference)
	sample_directory = args.sample_dir

	subfolders = []
	for f in list_dir_abs(sample_directory):
		if os.path.isdir(f):
			subfolders.append(f)
	result = process_sample_folders(subfolders)

if __name__ == "__main__":
	main()