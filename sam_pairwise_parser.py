"""
Parser for some info I need from sam files.

I need:

1. reference fasta file parsing. get species names from IDs (build dict)
2. dataset name (from dataset_names.txt)
3. parameters used (@PG header)
4. % reads mapped (samtools flagstat output / total # of reads)

example usage:
python /home/kchan/scripts_thesis/sam_pairwise_parser.py -r /home/kchan/thesis/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta -d /home/kchan/thesis/processed/minimap2 -o /home/kchan/thesis/processed/minimap2 -f minimap2_alignment_summary.txt
"""

from collections import defaultdict
from Bio import SeqIO
import argparse
import os
import re

class ReferenceRecord:
	def __init__(self, id, taxonomy, sequence):
		self.id = id
		self.taxonomy = taxonomy
		self.sequence = sequence

class SamFile():
	def __init__(self):
		self.params = None
		self.avg_ref_aligned = {}

class OutputData():
	def __init__(self, total_percent_mapped, sam_file):
		self.total_percent_mapped = total_percent_mapped
		self.sam_file = sam_file

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

def get_cigar_len(cigar):
	length, curr_nums = 0, ""
	for sub in cigar:
		if sub.isdigit():
			curr_nums += sub
		elif sub == "I" or sub == "D":
			curr_nums = ""
		else:
			length += int(curr_nums)
			curr_nums = ""
	return length

def parse_sam_file(sam_file, ref_seqs):
	aligned_len_map = defaultdict(list)
	sam_file_obj = SamFile()
	with open(sam_file, "r") as infile:
		for line in infile:
			if line.startswith("@HD") or line.startswith("@SQ") or line.startswith("@RG") or line.startswith("@CO"):
				continue
			elif line.startswith("@PG"):
				pg_header = line.split("\t")
				for col in pg_header:
					if col.startswith("CL:"):
						sam_file_obj.params = col.rstrip("\n")
			else:
				data = line.split("\t")
				ref_id, cigar_string = data[2], data[5]
				cigar_len = get_cigar_len(cigar_string)
				percent_aligned = float(cigar_len) / len(ref_seqs[ref_id].sequence)
				aligned_len_map[ref_id].append(percent_aligned)
		for k in aligned_len_map.keys():
			sam_file_obj.avg_ref_aligned[k] = round(float(sum(aligned_len_map[k])) / len(aligned_len_map[k]), 2)
	return sam_file_obj

def process_sample_folders(folders, ref_seqs):
	ret = []
	for sample in folders:
		print "...sample %s" % sample
		for f in list_dir_abs(sample):
			if f.endswith(".sam"):
				txt = os.path.basename(f).replace(".sam", ".txt")
				flagstats_file = os.path.join(sample, "flagstats", txt)
				total_percent_mapped = get_total_percent_mapped(flagstats_file)
				sam_file = parse_sam_file(f, ref_seqs)
				output_data = OutputData(total_percent_mapped, sam_file)
				ret.append(output_data)
	return ret

def get_args():
	parser = argparse.ArgumentParser(description = "Script for parsing FASTA reference files and SAM files for some stuff I need. " +
		"Outputs a summary table.")
	parser.add_argument("-r", "--reference-fasta", required = True, help = "Path to the fasta file containing all reference sequences.")
	parser.add_argument("-d", "--sample-dir", required = True, help = "Directory containing folders for each sample.")
	parser.add_argument("-o", "--output-dir", default = ".", help = "Output directory. Default: current working directory.")
	parser.add_argument("-f", "--file-output", default = "", help = "Name of the output file. If not specified, default is 'alignment_summary.txt'.")
	args = parser.parse_args()
	return args

def main():
	args = get_args()
	fasta_reference = args.reference_fasta
	print "fetching reference sequences..."
	ref_seqs = build_ref_seq_map(fasta_reference)
	print "...done"
	sample_directory = args.sample_dir
	file_output = args.file_output if args.file_output else "alignment_summary.txt"

	subfolders = []
	print "fetching sample sam files..."
	for f in list_dir_abs(sample_directory):
		if os.path.isdir(f):
			subfolders.append(f)
	print "...done"
	print "processing sam files..."
	all_samples = process_sample_folders(subfolders, ref_seqs)
	print "...done"

	print "writing output file..."
	outpath = os.path.join(args.output_dir, file_output)
	outfile = open(outpath, "w")
	outfile.write("params_used\tpercent_total_mapped\treference_id\treference_organism\tavg_len_aligned\n")
	for sample in all_samples:
		params_used = sample.sam_file.params
		percent_total_mapped = sample.total_percent_mapped
		keys = sample.sam_file.avg_ref_aligned.keys() # not saved to table
		reference_id = ";".join(keys)
		reference_organism = ";".join(ref_seqs[k].taxonomy for k in keys)
		avg_len_aligned = ";".join("%s;%s" % (k, v) for k, v in sample.sam_file.avg_ref_aligned.items())
		outfile.write(str(params_used) + "\t" + str(percent_total_mapped) + "\t" + str(reference_id) + "\t" + str(reference_organism) + "\t" + str(avg_len_aligned) + "\n")
	print "...done. goodbye"

if __name__ == "__main__":
	main()
