"""
Parser for some info I need from sam files.

I need:

1. reference fasta file parsing. get species names from IDs (build dict)
2. dataset name (from dataset_names.txt)
3. parameters used (@PG header)
4. % reads mapped (samtools flagstat output / total # of reads)

example usage:
time python /home/kchan/scripts_thesis/sam_pairwise_parser.py -r /home/kchan/thesis/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta -d /home/kchan/thesis/processed/minimap2 -o /home/kchan/thesis/processed/minimap2 -f minimap2_alignment_summary.txt
"""

from collections import defaultdict, OrderedDict
from Bio import SeqIO
import argparse
import sys
import os

class ReferenceRecord:
	def __init__(self, id, taxonomy, sequence):
		self.id = id
		self.taxonomy = taxonomy
		self.sequence = sequence

class SamFile():
	def __init__(self):
		self.params = None
		self.avg_ref_aligned = {}
		self.num_unique_reads = 0
		self.headers = []
		self.top_hits = OrderedDict()

	def get_total_percent_mapped(self):
		if self.num_unique_reads == 0:
			print "encountered a file with no unique reads, ignoring..."
			return 0
		return round(float(len(self.top_hits)) / self.num_unique_reads * 100, 2)

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

# not used...as of now
# def get_total_percent_mapped(flagstats):
# 	percentage = ""
# 	f = open(flagstats, "r").readlines()
# 	match = re.match(r"(.)*(\(.*\))(.)*", f[4])
# 	if match:
# 		percentage = match.group(2)
# 		percentage = percentage.split(":")[0][1:].rstrip()
# 	return percentage if percentage != "N/A" else "0.0%"

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
		# TODO: remove this after regenerating SAMs with headers inside...
		filename = os.path.splitext(os.path.basename(sam_file))[0]
		sam_file_obj.params = filename

		for line in infile:
			line = line.strip()
			if line.startswith("@"):
				if line.startswith("@PG"):
					pg_header = line.split("\t")
					for col in pg_header:
						if col.startswith("CL:"):
							sam_file_obj.params = col.rstrip("\n")
				sam_file_obj.headers.append(line)
			# if line.startswith("@HD") or line.startswith("@SQ") or line.startswith("@RG") or line.startswith("@CO"):
			# 	sam_file_obj.headers.append(line)
			# elif line.startswith("@PG"):
			# 	pg_header = line.split("\t")
			# 	for col in pg_header:
			# 		if col.startswith("CL:"):
			# 			sam_file_obj.params = col.rstrip("\n")
			# 	sam_file_obj.headers.append(line)
			else:
				data = line.split("\t")
				try:
					read_name, mapq = data[0], data[4]
				except:
					print "FILENAME  THAT FAILED: %s" % sam_file
					sys.exit(1)

				if read_name not in sam_file_obj.top_hits:
					if int(mapq) > 0: 
						sam_file_obj.top_hits[read_name] = data
					sam_file_obj.num_unique_reads += 1
				else:
					stored_mapq = sam_file_obj.top_hits[read_name][4]
					if int(mapq) > int(stored_mapq):
						sam_file_obj.top_hits[read_name] = data

	for top_hit_line in sam_file_obj.top_hits.values():
		ref_id, cigar_string = top_hit_line[2], top_hit_line[5]
		cigar_len = get_cigar_len(cigar_string)
		try:
			percent_aligned = float(cigar_len) / len(ref_seqs[ref_id].sequence)
			aligned_len_map[ref_id].append(percent_aligned)
		except KeyError:
			print "ERROR: reference id %s found in SAM file but not in the reference FASTA. Did you pass in the correct reference file?" % ref_id
			sys.exit(1)

	for k in aligned_len_map.keys():
		sam_file_obj.avg_ref_aligned[k] = round(float(sum(aligned_len_map[k])) / len(aligned_len_map[k]), 2)

	return sam_file_obj

def write_top_hits_sam(filepath, sam_file_obj, overwrite):
	outfile_name, outfile_ext = os.path.splitext(os.path.basename(filepath))
	print "...writing sam file..."
	outpath = os.path.dirname(filepath)
	if not overwrite:
		outfile_name = "%s_top_hits" % outfile_name
	outfile = open(os.path.join(outpath, outfile_name + outfile_ext), "w")
	outfile.write("\n".join(sam_file_obj.headers) + "\n")
	for line in sam_file_obj.top_hits.values():
		outfile.write("\t".join(line) + "\n")
	outfile.close()

def process_sample_folders(folders, ref_seqs, overwrite):
	ret = []
	for sample in folders:
		print "...sample %s" % sample
		for f in list_dir_abs(sample):
			if f.endswith(".sam"):
				sam_file = parse_sam_file(f, ref_seqs)
				write_top_hits_sam(f, sam_file, overwrite)
				ret.append(sam_file)
	return ret

def get_args():
	parser = argparse.ArgumentParser(description = "Script for parsing FASTA reference files and SAM files for some stuff I need. " +
		"Outputs a summary table.")
	parser.add_argument("-r", "--reference-fasta", required = True, help = "Path to the fasta file containing all reference sequences.")
	parser.add_argument("-d", "--sample-dir", required = True, help = "Directory containing folders for each sample.")
	parser.add_argument("-o", "--output-dir", default = ".", help = "Output directory. Default: current working directory.")
	parser.add_argument("-f", "--file-output", default = "", help = "Name of the output file. If not specified, default is 'alignment_summary.txt'.")
	parser.add_argument("-w", "--overwrite", action = "store_true", help = "Specify to overwrite the given sam file with a new sam file containing ONLY the top hits. If not specified, the output sam file will be named FILENAME_top_hits.sam.")
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
	overwrite = args.overwrite

	subfolders = []
	print "fetching sample sam files..."
	for f in list_dir_abs(sample_directory):
		if os.path.isdir(f):
			subfolders.append(f)
	print "...done"
	print "processing sam files..."
	all_sam_files = process_sample_folders(subfolders, ref_seqs, overwrite)
	print "...done"

	print "writing output table..."
	outpath = os.path.join(args.output_dir, file_output)
	outfile = open(outpath, "w")
	outfile.write("params_used\tpercent_total_mapped\treference_id\treference_organism\tavg_len_aligned\n")
	for sam_file in all_sam_files:
		params_used = sam_file.params
		percent_total_mapped = sam_file.get_total_percent_mapped()
		keys = sam_file.avg_ref_aligned.keys() # not saved to table
		reference_id = ";".join(keys)
		reference_organism = ";".join(ref_seqs[k].taxonomy for k in keys)
		avg_len_aligned = ";".join("%s;%s" % (k, v) for k, v in sam_file.avg_ref_aligned.items())
		outfile.write(str(params_used) + "\t" + str(percent_total_mapped) + "\t" + str(reference_id) + "\t" + str(reference_organism) + "\t" + str(avg_len_aligned) + "\n")
	print "...done. goodbye"

if __name__ == "__main__":
	main()
