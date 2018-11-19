"""
Parser for some info I need from sam files.

I need:

1. reference fasta file parsing. get species names from IDs (build dict)
2. dataset name (from dataset_names.txt)
3. parameters used (@PG header)
4. % reads mapped (samtools flagstat output / total # of reads)

example usage:
time python /home/kchan/scripts_thesis/sam_pairwise_parser.py -r /home/kchan/thesis/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta -d /home/kchan/thesis/processed/minimap2 -o /home/kchan/thesis/processed/minimap2 -f minimap2_alignment_summary.txt

IMPORTANT: include --write-top-hits if you want to write the top hits file for each sam!
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
		self.dataset_name = ""
		self.avg_ref_aligned = {}
		self.count_reads_geq_50 = defaultdict(int)
		self.count_reads_geq_90 = defaultdict(int)
		self.percent_reads_geq_50 = defaultdict(float)
		self.percent_reads_geq_90 = defaultdict(float)
		self.num_unique_reads = 0
		self.headers = []
		self.reads_aligned_per_ref = defaultdict(int)
		self.top_hits = OrderedDict()
		self.top_mapq = OrderedDict()

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
		elif sub == "I" or sub == "D" or sub == "H" or sub == "S":
			curr_nums = ""
		else:
			length += int(curr_nums)
			curr_nums = ""
	return length

def parse_dataset_name(params_line):
	dataset_path = params_line.split(" ")[-1]
	dataset_path = os.path.basename(dataset_path)
	return os.path.splitext(dataset_path)[0]

def parse_sam_file(sam_file, ref_seqs):
	aligned_len_map = defaultdict(list) # key: reference ID, values: list of percent length of read mapped to reference
	sam_file_obj = SamFile()
	with open(sam_file, "r") as infile:
		for line in infile:
			line = line.strip()
			if line.startswith("@"):
				if line.startswith("@PG"):
					pg_header = line.split("\t")
					for col in pg_header:
						if col.startswith("CL:"):
							col_stripped = col.rstrip("\n")
							sam_file_obj.dataset_name = parse_dataset_name(col_stripped)
							sam_file_obj.params = col_stripped
				sam_file_obj.headers.append(line)
			else:
				data = line.split("\t")
				try:
					read_name, ref_id, mapq = data[0], data[2], data[4]
				except:
					print "WARNING: had a problem with %s, skipping..." % sam_file
					continue

				if read_name not in sam_file_obj.top_hits:
					if int(mapq) > 0: 
						sam_file_obj.top_hits[read_name] = data
					sam_file_obj.num_unique_reads += 1
				else:
					stored_mapq = sam_file_obj.top_hits[read_name][4]
					if int(mapq) > int(stored_mapq):
						sam_file_obj.top_hits[read_name] = data

	for read_name, top_hit_line in sam_file_obj.top_hits.items():
		ref_id, mapq, cigar_string = top_hit_line[2], top_hit_line[4], top_hit_line[5]
		sam_file_obj.top_mapq[ref_id] = mapq
		cigar_len = get_cigar_len(cigar_string)
		try:
			ref_length = len(ref_seqs[ref_id].sequence)
			percent_aligned = float(cigar_len) / ref_length * 100
			if percent_aligned == 100.0:
				print "read %s aligned 100% of reference %s!" % (read_name, ref_id)
			sam_file_obj.count_reads_geq_50[ref_id] += 1 if percent_aligned >= 50.0 else 0
			sam_file_obj.count_reads_geq_90[ref_id] += 1 if percent_aligned >= 90.0 else 0
			sam_file_obj.reads_aligned_per_ref[ref_id] += 1
			aligned_len_map[ref_id].append(percent_aligned)
		except KeyError:
			print "ERROR: reference id %s found in SAM file but not in the reference FASTA. Did you pass in the correct reference file?" % ref_id
			sys.exit(1)
	for ref_id in sam_file_obj.count_reads_geq_50.keys():
		# print "ref id %s, num reads over 50 %f, num reads over 90 %f, len ref seq %d, num aligned to this ref is %d" % (ref_id, percent_reads_geq_50[ref_id], percent_reads_geq_90[ref_id], len(ref_seqs[ref_id].sequence), len(aligned_len_map[ref_id]))
		sam_file_obj.percent_reads_geq_50[ref_id] = round(float(sam_file_obj.count_reads_geq_50[ref_id]) / sam_file_obj.reads_aligned_per_ref[ref_id] * 100, 2)
	for ref_id in sam_file_obj.count_reads_geq_90.keys():
		sam_file_obj.percent_reads_geq_90[ref_id] = round(float(sam_file_obj.count_reads_geq_90[ref_id]) / sam_file_obj.reads_aligned_per_ref[ref_id] * 100, 2)
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

def process_sample_folders(folders, ref_seqs, overwrite, write_top_hits):
	ret = []
	for sample in folders:
		print "...sample %s" % sample
		for f in list_dir_abs(sample):
			if f.endswith(".sam") and not f.endswith("_top_hits.sam"):
				sam_file = parse_sam_file(f, ref_seqs)
				if write_top_hits:
					write_top_hits_sam(f, sam_file, overwrite) # comment this out for stop writing output sams
				ret.append(sam_file)
	return ret

def get_args():
	parser = argparse.ArgumentParser(description = "Script for parsing FASTA reference files and SAM files for some stuff I need. " +
		"Outputs a summary table.")
	parser.add_argument("-r", "--reference-fasta", required = True, help = "Path to the fasta file containing all reference sequences.")
	parser.add_argument("-d", "--sample-dir", required = True, help = "Directory containing folders for each sample.")
	parser.add_argument("-o", "--output-dir", default = ".", help = "Path to the output directory. [current working directory]")
	parser.add_argument("-f", "--file-output", default = "", help = "Name of the output file. [alignment_summary.txt]")
	parser.add_argument("-w", "--overwrite", action = "store_true", help = "Specify to overwrite the given sam file with a new sam file containing ONLY the top hits. [FILENAME_top_hits.sam]")
	parser.add_argument("--write-top-hits", action = "store_true", help = "Specify to write a sam file containing top hits based on max mapping quality per read. [false]")
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
	out_name, out_ext = os.path.splitext(file_output)
	overwrite = args.overwrite
	write_top_hits = args.write_top_hits

	subfolders = []
	print "fetching sample sam files..."
	for f in list_dir_abs(sample_directory):
		if os.path.isdir(f):
			subfolders.append(f)
	print "...done"
	print "processing sam files..."
	all_sam_files = process_sample_folders(subfolders, ref_seqs, overwrite, write_top_hits)
	print "...done"

	print "writing output tables..."
	outpath_table_full = os.path.join(args.output_dir, file_output)
	out_table_full = open(outpath_table_full, "w")
	out_table_full.write("params_used\tnum_unique_reads\tpercent_total_mapped\tnum_aligned_ref_total\tnum_aligned_ref_50\tpercent_50\tnum_aligned_ref_90\tpercent_90\treference_id\treference_organism\tref_id_mapq\tavg_len_aligned\n")
	
	for sam_file in all_sam_files:
		params_used = sam_file.params
		num_unique_reads = sam_file.num_unique_reads
		percent_total_mapped = sam_file.get_total_percent_mapped()
		num_aligned_ref_total = ";".join("%s;%d" % (k, v) for k, v in sam_file.reads_aligned_per_ref.items())
		num_aligned_ref_50 = ";".join("%s;%d" % (k, v) for k, v in sam_file.count_reads_geq_50.items())
		percent_50 = ";".join("%s;%s" % (k, v) for k, v in sam_file.percent_reads_geq_50.items())
		num_aligned_ref_90 = ";".join("%s;%d" % (k, v) for k, v in sam_file.count_reads_geq_90.items())
		percent_90 = ";".join("%s;%s" % (k, v) for k, v in sam_file.percent_reads_geq_90.items())
		keys = sam_file.avg_ref_aligned.keys() # not saved to table
		reference_id = ";".join(keys)
		reference_organism = ";".join(ref_seqs[k].taxonomy for k in keys)
		ref_id_mapq = ";".join("%s;%s" % (k, v) for k, v in sam_file.top_mapq.items())
		avg_len_aligned = ";".join("%s;%s" % (k, v) for k, v in sam_file.avg_ref_aligned.items())
		
		out_line = "\t".join(map(str, [params_used, num_unique_reads, percent_total_mapped, num_aligned_ref_total, num_aligned_ref_50, percent_50, num_aligned_ref_90, percent_90, reference_id, reference_organism, ref_id_mapq, avg_len_aligned]))
		out_table_full.write(out_line + "\n")
	print "...done. goodbye"

if __name__ == "__main__":
	main()
