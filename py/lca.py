"""
Script to rank alignments based on LCA and proportion of reads mapped over 80% of the reference.

How it works:
1. Query JGI taxonomy server to get full taxonomies of the organisms in the reference fasta.
2. Do 1. again, but with the reported organism in each dataset.
3. Calculate the LCA between 1. and 2., and use this as the comparison point.
4. For each reference aligned to, get the LCA between it and 3. Multiply by # reads aligned >= 80% 
5. Repeat for all parameter sets, for all datasets.
6. Pick the param set with the most min. number of scores.

USAGE:
time python /home/kchan/scripts_thesis/py/lca.py -r /home/kchan/thesis/references/recA/fungene_9.5.1_recA_nucleotide_uclust99.fasta -s /home/kchan/thesis/processed/minimap2 -t /home/kchan/thesis/dataset_to_taxon.txt
"""

from __future__ import division
import sys
import os
import requests
import time
from collections import defaultdict
from argparse import ArgumentParser
from util.file_utils import list_dir_abs
from sam_pairwise_parser import build_ref_seq_map, get_cigar_len
from model.sam_parser_classes import SamFile

def get_args():
	parser = ArgumentParser(description = "Script to rank alignments based on LCA and proportion of reads mapped over 80% of the reference.")
	parser.add_argument("-r", "--reference-fasta", required = True, help = "Path to the fasta file containing all reference sequences.")
	parser.add_argument("-s", "--sample-dir", required = True, help = "Directory containing folders for each sample.")
	parser.add_argument("-t", "--dataset-taxonomies", required = True, help = "Tab delimited file matching dataset names to their reported taxonomies.")
	parser.add_argument("-o", "--output-dir", default = ".", help = "Path to the output directory. [current working directory]")
	parser.add_argument("-f", "--output-file", default = "lca", help = "Output file name. Extension will be '.txt'. [lca]")	
	parser.add_argument("-m", "--marker-gene", required = True, help = "Name of the marker gene being evaluated.")
	args = parser.parse_args()
	return args

def validate_query(entry):
	illegals = [" ", "<", ">", "#", "%", "{", "}", "|", "\"", "/", "^", "~", "[", "]", "`", "-", "'", '"', "."]
	for illegal in illegals:
		entry = entry.replace(illegal, "_")
	return entry

def query_jgi_error(entry):
	"""
	verifies whether or not the NCBI contains an entry for query 'entry'.

	@return True if the NCBI doesn't contain 'entry', false otherwise
	"""
	url = "http://taxonomy.jgi-psf.org/name"
	r = requests.get(os.path.join(url, entry))
	response = r.json()
	return "error" in response[entry].keys()

def query_jgi_taxa(entry):
	"""
	queries the jgi taxonomy server and returns the full lineage (as a list) for a given taxa string.

	@return full lineage of 'entry'
	"""
	url = "http://taxonomy.jgi-psf.org/name/sc"
	entry = validate_query(entry)
	if query_jgi_error(entry):
		print "[query_jgi_taxa] WARNING: query sequence %s. assigning an empty LCA." % (entry)
		return []
	r = requests.get(os.path.join(url, entry))
	try:
		response = r.json()
		error = response.values()[0]
		print "[query_jgi_taxa] ERROR: query with %s, %s failed with response '%s'. exiting." % (entry, error)
		sys.exit(1)
	except ValueError:
		response = r.text
	tax_list = response.split(";")
	no_rank_level = 1
	for i in range(len(tax_list)):
		level = tax_list[i]
		if len(level.split(":")) == 1:
			rank = "nr%d" % no_rank_level
			print "[query_jgi_taxa] WARNING: taxonomic rank not found for level '%s' assigning a rank of '%s'." % (level, rank)
			tax_list[i] = rank + ":" + level
			no_rank_level += 1
	return tax_list

def query_jgi_lca(entry_1, entry_2):
	"""
	queries the jgi taxonomy server and retrieves the LCA between two taxa strings.
	endpoints:
	/ancestor/entry_1,entry_2,entry_n: get the LCA between two or more taxonomy strings
	/sc/entry_1: 					   get the full taxa string as a semicolon delimited string

	possible to combine the two above (i.e. /ancestor/sc/...)

	@return a list containing the full lineage of the LCA
	"""
	url = "http://taxonomy.jgi-psf.org/name/ancestor/sc"

	entry_1 = validate_query(entry_1)
	entry_2 = validate_query(entry_2)
	if query_jgi_error(entry_1) or query_jgi_error(entry_2):
		print "[query_jgi_lca] WARNING: reference sequence %s or %s not found in NCBI. assigning an empty LCA." % (entry_1, entry_2)
		return []

	query = ",".join([entry_1, entry_2])
	r = requests.get(os.path.join(url, query))
	try:
		response = r.json()
		error = response.values()[0]
		print "[query_jgi_lca] ERROR: query with %s, %s failed with response '%s'. exiting." % (entry_1, entry_2, error)
		sys.exit(1)
	except ValueError:
		response = r.text
	tax_list = response.split(";")
	no_rank_level = 1
	for i in range(len(tax_list)):
		level = tax_list[i]
		if len(level.split(":")) == 1:
			rank = "nr%d" % no_rank_level
			print "[query_jgi_lca] WARNING: taxonomic rank not found for level '%s' assigning a rank of '%s'." % (level, rank)
			tax_list[i] = rank + ":" + level
			no_rank_level += 1
	return tax_list

def get_longest_lca(dataset_to_lca):
	"""
	returns the deepest lca between each dataset and the reference sequences.
	# TODO: what if multiple references are of the same length?

	@return dict, where keys = dataset id and values = deepest lca
	"""
	optimal_placements = {}
	for _id in dataset_to_lca:
		depth = -1
		optimal_placement = []
		for ref_id in dataset_to_lca[_id]:
			if len(dataset_to_lca[_id][ref_id]) > depth:
				optimal_placement = dataset_to_lca[_id][ref_id]
				depth = len(optimal_placement)
		optimal_placements[_id] = optimal_placement
	return optimal_placements

def parse_sam_file(sam_file, ref_seqs):
	"""
	parses a sam file and retrieves information specific to:
		1. parameters used to generate the sam file
		2. the dataset name (unique id)
		3. each reference aligned to, and the length of the alignment as a percentage

	@return aligned_len_map a dictionary mapping reference ids to their alignments as percentages
	@return sam_file_obj an object representation of the information from the sam file mentioned above
	"""
	aligned_len_map = defaultdict(list) # key: reference ID, values: list of percent length of read mapped to reference
	sam_file_obj = SamFile()
	with open(sam_file, "r") as infile:
		for line in infile:
			line = line.strip()
			if line.startswith("@"):
				if line.startswith("@PG"):
					pg_header = line.split("\t")
					for col in pg_header:
						if col.startswith("PN:"):
							col_stripped = col.rstrip("\n")
							sam_file_obj.software = col_stripped[3:]
						elif col.startswith("CL:"):
							col_stripped = col.rstrip("\n")
							sam_file_obj.dataset_name = sam_file_obj.parse_dataset_name(col_stripped)
							sam_file_obj.params = col_stripped
				sam_file_obj.headers.append(line)
			else:
				data = line.split("\t")
				try:
					ref_id, cigar_string = data[2], data[5]
				except:
					print "[parse_sam_file] WARNING: had a problem with %s, skipping..." % sam_file
					continue
				try:
					ref_length = len(ref_seqs[ref_id].sequence)
					cigar_len = get_cigar_len(cigar_string)
					percent_aligned = float(cigar_len) / ref_length * 100
					aligned_len_map[ref_id].append(percent_aligned)
				except KeyError:
					print "[parse_sam_file] ERROR: reference id %s found in SAM file but not in the reference FASTA. Did you pass in the correct reference file?" % ref_id
					sys.exit(1)
	return aligned_len_map, sam_file_obj

def taxonomic_distance(optimal_lineage, query_lineage):
	"""
	calculates the taxonomic distance between the optimal and query lineage,
	where optimal_lineage is the starting point of comparison. lineages
	are represented as lists, where each element is a taxonomic rank. the 
	highest taxonomic rank is the 0th element, second highest is 1st element, etc.

	@return the taxonomic distance between lineage_1 and lineage_2
	"""
	if not optimal_lineage or not query_lineage:
		# empty lineages == couldn't assign == return maximum distance
		return -1
	distance = 0
	n, m = len(optimal_lineage), len(query_lineage)
	for i in range(min(n, m)):
		if optimal_lineage[i] != query_lineage[i]:
			distance += 1
	return distance

def reads_over_threshold(reads, threshold):
	"""
	counts the number of reads in a list of reads over the given threshold

	@param reads a list of reads aligned to a reference as a percentage
	@param threshold the minimum percent aligned
	@return the number of reads aligning to a reference >= threshold 
	"""
	return sum(1 for read_p in reads if read_p >= threshold)

def calculate_distance(full_ref_lineage, reads_aligned_p, optimal_placement):
	"""
	calculate the distance between a reference aligned to and the optimal placement 
	for a given dataset. distance is measured as:

		taxonomic_distance(ref_aligned_to, optimal_placement) * number of reads aligned >= 80%

	@return a weighted distance between the reference aligned to and the optimal placment
	"""
	taxa_distance = taxonomic_distance(optimal_placement, full_ref_lineage)
	reads = reads_over_threshold(reads_aligned_p, 80)
	return taxa_distance * reads

def main():
	args = get_args()
	dataset_to_lca = defaultdict(dict)
	print "[main] INFO: gathering reference taxonomies from reference file..."
	ref_seq_map = build_ref_seq_map(args.reference_fasta)
	print "[main] INFO: getting LCAs between dataset organism and all references..."
	with open(args.dataset_taxonomies, "r") as f:
		for line in f:
			dataset_id, taxonomy = line.strip().split("\t")
			start_time = time.time()
			for ref_record in ref_seq_map.values():
				dataset_to_lca[dataset_id][ref_record.id] = query_jgi_lca(taxonomy, ref_record.taxonomy)
			print "[main] INFO: done for dataset %s. time elapsed: %ds" % (dataset_id, time.time() - start_time)
	print "[main] INFO: getting optimal placements (deepest LCA) for each dataset..."
	optimal_placements = get_longest_lca(dataset_to_lca)
	print "[main] INFO: running through each sam file for each dataset..."
	all_sam_files = []
	for item in list_dir_abs(args.sample_dir):
		if os.path.isdir(item):
			for sam_file in list_dir_abs(item):
				if sam_file.endswith("_top_hits.sam"): # TODO: may need to make this more flexible
					all_sam_files.append(parse_sam_file(sam_file, ref_seq_map))
	print "[main] INFO: calculating distances between each reference aligned to and its lca..."
	outlines = []
	for v in all_sam_files:
		# get full lineages of each reference aligned to 
		# calculate distance between that and lca
		# multiply by number of reads >= 80% (aligned_len_map in sam_obj)
		aligned_len_map, sam_obj = v[0], v[1]
		for ref_id in aligned_len_map:
			full_ref_lineage = query_jgi_taxa(ref_seq_map[ref_id].taxonomy)
			# print "###########################################################"
			# print "[main] INFO: full ref aligned to is lineage is: %s" % full_ref_lineage
			# print "[main] INFO: optimal placement is: %s" % optimal_placements[sam_obj.dataset_name]
			# print optimal_placements.keys()
			
			cumulative_distance = calculate_distance(full_ref_lineage, aligned_len_map[ref_id], optimal_placements[sam_obj.dataset_name])
			taxa_distance = taxonomic_distance(full_ref_lineage, optimal_placements[sam_obj.dataset_name])

			# print "[main] INFO: DISTANCE IS: %d" % distance
			# print "###########################################################"
			taxa = ref_seq_map[ref_id].taxonomy
			num_reads_aligned_total = len(aligned_len_map[ref_id])
			mean_alignment_length = sum(aligned_len_map[ref_id]) / num_reads_aligned_total
			num_reads_aligned_p80 = sum(1 for p in aligned_len_map[ref_id] if p >= 80)
			reads_aligned_p80 = num_reads_aligned_p80 / num_reads_aligned_total * 100
			line = [args.marker_gene, sam_obj.dataset_name, ref_id, sam_obj.software, taxa, str(mean_alignment_length), str(num_reads_aligned_total), str(num_reads_aligned_p80), str(reads_aligned_p80), str(taxa_distance), str(cumulative_distance)]
			outlines.append(line)

	out_header = "marker_gene\tdataset_id\tref_id\tsoftware\ttaxa\tmean_alignment_length\tnum_reads_aligned_total\tnum_reads_aligned_p80\treads_aligned_p80\ttaxa_distance\tcumulative_distance"	
	outfile = open(os.path.join(args.output_dir, args.output_file + ".txt"), "w")
	outfile.write(out_header + "\n")
	for line in outlines:
		outfile.write("\t".join(line) + "\n")
	outfile.close()
	# for dataset_id in optimal_placements:
	# 	outname = dataset_id + ".txt"
	# 	outfile = open(os.path.join(args.output_dir, outname), "w")
	# 	outfile.write(out_header + "\n")
	# 	for line in outlines:
	# 		if line[0] == dataset_id:
	# 			outfile.write("\t".join(line) + "\n")
	# 	outfile.close()
	print "[main] INFO: finished script. goodbye."

if __name__ == "__main__":
	main()
