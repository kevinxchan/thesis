"""
Script to rank alignments based on LCA and proportion of reads mapped over 80% of the reference.

How it works:
1. Query JGI taxonomy server to get full taxonomies of the organisms in the reference fasta.
2. Do 1. again, but with the reported organism in each dataset.
3. Calculate the LCA between 1. and 2., and use this as the comparison point.
4. For each reference aligned to, get the LCA between it and 3. Multiply by # reads aligned > 80% 
5. Repeat for all parameter sets, for all datasets.
6. Pick the param set with the most min. number of scores.

USAGE:
python /home/kchan/scripts_thesis/py/lca.py -r /home/kchan/thesis/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta -s /home/kchan/thesis/processed/minimap2 -t /home/kchan/thesis/dataset_to_taxon.txt
"""

import sys
import os
import requests
import time
from collections import defaultdict, OrderedDict
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
		print "WARNING: reference sequence %s or %s not found in NCBI. assigning an empty LCA." % (entry_1, entry_2)
		return []

	query = ",".join([entry_1, entry_2])
	r = requests.get(os.path.join(url, query))
	try:
		response = r.json()
		error = response.values()[0]
		print "ERROR: query with %s, %s failed with response '%s'. exiting." % (entry_1, entry_2, error)
		sys.exit(1)
	except ValueError:
		response = r.text
	tax_list = response.split(";")
	no_rank_level = 1
	for i in range(len(tax_list)):
		level = tax_list[i]
		if len(level.split(":")) == 1:
			rank = "nr%d" % no_rank_level
			print "WARNING: taxonomic rank not found for level '%s' assigning a rank of '%s'." % (level, rank)
			tax_list[i] = rank + ":" + level
			no_rank_level += 1
	return tax_list

def get_longest_lca(dataset_to_lca):
	"""
	returns the deepest lca between each dataset and the reference sequences.
	# TODO: what if multiple references are of the same length?

	@return deepest lca 
	"""
	optimal_placements = {}
	for _id in dataset_to_lca:
		depth = -1
		optimal_placement = []
		for ref_id in dataset_to_lca[_id]:
			if len(dataset_to_lca[_id][ref_id]) > depth:
				print _id, ref_id
				print dataset_to_lca[_id][ref_id]
				optimal_placement = dataset_to_lca[_id][ref_id]
				depth = len(optimal_placement)
		optimal_placements[_id] = optimal_placement
	return optimal_placements

def main():
	args = get_args()
	dataset_to_lca = defaultdict(dict)
	print "gathering reference taxonomies from reference file..."
	ref_seq_map = build_ref_seq_map(args.reference_fasta)
	print "getting LCAs between dataset organism and all references..."
	with open(args.dataset_taxonomies, "r") as f:
		for line in f:
			dataset_id, taxonomy = line.strip().split("\t")
			start_time = time.time()
			for ref_record in ref_seq_map.values()[0:10]:
				dataset_to_lca[dataset_id][ref_record.id] = query_jgi_lca(taxonomy, ref_record.taxonomy)
			print "done for dataset %s. time elapsed: %ds" % (dataset_id, time.time() - start_time)
	print "getting optimal placements (deepest LCA) for each dataset..."
	optimal_placements = get_longest_lca(dataset_to_lca)
	
	
if __name__ == "__main__":
	main()
