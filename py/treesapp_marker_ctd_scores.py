"""
script to gather CTD per marker genes, so a scatter plot of CTD vs marker gene can be made.

formula for CTD:
let d = distance from the optimal placement of a dataset to the marker gene M
let p = fraction of reads in a dataset aligning to M

then for M, and a dataset, CTD = sum d * p. therefore, a single data point represents a single dataset.
"""

import os
import logging
import requests
from collections import defaultdict
from argparse import ArgumentParser
from util.file_utils import list_dir_abs

def get_args():
	parser = ArgumentParser(description = "Script to gather CTD per marker genes, so a scatter plot of CTD vs marker gene can be made.")
	parser.add_argument("-s", "--sample-dir", required = True, help = "Directory containing folders for each sample.")
	parser.add_argument("-t", "--dataset-taxonomies", required = True, help = "Tab delimited file matching dataset names to their reported taxonomies.")
	parser.add_argument("-d", "--data-dir", required = True, help = "Directory containing the tax ids for each marker gene. Should be in the $TreeSAPP directory.")	
	parser.add_argument("-m", "--marker-genes", required = True, help = "Comma delimited list of marker genes tested.")	
	parser.add_argument("-o", "--output-dir", default = ".", help = "Path to the output directory.")
	args = parser.parse_args()
	return args

def read_lineages(lineage_file):
	d = {}
	with open(lineage_file, "r") as f:
		for line in f:
			_id, lineage = line.strip().split("\t")
			d[_id] = lineage
	return d

def read_placements(placements_file):
	d = defaultdict(dict)
	with open(placements_file, "r") as f:
		for line in f:
			_id, marker, lineage = line.strip().split("\t")
			d[_id][marker] = lineage
	return d

def read_sample_dir(sample_dir):
	d = defaultdict(dict)
	for item in list_dir_abs(sample_dir):
		if os.path.isdir(item):
			try:
				with open(os.path.join(item, "final_outputs", "marker_contig_map.tsv"), "r") as infile:
					for i, line in enumerate(infile):
						if i == 0:
							# header, skip
							continue
						else:
							line = line.strip().split("\t")
							sample, marker, confident_taxa = line[0], line[2], line[5]
							d[sample][marker] = confident_taxa
			except IOError:
				logging.warning("directory {} doesn't appear to be an output from TreeSAPP, skipping...".format(item))
	return d

def get_dataset_lineages(dataset_to_taxon):
	ret = {}
	with open(dataset_to_taxon, "r") as infile:
		for line in infile:
			dataset_accession, taxa_string = line.split("\t")
			full_taxa_string = query_jgi_taxa(taxa_string)
			ret[dataset_accession] = full_taxa_string
	return ret

def write_full_lineages(full_dataset_lineage, outdir):
	outfile = open(os.path.join(outdir, "full_dataset_lineages.txt"), "w")
	for _id, lineage in full_dataset_lineage.items():
		outfile.write("{}\t{}\n".format(_id, lineage))

def write_optimal_placements(placements, outdir):
	outfile = open(os.path.join(outdir, "optimal_placements.txt"), "w")
	for _id in placements:
		for marker in placements[_id]:
			outfile.write("{}\t{}\t{}\n".format(_id, marker, placements[_id][marker]))

def query_jgi_taxa(query):
	url = "http://taxonomy.jgi-psf.org/name/sc"
	r = requests.get(os.path.join(url, query))
	return clean_lineage(r.text) 

def clean_lineage(jgi_lineage):
	cleaned_lineage = []
	ranks = jgi_lineage.split(";")
	for rank in ranks:
		if ":" in rank:
			rank = rank.split(":")[1]
			cleaned_lineage.append(rank)
	return "; ".join(cleaned_lineage)

def get_optimal_placements(full_dataset_lineage, marker_file):
	depth = -1
	with open(marker_file, "r") as infile: 
		for line in infile:
			line = line.strip().split("\t")
			reference_lineage = line[2]
			curr_depth, lca = get_deepest_lca(full_dataset_lineage, reference_lineage)
			if depth == -1 or curr_depth > depth:
				optimal_placement = lca
				depth = curr_depth
	return optimal_placement

def get_deepest_lca(dataset_lineage, reference_lineage):
	depth, lca = 0, ""
	l1, l2 = dataset_lineage.split("; "), reference_lineage.split("; ")
	for i in range(min(len(l1), len(l2))):
		if l1[i] == l2[i]:
			depth += 1
			lca = l1[:i+1]
	return depth, "; ".join(lca)

def main():
	logging.basicConfig(level = logging.DEBUG)
	logging.info("\n###### STARTING SCRIPT ######\n")
	args = get_args()
	cached_lineages = os.path.join(args.output_dir, "full_dataset_lineages.txt")

	# read in full lineages for each dataset
	if os.path.isfile(cached_lineages):
		logging.info("found file mapping dataset ids to full lineages in output directory, reading...")
		full_dataset_lineage = read_lineages(cached_lineages)
	else:
		logging.info("gathering full lineages for each dataset...")
		full_dataset_lineage = get_dataset_lineages(args.dataset_taxonomies)
		write_full_lineages(full_dataset_lineage, args.output_dir)
	logging.info("done.")

	# gather optimal placements for each dataset, and each marker gene
	cached_optimal_placements = os.path.join(args.output_dir, "optimal_placements.txt")
	if os.path.isfile(cached_optimal_placements):
		logging.info("found file mapping dataset ids to optimal placements for each marker, reading...")
		dataset_optimal_placement = read_placements(cached_optimal_placements)
	else:
		markers = args.marker_genes.split(",")
		dataset_optimal_placement = defaultdict(dict)
		for dataset in full_dataset_lineage:
			logging.info("getting optimal placements for dataset {}...".format(dataset))
			for marker in markers:
				path_to_marker = os.path.join(args.data_dir, "tax_ids_{}.txt".format(marker))
				dataset_optimal_placement[dataset][marker] = get_optimal_placements(full_dataset_lineage[dataset], path_to_marker)
				logging.info("...done for marker gene: {}.".format(marker))
		logging.info("writing optimal placements to disk...")
		write_optimal_placements(dataset_optimal_placement, args.output_dir)
	logging.info("done.")

	# parse info from marker contig map files
	args.sample_dir = os.path.abspath(args.sample_dir)
	marker_contig_map = read_sample_dir(args.sample_dir)
	print(marker_contig_map)
	
	logging.info("\n###### DONE. GOODBYE ######\n")
if __name__ == "__main__":
	main()
