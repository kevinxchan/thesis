"""
script to classify mock community results as TP, FP, TN and FN, and calculate MCC.

formula:
MCC = (TP * TN) - (FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

classifications:
TP = read assigned the correct marker gene (overlap) and maps to the correct organism (within tax distance of 2)
FP = read assigned the incorrect taxonomy (tax distance > 2), or incorrect marker gene
TN = read which was not aligned in TreeSAPP, nor minimap2
FN = read which was not aligned in TreeSAPP, but was in minimap2
"""

import os
import logging
import sys
from collections import defaultdict
from argparse import ArgumentParser
from util.file_utils import list_dir_abs
from util.jgi_utils import query_jgi_taxa
from model.mock_analysis_classes import *

def get_args():
	parser = ArgumentParser(description = "Script to calculate MCC using a mock community and TreeSAPP.")
	parser.add_argument("-f", "--minimap2-file", required = False, help = "Path to the minimap2 output in PAF format.")
	parser.add_argument("-r", "--reference-dir", required = True, help = "Directory containing folders for each reference organism.")
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
	d = defaultdict(lambda: defaultdict(list))
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
							d[sample][marker].append(confident_taxa)
			except IOError:
				logging.warning("directory {} doesn't appear to be an output from TreeSAPP, skipping...".format(item))
	return d

def read_parsed_gffs(parsed_gffs_file):
	d = defaultdict(list)
	with open(parsed_gffs_file, "r") as f:
		for line in f:
			hit_name, start, stop, strand, ref_gene, ref_name, full_lineage = line.strip.split("\t")
			gff_line = GFFLine(hit_name, start, stop, strand, ref_gene, ref_name, full_lineage)
			d[hit_name].append(gff_line)
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

def write_final_scores(scores, outdir):
	outfile = open(os.path.join(outdir, "ctd_scores.txt"), "w")
	outfile.write("Dataset\tMarker\tCTD_score\n")
	for score_line in scores:
		outfile.write("\t".join(map(str, score_line)) + "\n")

def write_parsed_gffs(parsed_gffs_file, output_dir):
	outpath = os.path.join(output_dir, "parsed_gffs.txt")
	outfile = open(outpath, "w")
	for v in parsed_gffs_file.values():
		for gff_line in v:
			outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
				gff_line.hit_name, gff_line.start, gff_line.stop, gff_line.strand, gff_line.ref_gene, 
				gff_line.ref_name, gff_line.full_lineage))

def get_optimal_placements(full_dataset_lineage, marker_file):
	depth = -1
	with open(marker_file, "r") as infile: 
		for line in infile:
			line = line.strip().split("\t")
			reference_lineage = "Root; " + line[2]
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

def check_for_marker(description):
	"""
	parse through GFF lines for markers of interest. this is hardcoded!
	"""
	desc_list = description.split(";")
	for item in desc_list:
		item_lower = item.lower()
		if "reca" in item_lower or "recombinase a" in item_lower:
			return ("recA", True)
		elif "rpob" in item_lower or "rna polymerase subunit beta" in item_lower:
			return ("rpoB", True)
		elif "ribosomal protein s8" in item_lower:
			return ("rps8", True)
		elif "ribosomal protein l10" in item_lower:
			return ("rL10P", True)
		elif "pyrg" in item_lower or "ctp synthase" in item_lower:
			return ("pyrG", True)
	return ("", False)

def parse_gffs(reference_path, full_dataset_lineage):
	d = defaultdict(list)
	for reference in full_dataset_lineage:
		ref_org_path = os.path.join(reference_path, reference)
		for f in list_dir_abs(ref_org_path):
			if f.endswith(".gff"):
				with open(f, "r") as infile:
					for line in infile:
						line = line.strip()
						if line.startswith("##"):
							if line.startswith("##FASTA"):
								break
							else:
								continue
						else:
							line = line.split("\t")
							hit_name, start, stop, strand, desc = line[0], int(line[3]), int(line[4]), line[6], line[8]
							ref_gene, contains_marker = check_for_marker(desc)
							if contains_marker:
								gff_line = GFFLine(hit_name, start-1, stop-1, strand, ref_gene, reference, full_dataset_lineage[reference])
								d[hit_name].append(gff_line)
	return d

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
		logging.info("found file mapping dataset ids to optimal placements for each marker in output directory, reading...")
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

	# parse through all GFFs and collect hits to marker genes of interest
	cached_parsed_gffs = os.path.join(args.output_dir, "parsed_gffs.txt")
	if os.path.isfile(cached_parsed_gffs):
		logging.info("found file containing parsed GFF files with marker genes in output directory, reading...")
		reference_gff_map = read_parsed_gffs(cached_parsed_gffs)
	else:
		logging.info("parsing GFF files for marker genes...")
		reference_gff_map = parse_gffs(args.reference_dir, full_dataset_lineage)
		write_parsed_gffs(reference_gff_map, args.output_dir)
	# parse info from marker contig map files
	# logging.info("parsing through marker contig map and calculating MCC...")
	# args.sample_dir = os.path.abspath(args.sample_dir)
	# marker_contig_map = read_sample_dir(args.sample_dir)
	logging.info("done.")

	logging.info("FINAL OUTPUT:\n")
	
	logging.info("\n###### DONE. GOODBYE ######\n")
if __name__ == "__main__":
	main()
