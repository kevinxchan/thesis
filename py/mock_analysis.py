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
from math import sqrt
from collections import defaultdict, Counter
from argparse import ArgumentParser
from Bio import SeqIO
from util.file_utils import list_dir_abs
from util.jgi_utils import query_jgi_taxa
from model.mock_analysis_classes import *

def get_args():
	parser = ArgumentParser(description = "Script to calculate MCC using a mock community and TreeSAPP.")
	parser.add_argument("-f", "--minimap2-file", required = True, help = "Path to the minimap2 output in PAF format.")
	parser.add_argument("-r", "--reference-dir", required = True, help = "Directory containing folders for each reference organism.")
	parser.add_argument("-t", "--dataset-taxonomies", required = True, help = "Tab delimited file matching dataset names to their reported taxonomies.")
	parser.add_argument("-d", "--data-dir", required = True, help = "Directory containing the tax ids for each marker gene. Should be in the $TreeSAPP directory.")	
	parser.add_argument("-m", "--marker-genes", required = True, help = "Comma delimited list of marker genes tested.")	
	parser.add_argument("--raw-data", required = True, help = "Path to raw FASTQ mock community data.")
	parser.add_argument("--treesapp-hits", required = True, help = "Path to the TreeSAPP marker_contig_map.tsv file.")
	parser.add_argument("--max-taxa-distance", required = True, help = "Maximum taxonomic distance allowed before counting as a false positive.")
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
			hit_name, start, stop, strand, ref_gene, ref_name, optimal_lineage = line.strip().split("\t")
			gff_line = GFFLine(hit_name, int(start), int(stop), strand, ref_gene, ref_name, optimal_lineage)
			d[hit_name].append(gff_line)
	return d

def read_marker_overlaps(marker_gene_overlaps):
	d = defaultdict(list)
	with open(marker_gene_overlaps, "r") as f:
		for line in f:
			query_name, gene, reference_name, optimal_placement = line.strip().split("\t")
			d[query_name].append(Overlap(query_name, gene, reference_name, optimal_placement))
	return d

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def get_dataset_lineages(dataset_to_taxon):
	ret = {}
	with open(dataset_to_taxon, "r") as infile:
		for line in infile:
			dataset_accession, taxa_string = line.split("\t")
			full_taxa_string = query_jgi_taxa(taxa_string)
			ret[dataset_accession] = full_taxa_string
	return ret

def read_marker_contig_map(marker_contig_map):
	# TODO: this will need to change when marker_contig_map.tsv is properly formatted
	d = {}
	with open(marker_contig_map, "r") as infile:
		for i, line in enumerate(infile):
			if i == 0:
				continue
			line = line.strip().split("\t")
			query_name = line[1]
			marker = line[2]
			confident_taxa = line[5]
			mcm = MarkerContigMap(query_name, marker, confident_taxa)
			d[query_name] = mcm
	return d

def write_full_lineages(full_dataset_lineage, outdir):
	outfile = open(os.path.join(outdir, "full_dataset_lineages.txt"), "w")
	for _id, lineage in full_dataset_lineage.items():
		outfile.write("{}\t{}\n".format(_id, lineage))

def write_optimal_placements(placements, outdir):
	outfile = open(os.path.join(outdir, "optimal_placements.txt"), "w")
	for _id in placements:
		for marker in placements[_id]:
			outfile.write("{}\t{}\t{}\n".format(_id, marker, placements[_id][marker]))

def write_marker_gene_overlaps(marker_genes, outdir):
	outfile = open(os.path.join(outdir, "marker_gene_overlaps.txt"), "w")
	for overlap_obj_list in marker_genes.values():
		for overlap_obj in overlap_obj_list:
			outfile.write("{}\t{}\t{}\t{}\n".format(overlap_obj.query_name, overlap_obj.gene, overlap_obj.reference_name, overlap_obj.optimal_placement))

def write_parsed_gffs(parsed_gffs_file, output_dir):
	outpath = os.path.join(output_dir, "parsed_gffs.txt")
	outfile = open(outpath, "w")
	for v in parsed_gffs_file.values():
		for gff_line in v:
			outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
				gff_line.hit_name, gff_line.start, gff_line.stop, gff_line.strand, gff_line.ref_gene, 
				gff_line.ref_name, gff_line.optimal_lineage))

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

def get_read_names(mock_fastq):
	s = set()
	file = open(mock_fastq, "r")
	for name, _, _ in readfq(file):
		s.add(name)
	return s

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

def parse_gffs(reference_path, dataset_optimal_placement):
	d = defaultdict(list)
	for reference in dataset_optimal_placement:
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
								gff_line = GFFLine(hit_name, start-1, stop-1, strand, ref_gene, reference, dataset_optimal_placement[reference][ref_gene])
								d[hit_name].append(gff_line)
	return d

# def parse_paf(paf_file):
#     """
#     Parse a Pairwise mApping Format (PAF) file, storing alignment information (e.g. read name, positions)
#     for reads that were mapped. Filters reads by maximum observed mapping quality. Assumes alignments
#     are sorted by read name (hence multiple alignments for a single read are successive).
#     :param paf_file: Path to a PAF file
#     :return: A dictionary mapping reference packages a list of PAF objects
#     """
#     hit_name_to_read = defaultdict(list)
#     with open(paf_file, "r") as infile:
#         prev_qname = ""
#         for line in infile:
#             data = line.split("\t")
#             paf_obj = PAFObj(data[0], int(data[1]), int(data[2]), int(data[3]), data[4], data[5], int(data[6]), int(data[7]),
#                              int(data[8]), int(data[9]), int(data[10]), int(data[11]))

#             hit_name = data[5]
#             if prev_qname != data[0]:
#                 if 0 < int(data[11]) < 255:
#                     hit_name_to_read[hit_name].append(paf_obj)
#             else:
#                 # filter reads by maximum mapping quality
#                 if len(hit_name_to_read[hit_name]) > 0:
#                     stored_mapq = hit_name_to_read[hit_name][-1].mapq
#                     if int(data[11]) > stored_mapq:
#                         hit_name_to_read[hit_name][-1] = paf_obj
#             prev_qname = data[0]
#     return hit_name_to_read

def parse_paf(paf_file):
    """
    Parse a Pairwise mApping Format (PAF) file, storing alignment information (e.g. read name, positions)
    for reads that were mapped. Filters reads by maximum observed mapping quality. Assumes alignments
    are sorted by read name (hence multiple alignments for a single read are successive).
    :param paf_file: Path to a PAF file
    :return: A dictionary of reference package names mapping to read names, mapping to a list of PAF objects
    """
    hit_name_to_read = defaultdict(list)
    parser_dict = defaultdict(set)
    missing_mapq_value = 255
    unmapped_mapq_value = 0
    with open(paf_file, "r") as infile:
        prev_qname = ""
        for line in infile:
            data = line.strip().split("\t")
            qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, n_match_bases, n_total_bases, mapq = data[:12]
            qlen = int(qlen)
            qstart = int(qstart)
            qend = int(qend)
            tlen = int(tlen)
            tstart = int(tstart)
            tend = int(tend)
            n_match_bases = int(n_match_bases)
            n_total_bases = int(n_total_bases)
            mapq = int(mapq)

            # upfront filter for unacceptable mapping qualities
            if mapq == missing_mapq_value or mapq == unmapped_mapq_value:
                continue
            else:
                paf_obj = PAFObj(qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, n_match_bases,
                                 n_total_bases, mapq)
                if prev_qname != qname:
                    parser_dict[qname].add(paf_obj)
                else:
                    # filter reads by maximum mapping quality if they overlap
                    is_curr_overlapping = False
                    for stored_paf_obj in parser_dict[qname].copy():
                        if stored_paf_obj.is_overlapping(qstart, qend):
                            is_curr_overlapping = True
                            stored_mapq = stored_paf_obj.mapq
                            if mapq > stored_mapq:
                                parser_dict[qname].remove(stored_paf_obj)
                                parser_dict[qname].add(paf_obj)
                    
                    if not is_curr_overlapping:
                        parser_dict[qname].add(paf_obj)

                prev_qname = qname

    for qname in parser_dict:
    	for stored_paf_obj in parser_dict[qname]:
    		hit_name_to_read[stored_paf_obj.tname].append(stored_paf_obj)
    return hit_name_to_read

def find_overlaps(minimap2_hits, gff_hits):
	marker_genes = defaultdict(list)
	for hit_name, alignments in minimap2_hits.items():
		for single_alignment in alignments:
			t_start, t_stop = single_alignment.tstart, single_alignment.tend
			for marker_hit in gff_hits[hit_name]:
				ref_start, ref_stop = marker_hit.start, marker_hit.stop
				if (t_start < ref_stop and t_stop > ref_start) or (ref_start < t_stop and ref_stop > t_start):
					marker_genes[single_alignment.qname].append(Overlap(single_alignment.qname, marker_hit.ref_gene, marker_hit.ref_name, marker_hit.optimal_lineage))
	return marker_genes

def within_taxa_distance(query_taxa, optimal_taxa, max_dist):
	query = query_taxa.split("; ")
	optimal = optimal_taxa.split("; ")
	distance = 0
	for i in range(max(len(query), len(optimal))):
		if i >= len(query) or i >= len(optimal):
			distance += 1
		else:
			if query[i] != optimal[i]:
				distance += 1
	return distance <= max_dist

def calculate_positives(marker_contig_map, marker_genes, max_taxa_dist):
	tp = fp = 0
	for mcm_obj in marker_contig_map.values():
		query_name = mcm_obj.query_name.split("_")[0] # remove start/stop coordinates
		if query_name not in marker_genes:
			fp += 1
		else:
			have_true_pos = False
			true_marker_list = marker_genes[query_name]
			for true_marker in true_marker_list:
				if true_marker.gene == mcm_obj.marker and within_taxa_distance(mcm_obj.confident_taxonomy, true_marker.optimal_placement, max_taxa_dist):
					tp += 1
					have_true_pos = True
					break
			if not have_true_pos:
				fp += 1
	return (tp, fp)

def calculate_negatives(marker_contig_map, marker_genes, all_read_names_set, num_true_positives, num_false_positives):
	true_marker_genes_count = sum(len(v) for v in marker_genes.values()) # total number of identified marker genes
	fn = true_marker_genes_count - num_true_positives
	tn = len(all_read_names_set) - (fn + num_true_positives + num_false_positives)
	return (tn, fn)

def calculate_MCC(tp, fp, tn, fn):
	numerator = (tp * tn) - (fp * fn)
	denominator = sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
	try:
		return numerator / denominator
	except ZeroDivisionError:
		logging.debug("MCC calculation resulted in division by 0!")
		print("TP: {}".format(tp))
		print("FP: {}".format(fp))
		print("TN: {}".format(tn))
		print("FN: {}".format(fn))
		sys.exit(1)

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
		hit_name_gff_map = read_parsed_gffs(cached_parsed_gffs)
	else:
		logging.info("parsing GFF files for marker genes...")
		hit_name_gff_map = parse_gffs(args.reference_dir, dataset_optimal_placement)
		write_parsed_gffs(hit_name_gff_map, args.output_dir)
	
	cached_marker_gene_hits = os.path.join(args.output_dir, "marker_gene_overlaps.txt")
	if os.path.isfile(cached_marker_gene_hits):
		logging.info("found file with reads mapping to marker gene predictions in output directory, reading...")
		marker_genes = read_marker_overlaps(cached_marker_gene_hits)
		logging.info("done.")
	else:
		# parse through minimap2 output, find overlaps
		logging.info("parsing through minimap2 file, saving top hits...")
		minimap2_file = os.path.abspath(args.minimap2_file)
		hit_name_to_read = parse_paf(minimap2_file)
		logging.info("done.")

		# find overlaps between minimap2 hits and predicted genes in GFF files
		logging.info("finding overlaps between minimap2 output and predicted marker genes in GFF files...")
		marker_genes = find_overlaps(hit_name_to_read, hit_name_gff_map)
		write_marker_gene_overlaps(marker_genes, args.output_dir)
		logging.info("done.")

	logging.info("parsing through marker contig map...")
	marker_contig_map = os.path.abspath(args.treesapp_hits)
	marker_contig_map = read_marker_contig_map(marker_contig_map)
	logging.info("done.")

	logging.info("parsing through fastq mock community and extracting read names...")
	mock_fastq = os.path.abspath(args.raw_data)
	all_read_names_set = get_read_names(mock_fastq)
	logging.info("done.")
	
	max_taxa_dist = int(args.max_taxa_distance)
	logging.info("counting true positives and false positives...")
	tp, fp = calculate_positives(marker_contig_map, marker_genes, max_taxa_dist)
	logging.info("done.")

	logging.info("counting true negatives and false negatives...")
	tn, fn = calculate_negatives(marker_contig_map, marker_genes, all_read_names_set, tp, fp)
	logging.info("done.")

	logging.info("calculating MCC...")
	MCC = calculate_MCC(tp, fp, tn, fn)
	logging.info("done.")

	print()
	print("FINAL OUTPUT:")
	print("TRUE POSITIVES: {}".format(tp))
	print("FALSE POSITIVES: {}".format(fp))
	print("TRUE NEGATIVES: {}".format(tn))
	print("FALSE NEGATIVES: {}".format(fn))
	print()
	print("MCC: {}".format(MCC))
	
	logging.info("\n###### DONE. GOODBYE ######\n")
if __name__ == "__main__":
	main()
