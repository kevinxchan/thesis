from collections import defaultdict, OrderedDict
import os
import numpy as np
import matplotlib 
matplotlib.use('Agg')
from matplotlib import pyplot as plt

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

	def parse_dataset_name(self, params_line):
		dataset_path = params_line.split(" ")[-1]
		dataset_path = os.path.basename(dataset_path)
		return os.path.splitext(dataset_path)[0].replace(".fastq", "")

class Histogram():

	def __init__(self, bin_size = 10):
		self.bin_start = 0
		self.bin_end = 100
		self.bin_size = bin_size
		self.counts_per_bin = defaultdict(list)

	def add_to_bin(self, params, dataset_name, ref_id, percent_aligned):
		self.counts_per_bin[(params, dataset_name, ref_id)].append(percent_aligned)

	def plot_histograms(self, filename, ref_seqs):
		bins = np.arange(self.bin_start, self.bin_end + self.bin_size, self.bin_size)
		for params, dataset_name, ref_id in self.counts_per_bin.keys():
			plt.figure()
			plt.hist(self.counts_per_bin[(params, dataset_name, ref_id)], bins = bins, label = "%s (%s)" % (ref_seqs[ref_id].taxonomy, ref_id))
			plt.title("%s (reference organism = %s, ID = %s)" % (dataset_name, ref_seqs[ref_id].taxonomy, ref_id), fontsize = 8)
			plt.xlabel("Bins (bin size = %d)" % self.bin_size)
			plt.ylabel("Count")
			out_taxa = ref_seqs[ref_id].taxonomy.replace(" ", "").replace("/", "")
			plt.savefig("%s_%s_%s.png" % (filename, ref_id, out_taxa))	
			plt.close()
