from collections import defaultdict, OrderedDict

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

