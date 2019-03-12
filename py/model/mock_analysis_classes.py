
class GFFLine:

	def __init__(self, hit_name, start, stop, strand, ref_gene, ref_name, full_lineage):
		self.hit_name = hit_name
		self.start = start
		self.stop = stop
		self.strand = strand
		self.ref_gene = ref_gene
		self.ref_name = ref_name
		self.full_lineage = full_lineage
