
class GFFLine:

	def __init__(self, hit_name, start, stop, strand, ref_gene, ref_name, optimal_lineage):
		self.hit_name = hit_name
		self.start = start
		self.stop = stop
		self.strand = strand
		self.ref_gene = ref_gene
		self.ref_name = ref_name
		self.optimal_lineage = optimal_lineage

class PAFObj:

    def __init__(self, qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, n_match_bases, n_total_bases, mapq):
        self.qname = qname
        self.qlen = qlen
        self.qstart = qstart
        self.qend = qend
        self.strand = strand
        self.tname = tname
        self.tlen = tlen
        self.tstart = tstart
        self.tend = tend
        self.n_match_bases = n_match_bases
        self.n_total_bases = n_total_bases
        self.mapq = mapq

    def is_overlapping(self, other_qstart, other_qend):
    	self_l = (self.qstart < other_qend) and (self.qend > other_qstart)
    	self_r = (other_qstart < self.qend) and (other_qend > self.qstart)
    	return self_l or self_r

class Overlap:

	def __init__(self, query_name, gene, reference_name, optimal_placement):
		self.query_name = query_name
		self.gene = gene
		self.reference_name = reference_name
		self.optimal_placement = optimal_placement

	def __repr__(self):
		return "query name: {}\nmarker gene: {}\nreference name: {}\noptimal placement: {}".format(
			self.query_name, self.gene, self.reference_name, self.optimal_placement)

	def __str__(self):
		return "query name: {}\nmarker gene: {}\nreference name: {}\noptimal placement: {}".format(
			self.query_name, self.gene, self.reference_name, self.optimal_placement)

class MarkerContigMap:

	def __init__(self, query_name, marker, confident_taxonomy):
		self.query_name = query_name
		self.marker = marker
		self.confident_taxonomy = confident_taxonomy
