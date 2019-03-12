import requests
import os

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
	return "Root; " + "; ".join(cleaned_lineage)
