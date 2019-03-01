import os
from Bio import Entrez, SeqIO
from argparse import ArgumentParser

def get_args():
	parser = ArgumentParser(description = "Script to scrape the NCBI website and download reference sequences.")
	parser.add_argument("-s", "--sources", required = True, help = "Path to the tab delimited sources metadata.")
	parser.add_argument("-o", "--output-dir", default = ".", help = "Path to the output directory. [current working directory]")
	args = parser.parse_args()
	return args

def scrape_ncbi(accessions, species_name, outfile):
	handle = Entrez.efetch(db = "nuccore", id = ",".join(accessions), rettype = "fasta", retmode = "text")
	for record in SeqIO.parse(handle, "fasta"):
		header = "{} {}".format(record.id, species_name)
		outfile.write(">{}\n{}\n".format(header, record.seq))
	handle.close()

def main():
	args = get_args()

	with open(args.sources, "r") as infile:
		for i, line in enumerate(infile):
			if i != 0:
				data = line.split("\t")
				species_name = data[1].replace(" ", "_")

				if data[6] != "":
					if data[5] != "":
						species_name += "_{}".format(data[5])
					outpath = os.path.join(args.output_dir, "{}.fasta".format(species_name))
					outfile = open(outpath, "w")
					accession_list = data[6].split(" ")
					scrape_ncbi(accession_list, species_name, outfile)
					outfile.close()

main()