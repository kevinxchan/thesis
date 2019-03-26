import os
from argparse import ArgumentParser
from Bio import SeqIO
from file_utils import list_dir_abs

def get_args():
	parser = ArgumentParser(description = "Script to reformat fasta files so headers match the filename, with an increment.")
	parser.add_argument("-r", "--reference-dir", required = False, help = "Directory containing folders for each reference organism.")
	parser.add_argument("-o", "--output-dir", default = ".", help = "Path to the output directory.")
	args = parser.parse_args()
	return args

def main():
	args = get_args()
	for f in list_dir_abs(args.reference_dir):
		if f.endswith(".fasta"):
			filename = os.path.splitext(os.path.basename(f))[0]
			outfile = "{}_formatted.fasta".format(filename) if args.output_dir == "." else "{}.fasta".format(filename)
			outfile = open(os.path.join(args.output_dir, outfile), "w")
			counter = 1
			for record in SeqIO.parse(f, "fasta"):
				new_header = "{}_{}".format(filename, counter)
				outfile.write(">{}\n{}\n".format(new_header, record.seq.upper()))
				counter += 1
			outfile.close()

if __name__ == "__main__":
	main()