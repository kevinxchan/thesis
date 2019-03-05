#!/usr/bin/bash

# automation of running TreeSAPP on all datasets.
# args:
# 1. absolute path to working directory
# 2. dataset names file (single dataset name per line, stored in .txt file)
# 3. full path to the reference
# USAGE: 
# /home/kchan/scripts_thesis/shell/treesapp.sh /home/kchan/thesis/dataset_names.txt /home/kchan/TreeSAPP /home/kchan/thesis/raw_data /home/kchan/thesis/marker_genes.txt /home/kchan/ts_out
dataset_names=$1
treesapp_dir=$2
input_dir=$3
marker_genes=$4
out_dir=$5

THREADS=12

if [[ -z $1 || -z $2 || -z $3 || -z $4 || -z $5 ]]; then
	echo "argument error: pass in list of dataset names, TreeSAPP directory, input and output path"
	exit 1
fi

echo "dataset names: $1"
echo "TreeSAPP directory: $2"
echo "input directory: $3"
echo "marker gene used: $4"
echo "output directory: $5"

mkdir -p $out_dir
source activate py36

while read name; do
	while read -r marker_id marker_name; do
		mkdir -p "$out_dir"/"$name"/"$marker_name"

		$treesapp_dir/treesapp.py -i $input_dir/$name.fastq -t $marker_id \
		--lr --molecule dna -o $out_dir/$name/$marker_name -T $THREADS --overwrite
	done < $marker_genes
done < $dataset_names
