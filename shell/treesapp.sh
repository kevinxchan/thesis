#!/usr/bin/bash

# automation of running TreeSAPP on all datasets.
# args:
# 1. absolute path to working directory
# 2. dataset names file (single dataset name per line, stored in .txt file)
# 3. full path to the reference
# USAGE: 
# /home/kchan/scripts_thesis/shell/treesapp.sh /home/kchan/thesis/dataset_names.txt /home/kchan/TreeSAPP /home/kchan/thesis/raw_data /home/kchan/ts_out
dataset_names=$1
treesapp_dir=$2
input_dir=$3
out_dir=$4

THREADS=16

if [[ -z $1 || -z $2 || -z $3 || -z $4 ]]; then
	echo "argument error: pass in list of dataset names, TreeSAPP directory, input and output path"
	exit 1
fi

echo
echo "dataset names: $1"
echo "TreeSAPP directory: $2"
echo "input directory: $3"
echo "output directory: $4"
echo

mkdir -p $out_dir
source activate py36

while read name; do
	mkdir -p "$out_dir"/"$name"

	$treesapp_dir/treesapp.py -i $input_dir/$name.fastq \
	--lr --molecule dna -o $out_dir/$name -T $THREADS --overwrite
done < $dataset_names
