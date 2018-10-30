#!/usr/bin/bash

# automation of running graphmap.
# args:
# 1. absolute path to working directory
# 2. dataset names file (single dataset name per line, stored in .txt file)
# 3. full path to the reference
# USAGE: 
# /home/kchan/scripts_thesis/run_graphmap.sh /home/kchan/thesis /home/kchan/thesis/dataset_names.txt /home/kchan/thesis/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
WORK_DIR=$1
dataset_names=$2
reference=$3
THREADS=5

GRAPHMAP=$(which graphmap)

if [[ -z $1 || -z $2 || -z $3 ]]; then
	echo "argument error: pass in working directory, dataset names file and path to reference"
	exit 1
fi

cd $WORK_DIR

# build the index. will be stored in index_path
index_path=/home/kchan/thesis/indexes/graphmap
echo "checking for index file..."
if [ -z "$(ls $index_path)" ]; then
	echo "not found, building index..."
	$GRAPHMAP align -I -r $reference && mv $(basename "$reference".gmidx) $index_path 
else
	echo "found index file, skipping..."
fi

echo
echo "aligning reads with..."
echo

echo
echo "##################"
echo "# DEFAULT PARAMS #"
echo "##################"
echo

index_file=$index_path/"$(basename $reference.gmidx)"

while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/graphmap/$name

	$GRAPHMAP align -r $reference --index $index_file --reads $WORK_DIR/raw_data/$name.fastq.gz -o $WORK_DIR/processed/graphmap/$name/default.sam --extcigar --threads $THREADS

done < $dataset_names

echo
echo "############"
echo "# FINISHED #"
echo "############"
echo
echo "goodbye."



