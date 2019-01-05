#!/usr/bin/bash

# automation of running graphmap.
# args:
# 1. absolute path to working directory
# 2. dataset names file (single dataset name per line, stored in .txt file)
# 3. full path to the reference
# USAGE: 
# time /home/kchan/scripts_thesis/shell/run_graphmap.sh /home/kchan/thesis /home/kchan/thesis/dataset_names.txt /home/kchan/thesis/references/rpoB/fungene_9.6_rpoB_nucleotide_uclust99.fasta rpoB &> graphmap_rpoB_log.log
WORK_DIR=$1
dataset_names=$2
reference=$3
marker_gene=$4
THREADS=8
CURR_DIR=$PWD

GRAPHMAP=$(which graphmap)

if [[ -z $1 || -z $2 || -z $3 || -z $4 ]]; then
	echo "argument error: pass in working directory, dataset names file, path to reference and marker gene name"
	exit 1
else
	echo "starting script"
	echo "working directory: $WORK_DIR"
	echo "dataset names: $dataset_names"
	echo "reference file path: $reference"
	echo "marker gene name: $marker_gene"
fi

cd $WORK_DIR


echo
echo "building index..."
echo

# build the index. will be stored in index_path
index_path=/home/kchan/thesis/indexes/graphmap
$GRAPHMAP align -I -r $reference && mv "$reference".gmidx $index_path 
mv $index_path/$reference.gmidx $index_path/$marker_gene.gmidx
index_file=$index_path/"$marker_gene.gmidx"

echo
echo "aligning reads with..."
echo

echo
echo "##################"
echo "# DEFAULT PARAMS #"
echo "##################"
echo

while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/$marker_gene/graphmap/$name

	$GRAPHMAP align --threads $THREADS -r $reference --index $index_file --extcigar \
	--reads $WORK_DIR/raw_data/$name.fastq.gz \
	-o $WORK_DIR/processed/$marker_gene/graphmap/$name/default.sam
done < $dataset_names

cd $CURR_DIR

echo
echo "############"
echo "# FINISHED #"
echo "############"
echo
echo "goodbye."



