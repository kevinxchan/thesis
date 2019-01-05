#!/usr/bin/bash

# automation of running LAST.
# args:
# 1. absolute path to working directory
# 2. dataset names file (single dataset name per line, stored in .txt file)
# 3. full path to the reference
# 4. marker gene name
# USAGE: 
# time /home/kchan/scripts_thesis/shell/run_last.sh /home/kchan/thesis /home/kchan/thesis/dataset_names.txt /home/kchan/thesis/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta &> /home/kchan/thesis/processed/last/last_log.log
WORK_DIR=$1
dataset_names=$2
reference=$3
marker_gene=$4
THREADS=8
CURR_DIR=$PWD

LASTDB=$(which lastdb)
LAST_TRAIN=$(which last-train)
LASTAL=$(which lastal)
REFORMAT=$(which reformat.sh)

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
echo "###############"
echo "# CREATING DB #"
echo "###############"
echo

index_path=/home/kchan/thesis/indexes/last/
dbname=$(basename $reference)
dbname="${dbname%.*}"
index_file="$index_path/$dbname"

echo "running with params: lastdb -P 8 -u NEAR"
mkdir -p $index_path
$LASTDB -P 8 -u NEAR $index_path/$dbname $reference 

outdir="$WORK_DIR/processed/last/$dbname"
mkdir -p $outdir

# echo
# echo "#############################"
# echo "# FASTQ TO FASTA CONVERSION #"
# echo "#############################"
# echo

# while read name; do
# 	echo
# 	echo "FOR DATASET $name"
# 	echo

# 	$REFORMAT in=$WORK_DIR/raw_data/$name.fastq.gz out=$WORK_DIR/raw_data/$name.fasta.gz overwrite
# done < $dataset_names

# echo
# echo "###################"
# echo "# TRAINING PARAMS #"
# echo "###################"
# echo

# mkdir -p $outdir/training_params

# while read name; do
# 	echo
# 	echo "FOR DATASET $name"
# 	echo

# 	$LAST_TRAIN -P8 -D100 $index_file $WORK_DIR/raw_data/$name.fasta.gz > $outdir/training_params/$name.par
# done < $dataset_names

echo
echo "aligning reads with..."
echo

echo
echo "########################"
echo "# PARAMS FROM TRAINING #"
echo "########################"
echo

mkdir -p $WORK_DIR/processed/$marker_gene/last

while read name; do
	echo
	echo "FOR DATASET $name"
	echo

	$LASTAL -P8 -f BlastTab+ -r6 -q18 -a21 -b9 $index_file $WORK_DIR/raw_data/$name.fasta.gz > $WORK_DIR/processed/$marker_gene/last/$name.txt
done < $dataset_names

cd $CURR_DIR

echo
echo "############"
echo "# FINISHED #"
echo "############"
echo
echo "goodbye."



