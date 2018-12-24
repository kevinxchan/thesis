#!/usr/bin/bash

# automation of running LAST.
# args:
# 1. absolute path to working directory
# 2. dataset names file (single dataset name per line, stored in .txt file)
# 3. full path to the reference
# USAGE: 
# time /home/kchan/scripts_thesis/shell/run_last.sh /home/kchan/thesis /home/kchan/thesis/dataset_names.txt /home/kchan/thesis/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta &> /home/kchan/thesis/processed/last/last_log.log
WORK_DIR=$1
dataset_names=$2
reference=$3
THREADS=8
CURR_DIR=$PWD

LASTDB=$(which lastdb)
LAST_TRAIN=$(which last-train)
LASTAL=$(which lastal)
REFORMAT=$(which reformat.sh)

if [[ -z $1 || -z $2 || -z $3 ]]; then
	echo "argument error: pass in working directory, dataset names file and path to reference"
	exit 1
fi

cd $WORK_DIR

echo
echo "###############"
echo "# CREATING DB #"
echo "###############"
echo

index_path=/home/kchan/thesis/indexes/last/YASS
dbname=$(basename $reference)
dbname="${dbname%.*}"
index_file="$index_path/$dbname"

echo "running with params: lastdb -P 8 -u YASS"
mkdir -p $index_path
$LASTDB -P 8 -u YASS $index_path/$dbname $reference 

outdir="$WORK_DIR/processed/last/$reference"
mkdir -p $outdir

echo
echo "#############################"
echo "# FASTQ TO FASTA CONVERSION #"
echo "#############################"
echo

while read name; do
	echo
	echo "FOR DATASET $name"
	echo

	$REFORMAT in=$WORK_DIR/raw_data/$name.fastq.gz out=$WORK_DIR/raw_data/$name.fasta.gz overwrite
done < $dataset_names

echo
echo "###################"
echo "# TRAINING PARAMS #"
echo "###################"
echo

mkdir -p $outdir/training_params

while read name; do
	echo
	echo "FOR DATASET $name"
	echo

	$LAST_TRAIN -P8 -D100 $index_file $WORK_DIR/raw_data/$name.fasta.gz > $outdir/training_params/$name.par
done < $dataset_names

echo
echo "aligning reads with..."
echo

echo
echo "########################"
echo "# PARAMS FROM TRAINING #"
echo "########################"
echo

while read name; do
	echo
	echo "FOR DATASET $name"
	echo

	$LASTAL -P8 -f BlastTab+ -p $outdir/training_params/$name.par $index_file $WORK_DIR/raw_data/$name.fasta.gz > $outdir/$name.txt
done < $dataset_names

cd $CURR_DIR

echo
echo "############"
echo "# FINISHED #"
echo "############"
echo
echo "goodbye."



