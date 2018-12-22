#!/usr/bin/bash

# automation of running snap-aligner.
# args:
# 1. absolute path to working directory
# 2. dataset names file (single dataset name per line, stored in .txt file)
# 3. full path to the reference
# USAGE: 
# time /home/kchan/scripts_thesis/shell/run_snap.sh /home/kchan/thesis /home/kchan/thesis/dataset_names.txt /home/kchan/thesis/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta &> snap_log.log
WORK_DIR=$1
dataset_names=$2
reference=$3
THREADS=8
CURR_DIR=$PWD

SNAP=$(which snap-aligner)

if [[ -z $1 || -z $2 || -z $3 ]]; then
	echo "argument error: pass in working directory, dataset names file and path to reference"
	exit 1
fi

cd $WORK_DIR
index_path=/home/kchan/thesis/indexes/snap
mkdir -p index_path

echo "building index with seed length 22..."
$SNAP index $reference $index_path/s22 -s 22
echo "building index with seed length 32 (max length)..."
$SNAP index $reference $index_path/s32 -keysize 8 -s 32

echo
echo "aligning reads with..."
echo

echo
echo "#########"
echo "# -s 22 #"
echo "#########"
echo

while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/snap/$name

	$SNAP single $index_path/s22 $WORK_DIR/raw_data/$name.fastq.gz -t $THREADS -= \
	-o $WORK_DIR/processed/snap/$name/s22.sam
done < $dataset_names

echo
echo "###############"
echo "# -s 22 -d 30 #"
echo "###############"
echo

while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/snap/$name

	$SNAP single $index_path/s22 $WORK_DIR/raw_data/$name.fastq.gz -t $THREADS -= -d 30 \
	-o $WORK_DIR/processed/snap/$name/s22_d30.sam
done < $dataset_names

echo
echo "#########"
echo "# -s 32 #"
echo "#########"
echo

while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/snap/$name

	$SNAP single $index_path/s32 $WORK_DIR/raw_data/$name.fastq.gz -t $THREADS -= \
	-o $WORK_DIR/processed/snap/$name/s32.sam
done < $dataset_names

echo
echo "###############"
echo "# -s 32 -d 30 #"
echo "###############"
echo

while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/snap/$name

	$SNAP single $index_path/s32 $WORK_DIR/raw_data/$name.fastq.gz -t $THREADS -= \
	-o $WORK_DIR/processed/snap/$name/s32_d30.sam
done < $dataset_names

cd $CURR_DIR

echo
echo "############"
echo "# FINISHED #"
echo "############"
echo
echo "goodbye."

