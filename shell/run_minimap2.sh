#!/usr/bin/bash

##############################################
# MAKE SURE YOU MAKE dataset_names.txt FIRST # 
##############################################
# automation of running minimap2. includes building indexes and running the program
# usage: time /home/kchan/scripts_thesis/shell/run_minimap2.sh /home/kchan/thesis /home/kchan/thesis/dataset_names.txt /home/kchan/thesis/references/rpoB/fungene_9.6_rpoB_nucleotide_uclust99.fasta rpoB
WORK_DIR=$1
dataset_names=$2
reference=$3
marker_gene=$4
THREADS=5
CURR_DIR=$PWD

MINIMAP2=$(which minimap2)

cd $WORK_DIR
if [[ -z $1 || -z $2 || -z $3 || -z $4 ]]; then
	echo "argument error: pass in working directory, dataset names file, full path to reference, and marker gene name"
	exit 1
else
	echo "starting script"
	echo "working directory: $WORK_DIR"
	echo "dataset names: $dataset_names"
	echo "reference file path: $reference"
	echo "marker gene name: $marker_gene"
fi

echo
echo "building index..."
echo

$MINIMAP2 -d $WORK_DIR/indexes/minimap2/$marker_gene.mmi $reference

############
# MINIMAP2 #
############
cd $WORK_DIR

echo
echo "aligning reads with..."
echo

######### TEST COMBO, DIFF PARAMS #########
echo
echo "#################"
echo "# DEFAULT INDEX #"
echo "#################"
echo

while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/$marker_gene/minimap2/$name

	minimap2 -t $THREADS -ax map-ont --eqx $WORK_DIR/indexes/minimap2/$marker_gene.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/$marker_gene/minimap2/$name/default_map_ont.sam

done < $dataset_names

cd $CURR_DIR

echo
echo "############"
echo "# FINISHED #"
echo "############"
echo
echo "goodbye."





















