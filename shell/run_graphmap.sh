#!/usr/bin/bash

# automation of running graphmap.
# args:
# 1. absolute path to working directory
# 2. dataset names file (single dataset name per line, stored in .txt file)
# 3. full path to the reference
# USAGE: 
# /home/kchan/scripts_thesis/shell/run_graphmap.sh /home/kchan/thesis /home/kchan/thesis/dataset_names.txt /home/kchan/thesis/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
WORK_DIR=$1
dataset_names=$2
reference=$3
THREADS=8
CURR_DIR=$PWD

GRAPHMAP=$(which graphmap)
SAMTOOLS=$(which samtools)

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

index_file=$index_path/"$(basename $reference.gmidx)"

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
	mkdir -p $WORK_DIR/processed/graphmap/$name

	$GRAPHMAP align --threads $THREADS -r $reference --index $index_file --extcigar \
	--reads $WORK_DIR/raw_data/$name.fastq.gz \
	-o $WORK_DIR/processed/graphmap/$name/default.sam
done < $dataset_names

echo
echo "############################"
echo "# DEFAULTS; E-VALUE < 1e-5 #"
echo "############################"
echo

while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/graphmap/$name

	$GRAPHMAP align --evalue 1e-5 --threads $THREADS -r $reference --index $index_file --extcigar \
	--reads $WORK_DIR/raw_data/$name.fastq.gz \
	-o $WORK_DIR/processed/graphmap/$name/default_eval1e-5.sam
done < $dataset_names

echo
echo "##########################"
echo "# SG ALG; E-VALUE < 1e-5 #"
echo "##########################"
echo

while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/graphmap/$name

	$GRAPHMAP align --alg sg --evalue 1e-5 --threads $THREADS -r $reference --index $index_file --extcigar \
	--reads $WORK_DIR/raw_data/$name.fastq.gz \
	-o $WORK_DIR/processed/graphmap/$name/sg_eval1e-5.sam
done < $dataset_names

echo
echo "##########################"
echo "# SGGOTOH ALG; E-VALUE < 1e-5 #"
echo "##########################"
echo

while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/graphmap/$name

	$GRAPHMAP align --alg sg --evalue 1e-5 --threads $THREADS -r $reference --index $index_file --extcigar \
	--reads $WORK_DIR/raw_data/$name.fastq.gz \
	-o $WORK_DIR/processed/graphmap/$name/sggotoh_eval1e-5.sam
done < $dataset_names

echo
echo "###################################"
echo "# ANCHORGOTOH ALG; E-VALUE < 1e-5 #"
echo "###################################"
echo

while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/graphmap/$name

	$GRAPHMAP align --alg sg --evalue 1e-5 --threads $THREADS -r $reference --index $index_file --extcigar \
	--reads $WORK_DIR/raw_data/$name.fastq.gz \
	-o $WORK_DIR/processed/graphmap/$name/anchorgotoh_eval1e-5.sam
done < $dataset_names

cd $WORK_DIR
echo
echo "getting flagstats for each sample..."
echo
while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/graphmap/$name/flagstats
	cd $WORK_DIR/processed/graphmap/$name

	for f in *.sam; do
		filename=${f//.sam}
		samtools flagstat $WORK_DIR/processed/graphmap/$name/$filename.sam > $WORK_DIR/processed/graphmap/$name/flagstats/$filename.txt
	done

done < $dataset_names


echo
echo "filtering sam files, generating bams..."
echo

cd $WORK_DIR
while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/graphmap/$name/sorted_bams
	cd $WORK_DIR/processed/graphmap/$name
	for f in *.sam; do
		filename=${f//.sam}
		out=${f/.sam/_with_md.bam}
		samtools view -h -F 4 $filename.sam > tmp.sam && mv tmp.sam $filename.sam
		samtools calmd -S "$filename".sam $reference | samtools view -bS - > $out
		mv $out $WORK_DIR/processed/graphmap/$name/sorted_bams
	done
done < $dataset_names

echo
echo "sorting and indexing bams..."
echo

cd $WORK_DIR
while read name; do
	echo "FOR DATASET $name"
	cd $WORK_DIR/processed/graphmap/$name/sorted_bams
	for f in *_with_md.bam; do
		filename=${f/.bam/_sorted.bam}
		samtools sort --threads $THREADS -o $filename $f && samtools index $filename && rm $f
	done
done < $dataset_names

cd $CURR_DIR

echo
echo "############"
echo "# FINISHED #"
echo "############"
echo
echo "goodbye."



