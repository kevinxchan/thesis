#!/usr/bin/bash

##############################################
# MAKE SURE YOU MAKE dataset_names.txt FIRST # 
##############################################
# automation of running minimap2. includes building indexes and running the program
# usage: time /home/kchan/scripts_thesis/shell/run_minimap2.sh /home/kchan/thesis /home/kchan/thesis/dataset_names.txt /home/kchan/thesis/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
WORK_DIR=$1
dataset_names=$2
reference=$3
THREADS=5
CURR_DIR=$PWD

cd $WORK_DIR
if [[ -z $1 || -z $2 || -z $3 ]]; then
	echo "argument error: pass in working directory, dataset names file, and full path to reference"
	exit 1
fi

# build index for recA
# echo
# echo "building indexes for recA..."
# echo

# minimap2 -d indexes/minimap2/recA_99_default.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta

# echo
# echo "adjusting k-mer sizes..."
# echo

# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k11.mmi -k 11 references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k13.mmi -k 13 references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k17.mmi -k 17 references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k19.mmi -k 19 references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k21.mmi -k 21 references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k28.mmi -k 28 $WORK_DIR/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_H_k19.mmi -Hk19 $WORK_DIR/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta

# echo
# echo "adjusting window sizes..."
# echo 

# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_default_w5.mmi -w 5 references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_default_w30.mmi -w 30 $WORK_DIR/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k11_w3.mmi -k 11 -w 3 references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k13_w4.mmi -k 13 -w 4 references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k17_w5.mmi -k 17 -w 5 references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k19_w6.mmi -k 19 -w 6 references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k21_w7.mmi -k 21 -w 7 references/fungene_9.5.1_recA_nucleotide_uclust99.fasta

############
# MINIMAP2 #
############
cd $WORK_DIR

echo
echo "aligning reads with..."
echo

######### TEST COMBO, DIFF PARAMS #########
# echo
# echo "#################"
# echo "# DEFAULT INDEX #"
# echo "#################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p $WORK_DIR/processed/minimap2/$name

# 	# minimap2 -t $THREADS -a --eqx --MD $WORK_DIR/indexes/minimap2/recA_99_default.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/default.sam
# 	minimap2 -t $THREADS -ax map-ont --eqx $WORK_DIR/indexes/minimap2/recA_99_default.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/default_map_ont.sam

# done < $dataset_names

# echo
# echo "###########################"
# echo "# -H -k 19 (map-pb) INDEX #"
# echo "###########################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -ax map-pb --eqx $WORK_DIR/indexes/minimap2/recA_99_H_k19.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/default_map_pb.sam

# done < $dataset_names

# echo
# echo "###############"
# echo "# -k 11 INDEX #"
# echo "###############"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p $WORK_DIR/processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 11 --eqx $WORK_DIR/indexes/minimap2/recA_99_k11.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k11.sam
# 	minimap2 -t $THREADS -a -k 11 --eqx -I 500000 $WORK_DIR/indexes/minimap2/recA_99_k11.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k11_I500k.sam	
# 	minimap2 -t $THREADS -a -k 11 --eqx -I 100000 $WORK_DIR/indexes/minimap2/recA_99_k11.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k11_I100k.sam	

# done < $dataset_names

# echo
# echo "###############"
# echo "# -k 13 INDEX #"
# echo "###############"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p $WORK_DIR/processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 13 --eqx $WORK_DIR/indexes/minimap2/recA_99_k13.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k13.sam
# 	minimap2 -t $THREADS -a -k 13 --eqx -I 500000 $WORK_DIR/indexes/minimap2/recA_99_k13.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k13_I500k.sam	
# 	minimap2 -t $THREADS -a -k 13 --eqx -I 100000 $WORK_DIR/indexes/minimap2/recA_99_k13.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k13_I100k.sam	

# done < $dataset_names

# echo
# echo "###############"
# echo "# -k 17 INDEX #"
# echo "###############"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p $WORK_DIR/processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 17 --eqx $WORK_DIR/indexes/minimap2/recA_99_k17.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k17.sam
# 	minimap2 -t $THREADS -a -k 17 --eqx -I 500000 $WORK_DIR/indexes/minimap2/recA_99_k17.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k17_I500k.sam	
# 	minimap2 -t $THREADS -a -k 17 --eqx -I 100000 $WORK_DIR/indexes/minimap2/recA_99_k17.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k17_I100k.sam	

# done < $dataset_names

# echo
# echo "###############"
# echo "# -k 19 INDEX #"
# echo "###############"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p $WORK_DIR/processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 19 --eqx $WORK_DIR/indexes/minimap2/recA_99_k19.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k19.sam
# 	minimap2 -t $THREADS -a -k 19 --eqx -I 500000 $WORK_DIR/indexes/minimap2/recA_99_k19.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k19_I500k.sam	
# 	minimap2 -t $THREADS -a -k 19 --eqx -I 100000 $WORK_DIR/indexes/minimap2/recA_99_k19.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k19_I100k.sam	

# done < $dataset_names

# echo
# echo "###############"
# echo "# -k 21 INDEX #"
# echo "###############"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 21 --eqx $WORK_DIR/indexes/minimap2/recA_99_k21.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k21.sam
# 	minimap2 -t $THREADS -a -k 21 --eqx $WORK_DIR/indexes/minimap2/recA_99_k21.mmi -I 500000 $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k21_I500k.sam
# 	minimap2 -t $THREADS -a -k 21 --eqx $WORK_DIR/indexes/minimap2/recA_99_k21.mmi -I 100000 $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k21_I500k.sam

# done < $dataset_names

# echo
# echo "###############"
# echo "# -k 28 INDEX #"
# echo "###############"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 28 --eqx $WORK_DIR/indexes/minimap2/recA_99_k28.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k28.sam

# done < $dataset_names

echo
echo "###############"
echo "# -w 30 INDEX #"
echo "###############"
echo

while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/minimap2/$name

	minimap2 -t $THREADS -a -w 30 --eqx $WORK_DIR/indexes/minimap2/recA_99_default_w30.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/w30.sam

done < $dataset_names

# echo
# echo "####################"
# echo "# -k 11 -w 3 INDEX #"
# echo "####################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p $WORK_DIR/processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 11 -w 3 --eqx $WORK_DIR/indexes/minimap2/recA_99_k13_w3.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k11_w3.sam
# 	minimap2 -t $THREADS -a -k 11 -w 3 --eqx -I 500000 $WORK_DIR/indexes/minimap2/recA_99_k13_w3.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k11_w3_I500k.sam	
# 	minimap2 -t $THREADS -a -k 11 -w 3 --eqx -I 100000 $WORK_DIR/indexes/minimap2/recA_99_k13_w3.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k11_w3_I100k.sam	

# done < $dataset_names

# echo
# echo "####################"
# echo "# -k 13 -w 4 INDEX #"
# echo "####################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p $WORK_DIR/processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 13 -w 4 --eqx $WORK_DIR/indexes/minimap2/recA_99_k13_w4.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k13_w4.sam
# 	minimap2 -t $THREADS -a -k 13 -w 4 --eqx -I 500000 $WORK_DIR/indexes/minimap2/recA_99_k13_w4.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k13_w4_I500k.sam	
# 	minimap2 -t $THREADS -a -k 13 -w 4 --eqx -I 100000 $WORK_DIR/indexes/minimap2/recA_99_k13_w4.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k13_w4_I100k.sam	

# done < $dataset_names

# echo
# echo "####################"
# echo "# -k 17 -w 5 INDEX #"
# echo "####################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 17 -w 5 --eqx $WORK_DIR/indexes/minimap2/recA_99_k17_w5.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k17_w5.sam
# 	minimap2 -t $THREADS -a -k 17 -w 5 --eqx -I 500000 $WORK_DIR/indexes/minimap2/recA_99_k17_w5.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k17_w5_I500k.sam	
# 	minimap2 -t $THREADS -a -k 17 -w 5 --eqx -I 100000 $WORK_DIR/indexes/minimap2/recA_99_k17_w5.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k17_w5_I100k.sam	

# done < $dataset_names

# echo
# echo "####################"
# echo "# -k 19 -w 6 INDEX #"
# echo "####################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p $WORK_DIR/processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 19 -w 6 --eqx $WORK_DIR/indexes/minimap2/recA_99_k19_w6.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k19_w6.sam
# 	minimap2 -t $THREADS -a -k 19 -w 6 --eqx -I 500000 $WORK_DIR/indexes/minimap2/recA_99_k19_w6.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k19_w6_I500k.sam	
# 	minimap2 -t $THREADS -a -k 19 -w 6 --eqx -I 100000 $WORK_DIR/indexes/minimap2/recA_99_k19_w6.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k19_w6_I100k.sam	

# done < $dataset_names

# echo
# echo "####################"
# echo "# -k 21 -w 7 INDEX #"
# echo "####################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p $WORK_DIR/processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 21 -w 7 --eqx $WORK_DIR/indexes/minimap2/recA_99_k21_w7.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k21_w7.sam
# 	minimap2 -t $THREADS -a -k 21 -w 7 --eqx -I 500000 $WORK_DIR/indexes/minimap2/recA_99_k21_w7.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k21_w7_I500k.sam	
# 	minimap2 -t $THREADS -a -k 21 -w 7 --eqx -I 100000 $WORK_DIR/indexes/minimap2/recA_99_k21_w7.mmi $WORK_DIR/raw_data/$name.fastq.gz > $WORK_DIR/processed/minimap2/$name/k21_w7_I100k.sam	

# done < $dataset_names

#######################
# SAMTOOLS FORMATTING #
#######################
# cd $WORK_DIR

# echo
# echo "sorting sam files by read names..."
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p $WORK_DIR/processed/minimap2/$name/
# 	cd $WORK_DIR/processed/minimap2/$name

# 	for f in *.sam; do
# 		filename=${f//.sam}
# 		samtools view -b -h $filename.sam | samtools sort -n - | samtools view -h > tmp.sam && mv tmp.sam $filename.sam
# 	done

# done < $dataset_names

# cd $WORK_DIR
# echo -e "\ngetting stats for each alignment...\n"
# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p $WORK_DIR/processed/minimap2/$name/flagstats
# 	cd $WORK_DIR/processed/minimap2/$name

# 	for f in *.sam; do
# 		filename=${f//.sam}
# 		samtools flagstat $WORK_DIR/processed/minimap2/$name/$filename.sam > $WORK_DIR/processed/minimap2/$name/flagstats/$filename.txt
# 	done

# done < $dataset_names

# cd $WORK_DIR
# echo -e "\nbam conversion and adding MD tags...\n"
# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p $WORK_DIR/processed/minimap2/$name/sorted_bams
# 	cd $WORK_DIR/processed/minimap2/$name
# 	for f in *.sam; do
# 		filename=${f//.sam}
# 		out=${f/.sam/_with_md.bam}
# 		samtools view -h -F 4 $filename.sam > tmp.sam && mv tmp.sam $filename.sam
# 		samtools calmd -S $filename.sam $WORK_DIR/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta | samtools view -bS - > $out
# 		mv $out $WORK_DIR/processed/minimap2/$name/sorted_bams
# 	done
# done < $dataset_names

# cd $WORK_DIR
# echo -e "\nsorting and indexing...\n"
# while read name; do
# 	echo "FOR DATASET $name"
# 	cd $WORK_DIR/processed/minimap2/$name/sorted_bams
# 	for f in *_with_md.bam; do
# 		filename=${f/.bam/_sorted.bam}
# 		samtools sort --threads $THREADS -o $filename $f && samtools index $filename && rm $f
# 	done
# done < $dataset_names

# cd $WORK_DIR
# echo -e "\ngenerating alignments with sam2pairwise...\n"
# while read name; do
# 	echo "FOR DATASET $name"
# 	cd $WORK_DIR/processed/minimap2/$name
# 	for f in *_sorted.bam; do
# 		filename=${f/.bam/_alns.out}
# 		samtools view $f | sam2pairwise > $filename
# 	done
# done < $dataset_names

cd $CURR_DIR

echo
echo "############"
echo "# FINISHED #"
echo "############"
echo
echo "goodbye."





















