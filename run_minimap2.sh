#!/usr/bin/bash

##############################################
# MAKE SURE YOU MAKE dataset_names.txt FIRST # 
##############################################
# automation of running minimap2. includes building indexes and running the program
# usage: ./run_minimap2.sh /home/kchan/thesis /home/kchan/thesis/dataset_names.txt
WORK_DIR=$1
dataset_names=$2
cd $WORK_DIR
THREADS=5

if [[ -z $1 || -z $2 ]]; then
	echo "argument error: pass in working directory and dataset names file"
	exit 1
fi

# build index for recA
# echo -e "\nbuilding indexes for recA...\n"
# minimap2 -d indexes/minimap2/recA_99_default.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# echo -e "\nadjusting k-mer sizes...\n"
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k11.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k13.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k17.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k19.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k21.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta

# echo -e "\nadjusting window sizes...\n"
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_default_w5.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k11_w3.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k13_w4.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k17_w5.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k19_w6.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta
# minimap2 -d $WORK_DIR/indexes/minimap2/recA_99_k21_w7.mmi references/fungene_9.5.1_recA_nucleotide_uclust99.fasta

echo "\nchecking for dataset names text file...\n"
if [[ -n "$dataset_names" ]]; then
	echo "received dataset_names arg"
else
	echo "argument error: pass in full path to dataset_names.txt"
	exit 1
fi

############
# MINIMAP2 #
############
cd $WORK_DIR
echo -e "\nrunning minimap2...\n"

######### TEST COMBO, DIFF PARAMS #########
# echo
# echo "#################"
# echo "# DEFAULT INDEX #"
# echo "#################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a --eqx indexes/minimap2/recA_99_default.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/default.sam
# 	minimap2 -t $THREADS -a --eqx -I 500000 indexes/minimap2/recA_99_default.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/default_I500k.sam	
# 	minimap2 -t $THREADS -a --eqx -I 100000 indexes/minimap2/recA_99_default.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/default_I100k.sam	

# 	minimap2 -t $THREADS -ax map-ont --eqx indexes/minimap2/recA_99_default.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/default_map_ont.sam
# 	minimap2 -t $THREADS -ax map-ont --eqx -I 500000 indexes/minimap2/recA_99_default.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/default_map_ont_I500k.sam	
# 	minimap2 -t $THREADS -ax map-ont --eqx -I 100000 indexes/minimap2/recA_99_default.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/default_map_ont_I100k.sam	

# done < $dataset_names

# echo
# echo "###############"
# echo "# -k 11 INDEX #"
# echo "###############"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 11 --eqx indexes/minimap2/recA_99_k11.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k11.sam
# 	minimap2 -t $THREADS -a -k 11 --eqx -I 500000 indexes/minimap2/recA_99_k11.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k11_I500k.sam	
# 	minimap2 -t $THREADS -a -k 11 --eqx -I 100000 indexes/minimap2/recA_99_k11.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k11_I100k.sam	
	
# 	minimap2 -t $THREADS -ax map-ont -k 11 --eqx indexes/minimap2/recA_99_k11.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k11_map_ont.sam
# 	minimap2 -t $THREADS -ax map-ont -k 11 --eqx -I 500000 indexes/minimap2/recA_99_k11.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k11_map_ont_I500k.sam	
# 	minimap2 -t $THREADS -ax map-ont -k 11 --eqx -I 100000 indexes/minimap2/recA_99_k11.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k11_map_ont_I100k.sam	

# done < $dataset_names

# echo
# echo "###############"
# echo "# -k 13 INDEX #"
# echo "###############"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 13 --eqx indexes/minimap2/recA_99_k13.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k13.sam
# 	minimap2 -t $THREADS -a -k 13 --eqx -I 500000 indexes/minimap2/recA_99_k13.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k13_I500k.sam	
# 	minimap2 -t $THREADS -a -k 13 --eqx -I 100000 indexes/minimap2/recA_99_k13.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k13_I100k.sam	

# 	minimap2 -t $THREADS -ax map-ont -k 13 --eqx indexes/minimap2/recA_99_k13.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k13_map_ont.sam
# 	minimap2 -t $THREADS -ax map-ont -k 13 --eqx -I 500000 indexes/minimap2/recA_99_k13.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k13_map_ont_I500k.sam	
# 	minimap2 -t $THREADS -ax map-ont -k 13 --eqx -I 100000 indexes/minimap2/recA_99_k13.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k13_map_ont_I100k.sam	

# done < $dataset_names

# echo
# echo "###############"
# echo "# -k 17 INDEX #"
# echo "###############"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 17 --eqx indexes/minimap2/recA_99_k17.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k17.sam
# 	minimap2 -t $THREADS -a -k 17 --eqx -I 500000 indexes/minimap2/recA_99_k17.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k17_I500k.sam	
# 	minimap2 -t $THREADS -a -k 17 --eqx -I 100000 indexes/minimap2/recA_99_k17.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k17_I100k.sam	

# 	minimap2 -t $THREADS -ax map-ont -k 17 --eqx indexes/minimap2/recA_99_k17.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k17_map_ont.sam
# 	minimap2 -t $THREADS -ax map-ont -k 17 --eqx -I 500000 indexes/minimap2/recA_99_k17.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k17_map_ont_I500k.sam	
# 	minimap2 -t $THREADS -ax map-ont -k 17 --eqx -I 100000 indexes/minimap2/recA_99_k17.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k17_map_ont_I100k.sam	

# done < $dataset_names

# echo
# echo "###############"
# echo "# -k 19 INDEX #"
# echo "###############"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 19 --eqx indexes/minimap2/recA_99_k19.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k19.sam
# 	minimap2 -t $THREADS -a -k 19 --eqx -I 500000 indexes/minimap2/recA_99_k19.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k19_I500k.sam	
# 	minimap2 -t $THREADS -a -k 19 --eqx -I 100000 indexes/minimap2/recA_99_k19.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k19_I100k.sam	

# 	minimap2 -t $THREADS -ax map-ont -k 19 --eqx indexes/minimap2/recA_99_k19.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k19_map_ont.sam
# 	minimap2 -t $THREADS -ax map-ont -k 19 --eqx -I 500000 indexes/minimap2/recA_99_k19.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k19_map_ont_I500k.sam	
# 	minimap2 -t $THREADS -ax map-ont -k 19 --eqx -I 100000 indexes/minimap2/recA_99_k19.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k19_map_ont_I100k.sam	

# done < $dataset_names

# echo
# echo "###############"
# echo "# -k 21 INDEX #"
# echo "###############"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 21 --eqx indexes/minimap2/recA_99_k21.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k21.sam
# 	minimap2 -t $THREADS -ax map-ont -k 21 --eqx indexes/minimap2/recA_99_k21.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k21_map_ont.sam
# 	minimap2 -t $THREADS -ax map-ont -k 21 --eqx -I 500000 indexes/minimap2/recA_99_k21.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k21_map_ont_I500k.sam	
# 	minimap2 -t $THREADS -ax map-ont -k 21 --eqx -I 100000 indexes/minimap2/recA_99_k21.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k21_map_ont_I100k.sam	

# done < $dataset_names

# echo
# echo "####################"
# echo "# -k 11 -w 3 INDEX #"
# echo "####################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 11 -w 3 --eqx indexes/minimap2/recA_99_k13_w3.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k11_w3.sam
# 	minimap2 -t $THREADS -a -k 11 -w 3 --eqx -I 500000 indexes/minimap2/recA_99_k13_w3.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k11_w3_I500k.sam	
# 	minimap2 -t $THREADS -a -k 11 -w 3 --eqx -I 100000 indexes/minimap2/recA_99_k13_w3.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k11_w3_I100k.sam	

# 	minimap2 -t $THREADS -ax map-ont -k 11 -w 3 --eqx indexes/minimap2/recA_99_k13_w3.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k11_w3_map_ont.sam
# 	minimap2 -t $THREADS -ax map-ont -k 11 -w 3 --eqx -I 500000 indexes/minimap2/recA_99_k13_w3.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k11_w3_map_ont_I500k.sam	
# 	minimap2 -t $THREADS -ax map-ont -k 11 -w 3 --eqx -I 100000 indexes/minimap2/recA_99_k13_w3.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k11_w3_map_ont_I100k.sam	

# done < $dataset_names

# echo
# echo "####################"
# echo "# -k 13 -w 4 INDEX #"
# echo "####################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 13 -w 4 --eqx indexes/minimap2/recA_99_k13_w4.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k13_w4.sam
# 	minimap2 -t $THREADS -a -k 13 -w 4 --eqx -I 500000 indexes/minimap2/recA_99_k13_w4.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k13_w4_I500k.sam	
# 	minimap2 -t $THREADS -a -k 13 -w 4 --eqx -I 100000 indexes/minimap2/recA_99_k13_w4.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k13_w4_I100k.sam	

# 	minimap2 -t $THREADS -ax map-ont -k 13 -w 4 --eqx indexes/minimap2/recA_99_k13_w4.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k13_w4_map_ont.sam
# 	minimap2 -t $THREADS -ax map-ont -k 13 -w 4 --eqx -I 500000 indexes/minimap2/recA_99_k13_w4.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k13_w4_map_ont_I500k.sam	
# 	minimap2 -t $THREADS -ax map-ont -k 13 -w 4 --eqx -I 100000 indexes/minimap2/recA_99_k13_w4.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k13_w4_map_ont_I100k.sam	

# done < $dataset_names

# echo
# echo "####################"
# echo "# -k 17 -w 5 INDEX #"
# echo "####################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 17 -w 5 --eqx indexes/minimap2/recA_99_k17_w5.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k17_w5.sam
# 	minimap2 -t $THREADS -a -k 17 -w 5 --eqx -I 500000 indexes/minimap2/recA_99_k17_w5.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k17_w5_I500k.sam	
# 	minimap2 -t $THREADS -a -k 17 -w 5 --eqx -I 100000 indexes/minimap2/recA_99_k17_w5.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k17_w5_I100k.sam	

# 	minimap2 -t $THREADS -ax map-ont -k 17 -w 5 --eqx indexes/minimap2/recA_99_k17_w5.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k17_w5_map_ont.sam
# 	minimap2 -t $THREADS -ax map-ont -k 17 -w 5 --eqx -I 500000 indexes/minimap2/recA_99_k17_w5.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k17_w5_map_ont_I500k.sam	
# 	minimap2 -t $THREADS -ax map-ont -k 17 -w 5 --eqx -I 100000 indexes/minimap2/recA_99_k17_w5.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k17_w5_map_ont_I100k.sam	

# done < $dataset_names

# echo
# echo "####################"
# echo "# -k 19 -w 6 INDEX #"
# echo "####################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 19 -w 6 --eqx indexes/minimap2/recA_99_k19_w6.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k19_w6.sam
# 	minimap2 -t $THREADS -a -k 19 -w 6 --eqx -I 500000 indexes/minimap2/recA_99_k19_w6.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k19_w6_I500k.sam	
# 	minimap2 -t $THREADS -a -k 19 -w 6 --eqx -I 100000 indexes/minimap2/recA_99_k19_w6.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k19_w6_I100k.sam	

# 	minimap2 -t $THREADS -ax map-ont -k 19 -w 6 --eqx indexes/minimap2/recA_99_k19_w6.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k19_w6_map_ont.sam
# 	minimap2 -t $THREADS -ax map-ont -k 19 -w 6 --eqx -I 500000 indexes/minimap2/recA_99_k19_w6.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k19_w6_map_ont_I500k.sam	
# 	minimap2 -t $THREADS -ax map-ont -k 19 -w 6 --eqx -I 100000 indexes/minimap2/recA_99_k19_w6.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k19_w6_map_ont_I100k.sam	

# done < $dataset_names

# echo
# echo "####################"
# echo "# -k 21 -w 7 INDEX #"
# echo "####################"
# echo

# while read name; do
# 	echo "FOR DATASET $name"
# 	mkdir -p processed/minimap2/$name

# 	minimap2 -t $THREADS -a -k 21 -w 7 --eqx indexes/minimap2/recA_99_k21_w7.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k21_w7.sam
# 	minimap2 -t $THREADS -a -k 21 -w 7 --eqx -I 500000 indexes/minimap2/recA_99_k21_w7.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k21_w7_I500k.sam	
# 	minimap2 -t $THREADS -a -k 21 -w 7 --eqx -I 100000 indexes/minimap2/recA_99_k21_w7.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k21_w7_I100k.sam	

# 	minimap2 -t $THREADS -ax map-ont -k 21 -w 7 --eqx indexes/minimap2/recA_99_k21_w7.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k21_w7_map_ont.sam
# 	minimap2 -t $THREADS -ax map-ont -k 21 -w 7 --eqx -I 500000 indexes/minimap2/recA_99_k21_w7.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k21_w7_map_ont_I500k.sam	
# 	minimap2 -t $THREADS -ax map-ont -k 21 -w 7 --eqx -I 100000 indexes/minimap2/recA_99_k21_w7.mmi raw_data/$name.fastq.gz > processed/minimap2/$name/k21_w7_map_ont_I100k.sam	

# done < $dataset_names

#######################
# SAMTOOLS FORMATTING #
#######################
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

cd $WORK_DIR
echo -e "\nbam conversion and adding MD tags...\n"
while read name; do
	echo "FOR DATASET $name"
	mkdir -p $WORK_DIR/processed/minimap2/$name/sorted_bams
	cd $WORK_DIR/processed/minimap2/$name
	for f in *.sam; do
		filename=${f//.sam}
		out=${f/.sam/_with_md.bam}
		samtools view -h -F 4 $filename.sam > tmp.sam && mv tmp.sam $filename.sam
		samtools calmd -S $filename.sam $WORK_DIR/references/fungene_9.5.1_recA_nucleotide_uclust99.fasta | samtools view -bS - > $out
		mv $out $WORK_DIR/processed/minimap2/$name/sorted_bams
	done
done < $dataset_names

cd $WORK_DIR
echo -e "\nsorting and indexing...\n"
while read name; do
	echo "FOR DATASET $name"
	cd $WORK_DIR/processed/minimap2/$name/sorted_bams
	for f in *_with_md.bam; do
		filename=${f/.bam/_sorted.bam}
		samtools sort --threads $THREADS -o $filename $f && samtools index $filename && rm $f
	done
done < $dataset_names

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

cd $WORK_DIR
echo -e "\nfinished!\n"






















