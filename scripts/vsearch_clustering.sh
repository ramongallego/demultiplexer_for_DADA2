#!/usr/bin/env bash

## usage bash vsearch_clustering.sh <path_to_output_folder> <path_to_no_primers_folder> <USE_HASH> <LENR1> <LENR2>


# trim reads to desired length & and qc

OUTPUT_FOLDER=$1
NOPRIMERS_DIR=$2
HASH=$3
LENR1=$4
LENR2=$5

for file in "${NOPRIMERS_DIR}"/*R1.fastq; do
fwd_file=$(basename $file)
rev_file=$(echo $fwd_file | sed 's/R1.fastq$/R2.fastq/')

head "${NOPRIMERS_DIR}"/$fwd_file
head "${NOPRIMERS_DIR}"/$rev_file

vsearch --fastx_filter "${NOPRIMERS_DIR}"/$fwd_file --reverse "${NOPRIMERS_DIR}"/$rev_file --fastq_trunclen "${LENR1}" --fastq_maxns 1  --fastqout_rev "${OUTPUT_FOLDER}"/"${rev_file}" --output "${OUTPUT_FOLDER}"/"${fwd_file}" 

done
