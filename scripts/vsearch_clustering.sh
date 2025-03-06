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

merged_file=$(echo $fwd_file | sed 's/R1.fastq$/merged.fastq/')

## TRIM to length with cutadapt, remove Ns


cutadapt -j 0 \
 -u 0 -U 0  \
 -l "${LENR1}" -L "${LENR2}" --max-n 0 \
 -o "${OUTPUT_FOLDER}"/"${fwd_file}" -p "${OUTPUT_FOLDER}"/"${rev_file}" \
 "${NOPRIMERS_DIR}"/$fwd_file "${NOPRIMERS_DIR}"/$rev_file


vsearch --fastq_mergepairs "${OUTPUT_FOLDER}"/$fwd_file --reverse "${OUTPUT_FOLDER}"/$rev_file --fastqout "${OUTPUT_FOLDER}"/"${merged_file}" 

done
