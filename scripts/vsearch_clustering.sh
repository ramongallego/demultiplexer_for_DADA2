#!/usr/bin/env bash

## usage bash vsearch_clustering.sh <path_to_output_folder> <path_to_no_primers_folder> <USE_HASH> <LENR1> <LENR2>


# trim reads to desired length & and qc

OUTPUT_FOLDER=$1
NOPRIMERS_DIR=$2
HASH=$3
LENR1=$4
LENR2=$5

for file in "${NOPRIMERS_DIR}"/*R1.fastq; do

rev_file=$(echo $file | sed 's/R1.fastq$/R2.fastq/')

echo $file
echo $rev_file

done
