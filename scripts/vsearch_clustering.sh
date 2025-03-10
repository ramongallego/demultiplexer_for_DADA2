#!/usr/bin/env bash

## usage bash vsearch_clustering.sh <path_to_output_folder> <path_to_no_primers_folder> <USE_HASH> <LENR1> <LENR2>

MAIN_DIR="$(dirname "$0")"
SCRIPT_DIR="${MAIN_DIR}"/scripts
for file in "${SCRIPT_DIR}"/* ; do
	source "${file}"
done

# trim reads to desired length & and qc

OUTPUT_FOLDER=$1
NOPRIMERS_DIR=$2
HASH=$3
LENR1=$4
LENR2=$5

for file in "${NOPRIMERS_DIR}"/*R1.fastq; do

    R1_file=$(basename $file)
    R2_file=$(echo $R1_file | sed 's/R1.fastq$/R2.fastq/')

    merged_file=$(echo $R1_file | sed 's/R1.fastq$/merged.fasta/')
    unmerged_file=$(echo $R1_file | sed 's/R1.fastq$/un_merged.fasta/')
    

    ## TRIM to length with cutadapt, remove Ns


    cutadapt -j 0 \
        -u 0 -U 0  \
        -l "${LENR1}" -L "${LENR2}" --max-n 0 \
        -o "${OUTPUT_FOLDER}"/"${R1_file}" -p "${OUTPUT_FOLDER}"/"${R2_file}" \
        "${NOPRIMERS_DIR}"/$R1_file "${NOPRIMERS_DIR}"/$R2_file


    vsearch --fastq_mergepairs "${OUTPUT_FOLDER}"/$R1_file --reverse "${OUTPUT_FOLDER}"/$R2_file --fastaout "${OUTPUT_FOLDER}"/"${merged_file}" \
        --fastaout_notmerged_fwd "${OUTPUT_FOLDER}"/"${unmerged_file}"

done

for file in "${OUTPUT_FOLDER}"/*_Fwd.merged.fasta$; do

    # We should reverse the _Rev files and concatenate them after the Fwd ones, but do that after merging R1 and R2

    fwd_file=$(basename $file)
    rev_file=$(echo $fwd_file | sed 's/_Fwd.merged.fasta$/_Rev.merged.fasta/')

    derep_file=$(echo $fwd_file | sed 's/_Fwd.merged.fasta$/_derep.fasta/')

    centroids_file=$(echo $fwd_file | sed 's/_Fwd.merged.fasta$/_centroids.fasta/')

    revcom "${OUTPUT_FOLDER}"/$rev_file >> "${OUTPUT_FOLDER}"/$fwd_file

    vsearch --fastx_uniques "${OUTPUT_FOLDER}"/"${merged_file}" --sizeout --fastaout "${OUTPUT_FOLDER}"/"${derep_file}"

    vsearch --cluster_unoise "${OUTPUT_FOLDER}"/"${derep_file}" --centroids "${OUTPUT_FOLDER}"/"${centroids_file}" --sizein --sizeout --minsize 1



done
