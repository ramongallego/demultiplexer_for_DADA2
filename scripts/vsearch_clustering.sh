#!/usr/bin/env bash

## usage bash vsearch_clustering.sh <path_to_output_folder> <path_to_no_primers_folder> <USE_HASH> <LENR1> <LENR2>

MAIN_DIR="$(dirname "$0")"

source "${MAIN_DIR}"/revcom.sh


# trim reads to desired length & and qc

OUTPUT_FOLDER=$1
NOPRIMERS_DIR=$2
HASH=$3
LENR1=$4
LENR2=$5

OUTPUT_SUMMARY="${OUTPUT_FOLDER}"/vsearch_summary.csv
echo "Sample, step, nreads" > "${OUTPUT_SUMMARY}"
MIDFILES="${OUTPUT_FOLDER}"/midfiles
mkdir "${MIDFILES}"

for file in "${NOPRIMERS_DIR}"/*R1.fastq; do

    R1_file=$(basename $file)
    R2_file=$(echo $R1_file | sed 's/R1.fastq$/R2.fastq/')

    merged_file=$(echo $R1_file | sed 's/R1.fastq$/merged.fasta/')
    unmerged_file=$(echo $R1_file | sed 's/R1.fastq$/un_merged.fasta/')
    sample=$(echo $R1_file | sed 's/.R1.fastq$//')
    

    ## TRIM to length with cutadapt, remove Ns


    cutadapt -j 0 \
        -u 0 -U 0  \
        -l "${LENR1}" -L "${LENR2}" --max-n 0 \
        -o "${MIDFILES}"/"${R1_file}" -p "${MIDFILES}"/"${R2_file}" \
        "${NOPRIMERS_DIR}"/$R1_file "${NOPRIMERS_DIR}"/$R2_file  > "${MIDFILES}"/cutadapt_logqcontrol.txt

    num=$(grep "Pairs written (passing filters)" "${MIDFILES}"/cutadapt_logqcontrol.txt | awk '{print $5}' | tr -d ',')
    
    echo "${sample}, filtering, ${num}" >> "${OUTPUT_SUMMARY}"


    num=$(vsearch --fastq_mergepairs "${MIDFILES}"/$R1_file --reverse "${MIDFILES}"/$R2_file --fastaout "${MIDFILES}"/"${merged_file}" \
        --fastaout_notmerged_fwd "${MIDFILES}"/"${unmerged_file}" 2>&1 | grep "Merged ("  | awk '{print $1}')
        # )

     echo "${sample}, merging, ${num}" >> "${OUTPUT_SUMMARY}"   
    
done

for file in "${MIDFILES}"/*_Fwd.merged.fasta; do

    # We should reverse the _Rev files and concatenate them after the Fwd ones, but do that after merging R1 and R2

    fwd_file=$(basename $file)
    rev_file=$(echo $fwd_file | sed 's/_Fwd.merged.fasta$/_Rev.merged.fasta/')

    derep_file=$(echo $fwd_file | sed 's/_Fwd.merged.fasta$/_derep.fasta/')

    centroids_file=$(echo $fwd_file | sed 's/_Fwd.merged.fasta$/_centroids.fasta/')

    non_chimeras_file=$(echo $fwd_file | sed 's/_Fwd.merged.fasta$/_non_chimeras.fasta/')

    sample=$(echo $fwd_file | sed 's/_Fwd.merged.fasta$//')

    # reversing Rev reads and adding them at the end of Fwd file  
    revcom "${MIDFILES}"/$rev_file >> "${MIDFILES}"/$fwd_file
    # dereplicate
    vsearch --fastx_uniques "${MIDFILES}"/"${merged_file}" --sizeout --fastaout "${MIDFILES}"/"${derep_file}"
    # denoise 
    vsearch --cluster_unoise "${MIDFILES}"/"${derep_file}"  --sizein --sizeout --minsize 1 --centroids "${MIDFILES}"/"${centroids_file}"

        # calculate number of reads after denoising 

        denoised_reads=$(grep -oP '(?<=size=)[0-9]+' "${MIDFILES}"/"${centroids_file}" | awk '{sum+=$1} END {print sum}')

        echo "${sample}, denoising, ${denoised_reads}" >> "${OUTPUT_SUMMARY}"
    
    # chimera checking
    vsearch --uchime3_denovo "${MIDFILES}"/"${centroids_file}"  --sizein --sizeout --nonchimeras - | seqkit seq -w 0 > "${MIDFILES}"/"${non_chimeras_file}" 
    
        # calculate number of reads after chimera checking
        nonchim_reads=$(grep -oP '(?<=size=)[0-9]+' "${MIDFILES}"/"${non_chimeras_file}" | awk '{sum+=$1} END {print sum}')
        echo "${sample}, chimeras, ${nonchim_reads}" >> "${OUTPUT_SUMMARY}"

done

## launch Parsing rscript

Rscript "${MAIN_DIR}"/Parse_Abundances.R "${OUTPUT_FOLDER}" "${HASH}"