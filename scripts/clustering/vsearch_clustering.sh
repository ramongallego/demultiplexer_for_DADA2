#!/usr/bin/env bash

## usage bash vsearch_clustering.sh <path_to_output_folder> <path_to_no_primers_folder> <USE_HASH> <LENR1> <LENR2>

MAIN_DIR="$(dirname "$0")"

# trim reads to desired length & and qc

OUTPUT_FOLDER=$1
NOPRIMERS_DIR=$2
HASH=$3
LENR1=$4
LENR2=$5

OUTPUT_SUMMARY="${OUTPUT_FOLDER}"/vsearch_summary.csv
echo "Sample, step, nreads" > "${OUTPUT_SUMMARY}"
MERGING_SUMMARY="${OUTPUT_FOLDER}"/merge_summary.csv
# echo "Sample,Pairs,Merged,Unmerged,kmers,multiple,lowscore,tooshort,staggered"  > "${MERGING_SUMMARY}"
echo "Sample,Pairs,Merged,Unmerged"  > "${MERGING_SUMMARY}"

MIDFILES="${OUTPUT_FOLDER}"/midfiles
mkdir "${MIDFILES}"

for file in "${NOPRIMERS_DIR}"/*R1.fastq; do

    R1_file=$(basename "$file")
    R2_file="${R1_file/R1.fastq/R2.fastq}"
    sample="${R1_file/.R1.fastq/}"

    merged_file="${sample}.merged.fastq"
    unmerged_file="${sample}.un_merged.fasta"
    cutadapt_log="${MIDFILES}/${sample}_cutadapt_log.txt"
    merge_log="${MIDFILES}/${sample}_merge.log" 

    echo "Processing sample: ${sample}"
    

    ## TRIM to length with cutadapt, remove Ns


    cutadapt -j 0 \
        -u 0 -U 0  \
        -l "${LENR1}" -L "${LENR2}" --max-n 0 \
        -o "${MIDFILES}"/"${R1_file}" -p "${MIDFILES}"/"${R2_file}" \
        "${NOPRIMERS_DIR}"/$R1_file "${NOPRIMERS_DIR}"/$R2_file  > "${MIDFILES}"/cutadapt_logqcontrol.txt

    num=$(grep "Pairs written (passing filters)" "${MIDFILES}"/cutadapt_logqcontrol.txt | awk '{print $5}' | tr -d ',')
    
    echo "${sample}, filtering, ${num}" >> "${OUTPUT_SUMMARY}"
    
    ## merge pairs and stats
    ## vsearch is terrible at merging super long overlaps, we are using flash2 instead

    flash2 "${MIDFILES}"/$R1_file "${MIDFILES}"/$R2_file -c > "${MIDFILES}"/"${merged_file}" 2> "${merge_log}"

#    vsearch --fastq_mergepairs "${MIDFILES}"/$R1_file --reverse "${MIDFILES}"/$R2_file --fastaout "${MIDFILES}"/"${merged_file}" \
#     --fastq_maxdiffs 2 \
#     --fastaout_notmerged_fwd "${MIDFILES}"/"${unmerged_file}" 2> "${merge_log}"

        pairs=$(grep "Total pairs:" "${merge_log}" | awk '{print $NF}')
        merged=$(grep "Combined pairs:" "${merge_log}" | awk '{print $NF}')
        unmerged=$(grep "Uncombined pairs:" "${merge_log}" | awk '{print $NF}')
            
    # #Extract key stats
    #     pairs=$(grep -oP '^\d+(?=\s+Pairs)' "$merge_log" || echo 0)
    #     merged=$(grep -oP '^\d+(?=\s+Merged)' "$merge_log" || echo 0)
    #     unmerged=$(grep -oP '^\d+(?=\s+Not merged)' "$merge_log" || echo 0)

    # # Extract failure reasons, setting defaults if missing
    #     kmers=$(grep -oP '^\d+(?=\s+too few kmers)' "$merge_log" || echo 0)
    #     multiple=$(grep -oP '^\d+(?=\s+multiple potential alignments)' "$merge_log" || echo 0)
    #     lowscore=$(grep -oP '^\d+(?=\s+alignment score too low)' "$merge_log" || echo 0)
    #     tooshort=$(grep -oP '^\d+(?=\s+overlap too short)' "$merge_log" || echo 0)
    #     staggered=$(grep -oP '^\d+(?=\s+staggered read pairs)' "$merge_log" || echo 0)

    echo "${sample}, merging, ${merged}" >> "${OUTPUT_SUMMARY}"  

    # Append to merge summary file
    # echo "${sample},${pairs},${merged},${unmerged},${kmers},${multiple},${lowscore},${tooshort},${staggered}" >> "${MERGING_SUMMARY}"
    echo "${sample},${pairs},${merged},${unmerged}" >> "${MERGING_SUMMARY}"
done

# for file in "${MIDFILES}"/*_Fwd.merged.fasta; do
for file in "${MIDFILES}"/*_Fwd.merged.fastq; do


    # We should reverse the _Rev files and concatenate them after the Fwd ones, but do that after merging R1 and R2

    fwd_file=$(basename "$file")
    # sample="${fwd_file/_Fwd.merged.fasta/}"
        sample="${fwd_file/_Fwd.merged.fastq/}"


    rev_file="${sample}_Rev.merged.fastq"
    derep_file="${sample}_derep.fasta"
    centroids_file="${sample}_centroids.fasta"
    non_chimeras_file="${sample}_non_chimeras.fasta"

    echo "Processing sample: ${sample} (Dereplication, Denoising, Chimera checking)"

    # reversing Rev reads and adding them at the end of Fwd file  
    seqkit seq -t dna -r -p -w 0 "${MIDFILES}"/"${rev_file}" >> "${MIDFILES}"/"${fwd_file}"
    # dereplicate
    vsearch --fastx_uniques "${MIDFILES}"/"${fwd_file}" --sizeout --fastaout "${MIDFILES}"/"${derep_file}"
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