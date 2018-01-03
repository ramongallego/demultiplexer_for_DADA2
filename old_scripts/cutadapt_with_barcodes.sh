#Approach to split fastqs into fastqs by barcode

#Idea is to use the -p option in cutadapt to split them by anchored 5' ends

#Reads for sample 1 could start by its heavy or light strand, but the adapter
#will be the same

#Get the 1. and .2 subsetted
file.1.fastq="Lib-A_S1_L001_R1_001_sub.fastq"
file.2.fastq="Lib-A_S1_L001_R2_001_sub.fastq"
# Options used here explained
# -u 3 delete the fist 3 bases - a buffer in our case
# -g ^ We are after anchored 5' adapters


cutadapt -u 3 -q 10 -a ADAPTER_FWD --minimum-length 20 -o tmp.1.fastq -p tmp.2.fastq reads.1.fastq reads.2.fastq


cutadapt -a first=ACACAC -a second=ACAGCA -A ACGTACGT -A TGCATGCA -o trimmed-{name}.1.fastq.gz -p trimmed-{name}.2.fastq.gz input.1.fastq.gz input.2.fastq.gz

cutadapt -c 3 -g file:/Users/Moncho/trial_barcode_splitter/barcodes_F.fasta -o trimmed-{name}.1.fastq.gz -p trimmed-{name}.2.fastq.gz Lib-A_S1_L001_R1_001_sub.fastq Lib-A_S1_L001_R2_001_sub.fastq

No funciona : too many parameters

try first with only 2 adapters, on the 5

cutadapt -c 3 -g first=^ACACAC -g second=^ACAGCA -o trimmed-{name}.1.fastq.gz -p trimmed-{name}.2.fastq.gz Lib-A_S1_L001_R1_001_sub.fastq Lib-A_S1_L001_R2_001_sub.fastq

same error, try with file names the way cutadapt likes them

cp Lib-A_S1_L001_R1_001_sub.fastq Lib_A.1.fastq
cp Lib-A_S1_L001_R2_001_sub.fastq Lib_A.2.fastq

cutadapt -c 3 -g first=^ACACAC -g second=^ACAGCA -o trimmed-{name}.1.fastq.gz -p trimmed-{name}.2.fastq.gz Lib_A.1.fastq Lib_A.2.fastq

Same error - maybe is the -c3?
cutadapt -u 3 Lib-A_S1_L001_R1_001_sub.fastq -o Lib_A.1.fastq
cutadapt -u 3 Lib-A_S1_L001_R2_001_sub.fastq -o Lib_A.2.fastq
cutadapt -g first=^ACACAC -g second=^ACAGCA -o trimmed-{name}.1.fastq.gz -p trimmed-{name}.2.fastq.gz Lib_A.1.fastq Lib_A.2.fastq

Lo que funciona es :

cutadapt -g file:/Users/Moncho/trial_barcode_splitter/barcodes_F.fasta --no-trim  -o split_by_barcode/trimmed-{name}.1.fastq.gz -p split_by_barcode/trimmed-{name}.2.fastq.gz Lib_A.1.fastq Lib_A.2.fastq

FILE1=($(awk -F',' -v COLNUM=$COLNUM_FILE1 \
  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
  sort | uniq))

FILE2=($(awk -F',' -v COLNUM=$COLNUM_FILE2 \
  'NR>1 {print $COLNUM}' $SEQUENCING_METADATA |\
  sort | uniq ))

NFILE1="${#FILE1[@]}"
NFILE2="${#FILE2[@]}"
#make a subdirectory for all demultiplexed fastqs
DEMULT_DIR="${PARENT_DIR}"/demultiplexed
  mkdir "${DEMULT_DIR}"

for (( i=0; i < "${#FILE1[@]}"; i++ )); do
  CURRENT_FILE[0]="${FILE1[i]}"
	CURRENT_FILE[1]="${FILE2[i]}"
  for j in 0 1 ; do
    FILEPATH=$( find "${PARENT_DIR}" -name "${CURRENT_FILE[j]}"*  )
    if file "${FILEPATH}" | grep -q gzip ; then
      echo $(date +%Y-%m-%d\ %H:%M) "decompressing "${FILEPATH}""
      "${ZIPPER}" -d "${FILEPATH}"
      READ[j]=$( find "${PARENT_DIR}" -name "${CURRENT_FILE[j]}"*  )
    else
      READ[j]="${FILEPATH}"
    fi
  done
  READ1="${READ[0]}"
  READ2="${READ[1]}"
  OUTPUT_LIB="${DEMULT_DIR}"/Lib_"${i}"_Tag

  #run the cutadapt for DEMULTIPLEXING
cutadapt -j 4 -g file:"${FASTA_file}" --no-trim  -o "${OUTPUT_LIB}"-{name}.1.fastq -p "${OUTPUT_LIB}"-{name}.2.fastq.gz \
 "${READ1}" "${READ2}"
done
