#!/usr/bin/env bash

# Usage bash demultiplex_both_fastqs.sh banzai_params.sh
#This script is built using banzai (github.com/jimmyodonnell/banzai) as template

#We need to gather: Location of functions  and fastqs:
MAIN_DIR="$(dirname "$0")"
SCRIPT_DIR="${MAIN_DIR}"/scripts
for file in "${SCRIPT_DIR}"/* ; do
	source "${file}"
done

param_file=${1}

echo "Reading analysis parameters from:"
echo "${param_file}"
source "${param_file}"

#Now gather the info we need:
if [[ $( file "${SEQUENCING_METADATA}" ) == *"CRLF"* ]]; then

  echo "The file has CRLF endings. Let me fix that for you..."

  BASE="${SEQUENCING_METADATA%.*}"

  EXT="${SEQUENCING_METADATA##*.}"

  NEWLINES_FIXED="${BASE}"_fix."${EXT}"

  tr -d '\r' < "${SEQUENCING_METADATA}" > "${NEWLINES_FIXED}"

  echo "the old file was: ${SEQUENCING_METADATA}"

  echo "The new file is here:"

  echo "${NEWLINES_FIXED}"

else

  echo "The file passes test for CRLF. Everybody dance!"
  echo

fi

if [[ -s "${NEWLINES_FIXED}" ]]; then
	SEQUENCING_METADATA="${NEWLINES_FIXED}"
fi

#Create output directory
START_TIME=$(date +%Y%m%d_%H%M)
DEMULT_DIR="${OUTPUT_DIRECTORY}"/demultiplexed_"${START_TIME}"

if [[ -d "${OUTPUT_DIRECTORY}" ]]; then
  mkdir "${DEMULT_DIR}"
  echo "Output files would be in ${DEMULT_DIR}"
else
  mkdir "${OUTPUT_DIRECTORY}"
  mkdir "${DEMULT_DIR}"
  echo "Output files would be in ${DEMULT_DIR}"
fi
mkdir "${DEMULT_DIR}"/cleaned
################################################################################
# READ METADATA
################################################################################
# report metadata dimensions
METADATA_DIM=($( awk -F, 'END{print NR, NF}' "${SEQUENCING_METADATA}" ))
echo "Metadata has" "${METADATA_DIM[0]}" "rows and" "${METADATA_DIM[1]}" "columns including header."
N_SAMPLES=$( echo "${METADATA_DIM[0]}" - 1 | bc )
echo "Expecting" "${N_SAMPLES}" "samples total."
echo
##NOW WE HAVE LOADED THE SEQUENCING_METADATA - WE NEED IT
# Filnames
COLNUM_FILE1=$( get_colnum "${COLNAME_FILE1}" "${SEQUENCING_METADATA}")
COLNUM_FILE2=$( get_colnum "${COLNAME_FILE2}" "${SEQUENCING_METADATA}")

# Library names
COLNUM_ID1=$( get_colnum "${COLNAME_ID1_NAME}" "${SEQUENCING_METADATA}")
COLNUM_ID1_SEQ=$( get_colnum "${COLNAME_ID1_SEQ}" "${SEQUENCING_METADATA}")

# Secondary indices
COLNUM_ID2=$( get_colnum "${COLNAME_ID2_SEQ}" "${SEQUENCING_METADATA}")

# Secondary index sequence positions
COLNUM_ID2_START=$( get_colnum "${COLNAME_ID2_START}" "${SEQUENCING_METADATA}")

# Sample names
COLNUM_SAMPLE=$( get_colnum "${COLNAME_SAMPLE_ID}" "${SEQUENCING_METADATA}")

# Primers
COLNUM_PRIMER1=$( get_colnum "${COLNAME_PRIMER1}" "${SEQUENCING_METADATA}")
COLNUM_PRIMER2=$( get_colnum "${COLNAME_PRIMER2}" "${SEQUENCING_METADATA}")

#see if we make it work with both fastqs
################################################################################
# CHECK FILES
################################################################################
FILE1=($(awk -F',' -v COLNUM=$COLNUM_FILE1 \
  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
  sort | uniq))

FILE2=($(awk -F',' -v COLNUM=$COLNUM_FILE2 \
  'NR>1 {print $COLNUM}' $SEQUENCING_METADATA |\
  sort | uniq ))

NFILE1="${#FILE1[@]}"
NFILE2="${#FILE2[@]}"
if [ "${NFILE1}" != "${NFILE2}" ]; then
	echo "ERROR: Whoa! different number of forward and reverse files"
fi

if [[ -n "${FILE1}" && -n "${FILE2}" ]]; then
  echo 'Files read from metadata columns' "${COLNUM_FILE1}" 'and' "${COLNUM_FILE2}"
  echo 'File names:'
	for (( i=0; i < "${NFILE1}"; ++i)); do
		printf '%s\t%s\n' "${FILE1[i]}" "${FILE2[i]}"
	done
	echo
else
  echo 'ERROR:' 'At least one file is not valid'
  echo 'Looked in metadata columns' "${COLNUM_FILE1}" 'and' "${COLNUM_FILE2}"
  echo 'Aborting script'
  exit 1
fi
#here we play again
if [[ "${SECONDARY_INDEX}" == "YES" ]]; then

	ID2S=($(awk -F',' -v COLNUM=$COLNUM_ID2 \
	  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq))
	N_index_sequences="${#ID2S[@]}"
	ID2_LENGTH=${#ID2S[0]}
	ID2_START=($(awk -F',' -v COLNUM=$COLNUM_ID2_START \
	  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq))

	# check if number of indexes is greater than one:
	if [[ "${N_index_sequences}" -gt 1 ]]; then
		echo "Secondary indexes read from sequencing metadata (""${N_index_sequences}"" total)"
		echo
	else
	  echo
	  echo 'ERROR:' "${N_index_sequences}" 'index sequences found. There should probably be more than 1.'
	  echo
	  echo 'Aborting script.'
		exit
	fi

fi
echo "These are the secondary barcodes"
echo "${ID2S[@]}"
echo "that is, ${#ID2S[@]} unique barcodes"
echo "and they seem to be sorted alphabetically?"
echo "and they are this long "
echo "ID2_LENGTH  es ${ID2_LENGTH}"

##ANOTHER CHUNK ##
################################################################################
# Read in primers
################################################################################
PRIMER1=($(awk -F',' -v COLNUM=$COLNUM_PRIMER1 \
  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA |\
  sort | uniq ))

PRIMER2=($(awk -F',' -v COLNUM=$COLNUM_PRIMER2 \
  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA |\
  sort | uniq ))

if [[ -n "${PRIMER1}" && -n "${PRIMER2}" ]]; then
  echo 'Primers read from metadata columns' "${COLNUM_PRIMER1}" 'and' "${COLNUM_PRIMER2}"
  echo 'Primer sequences:' "${PRIMER1}" "${PRIMER2}"
	echo
else
  echo 'ERROR:' 'At least one primer is not valid'
  echo 'Looked in metadata columns' "${COLNUM_PRIMER1}" 'and' "${COLNUM_PRIMER2}"
  echo 'Aborting script'
  exit
fi

#######
#Unique samples are given by combining the primary and secondary indexes
######
ID_COMBO=$( awk -F',' -v COLNUM1=$COLNUM_ID1 -v COLNUM2=$COLNUM_ID2 \
'NR>1 {
  print ";ID1=" $COLNUM1 ";ID2=" $COLNUM2
}' "${SEQUENCING_METADATA}" )

SAMPLE_NAMES=($(awk -F',' -v COLNUM=$COLNUM_SAMPLE \
  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))

ID1_ALL=($(awk -F',' -v COLNUM=$COLNUM_ID1 \
  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))
ID1S=($(awk -F',' -v COLNUM=$COLNUM_ID1 \
  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}"  |\
		sort | uniq))
ID2_ALL=($(awk -F',' -v COLNUM=$COLNUM_ID2 \
  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))
ID2_ALL_RC=($( for i in "${ID2_ALL[@]}"; do revcom $i; done))

# write file for translating demultiplexed output to samples
SAMPLE_TRANS_FILE="${DEMULT_DIR}"/sample_trans.tmp
for (( i=0; i < "${#ID2_ALL[@]}"; i++ )); do
  printf "ID1=%s;ID2A=%s;ID2B=%s\t%s_%s\t%s\n" \
	"${ID1_ALL[i]}" "${ID2_ALL[i]}" "${ID2_ALL_RC[i]}" \
	"${ID1_ALL[i]}" "${ID2_ALL[i]}" \
	"${SAMPLE_NAMES[i]}" >> "${SAMPLE_TRANS_FILE}"
done
for (( i=0; i < "${#ID1S[@]}"; i++ )); do
  printf "File1:%s\tFile2:%s\tLib:%s\n" \
  "${FILE1[i]}" "${FILE2[i]}" "${ID1S[i]}"


done


#Create the fasta file of the barcodes

FASTA_file="$DEMULT_DIR"/barcodes.fasta
for (( i=0; i < "${#ID2S[@]}"; i++ )); do
  printf ">%s\n^NNN%s\n" \
	"${ID2S[i]}" "${ID2S[i]}" >> "${FASTA_file}"
done
primers_file="${DEMULT_DIR}"/pcr_primers.fasta

printf ">FWD\n${PRIMER1}\n>REV\n${PRIMER2}\n" > "${primers_file}"

source "${SCRIPT_DIR}"/functions/check_primers.sh "${primers_file}"

#Hooray it works
#Create a dir for all the demultiplexed files

# now we have to remove all hard-coded stuff and link it to
#banzai_params
# to get the .1 files trimmed and the .2 selected along

OUTPUT_SUMMARY="${DEMULT_DIR}/summary.csv"
printf "First_trim_R1,nReads_R1,First_trim_R2,nReads_R2,Second_trim_R1,nReads_StR1,Second_trim_R2,nReads_StR2\n" \
> "${OUTPUT_SUMMARY}"

################################################################################
# BEGIN LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
################################################################################

for (( i=0; i < "${#FILE1[@]}"; i++ )); do
  # Identify the forward and reverse fastq files.

  READ1="${PARENT_DIR}/${FILE1[i]}"
	READ2="${PARENT_DIR}/${FILE2[i]}"

  BASE1="${FILE1[i]%.*}"
  BASE2="${FILE2[i]%.*}"

  mkdir "${DEMULT_DIR}"/"${ID1S[i]}"





#DEMULT_DIR="/Users/Moncho/trial_barcode_splitter/demult"
#FASTA_file="/Users/Moncho/banzai_out_20171228_1748/barcodes.fasta"
#READ1="/Users/Moncho/Google_Drive/Run_Nov17/OA_COI/171122_sub/Lib-A_S1_L001_R1_001_sub.fastq"
#READ2="/Users/Moncho/Google_Drive/Run_Nov17/OA_COI/171122_sub/Lib-A_S1_L001_R2_001_sub.fastq"


#version from banzai, comment it out from now on


cutadapt -g file:"${FASTA_file}" -o "${DEMULT_DIR}"/${ID1S[i]}/${ID1S[i]}-{name}_round1.1.fastq -p "${DEMULT_DIR}"/${ID1S[i]}/${ID1S[i]}-{name}_round1.2.fastq \
 "${READ1}" "${READ2}" --quiet --discard-untrimmed

#This split each pair of fastqs into as many pairs of fastqs as barcodes are
#but only looking at them on the .1 file -> do the same on the other file, and keep
#the order of reads similar in both files

 for file in "${DEMULT_DIR}"/"${ID1S[i]}"/*round1.2.fastq; do
  # List all files being used:
  # file: .2.fastq; r1file, mid1,mid2, nof1, nof2, nor1,nor2
  r1file=$(echo ${file} | sed 's/.2.fastq/.1.fastq/g' ) # .1.fastq
  MID_OUTPUT1="${DEMULT_DIR}"/"${ID1S[i]}"_"${RIGHT_BARCODE}"_mid.1.fastq #double trimmed
  MID_OUTPUT2="${DEMULT_DIR}"/"${ID1S[i]}"_"${RIGHT_BARCODE}"_mid.2.fastq #double trimmed
	NEW_OUTPUT_Fwd_1="${DEMULT_DIR}"/"${ID1S[i]}"_"${RIGHT_BARCODE}"_Fwd.1.fastq
	NEW_OUTPUT_Fwd_2="${DEMULT_DIR}"/"${ID1S[i]}"_"${RIGHT_BARCODE}"_Fwd.2.fastq
	NEW_OUTPUT_Rev_1="${DEMULT_DIR}"/"${ID1S[i]}"_"${RIGHT_BARCODE}"_Rev.1.fastq
	NEW_OUTPUT_Rev_2="${DEMULT_DIR}"/"${ID1S[i]}"_"${RIGHT_BARCODE}"_Rev.2.fastq




  #NOw, shorter file names for output
  source "${SCRIPT_DIR}"/functions/strip_path.sh "${file}" file
  source "${SCRIPT_DIR}"/functions/strip_path.sh "${r1file}" r1file
  source "${SCRIPT_DIR}"/functions/strip_path.sh "${MID_OUTPUT1}" MID_OUTPUT1
  source "${SCRIPT_DIR}"/functions/strip_path.sh "${MID_OUTPUT2}" MID_OUTPUT2
	source "${SCRIPT_DIR}"/functions/strip_path.sh "${NEW_OUTPUT_Fwd_1}" NEW_OUTPUT_Fwd_1
	source "${SCRIPT_DIR}"/functions/strip_path.sh "${NEW_OUTPUT_Fwd_2}" NEW_OUTPUT_Fwd_2
	source "${SCRIPT_DIR}"/functions/strip_path.sh "${NEW_OUTPUT_Rev_1}" NEW_OUTPUT_Rev_1
	source "${SCRIPT_DIR}"/functions/strip_path.sh "${NEW_OUTPUT_Rev_2}" NEW_OUTPUT_Rev_2

  #path_to_delete="$(dirname "${file}" | sed 's_/_\\/_g')" # this turns
  # the path into \/ so it can be trimmed out
  #echo "${path_to_delete}"

  #short_file=$(echo ${file} | sed "s/${path_to_delete}//g")
  #short_r1_file=$()
  #short_NO1
  #short_NO2



  echo "${short_MID_OUTPUT2}"


  echo "${short_file}"
  nseq_file=$(cat "${file}" | wc -l)
  echo "${nseq_file} reads before retrimming"


  echo "the other half reads are (.1.)"
  echo "${short_r1file}"
  nseq_r1file=$(cat "${r1file}" |  wc -l)
  echo "with ${nseq_r1file} reads"

  RIGHT_BARCODE=$(echo ${file} | awk '/_round1/ {
    match($0, /_round1/); print substr($0, RSTART - 6, 6);
    }')

  #Now use that as an argunment for cutadapt
  echo "this is the right barcode"
  echo ${RIGHT_BARCODE}


  cutadapt -g ^NNN"${RIGHT_BARCODE}" -o "${MID_OUTPUT2}" \
  -p "${MID_OUTPUT1}" "${file}" "${r1file}" --quiet --discard-untrimmed


  nseq_s2r1file=$(cat "${MID_OUTPUT1}" |  wc -l)
  nseq_s2r2file=$(cat "${MID_OUTPUT2}" |  wc -l)
  echo "This is the mid R1 file"
  echo "${short_MID_OUTPUT1}"
  echo "and it has ${nseq_s2r1file} reads after trimming"
  echo "Hopefully the same number of lines as the mid R2 ${nseq_s2r2file}"

# Now remove the pcr primers
# This is an important point for libraries prepared by ligation:
# The i7 adapters can ligate on either Fwd or Rev primer end- so you have
# roughly half the sequence in one direction and half on the other direction
# What to do with them is up to you: we'll generate 4 fastqs per sample
#FWD.1
#FWD.2
#REV.1
#REV.2
#You can either choose one pair and discard 50% of your data,
#Add #rev to the header of those affected and leave them as they are
#Add #rev and RC those affected
# do the analysis twice
#
#First remove the primers from .1 and select those from the file and then we'll see
#will later assume that if the barcodes were fine, then the primer will be fine too

cutadapt -g file:"${primers_file}" --discard-untrimmed\
 -o "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_{name}_clean.1.fastq \
 -p "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_{name}_clean.2.fastq \
 "${MID_OUTPUT1}" "${MID_OUTPUT2}"

#Now remove the rev primer at the beggining of the .2 for those READS
#in which we found the FWD primer at the beggining of .1
cutadapt -g "${PRIMER2}" --discard-untrimmed \
-o "${NEW_OUTPUT_Fwd_2}" \
-p "${NEW_OUTPUT_Fwd_1}" \
"${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_clean.2.fastq \
"${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_clean.1.fastq

#Now do similarly for those in which we found rev at the beggining of .1

cutadapt -g "${PRIMER1}" --discard-untrimmed \
-o "${NEW_OUTPUT_Rev_2}" \
-p "${NEW_OUTPUT_Rev_1}" \
"${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_clean.2.fastq \
"${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_clean.1.fastq


nseq_NOF1=$(cat ${NEW_OUTPUT_Fwd_1} | wc -l)
nseq_NOR1=$(cat ${NEW_OUTPUT_Rev_1} | wc -l)


#print the summary information
  printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
  "${short_r1file}" "${nseq_r1file}" \
  "${short_file}" "${nseq_file}" \
  "${short_MID_OUTPUT1}" "${nseq_s2r1file}" \
  "${short_MID_OUTPUT2}" "${nseq_s2r2file}" \
	"${short_NEW_OUTPUT_Fwd_1}" "${nseq_NOF1}" \
  "${short_NEW_OUTPUT_Rev_1}" "${nseq_NOR1}" >> "${OUTPUT_SUMMARY}"
  # now clean the middle FILES
  rm "${file}"
  rm "${r1file}"
	rm "${MID_OUTPUT1}"
	rm "${MID_OUTPUT2}"
	rm "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_clean.2.fastq
	rm "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_clean.1.fastq
	rm "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_clean.2.fastq
	rm "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_clean.1.fastq



  done
  rm -r "${DEMULT_DIR}"/"${ID1S[i]}"





done

if "${SEARCH_ASVs}"=="YES"; then
	Rscript "${SCRIPT_DIR}"/r/dada2.r "${DEMULT_DIR}" "${SCRIPT_DIR}"
fi
