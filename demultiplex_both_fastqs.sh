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

# Check if the metadata file exists
if [[ -s "${SEQUENCING_METADATA}" ]]; then
	echo "Reading metadata from:"
	echo "${SEQUENCING_METADATA}"
else
	echo 'ERROR! Could not find metadata file. You specified the file path:'
	echo
	echo "${SEQUENCING_METADATA}"
	echo
	echo 'That file is empty or does not exist. Aborting script.'
	exit
fi


# Now fix line ends if needed

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
OUTPUT_DIR="${OUTPUT_DIRECTORY}"/demultiplexed_"${START_TIME}"

mkdir "${OUTPUT_DIR}"
echo "Output directory is ${OUTPUT_DIR}"
# copy metadata and parameters file to output directory
cp "${SEQUENCING_METADATA}" "${OUTPUT_DIR}"/metadata.csv
cp "${param_file}" "${OUTPUT_DIR}"/banzai_params.sh

# Write a log file
LOGFILE="${OUTPUT_DIR}"/logfile.txt
exec > >(tee "${LOGFILE}") 2>&1


mkdir "${OUTPUT_DIR}"/cleaned
DEMULT_DIR="${OUTPUT_DIR}"/demultiplexed
mkdir "${DEMULT_DIR}"
################################################################################
# READ METADATA
################################################################################
# report metadata dimensions
METADATA_DIM=($( awk -F, 'END{print NR, NF}' "${SEQUENCING_METADATA}" ))
echo "Metadata has" "${METADATA_DIM[0]}" "rows and" "${METADATA_DIM[1]}" "columns including header."
N_SAMPLES=$( echo "${METADATA_DIM[0]}" - 1 | bc )
echo "Expecting" "${N_SAMPLES}" "samples total."
echo
## NOW WE HAVE LOADED THE SEQUENCING_METADATA - WE NEED to find the columns specified
## in the params file. We should set up an alert & quit if a critical column is not found

# Filnames
COLNUM_FILE1=$( get_colnum "${COLNAME_FILE1}" "${SEQUENCING_METADATA}")
COLNUM_FILE2=$( get_colnum "${COLNAME_FILE2}" "${SEQUENCING_METADATA}")
# Pass check
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

# Run away from the script if any of the previous columns was not found

all_columns=( COLNUM_FILE1 COLNUM_FILE2 COLNUM_ID1 COLNUM_ID1_SEQ COLNUM_ID2 \
COLNUM_ID2_START COLNUM_SAMPLE COLNUM_PRIMER1 COLNUM_PRIMER2)

echo "Checking that all columns in metadata are there"

for column in "${all_columns[@]}" ; do

 if [ "${!column}" > 0 ]; then
	 echo "looking good, ${column}"
 else
  echo "Something went wrong with column name ${column}"
	echo "exiting script"
	exit
fi
done
echo "All columns passed test"


################################################################################
# ADDING TO PREVIOUS ANALYSIS?
################################################################################
if [[ "${ADD_TO_PREVIOUS}" = "YES" ]]; then
	echo "You chose to add this analysis to a previous one"
	if [[ -n "${FORMER_HASH}" && -n "${FORMER_ABUNDANCE}" ]]; then
		echo "Using hash database from ${FORMER_HASH}"
		echo "Using ASV table from ${FORMER_ABUNDANCE}"
	else
		echo "Uppss, at least one of these files is missing"
		echo " - A Hash / sequence conversion table"
		echo " - An Abundance dataset "
		echo "Set the path to these files in the params file"
		exit
	fi

	if [[ -s "${LOG_FILE}" ]] ; then
		echo "Adding merge information to ${LOG_FILE}"
	else
		echo "No logfile provided or found"
		echo "Starting a new merge logfile"
		LOG_FILE="${OUTPUT_DIR}"/database_log.csv
		echo "New file is ${LOG_FILE}"
	fi
fi


################################################################################
# CHECK FILES
################################################################################

#Check if we are redoing the analysis after demultiplexing
if [[ "${ALREADY_DEMULTIPLEXED}" != "YES" ]]; then


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
	  exit
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

#####
# Check that sample names are not repeated
#####
NSAMPLES="${#SAMPLE_NAMES[@]}"

# Now calculate the number of unique sample names
UNIQ_SAMPLES=( $(echo "${SAMPLE_NAMES[@]}" | tr ' ' '\n' | sort -u))
N_UNIQ_SAMPLES="${#UNIQ_SAMPLES[@]}"


if [[ "${NSAMPLES}" != "${N_UNIQ_SAMPLES}" ]]; then
	echo " At least one sample name is repeated "
	echo " I am not angry, just dissapointed. Exiting script"
	exit
fi


	ID1_ALL=($(awk -F',' -v COLNUM=$COLNUM_ID1 \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))
	ID1S=($(awk -F',' -v COLNUM=$COLNUM_ID1 \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}"  |\
			sort | uniq))
	ID2_ALL=($(awk -F',' -v COLNUM=$COLNUM_ID2 \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))
	ID2_ALL_RC=($( for i in "${ID2_ALL[@]}"; do revcom $i; done))

# write file for translating demultiplexed output to samples
	SAMPLE_TRANS_FILE="${OUTPUT_DIR}"/sample_trans.tmp
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

	Barcodes_file="$OUTPUT_DIR"/barcodes.fasta
	for (( i=0; i < "${#ID2S[@]}"; i++ )); do
	  printf ">%s\n^NNN%s\n" \
		"${ID2S[i]}" "${ID2S[i]}" >> "${Barcodes_file}"
	done

	primers_file="${OUTPUT_DIR}"/pcr_primers.fasta

	printf ">FWD\n${PRIMER1}\n>REV\n${PRIMER2}\n" > "${primers_file}"

	source "${SCRIPT_DIR}"/functions/check_primers.sh "${primers_file}"

#Hooray it works
#Create a dir for all the demultiplexed files

# now we have to remove all hard-coded stuff and link it to
#banzai_params
# to get the .1 files trimmed and the .2 selected along

	OUTPUT_SUMMARY="${OUTPUT_DIR}/summary.csv"
	printf "First_trim_R1,nReads_R1,First_trim_R2,nReads_R2,Second_trim_R1,nReads_StR1,Second_trim_R2,nReads_StR2,NEW_OUTPUT_Fwd_1,nreadsFwd,NEW_OUTPUT_Rev_1,nreadsRev\n" \
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

	  mkdir "${OUTPUT_DIR}"/"${ID1S[i]}"


		mkdir "${OUTPUT_DIR}"/cleaned/"${ID1S[i]}"

		echo "Working on Library $[i+1] out of ${#FILE1[@]}"

	##First cutdapt:
	#TODO: use only the number of barcodes used for this Library


	cutadapt -g file:"${Barcodes_file}" -o "${OUTPUT_DIR}"/${ID1S[i]}/${ID1S[i]}-{name}_round1.1.fastq -p "${OUTPUT_DIR}"/${ID1S[i]}/${ID1S[i]}-{name}_round1.2.fastq \
	 "${READ1}" "${READ2}" --quiet --discard-untrimmed


	#This split each pair of fastqs into as many pairs of fastqs as barcodes are
	#but only looking at them on the .1 file -> do the same on the other file, and keep
	#the order of reads similar in both files

		n_files=("${OUTPUT_DIR}"/"${ID1S[i]}"/*round1.2.fastq)

		i_count=0

	 for file in "${n_files[@]}"; do
# We loop through all .2 files
		i_count=$((i_count+1))
 #The barcode detected on the .1 is written in the name, so we now look
 #for that barcode at the beggining of the .2 read

		RIGHT_BARCODE=$(echo ${file} | awk '/_round1/ {
	    match($0, /_round1/); print substr($0, RSTART - 6, 6);
	    }')

			short_file=$(basename "${file}")

		r1file=$(echo ${file} | sed 's/.2.fastq/.1.fastq/g' )
		 	short_r1file=$(basename "${r1file}") # .1.fastq
	  MID_OUTPUT1="${OUTPUT_DIR}"/"${ID1S[i]}"_"${RIGHT_BARCODE}"_mid.1.fastq
			short_MID_OUTPUT1=$(basename "${MID_OUTPUT1}") #double trimmed
	  MID_OUTPUT2="${OUTPUT_DIR}"/"${ID1S[i]}"_"${RIGHT_BARCODE}"_mid.2.fastq #double trimmed
			short_MID_OUTPUT2=$(basename "${MID_OUTPUT2}")
		NEW_OUTPUT_Fwd_1="${DEMULT_DIR}"/"${ID1S[i]}"_"${RIGHT_BARCODE}"_Fwd.1.fastq
			short_NEW_OUTPUT_Fwd_1=$(basename "${NEW_OUTPUT_Fwd_1}")
		NEW_OUTPUT_Fwd_2="${DEMULT_DIR}"/"${ID1S[i]}"_"${RIGHT_BARCODE}"_Fwd.2.fastq
			short_NEW_OUTPUT_Fwd_2=$(basename "${NEW_OUTPUT_Fwd_2}")
		NEW_OUTPUT_Rev_1="${DEMULT_DIR}"/"${ID1S[i]}"_"${RIGHT_BARCODE}"_Rev.1.fastq
			short_NEW_OUTPUT_Rev_1=$(basename "${NEW_OUTPUT_Rev_1}")
		NEW_OUTPUT_Rev_2="${DEMULT_DIR}"/"${ID1S[i]}"_"${RIGHT_BARCODE}"_Rev.2.fastq
			short_NEW_OUTPUT_Rev_2=$(basename "${NEW_OUTPUT_Rev_2}")

		#New messages so it's easier to see the progress of the script

		echo -ne "Working on sample ${i_count} of ${#n_files[@]}"'\r'

		#echo "${short_file}"
	  nseq_file=$(cat "${file}" | wc -l)
	  #echo "${nseq_file} reads before retrimming"


	  #echo "the other half reads are (.1.)"
	  #echo "${short_r1file}"
	  nseq_r1file=$(cat "${r1file}" |  wc -l)
	  #echo "with ${nseq_r1file} reads"


	  #Now use that as an argunment for cutadapt
	  #echo "this is the right barcode"
	  #echo ${RIGHT_BARCODE}

# try to make cutadapt quieter
	  cutadapt -g ^NNN"${RIGHT_BARCODE}" -o "${MID_OUTPUT2}" \
	  -p "${MID_OUTPUT1}" "${file}" "${r1file}" --quiet --discard-untrimmed 2>> "${LOGFILE}"

	  nseq_s2r1file=$(cat "${MID_OUTPUT1}" |  wc -l)
	  nseq_s2r2file=$(cat "${MID_OUTPUT2}" |  wc -l)
	  #echo "This is the mid R1 file"
	  #echo "${short_MID_OUTPUT1}"
	  #echo "and it has ${nseq_s2r1file} reads after trimming"
	  #echo "Hopefully the same number of lines as the mid R2 ${nseq_s2r2file}"


	# Now remove the pcr primers
	# This is an important point for libraries prepared by ligation:
	# The i7 adapters can ligate on either Fwd or Rev primer end- so you have
	# roughly half the sequences in one direction and half on the other direction
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
	 -o "${OUTPUT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_{name}_clean.1.fastq \
	 -p "${OUTPUT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_{name}_clean.2.fastq \
	 "${MID_OUTPUT1}" "${MID_OUTPUT2}" --quiet 2>> "${LOGFILE}"




	#Now remove the rev primer at the beggining of the .2 for those READS
	#in which we found the FWD primer at the beggining of .1
	cutadapt -g "${PRIMER2}" --discard-untrimmed \
	-o "${NEW_OUTPUT_Fwd_2}" \
	-p "${NEW_OUTPUT_Fwd_1}" \
	"${OUTPUT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_clean.2.fastq \
	"${OUTPUT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_clean.1.fastq --quiet 2>> "${LOGFILE}"

	#Now do similarly for those in which we found rev at the beggining of .1

	cutadapt -g "${PRIMER1}" --discard-untrimmed \
	-o "${NEW_OUTPUT_Rev_2}" \
	-p "${NEW_OUTPUT_Rev_1}" \
	"${OUTPUT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_clean.2.fastq \
	"${OUTPUT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_clean.1.fastq --quiet 2>> "${LOGFILE}"


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
	  # now clean the middle FILES - checking first if they do exist
	  rm "${file}"
	  rm "${r1file}"
		if [[ -s "${MID_OUTPUT1}" ]]; then
			rm "${MID_OUTPUT1}"
			rm "${MID_OUTPUT2}"
		fi
		if [[ -s "${OUTPUT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_clean.1.fastq ]]; then
			rm "${OUTPUT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_clean.2.fastq
			rm "${OUTPUT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_clean.1.fastq
		fi
		if [[ -s "${OUTPUT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_clean.2.fastq ]]; then
			rm "${OUTPUT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_clean.2.fastq
			rm "${OUTPUT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_clean.1.fastq
		fi

	  done
	  rm -r "${OUTPUT_DIR}"/"${ID1S[i]}"

	done

	rm -rf "${OUTPUT_DIR}"/cleaned

else #In case you already demultiplexed your samples, then cp the files you need
	cp "${DEMULT_OUTPUT}"/sample_trans.tmp "${OUTPUT_DIR}"
	cp "${DEMULT_OUTPUT}"/barcodes.fasta "${OUTPUT_DIR}"
	cp "${DEMULT_OUTPUT}"/summary.csv "${OUTPUT_DIR}"
	cp "${DEMULT_OUTPUT}"/pcr_primers.fasta "${OUTPUT_DIR}"

	DEMULT_DIR="${DEMULT_OUTPUT}"/demultiplexed

fi #This finishes the control flow in case you already demultiplexed
# We are selecting a pair of fastq files so we can check the direction of the
# ASVs
FILE1=($(awk -F',' -v COLNUM=$COLNUM_FILE1 \
	'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
	sort | uniq))

FILE2=($(awk -F',' -v COLNUM=$COLNUM_FILE2 \
	'NR>1 {print $COLNUM}' $SEQUENCING_METADATA |\
	sort | uniq ))
READ1="${PARENT_DIR}/${FILE1[1]}"
READ2="${PARENT_DIR}/${FILE2[1]}"

if [[ "${SEARCH_ASVs}" = "YES" ]]; then
	echo "This is read1 ${READ1}"
	Rscript "${SCRIPT_DIR}"/r/dada2.r "${OUTPUT_DIR}" "${DEMULT_DIR}" "${SCRIPT_DIR}" "${USE_HASH}" "${READ1}" "${READ2}"\
	"${ADD_TO_PREVIOUS}" "${FORMER_HASH}" "${FORMER_ABUNDANCE}" "${LOG_FILE}"
fi
