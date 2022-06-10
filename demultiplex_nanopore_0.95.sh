#!/bin/bash
# Start by activating decona
# Usage bash demultiplex_both_fastqs.sh banzai_params.sh ...
#This script is built using banzai (github.com/jimmyodonnell/banzai) as template

# Can we overwrite any arguments?

#We need to gather: Location of functions  and fastqs:
MAIN_DIR="$(dirname "$0")"
Current_dir=$(pwd)
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
OUTPUT_DIR="/home/mk1b/test_ns/n_0.95"

mkdir "${OUTPUT_DIR}"
echo "Output directory is ${OUTPUT_DIR}"
# copy metadata and parameters file to output directory
cp "${SEQUENCING_METADATA}" "${OUTPUT_DIR}"/metadata.csv
cp "${param_file}" "${OUTPUT_DIR}"/banzai_params.sh

# Write a log file
LOGFILE="${OUTPUT_DIR}"/logfile.txt
exec > >(tee "${LOGFILE}") 2>&1

FINAL_DIR="${OUTPUT_DIR}"/noprimers
mkdir "${FINAL_DIR}"

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
#COLNUM_FILE2=$( get_colnum "${COLNAME_FILE2}" "${SEQUENCING_METADATA}")
# Pass check
# Library names
COLNUM_ID1=$( get_colnum "${COLNAME_ID1_NAME}" "${SEQUENCING_METADATA}")

COLNUM_ID1_SEQ=$( get_colnum "${COLNAME_ID1_SEQ}" "${SEQUENCING_METADATA}")

# Secondary indices
COLNUM_ID2=$( get_colnum "${COLNAME_ID2_SEQ}" "${SEQUENCING_METADATA}")


# Sample names
COLNUM_SAMPLE=$( get_colnum "${COLNAME_SAMPLE_ID}" "${SEQUENCING_METADATA}")

# Primers
COLNUM_PRIMER1=$( get_colnum "${COLNAME_PRIMER1}" "${SEQUENCING_METADATA}")
COLNUM_PRIMER2=$( get_colnum "${COLNAME_PRIMER2}" "${SEQUENCING_METADATA}")
COLNUM_LOCUS=$( get_colnum "${COLNAME_LOCUS}" "${SEQUENCING_METADATA}")

# Run away from the script if any of the previous columns was not found

all_columns=( COLNUM_FILE1 COLNUM_ID1 COLNUM_ID1_SEQ COLNUM_ID2 \
COLNUM_SAMPLE COLNUM_PRIMER1 COLNUM_PRIMER2 COLNUM_LOCUS)

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
# if [[ "${ADD_TO_PREVIOUS}" = "YES" ]]; then
# 	echo "You chose to add this analysis to a previous one"
# 	if [[ -n "${FORMER_HASH}" && -n "${FORMER_ABUNDANCE}" ]]; then
# 		echo "Using hash database from ${FORMER_HASH}"
# 		echo "Using ASV table from ${FORMER_ABUNDANCE}"
# 	else
# 		echo "Uppss, at least one of these files is missing"
# 		echo " - A Hash / sequence conversion table"
# 		echo " - An Abundance dataset "
# 		echo "Set the path to these files in the params file"
# 		exit
# 	fi
#
# 	if [[ -s "${LOG_FILE}" ]] ; then
# 		echo "Adding merge information to ${LOG_FILE}"
# 	else
# 		echo "No logfile provided or found"
# 		echo "Starting a new merge logfile"
# 		LOG_FILE="${OUTPUT_DIR}"/database_log.csv
# 		echo "New file is ${LOG_FILE}"
# 	fi
# fi


################################################################################
# CHECK FILES
################################################################################

#Check if we are redoing the analysis after demultiplexing
	FILE1=($(awk -F',' -v COLNUM=$COLNUM_FILE1 \
	  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq))

	# FILE2=($(awk -F',' -v COLNUM=$COLNUM_FILE2 \
	#   'NR>1 {print $COLNUM}' $SEQUENCING_METADATA |\
	#   sort | uniq ))

	NFILE1="${#FILE1[@]}"
	# NFILE2="${#FILE2[@]}"
	# if [ "${NFILE1}" != "${NFILE2}" ]; then
	# 	echo "ERROR: Whoa! different number of forward and reverse files"
	# fi

	# if [[ -n "${FILE1}" && -n "${FILE2}" ]]; then
	  echo 'Files read from metadata columns' "${COLNUM_FILE1}"
	  echo 'File names:'
		for (( i=0; i < "${NFILE1}"; ++i)); do
			printf '%s\t%s\n' "${FILE1[i]}"
		done
		echo
	# else
	#   echo 'ERROR:' 'At least one file is not valid'
	#   echo 'Looked in metadata columns' "${COLNUM_FILE1}"
	#   echo 'Aborting script'
	#   exit
	# fi
	#here we play again
	if [[ "${SECONDARY_INDEX}" == "YES" ]]; then

		ID2S=($(awk -F',' -v COLNUM=$COLNUM_ID2 \
		  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
		  sort | uniq))
		N_index_sequences="${#ID2S[@]}"
		ID2_LENGTH=${#ID2S[0]}
		# ID2_START=($(awk -F',' -v COLNUM=$COLNUM_ID2_START \
		#   'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
		#   sort | uniq))

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

	LOCI=($(awk -F',' -v COLNUM=$COLNUM_LOCUS \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq ))

	PRIMER2_RC=($( for i in "${PRIMER2[@]}"; do revcom $i; done))

	if [[ -n "${PRIMER1}" && -n "${PRIMER2}" ]]; then
	  echo 'Primers read from metadata columns' "${COLNUM_PRIMER1}" 'and' "${COLNUM_PRIMER2}"
	  echo 'Primer sequences:' "${PRIMER1}" "${PRIMER2}"
		echo
		echo 'Reverse Complemented ' "${PRIMER2}" 'into' "${PRIMER2_RC}"
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
	ID1S=($(awk -F',' -v COLNUM=$COLNUM_ID1_SEQ \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}"  |\
			sort | uniq))

	ID1S_many=($(awk -F',' -v COLNUM=$COLNUM_ID1_SEQ \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}"))

	echo "This is what I was printing"

	echo "${ID1S}"

	echo "Is this better?"

	echo "${ID1S_many}"

  echo "Which column is this"

  echo "${COLNUM_ID1_SEQ}"

  echo "Show me column 3 then"




	ID2_ALL=($(awk -F',' -v COLNUM=$COLNUM_ID2 \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))
	ID2_ALL_RC=($( for i in "${ID2_ALL[@]}"; do revcom $i; done))

  echo "${ID2_ALL_RC}"

# write file for translating demultiplexed output to samples
	SAMPLE_TRANS_FILE="${OUTPUT_DIR}"/sample_trans.tmp
	for (( i=0; i < "${#ID2_ALL[@]}"; i++ )); do
	  printf "ID1=%s;ID2A=%s;ID2B=%s\t%s_%s\t%s\n" \
		"${ID1_ALL[i]}" "${ID2_ALL[i]}" "${ID2_ALL_RC[i]}" \
		"${ID1S_many[i]}" "${ID2_ALL[i]}" \
		"${SAMPLE_NAMES[i]}" >> "${SAMPLE_TRANS_FILE}"
	done
	for (( i=0; i < "${#ID1S[@]}"; i++ )); do
	  printf "File1:%s\tFile2:%s\tLib:%s\n" \
	  "${FILE1[i]}" "${FILE2[i]}" "${ID1S[i]}"


	done
	# Create the BARCODE FILE for nanopore
	Barcodes_file="$OUTPUT_DIR"/barcodes.fasta


	for (( i=0; i < "${#ID2_ALL[@]}"; i++ )); do


	printf ">%s_%s\n%s"..."%s\n" \
		"${ID1_ALL[i]}" "${ID2_ALL[i]}" \
		"${ID1S_many[i]}" "${ID2_ALL_RC[i]}" >> "${Barcodes_file}"
	done

	head "${Barcodes_file}"

# #Create the fasta file of the barcodes
#
# 	Barcodes_file="$OUTPUT_DIR"/barcodes.fasta
# 	for (( i=0; i < "${#ID2S[@]}"; i++ )); do
# 	  printf ">%s\n^NNN%s\n" \
# 		"${ID2S[i]}" "${ID2S[i]}" >> "${Barcodes_file}"
# 	done

	primers_file="${OUTPUT_DIR}"/pcr_primers.fasta

	printf ">${LOCI}\n${PRIMER1}...${PRIMER2_RC}\n" > "${primers_file}"

	source "${SCRIPT_DIR}"/functions/check_primers.sh "${primers_file}"

#Hooray it works
#Create a dir for all the demultiplexed files

# now we have to remove all hard-coded stuff and link it to
#banzai_params
# to get the .1 files trimmed and the .2 selected along

	OUTPUT_SUMMARY="${OUTPUT_DIR}/summary.csv"
	printf "Sample,locus,demultiplexed,noprimers\n" \
	> "${OUTPUT_SUMMARY}"

################################################################################
# BEGIN LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
################################################################################

	for (( i=0; i < "${#FILE1[@]}"; i++ )); do
	  # Identify the forward and reverse fastq files.

	  READ1="${PARENT_DIR}/${FILE1[i]}"
		# READ2="${PARENT_DIR}/${FILE2[i]}"

	  BASE1="${FILE1[i]%.*}"
	  # BASE2="${FILE2[i]%.*}"

	  mkdir "${DEMULT_DIR}"/"${FILE1[i]}"
		mkdir "${FINAL_DIR}"/"${FILE1[i]}"



		echo "Working on Library $[i+1] out of ${#FILE1[@]}"

	##First cutdapt:


	# Look for the reverse primer and if found, reverse the output

	cutadapt -g "${PRIMER1}" -o "${READ1}".new.fastq --quiet --action=none --rc -e 0.3 "${READ1}"

  echo "Number of lines to begin with"

  wc -l "${READ1}"

  echo "Number of lines after"

  wc -l "${READ1}".new.fastq

  cat "${Barcodes_file}"

	cutadapt -g file:"${Barcodes_file}" -o "${DEMULT_DIR}"/${FILE1[i]}/${FILE1[i]}-{name}_round1.fastq \
	 "${READ1}".new.fastq --quiet --discard-untrimmed -e 0.2 --rc


	#This split each fastq into as fastqs as barcodes are
	#but only looking at them on the .1 file -> do the same on the other file, and keep
	#the order of reads similar in both files

		n_files=("${DEMULT_DIR}"/"${FILE1[i]}"/*round1.fastq)



	 ls "${DEMULT_DIR}"/"${FILE1[i]}"/

	 for file in "${n_files[@]}"; do

		 echo "working on file " "${file}"
		 echo ""


	 BASE_OUTPUT=$(basename "${file}" |  sed 's/_round1.fastq//g')

	 echo "${BASE_OUTPUT}"

	 mkdir "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"

	 	cutadapt -g file:"${primers_file}" --discard-untrimmed \
	 -o "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/"${BASE_OUTPUT}"_{name}.fastq \
	 "${file}"  --quiet -e 0.2 >> "${LOGFILE}"

	nseq_demult=$(cat "${file}" | wc -l)

	for file2 in "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/; do


	nseq_noprimer=$(cat ${file2} | wc -l)

	 printf "%s,%s,%s,%s\n" \
	 "${BASE_OUTPUT}" "${file2}" \
	  "${nseq_demult}" "${nseq_noprimer}" >> "${OUTPUT_SUMMARY}"
	done

cd "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"



decona -l "${MIN_LENGTH}" -m "${MAX_LENGTH}" -q 10 -c 0.80 -n 10 -M


done # End of loop across all samples within a file


	done

cd "${Current_dir}"

# if [[ "${SEARCH_ASVs}" = "YES" ]]; then
	echo "launching  gather decona"
# 	conda activate decona
#
	Rscript --vanilla "${SCRIPT_DIR}"/r/gather.decona.r "${OUTPUT_DIR}"
# fi
