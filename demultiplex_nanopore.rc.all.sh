#!/bin/bash
# Start by activating decona
# Usage bash demultiplex_both_fastqs.sh banzai_params.sh ...
#This script is built using banzai (github.com/jimmyodonnell/banzai) as template
# THIS VERSION OF the script reverse complements the strings so we can look at the right occurrences of the barcodes and primers

##### LOCATING SCRIPTS AND SOURCING

MAIN_DIR="$(dirname "$0")"
Current_dir=$(pwd)
SCRIPT_DIR="${MAIN_DIR}"/scripts
for file in "${SCRIPT_DIR}"/*.sh ; do
	echo "sourcing script ${file}"
	source "${file}"
done

###### LOCATING PARAMETERS FILE

param_file=${1}

echo "Reading analysis parameters from:"
echo "${param_file}"
source "${param_file}"

###### LOCATING METADATA

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

###### CREATING OUTPUT DIRECTORY

START_TIME=$(date +%Y%m%d_%H%M)
OUTPUT_DIR="${OUTPUT_DIRECTORY}"/demultiplexed_"${START_TIME}"

mkdir "${OUTPUT_DIR}"
echo "Output directory is ${OUTPUT_DIR}"
echo

# copy metadata and parameters file to output directory
cp "${SEQUENCING_METADATA}" "${OUTPUT_DIR}"/metadata.csv
SEQUENCING_METADATA="${OUTPUT_DIR}"/metadata.csv
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

echo " #### Reading metadata file"
echo ""

echo "Metadata has" "${METADATA_DIM[0]}" "rows and" "${METADATA_DIM[1]}" "columns including header."
N_SAMPLES=$( echo "${METADATA_DIM[0]}" - 1 | bc )
echo "Expecting" "${N_SAMPLES}" "samples total."
echo
## NOW WE HAVE LOADED THE SEQU ENCING_METADATA - WE NEED to find the columns specified
## in the params file. We should set up an alert & quit if a critical column is not found

# Filenames
COLNUM_FILE1=$( get_colnum "${COLNAME_FILE1}" "${SEQUENCING_METADATA}")



# PLATE names and SEQS
COLNUM_ID1=$( get_colnum "${COLNAME_ID1_NAME}" "${SEQUENCING_METADATA}")

COLNUM_ID1_SEQ=$( get_colnum "${COLNAME_ID1_SEQ}" "${SEQUENCING_METADATA}")

# WELL Names and SEQS
COLNUM_ID2_SEQ=$( get_colnum "${COLNAME_ID2_SEQ}" "${SEQUENCING_METADATA}")
COLNUM_ID2_WELL=$(get_colnum "${COLNAME_ID2_WELL}" "${SEQUENCING_METADATA}")

# Sample names
COLNUM_SAMPLE=$( get_colnum "${COLNAME_SAMPLE_ID}" "${SEQUENCING_METADATA}")

# Primers
COLNUM_PRIMER1=$( get_colnum "${COLNAME_PRIMER1}" "${SEQUENCING_METADATA}")
COLNUM_PRIMER2=$( get_colnum "${COLNAME_PRIMER2}" "${SEQUENCING_METADATA}")
COLNUM_LOCUS=$( get_colnum "${COLNAME_LOCUS}" "${SEQUENCING_METADATA}")

# Run away from the script if any of the previous columns was not found

all_columns=( COLNUM_FILE1 COLNUM_ID1 COLNUM_ID1_SEQ COLNUM_ID2_SEQ COLNUM_ID2_WELL \
COLNUM_SAMPLE COLNUM_PRIMER1 COLNUM_PRIMER2 COLNUM_LOCUS)

echo "Checking that all columns in metadata are there"

for column in "${all_columns[@]}" ; do

 if [ "${!column}" -gt 0 ]; then
	 echo "looking good, ${column}"
 else
  echo "Something went wrong with column name ${column}"
	echo "exiting script"
	exit
fi
done
echo "All columns passed test"

#############
# ADD Reverse complements of Well indices and PCR primers
#############

	 ## First rc the i7, then join it to the metadata

awk -F',' -v COLNUM="${COLNUM_ID2_SEQ}"    \
'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" | \
while read line ; do revcom $line ; done |\
sed '1s/^/RC_i7\n/'  > "${OUTPUT_DIR}"/rci7.txt

paste "${SEQUENCING_METADATA}" "${OUTPUT_DIR}"/rci7.txt -d ',' > "${OUTPUT_DIR}"/metadata.new.csv

   ## First rc the PCR PRIMER, then join it to the metadata


awk -F',' -v COLNUM="${COLNUM_PRIMER2}"    \
'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" | \
while read line ; do revcom $line ; done |\
sed '1s/^/RC_PCR_REV\n/'  > "${OUTPUT_DIR}"/rc_pcr.txt

paste "${OUTPUT_DIR}"/metadata.new.csv "${OUTPUT_DIR}"/rc_pcr.txt -d ',' > "${OUTPUT_DIR}"/metadata.new.new.csv

mv "${OUTPUT_DIR}"/metadata.new.new.csv "${SEQUENCING_METADATA}"

rm "${OUTPUT_DIR}"/metadata.new.csv


##LOCATE the COLUMN OF the RCs 
COLNUM_PRIMER2_RC=$( get_colnum "RC_PCR_REV" "${SEQUENCING_METADATA}")
COLNUM_ID2_RCSEQ=$( get_colnum "RC_i7" "${SEQUENCING_METADATA}")
## Update metadata DIM

METADATA_DIM=($( awk -F, 'END{print NR, NF}' "${SEQUENCING_METADATA}" ))

echo "Metadata has" "${METADATA_DIM[0]}" "rows and" "${METADATA_DIM[1]}" "columns including header."


################################################################################
# CHECK FILES
################################################################################


	FILE1=($(awk -F',' -v COLNUM=$COLNUM_FILE1 \
	  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq))



	NFILE1="${#FILE1[@]}"


	  echo 'Files read from metadata columns' "${COLNUM_FILE1}"
	  echo 'File names:'
		for (( i=0; i < "${NFILE1}"; ++i)); do
			printf '%s\t%s\n' "${FILE1[i]}"
		done
		echo

    echo "Primary Indices read from column ${COLNAME_ID1_SEQ}"

	ID1_ALL=($(awk -F',' -v COLNAME=$COLNUM_ID1 -v COLSEQ=$COLNUM_ID1_SEQ \
	  'NR>1 { print $COLNAME, $COLSEQ }' "${SEQUENCING_METADATA}" |\
			sort | uniq))
	NID1="${#ID1_ALL[@]}"

	echo "There are ${NID1} primary indices"

	awk -F',' -v COLNAME=$COLNUM_ID1 -v COLSEQ=$COLNUM_ID1_SEQ \
	  'NR>1 { print $COLNAME, $COLSEQ }' $SEQUENCING_METADATA |sort | uniq |\
		awk '{printf ">%s\n%s\n", $1, $2}' >> "${OUTPUT_DIR}"/barcodes_P5.fasta

cat "${OUTPUT_DIR}"/barcodes_P5.fasta



################################################################################
# Read in primers
################################################################################
# Modify this so in case there are multiple primer pairs, they are not artificially paired
# alphabetically

	PRIMER1=($(awk -F',' -v COLNUM=$COLNUM_PRIMER1 \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA ))

	PRIMER2=($(awk -F',' -v COLNUM=$COLNUM_PRIMER2 \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA ))

	LOCI=($(awk -F',' -v COLNUM=$COLNUM_LOCUS \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA ))

	PRIMER2_RC=($(awk -F',' -v COLNUM=$COLNUM_PRIMER2_RC \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA ))

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

##############################################
# Create a fasta with all Fwds to get all sequences in the same direction
############################################

awk -F','  -v FWD=$COLNUM_PRIMER1 -v LOCI=$COLNUM_LOCUS \
	  ' NR > 1 { print $LOCI , $FWD }' $SEQUENCING_METADATA | sort | uniq  |\
	  awk '{printf ">%s\n%s\n", $1, $2}' >> "${OUTPUT_DIR}"/direction.fasta


#######
#Unique samples are given by combining the primary and secondary indexes
######
	ID_COMBO=$( awk -F',' -v COLNUM1=$COLNUM_ID1 -v COLNUM2=$COLNUM_ID2_WELL \
	'NR>1 {
	  print "Plate=" $COLNUM1 ";Well=" $COLNUM2
	}' "${SEQUENCING_METADATA}" )

	SAMPLE_NAMES=($(awk -F',' -v COLNUM=$COLNUM_SAMPLE \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))



	OUTPUT_SUMMARY="${OUTPUT_DIR}/summary.csv"
	printf "Sample,locus,demultiplexed,noprimers\n" \
	> "${OUTPUT_SUMMARY}"

################################################################################
# BEGIN LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
################################################################################

	for (( i=0; i < "${#FILE1[@]}"; i++ )); do
	  # Identify the forward and reverse fastq files.

	  READ1="${PARENT_DIR}/${FILE1[i]}"

	  BASE1="${FILE1[i]%.*}"


	  mkdir "${DEMULT_DIR}"/"${FILE1[i]}"
		mkdir "${FINAL_DIR}"/"${FILE1[i]}"



		echo "Working on input file $[i+1] out of ${#FILE1[@]}"

	##First cutdapt:


	# Look for the fwd primer and if found on rc, reverse the output


	cutadapt -g file:"${OUTPUT_DIR}"/direction.fasta -o "${READ1}".new.fastq --quiet --action=none --rc -e 0.3 "${READ1}"

  echo "Number of lines to begin with"

  wc -l "${READ1}"

  echo "Number of lines after"

  wc -l "${READ1}".new.fastq



	cutadapt -g file:"${OUTPUT_DIR}"/barcodes_P5.fasta -o "${DEMULT_DIR}"/"${FILE1[i]}"/{name}_round1.fastq \
	 --quiet --untrimmed-output "${OUTPUT_DIR}"/${BASE1}_nop5.fastq -e 0.2 "${READ1}".new.fastq




		n_files=("${DEMULT_DIR}"/"${FILE1[i]}"/*round1.fastq)



	 ls "${DEMULT_DIR}"/"${FILE1[i]}"/

	 for file in "${n_files[@]}"; do

		 echo "working on file " "${file}"
		 echo ""


	 BASE_OUTPUT=$(basename "${file}" |  sed 's/_round1.fastq//g')

	 echo "${BASE_OUTPUT}"

	 mkdir "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"
	 mkdir "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"
	 ## subset in awk



#awk -F',' -v COLNAME="${COLNUM_ID1}" -v VALUE="${BASE_OUTPUT}" \
#-v COLSEQ2="${COLNUM_ID2_RCSEQ}" -v COLID2="${COLNUM_ID2_WELL}" \
#'NR>1 { if ($COLNAME == VALUE) {printf ">Well_%s\n%s\n", $COLID2, $COLSEQ2} }'  "${SEQUENCING_METADATA}" > "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/barcodes.p7.fasta

#
#	 	cutadapt -a file:"${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/barcodes.p7.fasta --untrimmed-output "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"_nop7.fastq \
#	 -o "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/"${BASE_OUTPUT}"_{name}.fastq \
#	 "${file}"  --quiet -e 0.2 >> "${LOGFILE}"

seqkit seq -r -p -t DNA "${file}" -o "${file}".rc.fastq

awk -F',' -v COLNAME="${COLNUM_ID1}" -v VALUE="${BASE_OUTPUT}" \
-v COLSEQ2="${COLNUM_ID2_SEQ}" -v COLID2="${COLNUM_ID2_WELL}" \
'NR>1 { if ($COLNAME == VALUE) {printf ">Well_%s\n%s\n", $COLID2, $COLSEQ2} }'  "${SEQUENCING_METADATA}" > "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/barcodes.p7.rc.fasta


cutadapt -g file:"${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/barcodes.p7.rc.fasta --untrimmed-output "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"_nop7.fastq \
	 -o "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/"${BASE_OUTPUT}"_{name}.fastq \
	 "${file}".rc.fastq  --quiet -e 0.2 >> "${LOGFILE}"


### NOW FIND HOW MANY PRIMERS PER Plate-Well combo



	for file2 in "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/*.fastq; do

   BASE_P7=$(echo "${file2}" |  sed 's/.fastq//g' | sed 's/.*Well_//g' )

   echo "In plate ${BASE_OUTPUT} and Well  ${BASE_P7}"

   awk -F ',' -v ID1="${COLNUM_ID1}" -v ID2="${COLNUM_ID2_WELL}" \
   -v VALUE1="${BASE_OUTPUT}" -v VALUE2="${BASE_P7}" \
   -v SEQF=$COLNUM_PRIMER1 -v SEQR=$COLNUM_PRIMER2_RC -v LOCI=$COLNUM_LOCUS \
   'NR>1 { if ($ID1 == VALUE1 && $ID2 == VALUE2) {printf ">Locus_%s\n%s...%s\n", $LOCI, $SEQF, $SEQR} }'  "${SEQUENCING_METADATA}" > "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/pcr.fasta

#### And the final cutadapt
#### We will do the reverse primer first
awk -F ',' -v ID1="${COLNUM_ID1}" -v ID2="${COLNUM_ID2_WELL}" \
   -v VALUE1="${BASE_OUTPUT}" -v VALUE2="${BASE_P7}" \
   -v SEQF=$COLNUM_PRIMER1 -v SEQR=$COLNUM_PRIMER2 -v LOCI=$COLNUM_LOCUS \
   'NR>1 { if ($ID1 == VALUE1 && $ID2 == VALUE2) {printf ">Locus_%s\n%s\n", $LOCI, $SEQR} }'  "${SEQUENCING_METADATA}" > "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/pcr.rev.fasta


	 	cutadapt -g file:"${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/pcr.rev.fasta --untrimmed-output "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"_nopcr.fastq \
	 -o "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/"${BASE_OUTPUT}"_Well_"${BASE_P7}"_{name}.fastq \
	 "${file2}"  --quiet -e 0.2 >> "${LOGFILE}"
	 
### Now we would need to find, for each reverse primer its correspondent forward primer


    for file3 in "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/"${BASE_OUTPUT}"_Well_"${BASE_P7}"*.fastq; do  
 
 seqkit seq -r -p -t DNA "${file3}" -o "${file3}".temp.fastq  
 
 rm "${file3}"
 
 LOCI_NOW=$(basename "${file3}" |  sed 's/.fastq//g' | sed 's/.*Locus_//g' )
 
 awk -F ',' -v ID1="${COLNUM_ID1}" -v ID2="${COLNUM_ID2_WELL}" -v LOCI="${COLNUM_LOCUS}"  \
   -v VALUE1="${BASE_OUTPUT}" -v VALUE2="${BASE_P7}" -v VALUE3="${LOCI_NOW}" \
   -v SEQF=$COLNUM_PRIMER1 -v SEQR=$COLNUM_PRIMER2 \
   'NR>1 { if ($ID1 == VALUE1 && $ID2 == VALUE2 && $ID3 == VALUE3) {printf ">Locus_%s\n%s\n", $LOCI, $SEQF} }'  "${SEQUENCING_METADATA}" > "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/pcr.fwd.fasta

	 	cutadapt -g file:"${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/pcr.fwd.fasta --untrimmed-output "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"_nopcrfwd.fastq \
	 -o "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/"${BASE_OUTPUT}"_Well_"${BASE_P7}"_{name}_final.fastq \
	 "${file3}".temp.fastq  --quiet -e 0.2 >> "${LOGFILE}"
	 
	 

    
    done
	nseq_noprimer=$(cat ${file2} | wc -l)

	 printf "%s,%s,%s,%s\n" \
	 "${BASE_OUTPUT}" "${file2}" \
	  "${nseq_demult}" "${nseq_noprimer}" >> "${OUTPUT_SUMMARY}"

	done # End of loop across all wells within a plate

cd "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"



decona -l "${MIN_LENGTH}" -m "${MAX_LENGTH}" -q 10 -c "${CLUSTER_SIM}" -n 10 -M

cd "${Current_dir}"

done # End of loop across all plates within a file


	done # End of loop across all initial fastq files



# if [[ "${SEARCH_ASVs}" = "YES" ]]; then
	echo "launching  gather decona"
# 	conda activate decona
#
	Rscript --vanilla "${SCRIPT_DIR}"/r/gather.decona.r "${OUTPUT_DIR}"
# fi
