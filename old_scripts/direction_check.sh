#!/usr/bin/env bash

# A shell script that checks that the sequences in the output
# of demultiplex_both_fastq.sh are in the right direction
# Execute this as bash direction_check.sh <output.folder>

# Output folder must include seqnames.txt and banzai_params.sh
# Original fastq files must remain in the same location as when the pipeline
# was run on the first place
MAIN_DIR="$(dirname "$0")"
SCRIPT_DIR="${MAIN_DIR}"/../scripts
out_folder=${1}
LOGFILE="${out_folder}"/logfile_direction_report.txt
exec > >(tee "${LOGFILE}") 2>&1

if [[ -s "${out_folder}" ]]; then
  echo
	echo "Analyzing direction of sequences from ${out_folder}"
  echo
else
	echo 'ERROR! Could not find output folder. You specified the file path:'
	echo
	echo "${out_folder}"
	echo
	echo 'That folder does not exist. Aborting script.'
	exit
fi


params_file=${1}/banzai_params.sh

if [[ -s "${params_file}" ]]; then
	echo "Parameters file from ${params_file}"
  echo
else
	echo 'ERROR! Could not find parameter file. You specified the file path:'
	echo
	echo "${out_folder}"
	echo
	echo 'That folder does contain the params file. Aborting script.'
	exit
fi

seqs=${1}/seqnames.txt

if [[ -s "${seqs}" ]]; then
	echo "We are checking sequences from ${seqs}"
  echo
else
	echo 'ERROR! Could not find seqnames file. You specified the file path:'
	echo
  echo "${out_folder}"
	echo
	echo 'That folder does contain the seqnames file. Aborting script.'
	exit
fi
# Source parameters file load metadata and obtain the first pair of fastq files

source "${params_file}"

# Get to the metadata - try both the original and if not the output
if [[ -s "${SEQUENCING_METADATA}" ]]; then
	echo "Reading metadata from:"
	echo "${SEQUENCING_METADATA}"
else
  SEQUENCING_METADATA=${1}/metadata.csv
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
fi
# Fixing metadata endings if needed

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

source "${SCRIPT_DIR}"/get_colnum.sh

# Get the fastqs
COLNUM_FILE1=$( get_colnum "${COLNAME_FILE1}" "${SEQUENCING_METADATA}")
COLNUM_FILE2=$( get_colnum "${COLNAME_FILE2}" "${SEQUENCING_METADATA}")

COLNUM_PRIMER1=$( get_colnum "${COLNAME_PRIMER1}" "${SEQUENCING_METADATA}")
COLNUM_PRIMER2=$( get_colnum "${COLNAME_PRIMER2}" "${SEQUENCING_METADATA}")


all_columns=( COLNUM_FILE1 COLNUM_FILE2 COLNUM_PRIMER1 COLNUM_PRIMER2)

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
#PCR PRIMERS
PRIMER1=($(awk -F',' -v COLNUM=$COLNUM_PRIMER1 \
  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA |\
  sort | uniq ))

PRIMER2=($(awk -F',' -v COLNUM=$COLNUM_PRIMER2 \
  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA |\
  sort | uniq ))
#FASTQ FILES
FILE1=($(awk -F',' -v COLNUM=$COLNUM_FILE1 \
  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
  sort | uniq))

FILE2=($(awk -F',' -v COLNUM=$COLNUM_FILE2 \
  'NR>1 {print $COLNUM}' $SEQUENCING_METADATA |\
  sort | uniq ))
  READ1="${PARENT_DIR}/${FILE1[0]}"
  READ2="${PARENT_DIR}/${FILE2[0]}"

echo "Looking for sequences in ${seqs} in the original files :"
echo
echo "${READ1}"
echo "${READ2}"
echo
echo "Using PCR Primers :"
echo " Fwd: ${PRIMER1}"
echo
echo " Rev: ${PRIMER2}"

## Now we can start with the process itself

F1_found="${out_folder}"/matchesF.txt
F2_found="${out_folder}"/matchesR.txt

#READ1=$(head -n 50000 "${READ1}")
#READ2=$(head -n 50000 "${READ2}")

COIF=$(echo "${PRIMER1}" | tail -1 | sed 's/[^ACGT]/\[A-Z]/g')
COIR=$(echo "${PRIMER2}" | tail -1 | sed 's/[^ACGT]/\[A-Z]/g')

echo "${COIF}"


cat "${seqs}" | cut -c 1-120 | xargs -I '{}' sed -n 's/{}.*$//p' "${READ1}" | awk 'length($0)>15' > "${F1_found}"

cat "${seqs}" | cut -c 1-120 | xargs -I '{}' sed -n 's/{}.*$//p' "${READ2}" | awk 'length($0)>15' > "${F2_found}"

#This looks in those files for the Fwd and Rev PCR primers - we only should get the fwd primer

COIF=$(echo "${PRIMER1}" | tail -1 | sed 's/[^ACGT]/\[A-Z]/g')
COIR=$(echo "${PRIMER2}" | tail -1 | sed 's/[^ACGT]/\[A-Z]/g')

echo "${COIF}"

matchesF1=$(grep "${COIF}" "${F1_found}" --count)
matchesRRC1=$(grep "${COIR}" "${F1_found}" --count)
matchesF2=$(grep "${COIF}" "${F2_found}" --count)
matchesRRC2=$(grep "${COIR}" "${F2_found}" --count)

echo "${matchesF1}" "${matchesRRC1}" "${matchesF2}" "${matchesRRC2}"

printf "Fwd_Primer\tfile.1\t%s\nRev_Primer\tfile.1\t%s\nFwd_Primer\tfile.2\t%s\nRev_Primer\tfile.2\t%s\n" \
 "${matchesF1}" "${matchesRRC1}" "${matchesF2}" "${matchesRRC2}" > "${out_folder}"/matches.txt


Rscript "${MAIN_DIR}"/direction_check.R "${out_folder}" "${MAIN_DIR}"
