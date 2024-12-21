#!/usr/bin/env bash
module load cutadapt/4.1

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
NOPRIMERS_DIR="${OUTPUT_DIR}"/noprimers
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
COLNUM_LOCUS=$(get_colnum "${COLNAME_LOCUS}" "${SEQUENCING_METADATA}")
# Run away from the script if any of the previous columns was not found

all_columns=( COLNUM_FILE1 COLNUM_FILE2 COLNUM_ID1 COLNUM_ID2 \
COLNUM_ID2_START COLNUM_SAMPLE COLNUM_PRIMER1 COLNUM_PRIMER2 COLNUM_LOCUS)
#TODO:I am not using colnumID2 START
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
## Modify it so it can work with different loci in the same run
## First get_the loci

  LOCUS=($(awk -F',' -v COLNUM=$COLNUM_LOCUS \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq ))
## For each locus, get the primers and add them to two fasta files: fwd and rev. 
## We should do this for each library, so there is no interference between libraries and projects sharing a run

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

# 
# 	ID1_ALL=($(awk -F',' -v COLNUM=$COLNUM_ID1 \
# 	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))
# 	ID1S=($(awk -F',' -v COLNUM=$COLNUM_ID1 \
# 	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}"  |\
# 			sort | uniq))
# 	ID2_ALL=($(awk -F',' -v COLNUM=$COLNUM_ID2 \
# 	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))
# 	ID2_ALL_RC=($( for i in "${ID2_ALL[@]}"; do revcom $i; done))
# 
# # write file for translating demultiplexed output to samples
	SAMPLE_TRANS_FILE="${OUTPUT_DIR}"/sample_trans.tmp
# 	for (( i=0; i < "${#ID2_ALL[@]}"; i++ )); do
# 	  
# 	done
# 	for (( i=0; i < "${#ID1S[@]}"; i++ )); do
# 	  printf "File1:%s\tFile2:%s\tLib:%s\n" \
# 	  "${FILE1[i]}" "${FILE2[i]}" "${ID1S[i]}"
# 
# 
# 	done

# 
# #Create the fasta file of the barcodes
# 
# 	Barcodes_file="$OUTPUT_DIR"/barcodes.fasta
# 	for (( i=0; i < "${#ID2S[@]}"; i++ )); do
# 	  printf ">%s\n^%s\n" \
# 		"${ID2S[i]}" "${ID2S[i]}" >> "${Barcodes_file}"
# 	done
# 
# 	primers_file="${OUTPUT_DIR}"/pcr_primers.fasta
# 
# 	printf ">FWD\n${PRIMER1}\n>REV\n${PRIMER2}\n" > "${primers_file}"
# 
# 	source "${SCRIPT_DIR}"/functions/check_primers.sh "${primers_file}"

#Hooray it works
#Create a dir for all the demultiplexed files

# now we have to remove all hard-coded stuff and link it to
#banzai_params
# to get the .1 files trimmed and the .2 selected along

	OUTPUT_SUMMARY="${OUTPUT_DIR}/summary.csv"
	printf "library_sample,loci,step,nReads\n" \
	> "${OUTPUT_SUMMARY}"

################################################################################
# BEGIN LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
################################################################################

	for (( i=0; i < "${#FILE1[@]}"; i++ )); do
	  # Identify the forward and reverse fastq files.

	  READ1="${PARENT_DIR}/${FILE1[i]}"
	  READ2="${PARENT_DIR}/${FILE2[i]}"
	  
	  # Subset here to use the subsetting related to the file, and not dependent on the order 
    # of lib names: do this for ID, barcodes and primers
    
    ID1S=$( awk -F',' -v COLNUM=$COLNUM_FILE1 -v VALUE=${FILE1[i]} -v ID1=$COLNUM_ID1 \
	    ' {if ($COLNUM == VALUE) { print  $ID1 }} ' $SEQUENCING_METADATA | uniq)
	
	   echo ${ID1S}
	   
	  # Barcodes
	  
	  Barcodes_file="$OUTPUT_DIR"/barcodes_"${ID1S}".fasta
	  
	  awk -F',' -v COLNUM=$COLNUM_FILE1 -v VALUE=${FILE1[i]} -v ADAP=$COLNUM_ID2 \
	    '{if ($COLNUM == VALUE) { printf ">%s\n^%s\n", $ADAP, $ADAP } }' $SEQUENCING_METADATA > "${Barcodes_file}"
	  
	  # Primers and loci
	  primers_file_R1="$OUTPUT_DIR"/primers_"${ID1S}"_R1.fasta
    primers_file_R2="$OUTPUT_DIR"/primers_"${ID1S}"_R2.fasta
 
    awk -F',' -v COLNUM=$COLNUM_FILE1 -v VALUE=${FILE1[i]} -v LOCUS=$COLNUM_LOCUS \
      -v FWD=$COLNUM_PRIMER1 -v REV=$COLNUM_PRIMER2 \
	    '{if ($COLNUM == VALUE) { print $LOCUS,$FWD,$REV } }' $SEQUENCING_METADATA | sort|uniq > "${OUTPUT_DIR}"/unique_input.txt
	    
	  awk -v lib="${primers_file_R1}" '{
      file = lib ;
      fwd_header = ">Locus_" $1 "_Fwd";
      fwd_sequence = $2;
      rev_header = ">Locus_" $1 "_Rev";
      rev_sequence = $3;
      print fwd_header "\n" fwd_sequence "\n" rev_header "\n" rev_sequence >> file }' "${OUTPUT_DIR}"/unique_input.txt
      
    awk -v lib="${primers_file_R2}" '{
      file = lib ;
      fwd_header = ">Locus_" $1 "_Rev";
      fwd_sequence = $3;
      rev_header = ">Locus_" $1 "_Fwd";
      rev_sequence = $2;
      print fwd_header "\n" fwd_sequence "\n" rev_header "\n" rev_sequence >> file }' "${OUTPUT_DIR}"/unique_input.txt
      
    # Sample map  
	  
	  awk -F',' -v COLNUM=$COLNUM_FILE1 -v VALUE=${FILE1[i]} -v ID1=$COLNUM_ID1 \
	    -v ID2=$COLNUM_ID2 -v SAMPLE_NAME=$COLNUM_SAMPLE \
	    ' {if ($COLNUM == VALUE) { printf  "ID1=%s;ID2=%s\t%s_%s\t%s\n", $ID1, $ID2, $ID1, $ID2, $SAMPLE_NAME }} ' $SEQUENCING_METADATA >> "${SAMPLE_TRANS_FILE}"

	  mkdir "${OUTPUT_DIR}"/"${ID1S}"


		mkdir "${OUTPUT_DIR}"/cleaned/"${ID1S}"

		echo "Working on Library $[i+1] out of ${#FILE1[@]}"

	##First cutdapt:
	#TODO: use only the number of barcodes used for this Library
	
 # Anchoring the adapters seems like the only option 
	
 # Only one round of cutadapt is needed for demultiplexing
 
 cutadapt -g "file:"${Barcodes_file}";min_overlap=8" \
	  -G "file:"${Barcodes_file}";min_overlap=8" \
  	-o "${OUTPUT_DIR}"/"${ID1S}"/"${ID1S}"_{name}.R1.fastq \
	  -p "${OUTPUT_DIR}"/"${ID1S}"/"${ID1S}"_{name}.R2.fastq \
	  "${READ1}" "${READ2}" --discard-untrimmed -j 0 -e 1 --pair-adapters > "${OUTPUT_DIR}"/cutadapt_logfile.txt
	  
	 ## Now process the logfile to get the summary info: 
	
  if grep -A 2 '^=== \(First\|Second\) read: Adapter' "${OUTPUT_DIR}"/cutadapt_logfile.txt > "${OUTPUT_DIR}"/temp_log.txt; then
        awk -v Library="${ID1S}" '
        /^=== (First|Second) read: Adapter/ { 
            split($0, a, " "); 
            read=a[2]; 
        }
        /^Sequence:/ { 
            split($0, a, " "); 
            adapter_name=a[2]; 
            gsub(/;$/, "", adapter_name); 
            times=a[length(a)-1]; 
            gsub(/ times$/, "", times); 
            print Library "_" adapter_name",all_loci,demult_" read "," times;
        }' "${OUTPUT_DIR}"/temp_log.txt >> "${OUTPUT_SUMMARY}"
  else
        echo "iteration $IDS,Error,Error,Error" >> "${OUTPUT_SUMMARY}"
  fi



	n_files=("${OUTPUT_DIR}"/"${ID1S}"/*.R2.fastq)
		
		
		i_count=0

    for r2file in "${n_files[@]}"; do
	 
        # We loop through all .2 files
    		i_count=$((i_count+1))
    
    		short_r2file=$(basename "${r2file}"| sed 's/.R2.fastq$//')
    
    		r1file=$(echo ${r2file} | sed 's/.R2.fastq$/.R1.fastq/g' )
    		short_r1file=$(basename "${r1file}"| sed 's/.R1.fastq$//') 
    	 
    
    		#New messages so it's easier to see the progress of the script
    
    		echo -ne "Working on sample ${i_count} of ${#n_files[@]}"'\r'
    
    	cutadapt -g file:"${primers_file_R1}" -G file:"${primers_file_R2}" --discard-untrimmed \
    	 -o "${OUTPUT_DIR}"/cleaned/"${ID1S}"/"${short_r1file}"_{name}.R1.fastq \
    	 -p "${OUTPUT_DIR}"/cleaned/"${ID1S}"/"${short_r2file}"_{name}.R2.fastq \
    	 -j 0 "${r1file}" "${r2file}" --pair-adapters 2 > "${OUTPUT_DIR}"/cutadapt_logfile.txt
    	  
    
    ## Now process the logfile to get the summary info: 
    	
     if grep -A 2 '^=== \(First\|Second\) read: Adapter' "${OUTPUT_DIR}"/cutadapt_logfile.txt > "${OUTPUT_DIR}"/temp_log.txt; then
            awk  Sample="${short_r1file}" '
            /^=== (First|Second) read: Adapter/ { 
                split($0, a, " "); 
                read=a[2]; 
            }
            /^Sequence:/ { 
                split($0, a, " "); 
                primer_name=a[2]; 
                gsub(/;$/, "", adapter_name); 
                times=a[length(a)-1]; 
                gsub(/ times$/, "", times); 
                print Sample "," primer_name", demult_" read "," times;
            }' "${OUTPUT_DIR}"/temp_log.txt >> "${OUTPUT_SUMMARY}"
        else
            echo "iteration $IDS,Error,Error,Error" >> "${OUTPUT_SUMMARY}"
        fi
        
        mv "${r1file}" "${r2file}" "${DEMULT_DIR}"
        mv "${OUTPUT_DIR}"/cleaned/${ID1S}/* "${NOPRIMERS_DIR}"
    
    done # This finishes the for loop for all demulted files, getting the primers out
	  rm -r "${OUTPUT_DIR}"/"${ID1S}"

	done # This finishes the for loop for all libraries

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

module rm cutadapt/4.1

if [[ "${SEARCH_ASVs}" = "YES" ]]; then
	echo "This is read1 ${READ1}"
	module load R/4.3.1
	Rscript "${SCRIPT_DIR}"/r/code_dada2_cluster.r "${OUTPUT_DIR}" "${DEMULT_DIR}" "${SCRIPT_DIR}" "${USE_HASH}" "${READ1}" "${READ2}"\
	"${ADD_TO_PREVIOUS}" "${FORMER_HASH}" "${FORMER_ABUNDANCE}" "${LOG_FILE}"
fi
