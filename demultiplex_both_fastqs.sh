#!/usr/bin/env bash
module load cutadapt/4.1

# Usage bash demultiplex_both_fastqs.sh banzai_params.sh
#This script is built using banzai (github.com/jimmyodonnell/banzai) as template

#We need to gather: Location of functions  and fastqs:
MAIN_DIR="$(dirname "$0")"
SCRIPT_DIR="${MAIN_DIR}"/scripts
for file in "${SCRIPT_DIR}"/*.sh ; do
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
	mkdir "${NOPRIMERS_DIR}"
	BARCODES_DIR="${OUTPUT_DIR}"/barcodes_and_primers
	mkdir "${BARCODES_DIR}"
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

	# Sample names
	COLNUM_SAMPLE=$( get_colnum "${COLNAME_SAMPLE_ID}" "${SEQUENCING_METADATA}")

	# Primers
	COLNUM_PRIMER1=$( get_colnum "${COLNAME_PRIMER1}" "${SEQUENCING_METADATA}")
	COLNUM_PRIMER2=$( get_colnum "${COLNAME_PRIMER2}" "${SEQUENCING_METADATA}")
	COLNUM_LOCUS=$(get_colnum "${COLNAME_LOCUS}" "${SEQUENCING_METADATA}")
	# Run away from the script if any of the previous columns was not found

	all_columns=( COLNUM_FILE1 COLNUM_FILE2 COLNUM_ID1 COLNUM_ID2 \
	COLNUM_SAMPLE COLNUM_PRIMER1 COLNUM_PRIMER2 COLNUM_LOCUS)
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



 # Read in primers
 ## First get_the loci

	LOCUS=($(awk -F',' -v COLNUM=$COLNUM_LOCUS \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq ))
 ## For each locus, get the primers. We would create a primer file per library

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

 #Unique samples are given by combining the primary and secondary indexes

	ID_COMBO=$( awk -F',' -v COLNUM1=$COLNUM_ID1 -v COLNUM2=$COLNUM_ID2 \
	'NR>1 {
	  print ";ID1=" $COLNUM1 ";ID2=" $COLNUM2
	}' "${SEQUENCING_METADATA}" )

	SAMPLE_NAMES=($(awk -F',' -v COLNUM=$COLNUM_SAMPLE \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))

 
 # Check that sample names are not repeated
 
	NSAMPLES="${#SAMPLE_NAMES[@]}"

 # Now calculate the number of unique sample names and check uniqueness
	UNIQ_SAMPLES=( $(echo "${SAMPLE_NAMES[@]}" | tr ' ' '\n' | sort -u))
	N_UNIQ_SAMPLES="${#UNIQ_SAMPLES[@]}"


	if [[ "${NSAMPLES}" != "${N_UNIQ_SAMPLES}" ]]; then
		echo " At least one sample name is repeated "
		echo " I am not angry, just dissapointed. Exiting script"
		exit
	fi

 
 # write file for translating demultiplexed output to samples
	SAMPLE_TRANS_FILE="${OUTPUT_DIR}"/sample_trans.tmp
 # write summary file for translating demultiplexed output to samples
	OUTPUT_SUMMARY="${OUTPUT_DIR}/summary.csv"
	printf "fastq_header,locus,step,nReads\n" \
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
	    '{if ($COLNUM == VALUE) { printf ">%s\n%s\n", $ADAP, $ADAP } }' $SEQUENCING_METADATA > "${Barcodes_file}"
	  
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
		
		# Only one round of cutadapt is needed for demultiplexing, use all cores available
	
		cutadapt -g "file:"${Barcodes_file}";min_overlap=8" \
			-G "file:"${Barcodes_file}";min_overlap=8" \
			-o "${OUTPUT_DIR}"/"${ID1S}"/"${ID1S}"_{name}.R1.fastq \
			-p "${OUTPUT_DIR}"/"${ID1S}"/"${ID1S}"_{name}.R2.fastq \
			"${READ1}" "${READ2}" --discard-untrimmed -j 0 -e 1 --pair-adapters > "${OUTPUT_DIR}"/cutadapt_logfiledemult.txt
	  
		if grep -A 2 '^=== First read: Adapter' "${OUTPUT_DIR}"/cutadapt_logfiledemult.txt > "${OUTPUT_DIR}"/temp_log.txt; then
					awk -v Library="${ID1S}" '
					/^=== First read: Adapter/ { 
					split($0, a, " "); 
					adapter_name=a[5];

					}
					/^Sequence:/ { 
						split($0, a, " "); 
						times=a[length(a)-1]; 
						gsub(/ times$/, "", times); 
						print Library "_" adapter_name",all_loci,demultiplexing," times;
					}' "${OUTPUT_DIR}"/temp_log.txt >> "${OUTPUT_SUMMARY}"
		else
					echo "Library $IDS,Error,Error,Error" >> "${OUTPUT_SUMMARY}"
		fi



		n_files=("${OUTPUT_DIR}"/"${ID1S}"/*.R2.fastq)
		
		
		i_count=0
	  #Remove primers
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
			-j 0 "${r1file}" "${r2file}" --pair-adapters > "${OUTPUT_DIR}"/cutadapt_logfile.txt
    	  
    
			## Now process the logfile to get the summary info: 
			
			if grep -A 2 '^=== First read: Adapter' "${OUTPUT_DIR}"/cutadapt_logfile.txt > "${OUTPUT_DIR}"/temp_log.txt; then
				awk  -v Sample="${short_r1file}" '
				/^=== First read: Adapter/ { 
					split($0, a, " "); 
					primer_name=a[5];

				}
				/^Sequence:/ { 
					split($0, a, " "); 
					
					gsub(/;$/, "", adapter_name); 
					times=a[length(a)-1]; 
					gsub(/ times$/, "", times); 
					print Sample "," primer_name",noprimer_" read "," times;
				}' "${OUTPUT_DIR}"/temp_log.txt >> "${OUTPUT_SUMMARY}"
			else
				echo "iteration $IDS,Error,Error,Error" >> "${OUTPUT_SUMMARY}"
			fi
			
			mv "${r1file}" "${r2file}" "${DEMULT_DIR}"
			mv "${OUTPUT_DIR}"/cleaned/${ID1S}/* "${NOPRIMERS_DIR}"
		
		done # This finishes the for loop for all demulted files, getting the primers out
		rm -r "${OUTPUT_DIR}"/"${ID1S}"
		mv "${Barcodes_file}" "${BARCODES_DIR}"
		mv "${primers_file_R1}" "${BARCODES_DIR}"
		mv "${primers_file_R2}" "${BARCODES_DIR}"
	done # This finishes the for loop for all libraries

	rm -rf "${OUTPUT_DIR}"/cleaned
	rm  "${OUTPUT_DIR}"/cutadapt_logfile.txt
	rm  "${OUTPUT_DIR}"/cutadapt_logfiledemult.txt
	rm  "${OUTPUT_DIR}"/temp_log.txt
	rm  "${OUTPUT_DIR}"/unique_input.txt

else #In case you already demultiplexed your samples, then cp the files you need
	cp "${DEMULT_OUTPUT}"/sample_trans.tmp "${OUTPUT_DIR}"
	cp "${DEMULT_OUTPUT}"/summary.csv "${OUTPUT_DIR}"
	
	NOPRIMERS_DIR="${DEMULT_OUTPUT}"/noprimers

fi #This finishes the control flow in case you already demultiplexed

if [[ "${SEARCH_ASVs}" = "YES" ]]; then
	module load R/4.3.1
	Rscript "${SCRIPT_DIR}"/r/code_dada2_cluster.r "${OUTPUT_DIR}" "${NOPRIMERS_DIR}" "${USE_HASH}"  "${LENR1}" "${LENR2}"
fi
if [[ "${SEARCH_Unoise}" = "YES" ]]; then
	echo "Entering SEARCH_Uniose block"
	bash "${SCRIPT_DIR}"/clustering/vsearch_clustering.sh "${OUTPUT_DIR}" "${NOPRIMERS_DIR}" "${USE_HASH}"  "${LENR1}" "${LENR2}"
	fi

if [[ "${SECONDARY_SWARM}" = "YES" ]]; then
    echo "Preparing data for swarm..."
    Rscript "${SCRIPT_DIR}"/clustering/prepare_for_swarm.R "${OUTPUT_DIR}" 
   
    # Check if Rscript was successful
    if [[ $? -eq 0 ]]; then
        echo "Launching swarm..."
        bash "${SCRIPT_DIR}"/clustering/swarm.sh "${OUTPUT_DIR}"/swarm_input 
    else
        echo "Error: R script failed. Swarm will not be launched." >&2
        exit 1  # Exit with error status
    fi
	# Check if swarm was successful
    if [[ $? -eq 0 ]]; then
        echo "Parsing swarm..."
        Rscript "${SCRIPT_DIR}"/clustering/parsing_swarm.R "${OUTPUT_DIR}"
    else
        echo "Error: swarm script failed. " >&2
        exit 1  # Exit with error status
    fi

fi

if [[ "${HOARD}" = "NO" ]]; then

	rm -r "${OUTPUT_DIR}"/demultiplexed
	rm -r "${OUTPUT_DIR}"/noprimers
	rm -r "${OUTPUT_DIR}"/midfiles
	if [[  -d "${OUTPUT_DIR}"/swarm_input ]]; then
	rm -r "${OUTPUT_DIR}"/swarm_input
	fi
fi