#!/usr/bin/env bash


################################################################################
# INPUT
################################################################################
# What is the file path to the directory containing all of the libraries/reads?

PARENT_DIR="${MAIN_DIR}"/data

# Where is the sequencing metadata file? (SEE FORMATTING GUIDELINES IN README!)
SEQUENCING_METADATA="${PARENT_DIR}"/metadata.csv



################################################################################
# OUTPUT
################################################################################
# This script will generate a directory (folder) containing the output of the script.
# Where do you want this new folder to go?
OUTPUT_DIRECTORY="${HOME}"/fastqs_demultiplexed_for_DADA2 #"${PARENT_DIR%/*}"


################################################################################
# METADATA DETAILS
################################################################################
# Specify columns for raw sequencing files:
COLNAME_FILE1="file1"
COLNAME_FILE2="file2"

# MUST be unique for each row!
COLNAME_SAMPLE_ID="sample_id"


# Your metadata must have a column corresponding to the subfolders containing the raw reads.
# In order to make this flexible across both multiple and single library preps, you must include this even if you only sequenced one library (sorry!).
COLNAME_ID1_NAME="pri_index_name"
COLNAME_ID1_SEQ="pri_index_seq"

COLNAME_INSERT_SIZE="insert_size"

LENGTH_FRAG="385"


# if "NO", provide the following values for PEAR:
minimum_overlap="10" # [10]
assembled_max="10000" # [1000]
assembled_min="50" # [50]

################################################################################
# QUALITY FILTERING
################################################################################
# Substantial quality filtering (e.g. trimming, minimum length, etc) is performed by PEAR during read merging.
# You may also want to exclude sequences containing more than a specified threshold of 'expected errors'
# This number is equal to the sum of the error probabilities.
# For more information on this parameter, Google the usearch help
Perform_Expected_Error_Filter="YES" # [YES|NO]
Max_Expected_Errors="0.5"


################################################################################
# DEMULTIPLEXING
################################################################################

# Do the reads contain index sequences which identifies their sample of origin?
SECONDARY_INDEX="YES"

# Specify the nucleotide sequences that differentiate multiplexed samples
# (sometimes, confusingly referred to as "tags" or "barcodes")
# these are the secondary index -- the primary index added with the sequencing adapters should not be in the sequence data
# You can grab these from the file specified above (SEQUENCING_METADATA) by specifying the column name of index sequences.
COLNAME_ID2_SEQ="sec_index_seq"

# How many nucleotides pad the 5' end of the tag sequence?
# TODO build in flexibility (this number is unused right now)
TAG_Ns="3"
SECONDARY_INDEX_START="4"
COLNAME_ID2_START="sec_index_start"

################################################################################
# PRIMER REMOVAL
################################################################################
# Specify the primers used to generate these amplicons.
# As with the multiplex indexes, Banzai will grab these from the file SEQUENCING_METADATA.
# You must indicate the column names of the forward and reverse primers
COLNAME_PRIMER1="primerF_seq"
COLNAME_PRIMER2="primerR_seq"

################################################################################
# USE HASH
################################################################################
# Should the sequence ID after dereplication be the output of a hash algorithm?

USE_HASH="YES"

################################################################################
# CLUSTER OTUs: USING DADA2
################################################################################

SEARCH_ASVs="YES"

## TODO: Add variables to control the behaviour of dada2



################################################################################
# REANALYSIS
################################################################################
# Would you like to pick up where a previous analysis left off?

# Have the reads already been paired?
ALREADY_PEARED="NO" # YES/NO
PEAR_OUTPUT='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141030_2020/1_merged.assembled.fastq.gz'

# Have the merged reads been quality filtered?
ALREADY_FILTERED="NO" # [YES|NO]
FILTERED_OUTPUT='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141030_2020/2_filtered_renamed.fasta'


# If using ASVs, have you already demultiplexed your reads into .1 and .2 pairs per sample.
# Point towards the output folder (must include files: sample_trans.tmp,
# barcodes.fasta, summary.csv and pcr_primers.fasta; and the folder /demultiplexed
# so the pipeline can cp all necessary files

ALREADY_DEMULTIPLEXED="NO"
DEMULT_OUTPUT=""


################################################################################
# CONTINUING ANALYSIS
################################################################################
#Would you like to add this analysis to a previous set of samples already processed?

# You should provide a csv file with all sequences and sh1 hashes, and a csv with the previous abundance data
# It will add the new sequences and hashes to the first file, and the new samples and their sequence abundance to
# the second file. You can choose to overwrite or not the input files with the new output

ADD_TO_PREVIOUS="NO"
FORMER_HASH="/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180213_2333/hash_key.csv"
FORMER_ABUNDANCE="/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180213_2333/ASV_table.csv"
LOG_FILE=""


################################################################################
# GENERAL SETTINGS
################################################################################
# Would you like to save every single intermediate file as we go? YES | NO
# recommendation: NO, unless testing or troubleshooting
HOARD="YES"

# Would you like to compress extraneous intermediate files once the analysis is finished? YES/NO
PERFORM_CLEANUP="YES"
