#!/usr/bin/env bash


################################################################################
# INPUT
################################################################################
# What is the file path to the directory containing all of the libraries/reads?

PARENT_DIR="/home/meg/rgallego/Projects/Porcupine/raw"

# Where is the sequencing metadata file? (SEE FORMATTING GUIDELINES IN README!)
SEQUENCING_METADATA="/home/meg/rgallego/Projects/Porcupine/metadata.csv"



################################################################################
# OUTPUT
################################################################################
# This script will generate a directory (folder) containing the output of the script.
# Where do you want this new folder to go?
OUTPUT_DIRECTORY="/home/meg/rgallego/Projects/Porcupine/pipeline_output"

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

################################################################################
# QUALITY FILTERING
################################################################################
# Substantial quality filtering (e.g. trimming, minimum length, etc) is performed by DADA2
# Which also affects the length of the resulting amplicons: 
# TODO: make it a parameter that you can pass onto the R script

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
COLNAME_ID2_NAME="sec_index_seq"

################################################################################
# PRIMER REMOVAL
################################################################################
# Specify the primers used to generate these amplicons.
# As with the multiplex indexes, Banzai will grab these from the file SEQUENCING_METADATA.
# You must indicate the column names of the forward and reverse primers
COLNAME_PRIMER1="primerF_seq"
COLNAME_PRIMER2="primerR_seq"
COLNAME_LOCUS="locus"
################################################################################
# USE HASH
################################################################################
# Should the sequence ID after dereplication be the output of a hash algorithm?

USE_HASH="YES"

################################################################################
# CLUSTER OTUs: USING DADA2, Unoise and secondary swarm
################################################################################

SEARCH_ASVs="NO"

SEARCH_Unoise="YES"

SECONDARY_SWARM="YES"

## TODO: Add variables to control the behaviour of dada2
LENR1=240

LENR2=140


################################################################################
# REANALYSIS
################################################################################
# Would you like to pick up where a previous analysis left off?
ALREADY_DEMULTIPLEXED="NO"
DEMULT_OUTPUT="/home/meg/rgallego/Projects/Porcupine/pipeline_output/demultiplexed_20250303_1409"

########################################################################
# GENERAL SETTINGS
################################################################################
# Would you like to save every single intermediate file as we go? YES | NO
# recommendation: NO, unless testing or troubleshooting
HOARD="YES"

# Would you like to compress extraneous intermediate files once the analysis is finished? YES/NO
PERFORM_CLEANUP="YES"
