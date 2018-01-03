#!/usr/bin/env bash
#usage: bash find_barcodes.sh folder barcodes_file
#The goal is to transform the input .1 and .2 reads from Miseq output
# into split .1 and .2 by sample - assuming the first 6 characters
# are ok for splitting reads into samples - ie they are the secondary tags

my_folder="${1}"
out_dir="${my_dir}"_split

my_reverse_primer= bash revcom.sh
echo 'number of primer matches of Fwd Primer on Read1'
cat Lib-A_S1_L001_R1_001_sub.fastq | grep "GG[A-Z]AC[A-Z]GG[A-Z]TGAAC[A-Z]GT[A-Z]TA[A-Z]CC[A-Z]CC" --color -B1 --count

echo 'number of primer matches of Fwd Primer (RC) on Read1'
cat Lib-A_S1_L001_R1_001_sub.fastq | grep "GG[A-Z]GG[A-Z]TA[A-Z]AC[A-Z]GTTCA[A-Z]CC[A-Z]GT[A-Z]CC" --color -B1 --count

echo 'number of primer matches of Rev Primer on Read1'
cat Lib-A_S1_L001_R1_001_sub.fastq | grep "TA[A-Z]AC[A-Z]TC[A-Z]GG[A-Z]TG[A-Z]CC[A-Z]AA[A-Z]AA[A-Z]CA" --color -B1 --count

echo 'number of primer matches of Rev Primer (RC) on Read1'
cat Lib-A_S1_L001_R1_001_sub.fastq | grep "TG[A-Z]TT[A-Z]TT[A-Z]GG[A-Z]CA[A-Z]CC[A-Z]GA[A-Z]GT[A-Z]TA" --color -B1 --count

echo 'number of primer matches of Fwd Primer on Read2'
cat Lib-A_S1_L001_R2_001_sub.fastq | grep "GG[A-Z]AC[A-Z]GG[A-Z]TGAAC[A-Z]GT[A-Z]TA[A-Z]CC[A-Z]CC" --color -B1 --count

echo 'number of primer matches of Fwd Primer (RC) on Read2'
cat Lib-A_S1_L001_R2_001_sub.fastq | grep "GG[A-Z]GG[A-Z]TA[A-Z]AC[A-Z]GTTCA[A-Z]CC[A-Z]GT[A-Z]CC" --color -B1 --count

echo 'number of primer matches of Rev Primer on Read2'
cat Lib-A_S1_L001_R2_001_sub.fastq | grep "TA[A-Z]AC[A-Z]TC[A-Z]GG[A-Z]TG[A-Z]CC[A-Z]AA[A-Z]AA[A-Z]CA" --color -B1 --count

echo 'number of primer matches of Rev Primer (RC) on READ2'
cat Lib-A_S1_L001_R2_001_sub.fastq | grep "TG[A-Z]TT[A-Z]TT[A-Z]GG[A-Z]CA[A-Z]CC[A-Z]GA[A-Z]GT[A-Z]TA" --color -B1 --count
