#!/usr/bin/env bash

#The goal is to check that F.1.fastq files contain
#only Fwd reads, and that  F.2 only reverse READS

# Step 1 - get only the first 100 characters of a F1.fastq file

DEMULTF1=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/demulF1.fastq
DEMULTF2=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/demulF2.fastq
DEMULTR1=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/demulR1.fastq
DEMULTR2=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/demulR2.fastq

FILTEREDF1=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/filteredF1.fastq
FILTEREDF2=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/filteredF2.fastq
FILTEREDR1=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/filteredR1.fastq
FILTEREDR2=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/filteredR2.fastq

ORIGINAL1=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/all_lib.1.fastq
ORIGINAL2=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/all_lib.2.fastq

COIF="GG[A-Z]AC[A-Z]GG[A-Z]TGAAC[A-Z]GT[A-Z]TA[A-Z]CC[A-Z]CC"
COIF_RC="GG[A-Z]GG[A-Z]TA[A-Z]AC[A-Z]GTTCA[A-Z]CC[A-Z]GT[A-Z]CC"
COIR="TA[A-Z]AC[A-Z]TC[A-Z]GG[A-Z]TG[A-Z]CC[A-Z]AA[A-Z]AA[A-Z]CA"
COIR_RC="TG[A-Z]TT[A-Z]TT[A-Z]GG[A-Z]CA[A-Z]CC[A-Z]GA[A-Z]GT[A-Z]TA"

F1_found=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/temp1.txt
F2_found=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/temp2.txt
R1_found=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/temp3.txt
R2_found=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/temp4.txt
#sed -n '2~4p' "${DEMULTF1}" | sort | uniq | xargs -i sed -n 's/{}.*$//p' "${ORIGINAL1}" > "${F1_found}"

#sed -n '2~4p' "${DEMULTF2}" | sort | uniq | xargs -i sed -n 's/{}.*$//p' "${ORIGINAL2}" > "${F2_found}"

#sed -n '2~4p' "${DEMULTR1}" | sort | uniq | xargs -i sed -n 's/{}.*$//p' "${ORIGINAL1}" > "${R1_found}"

#sed -n '2~4p' "${DEMULTR2}" | sort | uniq | xargs -i sed -n 's/{}.*$//p' "${ORIGINAL2}" > "${R2_found}"
head "${F1_found}"
array=("${F1_found}" "${F2_found}" "${R1_found}" "${R2_found}")
echo "${array[@]}"
output_file=/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180507_2201/output.txt
for file in ${array[@]}; do

matchesF=$(grep "${COIF}" $file --count)
matchesRC_F=$(grep "${COIF_RC}" $file --count)
matchesR=$(grep "${COIR}" $file --count)
matchesRC_R=$(grep "${COIR_RC}" $file --count)
printf "File1:%s\tmatchesF:%s\tmatchesRC_F:%s\tmatchesR:%s\tmatchesRC_R:%s\n" \
"$file" "${matchesF}" "${matchesRC_F}" "${matchesR}" "${matchesRC_R}" >>"${output_file}"
done
##By looking at this one, we know that all our demultiplexed files are in the right direction.
## Now we should look at the filtered files, the first step in the DADA2 pipeline
sed -n '2~4p' "${FILTEREDF1}" | sort | uniq | xargs -i sed -n 's/{}.*$//p' "${ORIGINAL1}" > "${F1_found}"

sed -n '2~4p' "${FILTEREDF2}" | sort | uniq | xargs -i sed -n 's/{}.*$//p' "${ORIGINAL2}" > "${F2_found}"

sed -n '2~4p' "${FILTEREDR1}" | sort | uniq | xargs -i sed -n 's/{}.*$//p' "${ORIGINAL1}" > "${R1_found}"

sed -n '2~4p' "${FILTEREDR2}" | sort | uniq | xargs -i sed -n 's/{}.*$//p' "${ORIGINAL2}" > "${R2_found}"

array=("${F1_found}" "${F2_found}" "${R1_found}" "${R2_found}")
echo "${array[@]}"

for file in ${array[@]}; do

matchesF=$(grep "${COIF}" $file --count)
matchesRC_F=$(grep "${COIF_RC}" $file --count)
matchesR=$(grep "${COIR}" $file --count)
matchesRC_R=$(grep "${COIR_RC}" $file --count)
printf "File1:%s\tmatchesF:%s\tmatchesRC_F:%s\tmatchesR:%s\tmatchesRC_R:%s\n" \
"$file" "${matchesF}" "${matchesRC_F}" "${matchesR}" "${matchesRC_R}" >>"${output_file}"
done
