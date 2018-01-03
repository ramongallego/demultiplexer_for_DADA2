#!/usr/bin/env bash

# Usage is bash find_pcr_primer.sh folder

folder="${1}"

COIF="GG[A-Z]AC[A-Z]GG[A-Z]TGAAC[A-Z]GT[A-Z]TA[A-Z]CC[A-Z]CC"
COIF_RC="GG[A-Z]GG[A-Z]TA[A-Z]AC[A-Z]GTTCA[A-Z]CC[A-Z]GT[A-Z]CC"
COIR="TA[A-Z]AC[A-Z]TC[A-Z]GG[A-Z]TG[A-Z]CC[A-Z]AA[A-Z]AA[A-Z]CA"
COIR_RC="TG[A-Z]TT[A-Z]TT[A-Z]GG[A-Z]CA[A-Z]CC[A-Z]GA[A-Z]GT[A-Z]TA"


for file in "${folder}"/*.1.fastq; do


READ1="${file}"

READ2=$(echo ${file} | sed 's/.1.fastq/.2.fastq/g' )

echo "number of primer matches of Fwd Primer on ${READ1}"
cat "${READ1}" | grep "${COIF}" --color -B1 --count

echo 'number of primer matches of Fwd Primer (RC) on Read1'
cat "${READ1}" | grep "${COIF_RC}" --color -B1 --count

echo 'number of primer matches of Rev Primer on Read1'
cat "${READ1}" | grep "${COIR}" --color -B1 --count

echo 'number of primer matches of Rev Primer (RC) on Read1'
cat "${READ1}" | grep "${COIR_RC}" --color -B1 --count

echo 'number of primer matches of Fwd Primer on Read2'
cat "${READ2}" | grep "${COIF}" --color -B1 --count

echo 'number of primer matches of Fwd Primer (RC) on Read2'
cat "${READ2}" | grep "${COIF_RC}" --color -B1 --count

echo 'number of primer matches of Rev Primer on Read2'
cat "${READ2}" | grep "${COIR}" --color -B1 --count

echo 'number of primer matches of Rev Primer (RC) on READ2'
cat "${READ2}" | grep "${COIR_RC}" --color -B1 --count
done
