#!/usr/bin/time bash
#Trials to remove primers in two ways: reversing reads to have everything
#in the same direction, or keep them as they are to improve transition errors

#try with one demultiplexed pair of fastqs

READ1="/Users/Moncho/Google_Drive/Run_Nov17/OA_COI/171122_sub/demultiplexed/Lib_G_TCGCAT.1.fastq"
READ2="/Users/Moncho/Google_Drive/Run_Nov17/OA_COI/171122_sub/demultiplexed/Lib_G_TCGCAT.2.fastq"
DEMULT_DIR="/Users/Moncho/Google_Drive/Run_Nov17/OA_COI/171122_sub/demultiplexed"
mkdir "${DEMULT_DIR}"/cleaned
primerF_seq="GGWACWGGWTGAACWGTWTAYCCYCC"
primerR_seq="TANACYTCNGGRTGNCCRAARAAYCA"
ID1S="Lib_G"
RIGHT_BARCODE="TCGCAT"
primers="${DEMULT_DIR}"/pcr_primers.fasta

 printf ">FWD\n${RIGHT_BARCODE}${primerF_seq}\n>REV\n${RIGHT_BARCODE}${primerR_seq}\n" > "${primers}"



mkdir "${DEMULT_DIR}"/cleaned/"${ID1S}"
head "${READ1}"


#find and remove the fwd primer on the .1 read
#cutadapt -g FWD="${primerF_seq}"  -G "${primerR_seq}" \
#cutadapt -g REV="${primerR_seq}"  -G "${primerF_seq}" \
# -o "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_{name}_clean.1.fastq \
# -p "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_{name}_clean.2.fastq \
# "${READ1}" "${READ2}"


#cutadapt -g REV="${primerR_seq}"  -G "${primerF_seq}" \
#cutadapt -g FWD="${primerF_seq}"  -G "${primerR_seq}" \
# -o "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_{name}_clean.1.fastq \
# -p "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_{name}_clean.2.fastq \
# "${READ1}" "${READ2}"

 cutadapt -g file:"${primers}" --discard-untrimmed\
  -o "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_{name}_clean.1.fastq \
  -p "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_{name}_clean.2.fastq \
  "${READ1}" "${READ2}"

#And Now get rid of those in which the rev primer has not been found -

 cutadapt -g FWD="${primerR_seq}" --discard-untrimmed \
 -o "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_cleaner.2.fastq \
 -p "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_cleaner.1.fastq \
 "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_clean.2.fastq \
 "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_FWD_clean.1.fastq

 cutadapt -g REV="${primerF_seq}" --discard-untrimmed \
 -o "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_cleaner.2.fastq \
 -p "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_cleaner.1.fastq \
 "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_clean.2.fastq \
 "${DEMULT_DIR}"/cleaned/${ID1S[i]}/${ID1S[i]}-"${RIGHT_BARCODE}"_REV_clean.1.fastq


#cutadapt -g file:"${FASTA_file}" --no-trim  -o "${DEMULT_DIR}"/${ID1S[i]}/${ID1S[i]}-{name}_round1.1.fastq -p "${DEMULT_DIR}"/${ID1S[i]}/${ID1S[i]}-{name}_round1.2.fastq \
# "${READ1}" "${READ2}" --quiet --discard-untrimmed
#-g REV="${primerR_seq}"
#-G REV= "${primerF_seq}"
