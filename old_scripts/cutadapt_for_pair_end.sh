#!/usr/bin/env bash
#Process both reads at the same time and keep them as
#separated files
# Maybe we should go barcode by barcode?

cutadapt -a Tag1=GATGAC  -A Tag1=GTCATC -o out_Tag1.1.fastq -p out_Tag1.2.fastq ~/Google_Drive/Run_Nov17/OA_COI/171122_sub/Lib-B_S2_L001_R1_001_sub.fastq ~/Google_Drive/Run_Nov17/OA_COI/171122_sub/Lib-B_S2_L001_R2_001_sub.fastq


cutadapt -g Tag1=GATGAC  -G Tag1=GTCATC -o out_Tag1_b.1.fastq -p out_Tag1_b.2.fastq ~/Google_Drive/Run_Nov17/OA_COI/171122_sub/Lib-B_S2_L001_R1_001_sub.fastq ~/Google_Drive/Run_Nov17/OA_COI/171122_sub/Lib-B_S2_L001_R2_001_sub.fastq
