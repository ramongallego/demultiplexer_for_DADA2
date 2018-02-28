#!/usr/bin/env bash

#Are the sequences selected by Dada2 real or modeled?
#I run a small test, but we can do it for all of them

#Extract the dada2 output, the sequences: a llist of sequences per Library

# run this from the folder where all derepF1s are, then it will look for the right
# file in the demultiplexed folder

derepfiles=(*)


for (( i=0; i < "${#derepfiles[@]}"; i++ )); do
  fastqfile="../demultiplexed/${derepfiles[i]%.*}_Fwd.1.fastq"
  n_matches=$(grep -f "${derepfiles[i]}" "${fastqfile}" -c)

  printf '%s\t%s\t%s\n' "${derepfiles[i]}" "${fastqfile}" \
   "${n_matches}" >> "output.summary.txt"


done



while read j; do echo $j |grep -o "matchingString"| wc -l;  done < Lib_A_ACACAC.txt

# Another thing to check:
# Does the bash pipeline reverse complement the sequences? or do I still have
# the Fwd_.1.fastqs in the Fwd directon?


#TEst: look for the first 7 sequences on one of the Fwd.1.fastqs
# use them as an pattern argument for grep, targetting the original Lib .1
# retrieve the 30 characters before the hit
# Check if it found the Fwd or reverse primer

for  i in demultiplexed/*_Fwd.1.fastq ; do
echo $(basename "${i}")




 ; done

demult.files.1F=(demultiplexed/*_Fwd.1.fastq)
