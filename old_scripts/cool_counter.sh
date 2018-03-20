for i in {0..15}; do
  echo -ne "$i"'\r'

 done
  echo

# Number of files with pattern in a foilder
for j in *R1_001_sub.fastq; do
  echo "${#(*R1_001_sub.fastq)[@]}"
done
n_files=(*R1_001_sub.fastq)

echo "${#n_files[@]}"
