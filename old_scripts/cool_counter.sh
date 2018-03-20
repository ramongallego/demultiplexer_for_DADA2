for i in {0..15}; do
  echo -ne "$i"'\r'

 done
  echo

# Number of files with pattern in a foilder
n_files=(*R1_001_sub.fastq)
i_count=0
for j in "${n_files[@]}"; do
  i_count=$((i_count+1))
  echo -ne "Working on ${i_count} of ${#n_files[@]}"'\r'
  sleep 0.5
  #echo "${j}"
  #echo "${#n_files[@]}"
done
echo
