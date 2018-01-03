DEMULT_DIR="/Users/Moncho/banzai_outputs/banzai_out_20171228_1903/demultiplexed"

FILEPATH=$( find "${DEMULT_DIR}" -name *".fastq"  )
FILEPATH2=$( ls "${DEMULT_DIR}/"*".fastq" )
FILEPATH3="${DEMULT_DIR}"/*.fastq
echo "using find -name"
echo "${FILEPATH}" | wc -l
echo "Can you access the 2nd element of a list?"
echo ${FILEPATH[1]}
echo

echo "using ls"
echo "${FILEPATH2}" | wc -l
echo "Can you access the 2nd element of a list?"
echo ${FILEPATH2[1]}
echo

echo "using just /"
echo "${FILEPATH3}"
echo "Can you access the 2nd element of a list?"
echo ${FILEPATH3[1]}


#This is the good way of doing it, skipping the ls
for file in "${DEMULT_DIR}"/*.fastq; do

#file in "${FILEPATH2}"; do
  echo "${file}"
  cat "${file}"  |\
  #Now find the sequence line of each read (starting from the
  #second line, every fourlines),
  sed -n '2~4p' |
  #And now extract the 6 characters of the barcode
  cut -c 4-9 |
  #Now sort | uniq  to get the counts
  sort | uniq -c | sort -r | head -n 5

done
