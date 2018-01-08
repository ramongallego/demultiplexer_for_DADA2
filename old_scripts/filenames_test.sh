file="/Users/Moncho/demultiplexer_for_DADA2/README.md"

source /Users/Moncho/demultiplexer_for_DADA2/scripts/strip_path.sh "${file}" 'file'

barcodes="~/Google_Drive/Run_Nov17/OA_COI/171122_sub/demultiplexed/pcr_primers.fasta"
cat "${barcodes}"


echo $short_file

basename ${file}

source /Users/Moncho/demultiplexer_for_DADA2/scripts/functions/check_primers.sh "${barcodes}"
