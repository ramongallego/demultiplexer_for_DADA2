#!/bin/bash

# Can we overwrite arguments
# Usage bash test.overwrite.sh banzai_params.sh ...
echo "looping around params"

for n in 0.8 0.85 0.9 0.95; do
  echo "${n}"
  mkdir n_"${n}"
  # overwrite the parameter in the parasms file
  sed "s/CLUSTER_SIM=\"0.8\"/CLUSTER_SIM=\"${n}\"/g" ~/demultiplexer_for_DADA2/banzai_params_for_nanopore.sh > banzai_params_"${n}".sh
echo "the first sed works "
  # overwrite the manin script so it uses the custom output dirname
  sed "s/OUTPUT_DIR=\"$\{OUTPUT_DIRECTORY\}/demultiplexed_$\{START_TIME\}\"/OUTPUT_DIR=\"~/test_ns/n_$\{n\}\"/g" ~/demultiplexer_for_DADA2/demultiplex_nanopore.sh > ~/demultiplexer_for_DADA2/demultiplex_nanopore_"${n}".sh

done
