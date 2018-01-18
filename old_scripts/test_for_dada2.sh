#!/usr/bin/env bash

SCRIPT_DIR="/Users/rgallego/demultiplexer_for_dada2/scripts"

DEMULT_DIR="~/fastqs_demultiplexed_for_DADA2/demultiplexed_20180112_0930"

TEST="YES"

if [[ "${TEST}" = "YES" ]]; then

Rscript "${SCRIPT_DIR}"/r/biodiversity.r "${DEMULT_DIR}" "${SCRIPT_DIR}"

fi
#R -e "rmarkdown::render("dada2.Rmd", params=list(folder="${DEMULT_DIR}"))"
