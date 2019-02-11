# demultiplexer_for_DADA2
If you want to use DADA2 with libraries created by adding adapters through ligation and with two levels of adapters - see banzai (github.com/jimmyodonnell/banzai) - you need to demultiplex samples before pairing both ends.

This script will look for your adapters and pcr primers on your reads, and return 4 fastq files per unique sample: Fwd.1, Fwd.2, Rev.1 and Rev.2.

## Dependencies

This script has been tested in MacOS X. It should work in Linux as well, but I think different `sed` versions will return different results.

The script relies mostly on R pckages. Besides that, you will need the latest version of cutadapt (https://cutadapt.readthedocs.io/en/stable/index.html).

### R packages

Checking whether you have the right R packages installed before running R is complicated - and you don't want to realise you don't have them by the time the pipeline gets there. So please check that you have these packages installed:

* tidyverse
* devtools
* dada2
* rmarkdown
* Biostrings
* digest
