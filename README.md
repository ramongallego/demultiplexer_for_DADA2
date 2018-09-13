# demultiplexer_for_DADA2
If you want to use DADA2 with libraries created by adding adapters through ligation and with two levels of adapters - see banzai (github.com/jimmyodonnell/banzai) - you need to demultiplex samples before pairing both ends.

This script will look for your adapters and pcr primers on your reads, and return 4 fastq files per unique sample: Fwd.1, Fwd.2, Rev.1 and Rev.2.
