#!/usr/bin/env Rscript
# This is an R script that deals with the metadata way more cleanly than awk
# Parameters it needs: metadata file, output folder
# it will generate 
# File1--- Folder1 (Barcode.p5.fasta)
#           |
#           |-For each barcodep5
#               |
#               | SubfolderA (Barcode.p7.fasta)
#                   |
#                   |-For each subfolder(PCR.primer.fasta)

### usage Rscript --vanilla scripts/r/generate.barcode.fasta.r <METADATA> <OUTPUTFOLDER>
library(tidyverse)
library(seqinr)
library(Biostrings)

arguments <- commandArgs(TRUE)

metadata <- read_csv(arguments[1])

output.folder <- arguments[2]

# get the i5 barcode.fasta

metadata %>% 
  distinct(plate_name.p5, barcode.p5) -> p5.barcodes

dir.create(output.folder)

write.fasta(sequences = as.list(p5.barcodes$barcode.p5),
            names = as.list(p5.barcodes$plate_name.p5),
            file.out = file.path(output.folder, "p5.barcodes.fasta"))

metadata %>% 
  group_by(plate_name.p5) %>% 
  nest() %>% 
  mutate(action = map2(plate_name.p5, data, function(.x,.y) {
    
    subfolder <- file.path(output.folder, .x)
    
    dir.create(subfolder)
    
    .y %>%  
      distinct(Well.p7, barcode.p7) ->
      mutate(revcom.p7 =reverseComplement(DNAString(barcode.p7)) )
    
    
    
  }))