# A script to parse all sequences and abundances after running decona
# Rscript parse.decona.r <OUTPUT_FOlder>
library(tidyverse)
library(insect)
library(digest)
library(seqinr)

# First gather the parameter

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

args[1] 

# list all original fastq files

original.file <- list.dirs(path = file.path(args[1], "noprimers"), recursive = F, full.names = T)

# list all plates

original.plates <- map(original.file, ~ list.dirs(path = file.path(.x), recursive = F, full.names = T))

# list all wells

original.wells <- map(original.plates, ~ list.dirs(path = file.path(.x), recursive = F, full.names = T))

# get the metadata 

# trans.map <-read_delim(file.path(args[1], "sample_trans.tmp"), 
#                        delim = ";", escape_double = FALSE, col_names = c("Name.p5", "barcode.p7"), 
#                        trim_ws = TRUE)

metadata <- read_csv(file.path(args[1], "metadata.csv"))

tibble (Full.paths =  flatten(original.wells)) %>% 
  unnest() %>% 
  separate(Full.paths, into = c(NA, "file"), remove = F, sep = "noprimers/") %>% 
  separate(file, into = c("original.fastq", "Plate", "Well"), sep = "/") %>% 
  filter (!Well %in% c("pcr","barcodes")) -> all.files
  

# How many did get  a medaka consensus

check.medaka <- function(folder){
  
  ifelse(file.exists(file.path(folder,  "multi-seq", "all_medaka_fastas.fasta" )),
     T,
     F)
  
  
}

# Work only on those who made it through

all.files %>% 
  mutate(worked = map_lgl(Full.paths, check.medaka)) %>% 
  filter (worked) -> passing.medaka 



# map(all.outputs, check.medaka) %>% 
#   set_names(nm = list.dirs(original.file, recursive = F, full.names = F)) %>% 
#   purrr::discard(is.logical) -> passing.medaka

# For each of them , get the medaka consensus, extract the nreads

read.consensus <- function(folder){
  temp <- insect::readFASTA(file.path(folder, "multi-seq", "all_medaka_fastas.fasta" ),
                    bin = F)
  tibble (old.names = names(temp),
          seq = temp) %>% 
    separate(old.names, into = c(NA, NA, "nReads")) %>% 
    mutate (nReads = str_remove(nReads, ".fasta"),
            nReads = as.numeric(nReads)) %>% 
    group_by(seq) %>% 
    summarise(nReads = sum(nReads)) %>% # If things converge during medaka, they will be kept as separate fastas, converge this way
    group_by(seq) %>%
    mutate(Hash = digest::sha1(seq)) -> temp
  
  # return(list ("ASV_table" = temp %>% ungroup() %>% select(Sample, Hash, nReads),
  #              "Hash_key" = temp %>% select(Hash, seq)))
    return(temp)
}
passing.medaka %>% 
  mutate(all.consensus = map (Full.paths, read.consensus)) -> passing.medaka

passing.medaka %>% 
  select(original.fastq, Well, all.consensus) %>% 
  unnest(all.consensus)-> all.together

all.together %>% 
  ungroup %>% 
  unite(original.fastq, Well, col = "Sample", sep = "_Plate_") %>% 
  select(Sample, Hash, nReads) %>% 
  write_csv(file = file.path(args[1], "ASV_table.csv"))

all.together %>% 
  ungroup() %>% 
  select(Hash, seq) %>% 
  distinct() -> hashes 
  write_csv(hashes, file.path(args[1], "hash_key.csv"))
  
  seqinr::write.fasta(sequences = as.list(hashes$seq),
                      names = as.list(hashes$Hash),
                      file.out = file.path(args[1], "hash_key.fasta"))
  
