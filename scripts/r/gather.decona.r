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

# get the sample_trans  

trans.map <-read_delim(file.path(args[1], "sample_trans.tmp"), 
                       delim = ";", escape_double = FALSE, col_names = c("Name.p5", "barcode.p7"), 
                       trim_ws = TRUE)
metadata <- read_csv(file.path(args[1], "metadata.csv"))

# get the number of good outputs
all.outputs <- list.dirs(original.file, recursive = F, full.names = T)

# How many did get  a medaka consensus

check.medaka <- function(folder){
  
  ifelse(file.exists(file.path(folder, "all", "multi-seq", "all_medaka_fastas.fasta" )),
     folder,
     F)
  
  
}

# Work only on those who made it through

map(all.outputs, check.medaka) %>% 
  set_names(nm = list.dirs(original.file, recursive = F, full.names = F)) %>% 
  purrr::discard(is.logical) -> passing.medaka

# For each of them , get the medaka consensus, extract the nreads

read.consensus <- function(folder){
  temp <- insect::readFASTA(file.path(folder, "all", "multi-seq", "all_medaka_fastas.fasta" ),
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
map (passing.medaka, read.consensus) %>% 
  bind_rows(.id = "Sample") -> all.together

all.together %>% 
  ungroup %>% 
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
  
