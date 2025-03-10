arguments <- commandArgs(TRUE)

params <- list(folder = arguments[1],
                 hash = arguments[2])

params 

library (eDNAfuns)
library (tidyverse)
library (digest)

samples <- tibble (files = list.files(params$folder, pattern = "*centroids.fasta"))
sample_trans <- read_table(file.path(params$folder, "sample_trans.tmp"),
                           col_names = c("IDS", "Key", "Sample")) |> 
                           select(-IDS)

samples |> 
  separate(files, into = c("Key", "locus"), sep = "_Locus_", remove = F) |> 
  mutate(locus = str_remove(locus, "_centroids.fasta")) |> 
  mutate(seqs = map(files, ~fasta_reader(file.path(params$folder, .x))))-> samples

samples |> 
  inner_join(sample_trans) |> 
  select(-Key, -files) -> samples

samples |>
  unnest(seqs) |> 
  separate(header,
           into = c(NA, "nReads"),
           sep = ";size=",
           convert = T) -> samples

if ( grepl ("yes", params$hash, ignore.case = TRUE)) {

  samples |> 
    ungroup() |> 
    distinct(seq) |> 
    rowwise() |>
    mutate (Hash = sha1(seq)) -> Hash_key
  
  write_csv(Hash_key, file = file.path(params$folder,"Hash_key.csv"))
  
  eDNAfuns::fasta_writer(df = Hash_key,
                         sequence = sequence, 
                         header = Hash,
                         file.out = file.path(params$folder, "Hash_key.fasta"))
  
  samples |> 
    inner_join(Hash_key) |> 
    select (Sample, locus, Hash, nReads) -> Abundance
  
  write_csv(Abundance, file = file.path(params$folder,"ASV_table.csv"))
    

}else {
  samples |> 
    select (Sample, locus,seq, nReads) |> 
    write_csv( file.path(params$folder, "ASV_table.csv"))
}