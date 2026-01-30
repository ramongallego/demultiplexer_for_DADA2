arguments <- commandArgs(TRUE)

params <- list(folder = arguments[1])

library (eDNAfuns)
library (tidyverse)

centroids.paths <- list.files(file.path(params$folder,"swarm_input","centroids"),
                              pattern = "centroids.fasta") 

map(centroids.paths, ~fasta_reader(file.path(params$folder,"swarm_input","centroids", .x))) -> seqs.centroids

seqs.centroids |>  
  set_names(nm= centroids.paths)-> centroids

centroids |>
  bind_rows(.id = "Sample") |> 
  separate(header, into = c("Hash", "nReads"), sep = ";size=|;", convert = T) |> 
  mutate(Sample = str_remove(Sample, ".centroids.fasta")) -> new_ASV

  new_ASV |>
   select(Sample, Hash, nReads) |>
   write_csv(file.path(params$folder, "new_ASV_after_swarm.csv"))

new_ASV |>
    ungroup() |>
   distinct ( Hash, seq) -> new_hash

new_hash |>
  write_csv(file.path(params$folder, "new_Hash_key_after_swarm.csv"))

new_hash |>
    fasta_writer(sequence = seq, 
                 header = Hash,
                 file.out = file.path(params$folder, "new_Hash_key_after_swarm.fasta"))