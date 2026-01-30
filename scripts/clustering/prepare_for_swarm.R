arguments <- commandArgs(TRUE)

params <- list(folder = arguments[1])

library (eDNAfuns)
library (tidyverse)

ASV_table <- read_csv(file.path(params$folder, "ASV_table.csv"))
Hash_key  <- read_csv(file.path(params$folder, "Hash_key.csv"))

dir.create(file.path(params$folder, "swarm_input"))

ASV_table |> 
    inner_join(Hash_key) |>
    group_by(Sample) |>
    unite(Hash, nReads, col = "header", sep = ";size=") |> 
    nest() |> 
    mutate (write = walk2(Sample, data, function (.x, .y){
    
            .y |> 
                eDNAfuns::fasta_writer(sequence = sequence,
                   header = header, 
                   file.out =file.path(params$folder,"swarm_input" ,paste0(.x, ".fasta")))
  }))
