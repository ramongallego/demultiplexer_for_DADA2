
arguments <- commandArgs(TRUE)



params <- list(folder =arguments[1],
 hash = arguments[3],
 fastqs = arguments[2],
 len1= as.numeric(arguments[4]),
 len2= as.numeric(arguments[5]))

params 

## ----setup, include=FALSE-----------------------------------------------------
#TODO: make sure the Rscript uses all allocated resources


## ----loading packages, echo=FALSE ,message=FALSE------------------------------

library (tidyverse)
library (dada2)
library (Biostrings)
library (digest)
library (insect)
library (tictoc)


sample.map <- read_delim(file.path(params$folder,"/sample_trans.tmp"),col_names = c("Full_Id", "fastq_header","Sample"),delim = "\t")
print("The sample map looks like this")
head (sample.map)

path1 <- params$fastqs



## ----listing files------------------------------------------------------------
files.noprimers <- tibble (files = list.files(path1, full.names = TRUE))

print("These are the demulted files")
files.noprimers

files.noprimers |> 
  mutate(locus = str_extract(files, "(?<=_Locus_)[^_]+"),
         direction = str_extract(files, "(Fwd|Rev)\\.R[12]"),
         fastq_header = str_extract(basename(files), "^[^_]+_[^_]+")) |>
  pivot_wider (names_from = "direction",
               values_from = "files") -> files.noprimers

# Keep only real filenames and complete cases
print("Now there should be four files per sample/locus")
files.noprimers


files.noprimers |>
 drop_na() |>
 inner_join(sample.map |> 
             select(fastq_header, Sample )) -> files.noprimers

print("We should have reduced the dataset to those matching the sample map and added the real sample names")
files.noprimers

## ----filter and trim----------------------------------------------------------
filt_path <- file.path(params$folder, "/filtered") # Place filtered files in filtered/ subdirectory

filt_function <- function(file1, file2){
  filt1s <- file.path(filt_path,basename(file1))
  filt2s <- file.path(filt_path,basename(file2))
        filterAndTrim(file1, filt1s, file2, filt2s, 
                      truncLen = c(params$len1,params$len2), 
                      maxN=0, maxEE=c(2,2),
                      truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE)  |> 
        as_tibble()
    }

## TODO implement futuremap for multicore usage - testing now


tic("filtering")
files.noprimers |> 
  mutate(filtF1s = file.path(filt_path,basename(Fwd.R1)),
         filtF2s = file.path(filt_path,basename(Fwd.R2)),
         filtR1s = file.path(filt_path,basename(Rev.R1)),
         filtR2s = file.path(filt_path,basename(Rev.R2)),
         outFs = map2_dfr (Fwd.R1, Fwd.R2, filt_function ),
         outRs = map2_dfr (Rev.R1, Rev.R2, filt_function)) -> files.noprimers

# discard those with fewer than 100 seqs passing either filter
toc()

files.noprimers |>
 filter (outFs$reads.out >100 & outRs$reads.out > 100) -> goodqs

 files.noprimers |> 
  anti_join(goodqs) -> discarded

rm(files.noprimers)

#### Learn 4 errors objects: these are a function of the NEXTSEQ run and not of 
#### the sample, so it does not make sense to calculate them once per row.
#### We need to point towards the filtered files

## ----learning errors, echo=T- this is so intensive I would rather use all cores on each error calculation-------------------------------------------------

tic("Learning errors")

errF1 <- learnErrors(goodqs$filtF1s, multithread=TRUE,verbose = 0, nbases = 100e6)
errF2 <- learnErrors(goodqs$filtF2s, multithread=TRUE,verbose = 0, nbases = 100e6)
errR1 <- learnErrors(goodqs$filtR1s, multithread=TRUE,verbose = 0, nbases = 100e6)
errR2 <- learnErrors(goodqs$filtR2s, multithread=TRUE,verbose = 0, nbases = 100e6)

toc()
# Write errors to csv to see if they matter at all
tosave <- list(errF1, errF2, errR1, errR2)

saveRDS(tosave, file = "all.errors.rds")

## ----dereplication, echo=F,message=FALSE--------------------------------------
tic("Dereplicating")
goodqs |> 
  mutate (across( 
                  starts_with("filt"),
                  ~  map(.x,derepFastq),                 # Transformation (identity in this example)
                  .names = "{gsub('filt', 'derep', .col)}" # Rename by replacing "filt" with "derep"
    )
  ) -> goodqs
toc()
## ----dadaing, message=FALSE---------------------------------------------------

tic("dada")
goodqs |> 
  mutate (
          dadaF1s = map(derepF1s, ~dada(.x, err = errF1, multithread = TRUE, verbose = F)),
          dadaF2s = map(derepF2s, ~dada(.x, err = errF2, multithread = TRUE, verbose = F)),
          dadaR1s = map(derepR1s, ~dada(.x, err = errR1, multithread = TRUE, verbose = F)),
          dadaR2s = map(derepR2s, ~dada(.x, err = errR2, multithread = TRUE, verbose = F))) -> goodqs

## TODO: do this only if you have a HOARD=yes From here onwards takes very little time
saveRDS(goodqs, file = "tosave.rds")
toc()

## ----merging pairs------------------------------------------------------------
tic("merging")
goodqs |>
  mutate(mergersF = pmap(.l = list(dadaF1s,derepF1s, dadaF2s,derepF2s), 
                         .f = mergePairs,
                         minOverlap = 10),
         mergersR = pmap(.l = list(dadaR1s,derepR1s, dadaR2s,derepR2s), 
                         .f = mergePairs,
                         minOverlap = 10)) -> goodqs
toc()

## merging F and R , and chimeras1

goodqs |>
  mutate (mergersR = map (mergersR, ~ .x |> 
                            mutate(sequence = insect::rc(sequence))),
          joined   = map2(mergersF, mergersR, ~ bind_rows(.x, .y) |> 
                          filter (accept) |> 
                          group_by(sequence, accept) |> 
                          summarise(across(everything(), sum))),
          
          nochim = map(joined, removeBimeraDenovo)) -> goodqs



## ----tidying and writing------------------------------------------------------
goodqs|> 
  select(Sample,locus,nochim) |>
  unnest(nochim) |> 
  select(Sample,locus, sequence, nReads = abundance ) -> Abundance_table

if ( grepl ("yes", params$hash, ignore.case = TRUE)) {

  Abundance_table |> 
  distinct(sequence) |> 
  mutate (Hash = map_chr(sequence, digest::sha1)) -> Hash_key

write_csv(Hash_key, file.path(params$folder, "Hash_key.csv"))
eDNAfuns::fasta_writer(df = Hash_key,
                       sequence = sequence, 
                       header = Hash,
                       file.out = file.path(params$folder, "Hash_key.fasta"))

Abundance_table |> 
  inner_join(Hash_key) |> 
  select(Sample, locus,Hash, nReads)  |> 
  write_csv( file.path(params$folder, "ASV_table.csv"))

}else {write_csv(Abundance_table, file.path(params$folder, "ASV_table.csv"))}

### Summary stats 
getN <- function(x) sum(getUniques(x))

goodqs |> 
  select(Sample, locus, outFs, outRs, starts_with("dada"), starts_with ("merg"),joined, nochim) |> 
  unnest(outFs) |>
  dplyr::rename(inputF= reads.in, filtered.F = reads.out) |>
  unnest(outRs) |>
  dplyr::rename(inputR= reads.in, filtered.R = reads.out) |> 
  mutate (across(c(starts_with("dada"), starts_with ("merg"),joined, nochim),
                  ~map_dbl(.x, getN))) -> summary.good


if(nrow(discarded)>=1){              
bind_rows(summary.good, discarded|>
                        unnest(outFs) |>
                        dplyr::rename(inputF= reads.in, filtered.F = reads.out) |>
                        unnest(outRs) |>
                        dplyr::rename(inputR= reads.in, filtered.R = reads.out)) |>
  write_csv(file.path(params$folder, "summary_dada2.csv"))

} else {write_csv(summary.good, file.path(params$folder, "summary_dada2.csv"))}

