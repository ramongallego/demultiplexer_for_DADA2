
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
         fastq_header = str_extract(basename(files), "^.+(?=_Locus_)")) |>
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
names(files.noprimers$Fwd.R1) <- names(files.noprimers$Fwd.R2) <- names(files.noprimers$Rev.R1) <- names(files.noprimers$Rev.R2) <- files.noprimers$Sample
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
         filtR2s = file.path(filt_path,basename(Rev.R2))) -> files.noprimers
outFs <- filt_function(files.noprimers$Fwd.R1,
                       files.noprimers$Fwd.R2 ) |>
                        mutate(Sample = files.noprimers$Sample)
outFs
outRs <- filt_function(files.noprimers$Rev.R1,
                        files.noprimers$Rev.R2) |>
                        mutate(Sample = files.noprimers$Sample)
outRs

toc()

# discard those with fewer than 100 seqs passing either filter

files.noprimers |> 
  inner_join(outFs |> 
             select(Sample, reads.outF = reads.out), by = "Sample") |>
  inner_join(outRs |>
             select(Sample, reads.outR = reads.out), by = "Sample") |>
  filter (reads.outF >100 & reads.outR > 100) -> goodqs


goodqs

#### Learn 4 errors objects: these are a function of the NEXTSEQ run and not of 
#### the sample, so it does not make sense to calculate them once per row.
#### We need to point towards the filtered files

## ----learning errors, echo=T- this is so intensive I would rather use all cores on each error calculation-------------------------------------------------
print ("Learning errors")
tic("Learning errors")

errF1 <- learnErrors(goodqs$filtF1s, multithread=TRUE,verbose = 0, nbases = 100e6)
errF2 <- learnErrors(goodqs$filtF2s, multithread=TRUE,verbose = 0, nbases = 100e6)
errR1 <- learnErrors(goodqs$filtR1s, multithread=TRUE,verbose = 0, nbases = 100e6)
errR2 <- learnErrors(goodqs$filtR2s, multithread=TRUE,verbose = 0, nbases = 100e6)

toc()
# Write errors to csv to see if they matter at all
tosave <- list(errF1, errF2, errR1, errR2)

saveRDS(tosave, file = file.path(params$folder,"all.errors.rds"))

## ----dereplication, echo=F,message=FALSE--------------------------------------
tic("Dereplicating")

names(goodqs$filtF1s) <- goodqs$Sample
names(goodqs$filtF2s) <- goodqs$Sample
names(goodqs$filtR1s) <- goodqs$Sample
names(goodqs$filtR2s) <- goodqs$Sample

derepF1s <- derepFastq(goodqs$filtF1s, verbose = 0)
derepF2s <- derepFastq(goodqs$filtF2s, verbose = 0)
derepR1s <- derepFastq(goodqs$filtR1s, verbose = 0)
derepR2s <- derepFastq(goodqs$filtR2s, verbose = 0)

toc()
## ----dadaing, message=FALSE---------------------------------------------------

tic("dada")
dadaF1s <- dada(derepF1s, err = errF1, multithread = TRUE)
dadaF2s <- dada(derepF2s, err = errF2, multithread = TRUE)
dadaR1s <- dada(derepR1s, err = errR1, multithread = TRUE)
dadaR2s <- dada(derepR2s, err = errR2, multithread = TRUE)
toc()


## TODO: do this only if you have a HOARD=yes From here onwards takes very little time
saveRDS(goodqs, file = file.path(params$folder,"tosave.rds"))
toc()

## ----merging pairs------------------------------------------------------------
tic("merging R1 and R2")
mergersF <- mergePairs(dadaF1s,derepF1s, dadaF2s,derepF2s,
                      minOverlap = 10)
mergersR <- mergePairs(dadaR1s,derepR1s, dadaR2s,derepR2s,
                      minOverlap = 10) 
toc()

## merging F and R , and chimeras1
tic("merging F and R")
mergersR <- map (mergersR, ~ .x |> 
                            mutate(sequence = insect::rc(sequence)))

joined   <- map2(mergersF, mergersR, ~ bind_rows(.x, .y) |> 
                          filter (accept) |> 
                          group_by(sequence, accept) |> 
                          summarise(across(everything(), sum),  .groups = "drop"))
joined |>
  write_rds(file.path(params$folder, "joined.rds"))
toc()

tic("removing chimeras")
nochim <- removeBimeraDenovo(joined, multithread= TRUE) 
toc()


## ----tidying and writing------------------------------------------------------

nochim |> 
  bind_rows(.id = "Sample") |>
  inner_join(goodqs |> 
             select(Sample, locus), by = "Sample") |>
  select(Sample, locus, sequence, nReads = abundance) -> Abundance_table 

Abundance_table |>
  write_rds(file.path(params$folder, "Abundance_table.rds"))
 

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

files.noprimers |> 
  select(Sample, locus) |>
  inner_join(outFs |> 
             select(Sample, reads.outF = reads.out), by = "Sample") |>
  inner_join(outRs |>
             select(Sample, reads.outR = reads.out), by = "Sample") -> summary_dada2

# Example named lists
derepF_counts   <- sapply(derepF1s, getN)
derepR_counts   <- sapply(derepR1s, getN)
dadaF_counts    <- sapply(dadaF1s, getN)
dadaR_counts    <- sapply(dadaR1s, getN)
mergedFs       <- sapply(mergersF, getN) 
mergedRs       <- sapply(mergersR, getN)
joined_counts <- sapply(joined, getN)
nochim_counts <- sapply(nochim, getN)

# Build a tibble
names_good <- summary_dada2$Sample
tibble(
  Sample = names_good,
  derepF = derepF_counts[names_good],
  derepR = derepR_counts[names_good],
  dadaF = dadaF_counts[names_good],
  dadaR = dadaR_counts[names_good], 
  mergedF = mergedFs[names_good],
  mergedR = mergedRs[names_good],
  joined = joined_counts[names_good],
  nochim = nochim_counts[names_good]
) %>%
 inner_join(summary_dada2,., by = "Sample") -> summary_dada2

write_csv(summary_dada2, file.path(params$folder, "summary_dada2.csv"))
