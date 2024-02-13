
arguments <- commandArgs(TRUE)



params <- list(folder =arguments[1],
 hash = arguments[4],
  cont = arguments[7:10], fastqs = arguments[2], original = arguments[5:6])

## ----setup, include=FALSE-----------------------------------------------------


## ----loading packages, echo=FALSE ,message=FALSE------------------------------

library (tidyverse)
library (stringr)
library (dada2)
library (Biostrings)
library (digest)

#fastq.folder="/Users/rgallego/fastqs_demultiplexed_for_DADA2/demultiplexed_20180108_1539"
sample.map <- read_delim(paste0(params$folder,"/sample_trans.tmp"),col_names = c("Full_Id", "fastq_header","Sample"),delim = "\t")
head (sample.map)

path1 <- params$fastqs
params$original[1]
Sys.setenv(READ1=params$original[1])
Sys.setenv(READ2=params$original[2])


## ----listing files------------------------------------------------------------

F1s <- sort(list.files(path1, pattern="_Fwd.1.fastq", full.names = TRUE))
F2s <- sort(list.files(path1, pattern="_Fwd.2.fastq", full.names = TRUE))
R1s <- sort(list.files(path1, pattern="_Rev.1.fastq", full.names = TRUE))
R2s <- sort(list.files(path1, pattern="_Rev.2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- str_replace(basename(F1s), "_Fwd.1.fastq","")

#Now use only those that reflect a real sample

good.sample.names<-sample.names[sample.names %in% sample.map$fastq_header]
# Introduce here the biological counterpart of the fastq file

real.sample.name <- sample.map[,3][match(good.sample.names,sample.map$fastq_header),1]

# Now subset the number of files to be looked at

F1sgood<-F1s[sample.names %in% sample.map$fastq_header]
F2sgood<-F2s[sample.names %in% sample.map$fastq_header]
R1sgood<-R1s[sample.names %in% sample.map$fastq_header]
R2sgood<-R2s[sample.names %in% sample.map$fastq_header]




## ----filter and trim----------------------------------------------------------
filt_path <- file.path(params$folder, "/filtered") # Place filtered files in filtered/ subdirectory
filtF1s <- file.path(filt_path, paste0(good.sample.names, "_F1_filt.fastq.gz"))
filtF2s <- file.path(filt_path, paste0(good.sample.names, "_F2_filt.fastq.gz"))
out_Fs <- filterAndTrim(F1sgood, filtF1s, F2sgood, filtF2s, truncLen=c(200,200),
                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

filtR1s <- file.path(filt_path, paste0(good.sample.names, "_R1_filt.fastq.gz"))
filtR2s <- file.path(filt_path, paste0(good.sample.names, "_R2_filt.fastq.gz"))
out_Rs <- filterAndTrim(R1sgood, filtR1s, R2sgood, filtR2s, truncLen=c(200,200),
                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE)


## ----as tibble----------------------------------------------------------------
out_Fs_tibble<-tibble(file=dimnames(out_Fs)[[1]], reads.out=out_Fs[,2])
out_Rs_tibble<-tibble(file=dimnames(out_Rs)[[1]], reads.out=out_Rs[,2])
out_Fs_tibble


## ----learning errors, echo=T--------------------------------------------------
errF1 <- learnErrors(filtF1s, multithread=TRUE,verbose = 0)
errF2 <- learnErrors(filtF2s, multithread=TRUE,verbose = 0)
errR1 <- learnErrors(filtR1s, multithread=TRUE,verbose = 0)
errR2 <- learnErrors(filtR2s, multithread=TRUE,verbose = 0)

# Write errors to csv to see if they matter at all
tosave <- list(errF1, errF2, errR1, errR2)

saveRDS(tosave, file = "all.errors.rds")







## ----dereplication, echo=F,message=FALSE--------------------------------------
derepF1s <- derepFastq(filtF1s, verbose=TRUE)
derepF2s <- derepFastq(filtF2s, verbose=TRUE)
derepR1s <- derepFastq(filtR1s, verbose=TRUE)
derepR2s <- derepFastq(filtR2s, verbose=TRUE)

# Name the derep-class objects by the sample names
rownames(out_Fs) <- names(derepF1s) <- names(derepF2s) <- real.sample.name$Sample

rownames(out_Rs) <- names(derepR1s) <- names(derepR2s) <- real.sample.name$Sample



## ----dadaing, message=FALSE---------------------------------------------------
dadaF1s <- dada(derepF1s, err = errF1, multithread = TRUE)
dadaF2s <- dada(derepF2s, err = errF2, multithread = TRUE)
dadaR1s <- dada(derepR1s, err = errR1, multithread = TRUE)
dadaR2s <- dada(derepR2s, err = errR2, multithread = TRUE)



## ----saving to disk-----------------------------------------------------------
tosave <- list(derepF1s, derepF2s, derepR1s, derepR2s, dadaF1s, dadaF2s, dadaR1s,dadaR2s)

   saveRDS(tosave, file = "tosave.rds")




## ----merging pairs------------------------------------------------------------

to.keep.F <- map_lgl(dadaF1s, ~(sum(.x$denoised) > 150)) # Which samples have more than 150 reads passing filters in the Fwd direction
to.keep.R <- map_lgl(dadaR1s, ~(sum(.x$denoised) > 150)) #  Which samples have more than 150 reads passing filters in the Rev direction

to.keep <- (to.keep.F + to.keep.R) == 2 # Keep only those that pass boths filters

mergersF <- mergePairs(dadaF1s[to.keep],
   derepF1s[to.keep],
   dadaF2s[to.keep],
   derepF2s[to.keep],
   verbose = 0)

#Run a for loop that adds the number of unique reads that went into each ASV

for (j in 1:length(mergersF)){

  dadaF1s[[j]]@.Data[[2]] %>% rownames_to_column(var="forward") %>% select("forward", "nunq") ->Fwd
  Fwd$forward<-as.integer(Fwd$forward)
  dadaF2s[[j]]@.Data[[2]] %>% rownames_to_column(var="reverse") %>% select("reverse", "nunq") ->Rev
  Rev$reverse<-as.integer(Rev$reverse)

  mergersF[[j]] <- left_join(mergersF[[j]],Fwd, by="forward") %>% left_join(Rev, by="reverse") %>% mutate(nunq=pmin(nunq.x,nunq.y)) %>% select(-nunq.x,-nunq.y)


}

mergersR <- mergePairs(dadaR1s[to.keep],
   derepR1s[to.keep],
   dadaR2s[to.keep],
   derepR2s[to.keep],
   verbose = 0)

for (j in 1:length(mergersR)){

  dadaR1s[[j]]@.Data[[2]] %>% rownames_to_column(var="forward") %>% select("forward", "nunq") ->Fwd
  Fwd$forward<-as.integer(Fwd$forward)
  dadaR2s[[j]]@.Data[[2]] %>% rownames_to_column(var="reverse") %>% select("reverse", "nunq") ->Rev
  Rev$reverse<-as.integer(Rev$reverse)

  mergersR[[j]] <- left_join(mergersR[[j]],Fwd, by="forward") %>% left_join(Rev, by="reverse") %>% mutate(nunq=pmin(nunq.x,nunq.y)) %>% select(-nunq.x,-nunq.y)

}



## ----merging F and R (1)------------------------------------------------------

seqtabF <- makeSequenceTable(mergersF)

dim(seqtabF)

table(nchar(getSequences(seqtabF)))

seqtabR <- makeSequenceTable(mergersR)

dim(seqtabR)

table(nchar(getSequences(seqtabR)))


## ----merging F and R (2)------------------------------------------------------
reversed_sequences<-as.character(reverseComplement(DNAStringSet(colnames(seqtabR))))

summary (colnames(seqtabF) %in% reversed_sequences)

summary (reversed_sequences %in% colnames(seqtabF))

colnames(seqtabR)<-reversed_sequences



## ----merging F and R (3)------------------------------------------------------

seqtab.R.df=as.data.frame(seqtabR)

final.seqtab<-data.frame(row.names = rownames(seqtabF))


for (i in 1:ncol(seqtabF)) {   #for each column of the first dataframe

  current_seq<-colnames(seqtabF)[i]

  if (current_seq %in% colnames(seqtab.R.df)) { # is that column present on the second df?

    final.seqtab[,current_seq]<-seqtabF[,i] + seqtab.R.df[,current_seq] #if yes, the new df has the sum of reads

    seqtab.R.df[,current_seq]<-NULL
    #seqtab.R2.df<- seqtab.R2.df[,-current_seq] #we delete the column from the second df to speed up next search


  } else {        # if the column is not present, then the new df is the value of the first df

    final.seqtab[,current_seq] <- seqtabF[,i]

  }



}
# Now cbind both dataset
final.seqtab <- as.matrix (cbind(final.seqtab,seqtab.R.df))



## ----RemovingChimeras, message=F----------------------------------------------

seqtab.nochim <- removeBimeraDenovo(final.seqtab, method="consensus", multithread=TRUE)

dim(seqtab.nochim)




## ----write colnames-----------------------------------------------------------
ifelse(length(colnames(seqtab.nochim))>500, to_write<- sample(colnames(seqtab.nochim),size = 500, replace = F), to_write<-colnames(seqtab.nochim))
#to_write<- sample(colnames(seqtab.nochim),size = 500, replace = F)
write_lines(to_write, path = "seqnames.txt")
getwd()


## ----tidying and writing------------------------------------------------------

seqtab.nochim.df=as.data.frame(seqtab.nochim)

ASV.file <- file.path(params$folder, "ASV_table.csv")

# Now decide if you want hashing or not

if (grepl ("yes", params$hash, ignore.case = TRUE)) {

  conv_file <-  file.path(params$folder,"hash_key.csv")

  conv_table <- tibble( Hash = "", Sequence ="")

  #hashes <- list(NULL)
  
  map_chr (colnames(seqtab.nochim.df), ~ digest(.x, algo = "sha1", serialize = F, skip = "auto")) -> Hashes
  conv_table <- tibble (Hash = Hashes,
                        Sequence = colnames(seqtab.nochim.df))
  
  colnames(seqtab.nochim.df) <- Hashes

  # for (i in 1:ncol(seqtab.nochim.df)) {   #for each column of the dataframe
  # 
  #   current_seq <-colnames(seqtab.nochim.df)[i]
  # 
  #   current_hash <- digest(current_seq,algo = "sha1",serialize = F,skip = "auto")
  # 
  #   hashes[[i]] <-  current_hash
  # 
  #   conv_table [i,]<-c(current_hash, current_seq)
  # 
  #   colnames(seqtab.nochim.df)[i] <- current_hash
  # 
  # }

  write_csv(conv_table, conv_file)


seqtab.nochim.df %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to="Hash", values_to = "nReads") %>%
  filter (nReads > 0) -> current_asv

write_csv(current_asv, ASV.file)

} else { #What do we do if you don't want hashes: two things - Change the header of the ASV table, write only one file

  seqtab.nochim.df %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to="Sequence", values_to = "nReads") %>%
  filter (nReads > 0) ->  current_asv
  write_csv(current_asv, ASV.file)
}




## ----output_summary-----------------------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- as.data.frame(cbind(out_Fs, out_Rs,
                             sapply(dadaF1s, getN), sapply(dadaR1s, getN),
                             sapply(mergersF, getN),sapply(mergersR, getN),
                             rowSums(seqtabF),rowSums(seqtabR),
                             rowSums(final.seqtab), rowSums(seqtab.nochim)))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input_F", "input_R", "filtered_F","filtered_R",
                     "denoised_F", "denoised_R",
                     "merged_F", "merged_R",
                     "tabled_F", "tabled_R",
                     "tabled_together", "nonchim")

write_csv(track, paste0(params$folder,"/dada2_summary.csv"))



## ----merging with previous datasets-------------------------------------------

if (grepl ("yes", params$cont[1], ignore.case = TRUE)) {

old_hash = read_csv(params$cont[2])

old_asv = read_csv(params$cont[3])

all_hash = union_all(old_hash, conv_table)

new.hashes = setdiff(conv_table,old_hash)

reads.new.hashes <- left_join(new.hashes, current_asv, by = c("Hash" = "Sequence"))

all_asv = bind_rows (old_asv, current_asv)

merging.output = tibble(Date = Sys.Date(), Old_db = params$cont[2], Old_ASV = params$cont[3],
                        Input_folder= params$folder, Nseq_old_db = nrow(old_hash), Nseq_new_input = nrow(conv_table),
                        Nseq_new_db = nrow(all_hash), N_newseq = nrow(new.hashes), Nreads_newseq = sum(reads.new.hashes$nReads) )
ifelse(file.exists(params$cont[4]),
       write_csv(merging.output, params$cont[4], append = T),
       write_csv(merging.output, params$cont[4]) )

write_csv(all_hash, paste0(params$folder, "/Hash_db_",Sys.Date(),".csv"))
write_csv(all_asv, paste0(params$folder, "/ASV_all_",Sys.Date(),".csv"))

}


