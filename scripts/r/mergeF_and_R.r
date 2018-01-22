seqtabF
seqtabR

head(seqtabF)
head(seqtabR)

seqtab.F.df
seqtab.R2.df <- seqtab.R.df
final.seqtab<-data.frame(row.names = rownames(seqtab.F.df))


for (i in 1:ncol(seqtab.F.df)) {   #for each column of the first dataframe
  
  current_seq<-colnames(seqtab.F.df)[i]
  
  if (current_seq %in% colnames(seqtab.R2.df)) { # is that column present on the second df?
    
    final.seqtab[,current_seq]<-seqtab.F.df[,i] + seqtab.R2.df[,current_seq] #if yes, the new df has the sum of reads
    
    seqtab.R2.df[,current_seq]<-NULL
    #seqtab.R2.df<- seqtab.R2.df[,-current_seq] #we delete the column from the second df to speed up next search
    
    print (ncol (seqtab.R2.df)) #check that is true
    
    
  } else {        # if the column is not present, then the new df is the value of the first df
    
    final.seqtab[,current_seq] <- seqtab.F.df[,i]
    
  }
  
  
  
}
# Now cbind both dataset
final.seqtab <- cbind(final.seqtab,seqtab.R2.df)

fst<-as.tibble(gather(final.seqtab, key=sequence, value = abundance, -sample))

test2<-getUniques(fst)
