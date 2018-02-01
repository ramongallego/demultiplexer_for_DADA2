seqtab.R.df=as.data.frame(seqtabR)

final.seqtab<-data.frame(row.names = rownames(seqtab.F.df))


for (i in 1:ncol(seqtabF)) {   #for each column of the first dataframe
  
  current_seq<-colnames(seqtabF)[i]
  
  if (current_seq %in% colnames(seqtab.R.df)) { # is that column present on the second df?
    
    final.seqtab[,current_seq]<-seqtabF[,i] + seqtab.R.df[,current_seq] #if yes, the new df has the sum of reads
    
    seqtab.R.df[,current_seq]<-NULL
    #seqtab.R2.df<- seqtab.R2.df[,-current_seq] #we delete the column from the second df to speed up next search
    
    print (ncol (seqtab.R.df)) #check that is true
    
    
  } else {        # if the column is not present, then the new df is the value of the first df
    
    final.seqtab[,current_seq] <- seqtabF[,i]
    
  }
  
  
  
}
# Now cbind both dataset
final.seqtab <- as.matrix (cbind(final.seqtab,seqtab.R.df))
#TODO: All the above should be in a function, and just call it from the main script
removeBimeraDenovo(test2)
