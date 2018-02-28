input <- readLines("/Users/Moncho/Google_Drive/banzai_out_20180227_1643/all_lib/derep.fasta")
output <- file("/Users/Moncho/Google_Drive/banzai_out_20180227_1643/all_lib/derep.csv","w")

currentSeq <- 0
newLine <- 0

for(i in 1:length(input)) {
  if(strtrim(input[i], 1) == ">") {
    if(currentSeq == 0) {
      writeLines(paste(input[i],"\t"), output, sep="")
      currentSeq <- currentSeq + 1
    } else {
      writeLines(paste("\n",input[i],"\t", sep =""), output, sep="")
    }
  } else {
    writeLines(paste(input[i]), output, sep="")
  }
}

close(output)
