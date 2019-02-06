input <- readLines("/Users/rgallego/fastqs_demultiplexed_for_DADA2/banzai_out_20180207_1704/all_lib/derep.fasta")
output <- file("/Users/rgallego/fastqs_demultiplexed_for_DADA2/banzai_out_20180207_1704/all_lib/derep.csv","w")

currentSeq <- 0
newLine <- 0

for(i in 1:length(input)) {
  if(strtrim(input[i], 1) == ">") {
    if(currentSeq == 0) {
      writeLines(paste(input[i],","), output, sep="")
      currentSeq <- currentSeq + 1
    } else {
      writeLines(paste("\n",input[i],",", sep =""), output, sep="")
    }
  } else {
    writeLines(paste(input[i]), output, sep="")
  }
}

close(output)
