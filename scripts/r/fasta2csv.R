input <- readLines("/Users/rgallego/Anni_data/banzai_output/derep.fasta")
output <- file("/Users/rgallego/Anni_data/banzai_output/derep.csv","w")
input_path="/Users/rgallego/Anni_data/banzai_output/derep.fasta"
output_path=basename(input_path)
path_1=dirname(input_path)
output_file=paste(path_1,sub(".fasta",".csv", basename(input_path)),sep = "/")
currentSeq <- 0
newLine <- 0

banzai_ouput4<-"/Users/Moncho/Google_Drive/banzai_out_20180301_13_27_test_updated/all_lib/derep.fasta"

fasta2csv<- function (input_path){
  input <- readLines(input_path)
  output_path <- paste(dirname(input_path),sub(".fasta",".csv", basename(input_path)),sep = "/")
  output <-file (output_path, "w")

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
  
}


####Version without a function
input <- readLines("/Users/Moncho/Downloads/derep.fasta")
output <- file("/Users/Moncho/Downloads/derep.csv","w")
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
