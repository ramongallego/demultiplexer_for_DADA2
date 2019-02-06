library(tidyverse)
library(seqinr)

input <- read_csv("/Users/ramongallego/GoogleDrive/Kelly_Lab/Projects/OA_eDNA/Data/Growing_database/Hash_Key_2018-09-13.csv")
output <- "/Users/ramongallego/GoogleDrive/Kelly_Lab/Projects/OA_eDNA/Data/Growing_database/Hash_Key_2018-09-13.fasta"

write.fasta (sequences = as.list(input$Sequence),
             names = as.list(input$Hash),
             file.out = output)
