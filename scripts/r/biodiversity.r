#!/usr/bin/env Rscript

msg <-"the script gets executed"

#library (tidyverse)
#library(devtools)
#devtools::install_github("benjjneb/dada2")
#library(dada2)
library(rmarkdown)


print (msg)


arguments <- commandArgs(TRUE)

fastq.folder<-arguments[1]

scripts.folder <-arguments[2]

#list.files()
setwd (scripts.folder)
getwd ()
render("r/biodiversity_report.Rmd", output_file = paste0(fastq.folder,"/biodiversity_report.html"), params = list(folder=fastq.folder))

#print (paste0("folder ",fastq.folder, " file ", arguments[2]))

#list.files(path=fastq.folder, pattern = "Rev.1.fastq")


#sample.map<-read_table(paste0(fastq.folder,"/sample_trans.tmp"),col_names = c("Full_Id", "fastq_header","Sample"))
