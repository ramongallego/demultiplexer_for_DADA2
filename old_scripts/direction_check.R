# An Rscript that bridges between the shell script and the Rmarkdown file

#!/usr/bin/env Rscript


library(rmarkdown)

arguments <- commandArgs(TRUE)

output_folder <- arguments[1]

script_folder <- arguments[2]


render(paste0(script_folder,"/direction_check.Rmd"), output_file = paste0(output_folder,"/direction_report.html"),
 params = list(output_folder = output_folder))
