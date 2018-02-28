Banzai_test_1=read_delim("/Users/Moncho/banzai_out_20180227_1612/all_lib/derep.csv",
                         delim = "\t", col_names = c("Hash", "Sequence"))  %>%
  separate(Hash, c("Hash","size"), sep = ";size=") %>%
  mutate (Hash = str_replace(Hash,pattern = ">SHA1=", ""))

Banzai_test_2=read_delim("/Users/Moncho/banzai_out_20180227_161153/all_lib/derep.csv",
                         delim = "\t", col_names = c("Hash", "Sequence"))  %>%
  separate(Hash, c("Hash","size"), sep = ";size=") %>%
  mutate (Hash = str_replace(Hash,pattern = ">SHA1=", ""))

Banzai_test_3=read_delim("/Users/Moncho/Google_Drive/banzai_out_20180227_1643/all_lib/derep.csv",
                         delim = "\t", col_names = c("Hash", "Sequence"))  %>%
  separate(Hash, c("Hash","size"), sep = ";size=") %>%
  mutate (Hash = str_replace(Hash,pattern = ">SHA1=", ""))

Demul_Tides %>% transmute (Seq1= DNAStringSet(Sequence))


JOIN_ASV_Hood_Canal<-bind_rows(Demul_Tides_ASV,Demul_Jul_Aug_ASV)
JOIN_Key_Hood_Canal<- Demul_Jul_Aug %>% bind_rows(Demul_Tides) %>% group_by(Sequence) %>% mutate (nReads = sum(size)) %>% mutate (Percentil = round(size/nReads,2))
write_csv(JOIN_ASV_Hood_Canal, "~/Merged_HC.csv", col_names = T)
