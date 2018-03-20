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

Banzai_test_4_after_update=read_delim("/Users/Moncho/Google_Drive/banzai_out_20180301_1327_test_updated/all_lib/derep.csv",
  delim = ",", col_names = c("Hash", "Sequence"))  %>%
  separate(Hash, c("Hash","size"), sep = ";size=") %>%
  mutate (Hash = str_replace(Hash,pattern = ">SHA1=", ""))

Banzai_test_Ryan=read_delim("/Users/Moncho/Downloads/derep.csv",
                         delim = "\t", col_names = c("Hash", "Sequence"))  %>%
  separate(Hash, c("Hash","size"), sep = ";size=") %>%
  mutate (Hash = str_replace(Hash,pattern = ">SHA1=", ""))

Banzai_test_5_linux=read_csv("/Users/Moncho/Google_Drive/banzai_out_20180228_1126/all_lib/derep.csv",
                                      col_names = c("Hash", "Sequence"))  %>%
  separate(Hash, c("Hash","size"), sep = ";size=") %>%
  mutate (Hash = str_replace(Hash,pattern = ">SHA1=", ""))


Array=list(Banzai_test_1,Banzai_test_2,Banzai_test_3,Banzai_test_Ryan, Banzai_test_4_after_update, Banzai_test_5_linux)
names(Array) = c("Banzai_test_1","Banzai_test_2","Banzai_test_3","Banzai_test_Ryan", "Banzai_test_4_after_update", "Banzai_test_5_linux")

Demul_Tides %>% transmute (Seq1= DNAStringSet(Sequence))

Output_report="/Users/Moncho/Test_dada2_joining/sharedhash_between_runs.csv"

JOIN_ASV_Hood_Canal<-bind_rows(Demul_Tides_ASV,Demul_Jul_Aug_ASV)
JOIN_Key_Hood_Canal<- Demul_Jul_Aug %>% bind_rows(Demul_Tides) %>% group_by(Sequence) %>% mutate (nReads = sum(size)) %>% mutate (Percentil = round(size/nReads,2))
write_csv(JOIN_ASV_Hood_Canal, "~/Merged_HC.csv", col_names = T)
