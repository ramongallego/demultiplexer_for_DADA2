#Check whether sequences are shared between runs

#What this script aims for is to load the DUP table - derep.csv or ASV table from three datasets and checks for shared
#sequences between datasets and pipelines

# I run all pipelines with the Hash option set to "YES" so I can compare both unique sequences and hashed sequences

# I have converted the derep.fasta files from banzai into tab files (oddly I named them ".csv", sorry about that)
# To do so, I used the script fasta2csv.r, also attached

library ('tidyverse')

#For each banzai output, I separate the column Hash into hash and number of reads, and then remove the characters not needed

Banzai_Tides=read_delim("/Users/Moncho/Google_Drive/Tide_COI/full/COI_Tides_Miseq_Full_run/full/banzai_out_20170921_1011/all_lib/derep.csv",
                        delim = "\t", col_names = c("Hash", "Sequence"))  %>%
  separate(Hash, c("Hash","size"), sep = ";size=") %>%
  mutate (Hash = str_replace(Hash,pattern = ">SHA1=", ""))

Banzai_EJP=read_delim("/Users/Moncho/Google_Drive/Kelly_Lab/Projects/HALO/Data/banzai_out_20180215_1555_Full_dataset_banzai/all_lib/derep.csv",
                      delim = "\t", col_names = c("Hash", "Sequence"))  %>%
  separate(Hash, c("Hash","size"), sep = ";size=") %>%
  mutate (Hash = str_replace(Hash,pattern = ">SHA1=", ""))

Banzai_Jul_Aug=read_delim("/Users/Moncho/banzai_outputs/banzai_out_0already_paired/all_lib/derep.csv",
                          delim = "\t", col_names = c("Hash", "Sequence")) %>%
  separate(Hash, c("Hash","size"), sep = ";size=") %>%
  mutate (Hash = str_replace(Hash,pattern = ">SHA1=", ""))



#For each dataset run with DADA2, I load the hash key and the ASV table - I'll probably should load the DUP table
# for the other datasets too, if I am to compare community structure as well

# Create a function that does the same operation in all datasets
merdf<- function (hash_key, ASV_table) {
  ASV_table %>% group_by(Hash=Sequence) %>%
    summarise(size =sum(nReads)) %>%
    arrange(desc(size)) %>%
    inner_join(hash_key, by="Hash") -> hash_key
  return (hash_key)
}

Demul_Tides=read_csv("/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180221_1427_Tides_Miseq/hash_key_tides.csv")
Demul_Tides_ASV=read_csv("/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180221_1427_Tides_Miseq/ASV_table.csv")

Demul_test = read_csv("/Users/rgallego/fastqs_demultiplexed_for_DADA2/demultiplexed_20180228_1549/hash_key.csv")
Demul_test_ASV = read_csv("/Users/rgallego/fastqs_demultiplexed_for_DADA2/demultiplexed_20180228_1549/ASV_table.csv")

Demul_Tides_ASV %>%
  group_by(Hash=Sequence) %>%
  summarise(size =sum(nReads)) %>%
  arrange(desc(size)) %>%
  inner_join(Demul_Tides, by="Hash") -> Demul_Tides


Demul_Tides %>% transmute (Seq1= DNAStringSet(Sequence))

Demul_test<-merdf(hash_key = Demul_test, ASV_table = Demul_test_ASV)

Demul_EJP=read_csv("/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180221_1101_EJP_full/hash_key.csv")
Demul_EJP_ASV=read_csv("/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180221_1101_EJP_full/ASV_table.csv")

Demul_EJP_ASV %>%
  group_by(Hash=Sequence) %>%
  summarise(size =sum(nReads)) %>%
  arrange(desc(size)) %>%
  inner_join(Demul_EJP, by="Hash") -> Demul_EJP

Demul_Jul_Aug=read_csv("/Users/Moncho/Test_dada2_joining/demultiplexed_20180222_1157_Run_July_on_top_of_Tides/hash_key.csv")
Demul_Jul_Aug_ASV=read_csv("/Users/Moncho/Test_dada2_joining/demultiplexed_20180222_1157_Run_July_on_top_of_Tides/ASV_table.csv")

Demul_Jul_Aug_ASV %>%
  group_by(Hash=Sequence) %>%
  summarise(size =sum(nReads)) %>%
  arrange(desc(size)) %>%
  inner_join(Demul_Jul_Aug, by="Hash") -> Demul_Jul_Aug

Demul_Jul_Aug_rc=read_csv("/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180228_1137/hash_key.csv")
Demul_Jul_Aug_ASV_rc=read_csv("/Users/Moncho/fastqs_demultiplexed_for_DADA2/demultiplexed_20180228_1137/ASV_table.csv")

Demul_Jul_Aug_ASV_rc %>%
  group_by(Hash=Sequence) %>%
  summarise(size =sum(nReads)) %>%
  arrange(desc(size)) %>%
  inner_join(Demul_Jul_Aug_rc, by="Hash") -> Demul_Jul_Aug_rc

# Now I put all datasets to compare in a list
Array=list(Banzai_Tides,Banzai_EJP,Banzai_Jul_Aug,Demul_Tides,Demul_EJP,Demul_Jul_Aug,Demul_Jul_Aug_rc)
names(Array)=c("Banzai_Tides","Banzai_EJP","Banzai_Jul_Aug","Demul_Tides","Demul_EJP","Demul_Jul_Aug", "Demul_Jul_Aug_rc")

Output_report="/Users/Moncho/Test_dada2_joining/sharedseqs_between_runs.txt"

output_object=tibble(Query_File=character(),
                     Query_Pipeline=character(),
                     Matching_File=character(),
                     Matching_pipeline=character(),
                     nQueries=integer() ,
                     nMatches=integer(),
                     nMatchesHash=integer())

line_now=NULL
k=1

# Now do a loop to perform all pairwise comparisons between datasets

for (i in 1:length(Array)) {

  pipeline_i=str_split(names(Array)[i], pattern="_", simplify = T)[1]
  for (j in (1:length(Array))[-i] ) {

  Array[[i]] %>% inner_join(Array[[j]], by ="Sequence") %>% nrow() ->nmatches

    summary(Array[[i]]$Sequence %in% Array[[j]]$Sequence)

  Array[[i]] %>% inner_join(Array[[j]], by ="Hash") %>% nrow() ->nmatchesHash

    pipeline_j=str_split(names(Array)[j], pattern="_", simplify = T)[1]

  line_now= c(names(Array)[i],pipeline_i,names(Array)[j],pipeline_j,nrow(Array[[i]]),nmatches,nmatchesHash)


  output_object[k,]<-line_now

  k=k+1
  }

}

Shared_seqs<-matrix(0,ncol = length(Array), nrow = length(Array), dimnames = list (names(Array), names(Array)))
Shared_hashes<-matrix(0,ncol = length(Array), nrow = length(Array), dimnames = list (names(Array), names(Array)))

Shared_seqs[with(output_object, cbind(Query_File,Matching_File)) ] <- with(output_object, nMatches)
Shared_hashes[with(output_object, cbind(Query_File,Matching_File)) ] <- with(output_object, nMatchesHash)

Shared_seqs[upper.tri(Shared_seqs)]<-0
Shared_hashes[upper.tri(Shared_hashes)]<-0

write_csv(output_object, Output_report, col_names = T)

JOIN_ASV_Hood_Canal<-bind_rows(Demul_Tides_ASV,Demul_Jul_Aug_ASV)
JOIN_Key_Hood_Canal<- Demul_Jul_Aug %>% bind_rows(Demul_Tides) %>% group_by(Sequence) %>% mutate (nReads = sum(size)) %>% mutate (Percentil = round(size/nReads,2))
write_csv(JOIN_ASV_Hood_Canal, "~/Merged_HC.csv", col_names = T)
