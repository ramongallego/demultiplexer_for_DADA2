
library(tidyverse)
library(vegan)
library(gridExtra)
library(forcats)
 
## ----load a function that does the ggplot2 ---- 
MDS.plot.tibble <- function (otu_tibble, otu_field, metadata, transformation="log", distance = "bray"){
  
  if(! transformation %in% c( "log" , "sqrt" ,"") ){
    stop("transformation must be either 'log', 'sqrt' or ''")
  }
  METHODS <- c("manhattan", "euclidean", "canberra", "bray", 
               "kulczynski", "gower", "morisita", "horn", "mountford", 
               "jaccard", "raup", "binomial", "chao", "altGower", "cao", 
               "mahalanobis")
  if (! distance %in% METHODS ){
    stop(paste0("distance must be one of the following: ",paste(METHODS, collapse = ", ")))
  }
  
  #Step 1: turn the tibble into a community matrix
  otu_tibble %>%
    spread(key = otu_field, value = nReads) %>% select (-sample) %>% data.matrix() -> matrix_1
  
  # 2 : NAs to 0s
  matrix_1 [is.na(matrix_1)] <- 0
  
  # 3 : Transform if needed
  if (transformation=="log"){
    matrix_1= log(matrix_1 +1)
  } else {
    if (transformation=="sqrt"){
      matrix_1= sqrt(matrix_1)
    }
  }                        
  
  # 4: MDS plot coordinates
  
  matrix_1  %>% 
    vegdist(method = distance) %>% monoMDS() -> MDS_step_a
  
  MDS_step_b<- as_tibble(cbind(MDS_step_a$points, metadata$sample)) %>% select (sample= V1, MDS1, MDS2)
  
  
  
  # 5 : ggplot - return it as the object so it can be played with later?
  
  ggplot (data=MDS_step_b, aes(x= MDS1, y= MDS2))+
    theme_bw()+
    theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.x = element_blank(),axis.text.y = element_blank(),
          legend.title=element_blank())+
    geom_point(aes(color=metadata$Site),size=5) +
    
    geom_polygon(aes(group=metadata$biol), alpha= 0.2) + coord_equal()
  
  
  
}


## Load datasets: ---- 
### Dataset1: EJP - Banzai ----
Banzai_EJP_metadata=read_csv("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Banzai_Halo/metadata.csv") %>%
  select(sample=sample_id, Site, Position, Distance, Month) %>% separate (col=sample, into = c("biol","rep"), remove = F)

Banzai_EJP_DUP_table= read_delim("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Banzai_Halo/derep.map",
                                 delim = "\t", col_names = c("DUP","ID", "counts")) 

Banzai_EJP_sample_trans= read_delim("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Banzai_Halo/sample_trans.tmp",
                                    delim = "\t", col_names = c("ID","void", "sample")) %>%
  mutate (sample = str_replace(sample,pattern = "sample=", ""))

Banzai_EJP_dup_to_otus<-read_csv("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Banzai_Halo/dups_to_otus.csv")
## Now map reads with sample map

Banzai_EJP_DUP_table %>% inner_join(Banzai_EJP_sample_trans, by = "ID") %>% select(DUP,sample,counts) -> Banzai_EJP_DUP_table ; rm(Banzai_EJP_sample_trans)

# Now collapse them into otus

#1 dups to otus

Banzai_EJP_DUP_table %>%
  inner_join(Banzai_EJP_dup_to_otus, by = c("DUP"="Query" )) %>%
  group_by(Swarm_OTU=Match, sample) %>%
  summarise(nReads = sum(counts)) -> Banzai_EJP_Swarm_otu_table ; rm (Banzai_EJP, Banzai_EJP_DUP_table)

#Remove positive controls if any

Banzai_EJP_Swarm_otu_table <- filter (Banzai_EJP_Swarm_otu_table, !str_detect(sample, "Ostrich"))
Banzai_EJP_metadata <- filter (Banzai_EJP_metadata, !str_detect(sample, "Ostrich"))

###Dataset2: EJP - DADA2 (Use the same metadata as EJP-Banzai) ----

DADA2_EJP_ASV <- read_csv("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Dada2_Halo/ASV_table.csv") %>% filter (!str_detect(sample,"Ostrich"))


### Dataset3: Tides - Banzai ----
Banzai_Tides_metadata=read_csv("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Banzai_Tides/metadata.csv") %>%
  select(sample=sample_name, Site) %>% separate (col=sample, into = c("biol","rep"), remove = F)

Banzai_Tides_DUP_table= read_delim("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Banzai_Tides/derep.map",
                                 delim = "\t", col_names = c("DUP","ID", "counts")) 

Banzai_Tides_sample_trans= read_delim("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Banzai_Tides/sample_trans.tmp",
                                    delim = "\t", col_names = c("ID","void", "sample")) %>%
  mutate (sample = str_replace(sample,pattern = "sample=", ""))

Banzai_Tides_dup_to_otus<-read_csv("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Banzai_Tides/dups_to_otus.csv")  %>%
  separate(Query, c("Query","size"), sep = ";size=") %>% separate(Match, c("Match","size2"), sep = ";size=") 

## Now map reads with sample map

Banzai_Tides_DUP_table %>% inner_join(Banzai_Tides_sample_trans, by = "ID") %>% select(DUP,sample,counts) -> Banzai_Tides_DUP_table ; rm(Banzai_Tides_sample_trans)

# Now collapse them into otus

#1 dups to otus

Banzai_Tides_DUP_table %>%
  inner_join(Banzai_Tides_dup_to_otus, by = c("DUP"="Query" )) %>%
  group_by(Swarm_OTU = Match, sample) %>%
  summarise(nReads = sum(counts)) -> Banzai_Tides_Swarm_otu_table ; rm (Banzai_Tides, Banzai_Tides_DUP_table)

#Remove positive controls if any

Banzai_Tides_Swarm_otu_table <- filter (Banzai_Tides_Swarm_otu_table, !str_detect(sample, "Ostrich"))
Banzai_Tides_metadata <- filter (Banzai_Tides_metadata, !str_detect(sample, "Ostrich"))

###Dataset4: Tides - DADA2 (Use the same metadata as Tides-Banzai) ----

DADA2_Tides_ASV <- read_csv("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Dada2_Tides/ASV_table.csv") %>% filter (!str_detect(sample,"Ostrich"))

### Dataset5: Jul_Aug - Banzai ----
Banzai_Jul_Aug_metadata=read_csv("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Banzai_Jul_Aug/metadata.csv") %>%
  select(sample=sample_id, Site, Month) %>% separate (col=sample, into = c("biol","rep"), sep = "[.]" , remove = F)

Banzai_Jul_Aug_DUP_table= read_delim("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Banzai_Jul_Aug/derep.map",
                                   delim = "\t", col_names = c("DUP","ID", "counts")) 

Banzai_Jul_Aug_sample_trans= read_delim("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Banzai_Jul_Aug/sample_trans.tmp",
                                      delim = "\t", col_names = c("ID","void", "sample")) %>%
  mutate (sample = str_replace(sample,pattern = "sample=", ""))

Banzai_Jul_Aug_dup_to_otus<-read_csv("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Banzai_Jul_Aug/dups_to_otus.csv")  %>%
  separate(Query, c("Query","size"), sep = ";size=") %>% separate(Match, c("Match","size2"), sep = ";size=") 

## Now map reads with sample map

Banzai_Jul_Aug_DUP_table %>% inner_join(Banzai_Jul_Aug_sample_trans, by = "ID") %>% select(DUP,sample,counts) -> Banzai_Jul_Aug_DUP_table ; rm(Banzai_Jul_Aug_sample_trans)

# Now collapse them into otus

#1 dups to otus

Banzai_Jul_Aug_DUP_table %>%
  inner_join(Banzai_Jul_Aug_dup_to_otus, by = c("DUP"="Query" )) %>%
  group_by(Swarm_OTU = Match, sample) %>%
  summarise(nReads = sum(counts)) -> Banzai_Jul_Aug_Swarm_otu_table ; rm ( Banzai_Jul_Aug_DUP_table)

#Remove positive controls if any

Banzai_Jul_Aug_Swarm_otu_table <- filter (Banzai_Jul_Aug_Swarm_otu_table, !str_detect(sample, "Ostrich"))
Banzai_Jul_Aug_metadata <- filter (Banzai_Jul_Aug_metadata, !str_detect(sample, "Ostrich"))

###Dataset6: Jul_Aug - DADA2 (Use the same metadata as Jul_Aug-Banzai) ----

DADA2_Jul_Aug_ASV <- read_csv("/Users/rgallego/Google_Drive/comparing_dada2_banzai/Dada2_Jul_Aug/ASV_table.csv") %>% filter (!str_detect(sample,"Ostrich"))


# Calculate unvariate stats -----
 # Dataset1:
Banzai_EJP_Swarm_otu_table %>%
  group_by(sample) %>% mutate("Number_of_OTUs"  = n_distinct(Swarm_OTU),
                              "Number_of_Reads" = sum(nReads),
                              "Margalef_index" = round((Number_of_OTUs-1)/log(Number_of_Reads),2)) %>%                            
  group_by(sample,Swarm_OTU) %>% mutate ("pi" = nReads/Number_of_Reads) %>%
  group_by(sample) %>% mutate("Shannon_diversity" = -sum(pi * log(pi)),
                              "Pielou_equity" = Shannon_diversity/log(Number_of_OTUs)) %>%
  group_by(sample) %>% select (-Swarm_OTU, - nReads, - pi) %>% summarise_all(mean) -> Banzai_EJP_univariate

# Dataset 2:
DADA2_EJP_ASV %>%
  group_by(sample) %>% mutate("Number_of_OTUs"  = n_distinct(Sequence),
                              "Number_of_Reads" = sum(nReads),
                              "Margalef_index" = round((Number_of_OTUs-1)/log(Number_of_Reads),2)) %>%                            
  group_by(sample,Sequence) %>% mutate ("pi" = nReads/Number_of_Reads) %>%
  group_by(sample) %>% mutate("Shannon_diversity" = -sum(pi * log(pi)),
                              "Pielou_equity" = Shannon_diversity/log(Number_of_OTUs)) %>%
  group_by(sample) %>% select (-Sequence, - nReads, - pi) %>% summarise_all(mean) -> DADA2_EJP_univariate

# Dataset3:
Banzai_Tides_Swarm_otu_table %>%
  group_by(sample) %>% mutate("Number_of_OTUs"  = n_distinct(Swarm_OTU),
                              "Number_of_Reads" = sum(nReads),
                              "Margalef_index" = round((Number_of_OTUs-1)/log(Number_of_Reads),2)) %>%                            
  group_by(sample,Swarm_OTU) %>% mutate ("pi" = nReads/Number_of_Reads) %>%
  group_by(sample) %>% mutate("Shannon_diversity" = -sum(pi * log(pi)),
                              "Pielou_equity" = Shannon_diversity/log(Number_of_OTUs)) %>%
  group_by(sample) %>% select (-Swarm_OTU, - nReads, - pi) %>% summarise_all(mean) -> Banzai_Tides_univariate

# Dataset 4:
DADA2_Tides_ASV %>%
  group_by(sample) %>% mutate("Number_of_OTUs"  = n_distinct(Sequence),
                              "Number_of_Reads" = sum(nReads),
                              "Margalef_index" = round((Number_of_OTUs-1)/log(Number_of_Reads),2)) %>%                            
  group_by(sample,Sequence) %>% mutate ("pi" = nReads/Number_of_Reads) %>%
  group_by(sample) %>% mutate("Shannon_diversity" = -sum(pi * log(pi)),
                              "Pielou_equity" = Shannon_diversity/log(Number_of_OTUs)) %>%
  group_by(sample) %>% select (-Sequence, - nReads, - pi) %>% summarise_all(mean) -> DADA2_Tides_univariate

# Dataset5:
Banzai_Jul_Aug_Swarm_otu_table %>%
  group_by(sample) %>% mutate("Number_of_OTUs"  = n_distinct(Swarm_OTU),
                              "Number_of_Reads" = sum(nReads),
                              "Margalef_index" = round((Number_of_OTUs-1)/log(Number_of_Reads),2)) %>%                            
  group_by(sample,Swarm_OTU) %>% mutate ("pi" = nReads/Number_of_Reads) %>%
  group_by(sample) %>% mutate("Shannon_diversity" = -sum(pi * log(pi)),
                              "Pielou_equity" = Shannon_diversity/log(Number_of_OTUs)) %>%
  group_by(sample) %>% select (-Swarm_OTU, - nReads, - pi) %>% summarise_all(mean) -> Banzai_Jul_Aug_univariate

# Dataset 6:
DADA2_Jul_Aug_ASV %>%
  group_by(sample) %>% mutate("Number_of_OTUs"  = n_distinct(Sequence),
                              "Number_of_Reads" = sum(nReads),
                              "Margalef_index" = round((Number_of_OTUs-1)/log(Number_of_Reads),2)) %>%                            
  group_by(sample,Sequence) %>% mutate ("pi" = nReads/Number_of_Reads) %>%
  group_by(sample) %>% mutate("Shannon_diversity" = -sum(pi * log(pi)),
                              "Pielou_equity" = Shannon_diversity/log(Number_of_OTUs)) %>%
  group_by(sample) %>% select (-Sequence, - nReads, - pi) %>% summarise_all(mean) -> DADA2_Jul_Aug_univariate



# Calculate MDS plots of bc distances on raw, sqroot and log transformed data ----

#
MDS_EJP_Banzai = MDS.plot.tibble(otu_tibble = Banzai_EJP_Swarm_otu_table,
                                 otu_field = "Swarm_OTU",
                                 metadata = Banzai_EJP_metadata,
                                 transformation = "sqrt")
MDS_EJP_DADA2 = MDS.plot.tibble(otu_tibble = DADA2_EJP_ASV,
                                otu_field = "Sequence",
                                metadata = Banzai_EJP_metadata,
                                transformation = "sqrt")
MDS_Tides_Banzai = MDS.plot.tibble(otu_tibble = Banzai_Tides_Swarm_otu_table,
                                 otu_field = "Swarm_OTU",
                                 metadata = Banzai_Tides_metadata,
                                 transformation = "")
MDS_Tides_DADA2 = MDS.plot.tibble(otu_tibble = DADA2_Tides_ASV,
                                otu_field = "Sequence",
                                metadata = Banzai_Tides_metadata,
                                transformation = "")
MDS_Jul_Aug_Banzai = MDS.plot.tibble(otu_tibble = Banzai_Jul_Aug_Swarm_otu_table,
                                   otu_field = "Swarm_OTU",
                                   metadata = Banzai_Jul_Aug_metadata,
                                   transformation = "")
MDS_Jul_Aug_DADA2 = MDS.plot.tibble(otu_tibble = DADA2_Jul_Aug_ASV,
                                  otu_field = "Sequence",
                                  metadata = Banzai_Jul_Aug_metadata,
                                  transformation = "")


# print the plots
grid.arrange(MDS_EJP_Banzai,MDS_EJP_DADA2, nrow=1)
grid.arrange(MDS_Tides_Banzai,MDS_Tides_DADA2, nrow=1)
# Plot a few univariate stats
Univar_to_plot= bind_rows(Banzai_EJP_univariate,DADA2_EJP_univariate, .id = "df") %>% mutate (df = fct_recode(df,"Banzai"="1", "DADA2" ="2")) %>%
  inner_join(Banzai_EJP_metadata[,c("sample", "biol", "Site")], by = "sample") 

Pielou_plot <- ggplot(data = Univar_to_plot,aes(x = biol)) + 
  geom_boxplot(aes(y=Pielou_equity, fill= Site))+
  facet_grid(.~df)

Shannon_plot <- ggplot(data = Univar_to_plot,aes(x = biol)) + 
  geom_boxplot(aes(y=Shannon_diversity, fill= Site))+
  facet_grid(.~df)

grid.arrange(Pielou_plot,Shannon_plot, nrow=2)
