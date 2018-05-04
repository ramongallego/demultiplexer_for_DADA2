#A function to do from tibble to MDS plot

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