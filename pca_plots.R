# Plotting PCA of the signature space in a standardized way

library(tidyverse)
library(TDA)

target_name <- c('coke_bottle','cup_open','cup_capped','pipe_open','pipe_capped')

for(target in target_name ){
  data<-read_csv(paste0(target,'_echos.csv'),col_names=FALSE)

  pcadata <- data %>% 
    mutate(across(everything(),as.complex)) %>% 
    mutate(across(everything(),abs)) %>% 
    prcomp()

  pcadata$x %>% 
    as_tibble() %>% 
    mutate(angle=row_number()) %>% 
    ggplot(aes(PC1,PC2,color=angle)) + 
    geom_path()

  ggsave(paste0(target,'_signature_pca.png'))
  
#  data %>% 
#    mutate(across(everything(),as.complex)) %>% 
#    mutate(across(everything(),abs)) %>%
#    ripsDiag(maxdimension = 2,
#             maxscale = 1.0)
}

