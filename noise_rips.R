# Persistence diagrams for the varying noise experiment
library(tidyverse)
library(TDA)

noisedata <- tibble(noiselevel=c(), 
                    md=c())

for(file in dir('.')){
  if( grepl('pipe_open_echos_',file)){
    noiselevel=as.numeric(substring(file,first=17,last=nchar(file)-4));
    print(noiselevel);
    
    data<-read_csv(file,col_names=FALSE);
    
    pcadata <- data %>% 
      mutate(across(everything(),as.complex)) %>% 
      mutate(across(everything(),abs)) %>% 
      prcomp()
    
    p <- pcadata$x %>% 
      as_tibble() %>% 
      mutate(angle=row_number()) %>% 
      ggplot(aes(PC1,PC2,color=angle)) + 
      geom_path()
    
    print(p)
    
#    diag <- data %>% 
#        mutate(across(everything(),as.complex)) %>% 
#        mutate(across(everything(),abs)) %>%
#        ripsDiag(maxdimension = 2,
#                 maxscale = 0.2)
#
#    diagram <- tibble(dimension=diag$diagram[,1],
#                      birth=diag$diagram[,2],
#                      death=diag$diagram[,3])
#  
#    noisedata <- noisedata %>% 
#       bind_rows(diagram %>% 
#          filter(dimension == 1) %>%
#        summarize(md = max(death)) %>%
#      mutate(noiselevel=noiselevel) %>% 
#      first())
  }
}

noisedata %>% 
  ggplot(aes(x=noiselevel,y=md)) +
  geom_point()
  
