# Persistence diagrams for the varying noise experiment
library(tidyverse)
library(TDA)

noisedata <- tibble(noiselevel=c(), 
                    md=c())

for(file in dir('.')){
  if( grepl('pipe_open_echos_',file)){
    noiselevel=as.numeric(substring(file,first=17,last=nchar(file)-4));
    print(noiselevel);
    
    data<-read_csv(file,col_names=FALSE) %>% 
      mutate(across(everything(),as.complex)) %>% 
      mutate(across(everything(),abs)) 
  
    sliding <- cbind(data,
                     data[c(25:360,1:24),],
                     data[c(4:360,1:3),])/sqrt(3)

    
    pcadata <- sliding %>% 
      prcomp()
    
    p <- pcadata$x %>% 
      as_tibble() %>% 
      mutate(angle=row_number()) %>% 
      ggplot(aes(PC1,PC2,color=angle)) + 
      geom_path()
    
    print(p)
    
    diag <- sliding %>% 
        ripsDiag(maxdimension = 1,
                 maxscale = 0.7)

    diagram <- tibble(dimension=diag$diagram[,1],
                      birth=diag$diagram[,2],
                      death=diag$diagram[,3]) %>%
      mutate(persistence=death-birth)
  
    noisedata <- noisedata %>% 
       bind_rows(diagram %>% 
          filter(dimension == 1) %>%
        summarize(mp = max(persistence)) %>%
      mutate(noiselevel=noiselevel) %>% 
      first())
  }
}

noisedata %>% 
  filter(noiselevel < 0.025) %>%
  ggplot(aes(x=noiselevel,y=mp)) +
  geom_line() + 
  xlab('Noise level') +
  ylab('Maximum persistence')
  
diagram %>% 
  mutate(dimension=as.factor(dimension)) %>%
  ggplot(aes(x=birth,y=death,color=dimension)) + 
  geom_point()

