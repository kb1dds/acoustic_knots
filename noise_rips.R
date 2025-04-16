# Persistence diagrams for the varying noise experiment
library(tidyverse)
library(TDA)

noisedata <- tibble(noiselevel=c(), 
                    death=c(),
                    persistence=c())

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
    
    q <- diagram %>% 
      mutate(dimension=as.factor(dimension)) %>%
      ggplot(aes(x=birth,y=death,color=dimension)) + 
      geom_point() + 
      xlim(0,0.3) + ylim(0,0.3) +
      scale_color_manual(values=c('0'='black','1'='red')) + 
      geom_abline(slope=1,intercept=0)
    
    print(q)
  
    noisedata <- noisedata %>% 
       bind_rows(diagram %>% 
                   filter(dimension==1) %>% 
                   arrange(desc(persistence)) %>% 
                   first() %>% 
                   select(death,persistence) %>%
                   mutate(noiselevel=noiselevel))
  }
}

noisedata %>% 
  filter(noiselevel < 0.025) %>%
  mutate(SNR=10*log10(0.7/(noiselevel*10))) %>%
  pivot_longer(cols=c('death','persistence'),names_to='parameter',values_to='value') %>%
  ggplot(aes(x=SNR,y=value,linetype=parameter)) +
  geom_line() + 
  ylim(0,0.3) +
  xlab('SNR (dB)') +
  ylab('Parameter')

