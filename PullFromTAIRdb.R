library(tidyverse)
library(magrittr)
library(tibble)
library(org.At.tair.db)

path2tair <- as.list(org.At.tairPATH2TAIR) #123
t2p <- lapply(names(path2tair), function(x){
  tibble(
    pathwayID = x,
    Gene = path2tair[x] %>% unlist()
  )
}) %>% 
  bind_rows() %>% 
  distinct(pathwayID,Gene) %>% 
  split(f = .$Gene) %>% 
  lapply(function(x){
    x[["pathwayID"]]
  }) #3450



x <- org.At.tairPATH
mapped_gene <-  mappedkeys(x)
tair2path <- as.list(x[mapped_gene]) #3450
