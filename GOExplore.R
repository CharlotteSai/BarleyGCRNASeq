library(tidyverse)
library(magrittr)
library(GO.db)

# ABA 
# Phosphotransferase, phosphorylation and kinase
photo_tag <- c("GO:0016301", "GO:0016773", "GO:0016310",
               "GO:0006468","GO:0004674")
Term(photo_tag)
lapply(photo_tag, function(t){
  go2Gene[[t]] %>% 
    intersect(DEs$ABAde$GeneID) %>% 
    length()
})

# Hydrogen peroxidemetabolic
