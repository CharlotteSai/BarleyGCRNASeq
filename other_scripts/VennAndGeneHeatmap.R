library(tidyverse)
library(magrittr)

# inticial function for extract R object
extractorRData <- function(file, object){
  E <- new.env()
  load(file=file, envir=E)
  return(get(object, envir=E, inherits=F))
}

ABAgenes <- extractorRData("./DEGs_removeKont_S2.RData", "ABAgenes")
GABAgenes <- extractorRData("./DEGs_removeKont_S2.RData", "GABAgenes")


ABAtags <- ABAgenes %>% 
  mutate(ABA = ifelse(adj.p <= 0.05 & abs(logFC) >= 1, 
                      1,0)) %>% 
  dplyr::select(GeneID,ABA,logFC) %>% 
  data.table::setnames("logFC","ABA_logFC") 
GABAtags <- GABAgenes %>% 
  mutate(GABA = ifelse(adj.p <= 0.05 & abs(logFC) >= 1, 
                       1,0)) %>% 
  dplyr::select(GeneID,GABA,logFC) %>% 
  data.table::setnames("logFC","GABA_logFC")

# vennDiagram 
merge(ABAtags,
      GABAtags,
      by = "GeneID") %>%
  dplyr::select(-contains("logFC")) %>% 
  column_to_rownames("GeneID") %>%
  vennDiagram(circle.col=c("turquoise", "salmon"))

# common genes
patternDE <- merge(ABAtags,
                   GABAtags,
                   by = "GeneID") %>% 
  mutate(ABA = ifelse(ABA == 1 & ABA_logFC < 0,
                      -1, ABA),
         GABA = ifelse(GABA == 1 & GABA_logFC < 0,
                       -1, GABA)) %>% 
  column_to_rownames("GeneID") %>% 
  filter(abs(ABA) == 1 | abs(GABA) == 1)


# same direction
shareSame <- patternDE %>% 
  filter((ABA == 1 & GABA == 1) | (ABA == -1 & GABA == -1)) %>% 
  row.names() # 428

# opposite direction
shareOps <- patternDE %>% 
  filter((ABA == 1 & GABA == -1) | (ABA == -1 & GABA == 1)) %>% 
  row.names() # 111


DEheat <- patternDE %>% 
  mutate(ABA_logFC = ifelse(ABA == 0, 0, ABA_logFC),
         GABA_logFC = ifelse(GABA == 0, 0, GABA_logFC)) %>% 
  dplyr::select(-GABA, -ABA) 

paletteLength <- 100
myBreaks <- c(seq(min(DEheat), 0, 
                  length.out = ceiling(paletteLength/2)), 
              seq(max(DEheat)/paletteLength, max(DEheat), length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

DEheat %>% 
  pheatmap::pheatmap(cluster_cols = FALSE,
                     cluster_rows = TRUE,
                     show_rownames = FALSE,
                     color = myColor,
                     breaks = myBreaks) 


