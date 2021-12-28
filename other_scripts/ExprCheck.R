### This quick gene expression check is based on the DGE list ###
library(tidyverse)
library(magrittr)

ExpData <- countList %>% 
  cpm() %>% 
  as.data.frame()

ExpData <- log2(ExpData+1) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(AtHomo, by = c("GeneID"="ID")) %>%
  reshape2::melt(id.vars = c("GeneID", "Ath_besthit",
                             "Araport11_defline"),
                 variable.name = "Sample",
                 value.name = "Counts") %>% 
  mutate(treat = str_extract(Sample, "Kont|^ABA|^GABA"))

# find GLR3.3 (AT1G42540)
genes <- AtHomo %>% filter(Ath_besthit == "AT1G42540") %>% pull(ID)
intersect(DEs$GABAde$GeneID, genes)

ExpData %>% 
  filter(GeneID %in% genes# & Genotype != "gad1SALK") %>% 
  ) %>%
  ggplot(aes(y = Counts, 
             x = treat)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               position = position_dodge(1)) +
  theme_bw() +
  theme(text = element_text(size=16)) +
  labs(x = "",
       y = "log2(cpm+1)") +
  facet_wrap("GeneID", scales = "free_y") #-> geneExp

# export::graph2office(x = geneExp, type = "ppt",
#              file = "geneExp.pptx", width = 8,
#              aspectr = 1.5, append = TRUE)

