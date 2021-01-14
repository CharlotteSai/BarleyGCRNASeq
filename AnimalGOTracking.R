library(tidyverse)

original <- read_csv("./BarleyAnnos/gene2GO.csv")
animal <- read_csv("./SigGO.csv") %>% 
  filter(str_detect(Description,
                    "development|morphogenesis|cardioblast|cardiocyte|muscle|heart|animal|determination|cell|synapse"))
# overlap with original source
intersect(original$GO_terms,animal$GO_term) %>% length() #42
# total animal related GO
animal$GO_term %>% unique() %>% length() # 61

# this means more than half is coming from original source data