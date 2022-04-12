library(KEGGREST)
library(tidyverse)
library(magrittr)
library(tibble)
library(pathview)
# based on some info from https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html

# load Hv-At homolog genes and make trans id to gene id 
# (coming from processed result from Germany, miss 96+52 DEs will refill form blastp)
AtHomo <- read_csv("./BarleyAnnos/Barley_Arabidopsis_Homolog.csv") %>% 
  mutate(ID = str_remove(ID, ".\\d$")) # 29515 barley to 14464 arabidopsis

plyr::count(AtHomo, "Ath_besthit") %>% 
  arrange(-freq) %>% 
  head(20)

# load DE genes
# inticial function for extract R object
extractorRData <- function(file, object){
  E <- new.env()
  load(file=file, envir=E)
  return(get(object, envir=E, inherits=F))
}

DEs  <- extractorRData("./DEGs_removeKont_S2.RData", "DEs")

# convert Hv des to At best hit
AtDEs <- DEs %>% 
  lapply(function(x){
    x %<>% 
      left_join(AtHomo, by = c("GeneID" = "ID")) %>% 
      na.omit() # remove Hv did not have At -> NA
    unlist(x$Ath_besthit)
  })

AtDEs %>% 
  lapply(function(de){
    length(unique(de))
  })
# $ABAde
# [1] 2042 -> 1983 -> 1675
# 
# $GABAde
# [1] 2743 -> 2633 -> 2042

# Pull all pathways for AT  
pathways_list <- keggList("pathway", "ath")
keggSummary <- tibble(
  keggID = sub("path:", "", names(pathways_list)),
  Description = str_remove(pathways_list, "- Arabidopsis thaliana \\(thale cress\\)")
)
# Pull all genes for each pathway
pathway_codes <- sub("path:", "", names(pathways_list)) 
pathway2gene <- sapply(pathway_codes,
                        function(pwid){
                          pw <- keggGet(pwid)
                          if (is.null(pw[[1]]$GENE)) return(NA)
                          pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                          pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), 
                                               function(x)x[1]))
                          return(pw2)
                        }
) # 136

gene2pathway <- lapply(names(pathway2gene), function(x){
  tibble(
    pathwayID = x,
    Gene = pathway2gene[x] %>% unlist()
  )
}) %>% 
  bind_rows() %>% 
  na.omit() %>% 
  distinct(pathwayID,Gene) %>% 
  split(f = .$Gene) %>% 
  lapply(function(x){
    x[["pathwayID"]]
  }) # 5118

#DE genes enrich kegg
#--------- DE gene summary ----------#
AtDEs %>% 
  lapply(function(de){
    tibble(n.de = length(de),
           n.de2kegg = intersect(de, names(gene2pathway)) %>% length())
  })
#------------------------------------#
# $ABAde
# # A tibble: 1 x 2
# n.de n.de2kegg
# <int>     <int>
#   1  1983       321
# 
# $GABAde
# # A tibble: 1 x 2
# n.de n.de2kegg
# <int>     <int>
#   1  2633       400


keggRes <- AtDEs %>% 
  lapply(function(des){
    DEgenes <- des %>% intersect(names(gene2pathway))
    deGenes2kegg <- gene2pathway[DEgenes] 
    notDEGenes2kegg <- gene2pathway[setdiff(names(gene2pathway), DEgenes)] 
    # genome wild protein coding gene except genes in the DEgenes as background
    nDE <- length(deGenes2kegg)
    nNotDE <- length(notDEGenes2kegg)
    
    NotdeKEGG <- unlist(notDEGenes2kegg) %>% 
      table %>%
      as.data.frame() %>%
      set_names(c("keggID", "notDECount"))
    
    deKEGG <- unlist(deGenes2kegg) %>% 
      table %>%
      as.data.frame() %>%
      set_names(c("keggID", "DECount")) %>%
      left_join(NotdeKEGG, by = "keggID") %>%
      as_tibble() %>%
      # filter(DECount > 1) %>% # filter de count more than 1 for each go terms
      mutate(notDECount = ifelse(is.na(notDECount),0,notDECount)) %>% # avoid NA as a given GO associated gene all in DE 
      arrange(keggID) %>%
      droplevels() 
    
    deKEGG %>%
      split(f = .$keggID) %>%
      # extract(127) %>%
      lapply(function(df){
        mat <- matrix(c(df$DECount[1], df$notDECount[1],
                        nDE - df$DECount[1], nNotDE - df$notDECount[1]),
                      nrow = 2) %>%
          set_colnames(c("Genes With Kegg", "Genes Without Kegg")) %>%
          set_rownames(c("Genes of Interest", "Control Genes"))
        ft <- fisher.test(mat)
        mutate(df,
               Expected = nDE * df$notDECount[1] / nNotDE,
               # number of de expected with proportion calculated in control set
               GeneRatio = paste0(DECount,"/", nDE),
               BgRatio = paste0(df$notDECount[1], "/", nNotDE),
               p = ft$p.value,
               adjP = p.adjust(p, "bonferroni"),
               FDR = p.adjust(p, "fdr")
               )
      }) %>%
      bind_rows() %>%
      left_join(keggSummary, by = "keggID") %>%
      filter(DECount > Expected) %>%
      dplyr::select(keggID, Description, DECount, Expected, everything())
  })

# prepare filtered data for heatmap
keggResfilter <- names(keggRes) %>% 
  lapply(function(x){
    keggRes[[x]] %>% 
      mutate(Group = str_remove(x, "de")) #%>% 
      # dplyr::select(-notDECount)
  }) %>% 
  bind_rows() %>% 
  group_by(keggID) %>% 
  filter(any(adjP <= 0.05)) %>% 
  arrange(keggID)

# kegg heatmap
keggHeat <- keggResfilter %>%
  mutate(
    p = ifelse(
      adjP <= 0.05, p, 1
    ),
    p = -log10(p)
  ) %>% 
  reshape2::dcast(keggID + Description ~ Group,
                  value.var = "p") %>%
  mutate_all(~replace(., is.na(.), 0))

mapdata <- keggHeat %>% 
  column_to_rownames("keggID") %>% 
  dplyr::select(ABA,GABA) 

label <- keggHeat$Description 
names(label) <- keggHeat$keggID

# myColor <- colorspace::sequential_hcl(30, "Purples 2") %>% rev()
myColor <- colorRampPalette(c("white", "#511B54"))(30)
# myColor <- colorRampPalette(c("white", "#2C73D2"))(30)

pheatmap::pheatmap(mat = mapdata,
                   cluster_cols = FALSE,
                   # labels_row = label,
                   # show_rownames = FALSE,
                   # breaks = myBreaks,
                   color = myColor
) #-> KEGGheat

# export
# export::graph2office(x = KEGGheat, type = "ppt",
#              file = "KEGG_Res.pptx", width = 6,
#              aspectr = 1.4, append = TRUE)

# final kegg Results
keggFinal <- keggResfilter%>% 
  bind_rows() %>% 
  filter(adjP <= 0.05) %>% 
  split(f = .$Group)

# shared pathway
shareKEGG <- intersect(keggFinal[["ABA"]]$keggID,
                       keggFinal[["GABA"]]$keggID)
# [1] "ath00010" "ath00260" "ath00500" "ath00908" "ath00940"
# [6] "ath04016"
ABASpecial <- setdiff(keggFinal[["ABA"]]$keggID,
                       keggFinal[["GABA"]]$keggID)

GABASpecial <- setdiff(keggFinal[["GABA"]]$keggID,
                        keggFinal[["ABA"]]$keggID)

# names(keggRes) %>%
#   lapply(function(x){
#     keggRes[[x]] %>%
#       mutate(Group = str_remove(x, "de"))
#   }) %>%
#   bind_rows() %>% 
#   write_csv("./output_results/KeggResGABA_ABA.csv")

# ----------- about MAPK ----------- #
ABA_mapk <- intersect(pathway2gene$ath04016, AtDEs$ABAde) #20
GABA_mapk <- intersect(pathway2gene$ath04016, AtDEs$GABAde) #18
Hv_mapkShare <- filter(AtHomo, 
                       Ath_besthit %in% 
                         intersect(ABA_mapk,GABA_mapk))
DEs$ABAde %>% filter(GeneID %in% Hv_mapkShare$ID) -> a
DEs$GABAde %>% filter(GeneID %in% Hv_mapkShare$ID) -> g
# overall 9 genes (3 shares)
filter(AtHomo, ID %in% c(a$GeneID,g$GeneID))
# 3 share gene
filter(AtHomo, ID %in% intersect(a$GeneID,g$GeneID)) -> shared_mapkDE
bind_rows(DEs) %>% 
  filter(GeneID %in% shared_mapkDE$ID) %>% 
  arrange(GeneID) %>% 
  left_join(shared_mapkDE, by = c("GeneID"="ID")) %>% 
  select(GeneID,logFC,expression,comparision,
         Ath_besthit,Araport11_defline) %>% 
  arrange(Ath_besthit) %>% 
  mutate(logFC = round(logFC, digits = 4)) %>% 
  reshape2::dcast(GeneID+Ath_besthit+
                    Araport11_defline~comparision,
                  value.var = "logFC") %>% 
  select(GeneID,contains("v.s."),
         Ath_besthit,Araport11_defline) %>% 
  inner_join(go, by = "GeneID") %>% 
  inner_join(goSummaries, 
             by = c("GO_terms"="id")) %>% 
  filter(GeneID == "HORVU.MOREX.r2.7HG0555000" |
           shortest_path >= 4) %>%
  # select(-contains("_path", -terminal_node)) %>% 
  mutate(GO_description = annotate::Term(as.character(GO_terms))) %>% 
  write_csv("mapkShareGene.csv")

# all DE in MAPK
AtHomo %>% 
  filter(Ath_besthit %in% ABA_mapk) %>% 
  filter(ID %in% c(DEs$ABAde$GeneID)) -> Hv_ABA_mapk
AtHomo %>% 
  filter(Ath_besthit %in% GABA_mapk) %>% 
  filter(ID %in% c(DEs$GABAde$GeneID)) -> Hv_GABA_mapk
Hv_mapk <- filter(AtHomo, Ath_besthit %in% 
                    c(ABA_mapk,GABA_mapk)) %>% 
  filter(ID %in% c(DEs$ABAde$GeneID,DEs$GABAde$GeneID))
# prepare for export table
DEs$ABAde %>% filter(GeneID %in% Hv_ABA_mapk$ID) %>% 
  left_join(Hv_ABA_mapk, by = c("GeneID"="ID")) %>% 
  select(GeneID,logFC,comparision,
         Ath_besthit,Araport11_defline) %>% 
  arrange(Ath_besthit) -> Hv_a
DEs$GABAde %>% filter(GeneID %in% Hv_GABA_mapk$ID) %>% 
  left_join(Hv_GABA_mapk, by = c("GeneID"="ID")) %>% 
  select(GeneID,logFC,comparision,
         Ath_besthit,Araport11_defline) %>% 
  arrange(Ath_besthit) -> Hv_g

# rbind(Hv_a,Hv_g) %>% 
#   mutate(logFC = round(logFC, digits = 4)) %>% 
#   arrange(comparision,Ath_besthit) %>% 
#   select(comparision, everything()) %>% 
#   write_csv("mapkGenes.csv")

# ----------- about Carbon fixation in photosynthetic organisms ----------- #
ABA_photo <- intersect(pathway2gene$ath00710, AtDEs$ABAde) #6
GABA_photo <- intersect(pathway2gene$ath00710, AtDEs$GABAde) #11

AtHomo %>% 
  filter(Ath_besthit %in% ABA_photo) %>% 
  filter(ID %in% c(DEs$ABAde$GeneID)) -> Hv_ABA_photo
AtHomo %>% 
  filter(Ath_besthit %in% GABA_photo) %>% 
  filter(ID %in% c(DEs$GABAde$GeneID)) -> Hv_GABA_photo
Hv_photo <- filter(AtHomo, Ath_besthit %in% 
                     c(ABA_photo,GABA_photo)) %>% 
  filter(ID %in% c(DEs$ABAde$GeneID,DEs$GABAde$GeneID))
# prepare for export table
DEs$ABAde %>% filter(GeneID %in% Hv_ABA_photo$ID) %>% 
  left_join(Hv_ABA_photo, by = c("GeneID"="ID")) %>% 
  select(GeneID,logFC,comparision,
         Ath_besthit,Araport11_defline) %>% 
  arrange(Ath_besthit) -> Hv_a
DEs$GABAde %>% filter(GeneID %in% Hv_GABA_photo$ID) %>% 
  left_join(Hv_GABA_photo, by = c("GeneID"="ID")) %>% 
  select(GeneID,logFC,comparision,
         Ath_besthit,Araport11_defline) %>% 
  arrange(Ath_besthit) -> Hv_g

# rbind(Hv_a,Hv_g) %>% 
#   mutate(logFC = round(logFC, digits = 4)) %>% 
#   arrange(comparision,Ath_besthit) %>% 
#   select(comparision, everything()) %>% 
#   write_csv("photoGenes.csv")

