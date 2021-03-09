library(KEGGREST)
library(tidyverse)
library(magrittr)
library(tibble)
# based on some info from https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html

# load Hv-At homolog genes and make trans id to gene id
AtHomo <- read_csv("BarleyAnnos/Barley_Arabidopsis_Homolog.csv") %>% 
  mutate(ID = str_remove(ID, ".\\d$")) # 29515 barley to 14464 arabidopsis

# load DE genes
load("~/Box/Bioinfo/BarleyGCRNASeq/DEGs_removeKont_S2.RData")

# convert Hv des to At best hit
AtDEs <- DEs %>% 
  lapply(function(x){
    x %<>% 
      left_join(AtHomo, by = c("GeneID" = "ID")) %>% 
      na.omit() # remove Hv did not have At -> NA
    unlist(x$Ath_besthit) 
      # unique() # this will lose some for the representation as the basic of 
      # the over-representation test is counting
  })

AtDEs %>% 
  lapply(function(de){
    length(unique(de))
  })
# $ABAde
# [1] 2042 -> 1983 -> 1676
# 
# $GABAde
# [1] 2743 -> 2633 -> 2043

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
  distinct(pathwayID,Gene) %>% 
  split(f = .$Gene) %>% 
  lapply(function(x){
    x[["pathwayID"]]
  }) # 5109

#DE genes enrich kegg
#--------- DE gene summary ----------#
AtDEs %>% 
  lapply(function(de){
    tibble(n.de = length(de),
           n.de2kegg = intersect(de, names(gene2pathway)) %>% length())
  })
#------------------------------------#

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
      filter(DECount > 1) %>% # filter de count more than 1 for each go terms
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
          set_colnames(c("Genes With GO Term", "Genes Without GO Term")) %>%
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
      mutate(Group = str_remove(x, "de")) %>% 
      dplyr::select(-notDECount)
  }) %>% 
  bind_rows() %>% 
  group_by(keggID) %>% 
  filter(any(FDR <= 0.05)) %>% 
  arrange(keggID)

# kegg heatmap
keggHeat <- keggResfilter %>%
  mutate(
    p = ifelse(
      adjP <= 0.05, p, 1
    ),
    p = -log10(p)
  ) 

keggHeat  %<>% 
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
pheatmap::pheatmap(mat = mapdata,
                   cluster_cols = FALSE,
                   # labels_row = label,
                   # show_rownames = FALSE,
                   # breaks = myBreaks,
                   color = myColor
) 

# final kegg Results
keggFinal <- keggResfilter%>% 
  bind_rows() %>% 
  filter(adjP <= 0.05) %>% 
  split(f = .$Group)


# # clusterprofiler results
# KEGGRes <- AtDEs %>%
#   lapply(function(des){
#     clusterProfiler:: enrichKEGG(des, organism = "ath",
#                pAdjustMethod = "fdr")
#   })
# 
# KEGGRes[["ABAde"]]@result %>% 
#   filter(p.adjust <= 0.05)
# KEGGRes[["GABAde"]]@result %>% 
#   filter(p.adjust <= 0.05)
