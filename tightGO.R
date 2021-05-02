#--------- DE gene summary ----------#
# load("~/Box/Bioinfo/BarleyGCRNASeq/DEGs_removeKont_S2.RData")
tightDEs %>% 
  lapply(function(de){
    tibble(n.de = length(de$GeneID),
           n.deWithGO = intersect(de$GeneID, names(gene2GO)) %>% length())
  })
#------------------------------------#


goRes <- tightDEs%>% 
  lapply(function(treatDE){
    treatDE %>% 
      split(f = .$expression) %>% 
      lapply(function(de){
        DEgenes <- de$GeneID %>% intersect(names(gene2GO))
        deGenes2GO <- gene2GO[DEgenes]
        notDEGenes2GO <- gene2GO[setdiff(names(gene2GO), DEgenes)] 
        # genome wild protein coding gene except genes in the DEgenes as background
        nDE <- length(deGenes2GO)
        nNotDE <- length(notDEGenes2GO)
        
        NotdeGO <- unlist(notDEGenes2GO) %>% 
          table %>%
          as.data.frame() %>%
          set_names(c("GO_term", "notDECount")) 
        
        deGO <- unlist(deGenes2GO) %>% 
          table %>%
          as.data.frame() %>%
          set_names(c("GO_term", "DECount")) %>%
          left_join(NotdeGO, by = "GO_term") %>%
          as_tibble() %>%
          filter(DECount > 1) %>% # filter de count more than 1 for each go terms
          mutate(notDECount = ifelse(is.na(notDECount),0,notDECount)) %>% # avoid NA as a given GO associated gene all in DE 
          arrange(GO_term) %>%
          droplevels() 
        
        
        deGO %>%
          split(f = .$GO_term) %>%
          # extract(127) %>%
          lapply(function(df){
            mat <- matrix(c(df$DECount[1], df$notDECount[1],
                            nDE - df$DECount[1], nNotDE - df$notDECount[1]),
                          nrow = 2) %>% 
              set_colnames(c("Genes With GO Term", "Genes Without GO Term")) %>% 
              set_rownames(c("Genes of Interest", "Control Genes"))
            ft <- fisher.test(mat)
            mutate(df, 
                   N = sum(mat[,"Genes With GO Term"]),
                   Expected = nDE * df$notDECount[1] / nNotDE,
                   # number of de expected with proportion calculated in control set
                   p = ft$p.value)
          }) %>%
          bind_rows() %>%
          mutate(adjP = p.adjust(p, "bonferroni"),
                 FDR = p.adjust(p, "fdr"),
                 Description = Term(as.character(GO_term)),
                 expression = unique(de$expression)) %>%
          arrange(p) %>%
          left_join(goSummaries, by = c("GO_term" = "id")) %>%
          filter(DECount > Expected) %>% # consider the fisher test is two-sided
          dplyr::select(GO_term, Description, DECount, Expected, N, everything())
      }) %>% 
      bind_rows()
  })