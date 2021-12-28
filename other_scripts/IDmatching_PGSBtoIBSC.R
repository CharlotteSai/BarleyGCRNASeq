library(tidyverse)
library(magrittr)
library(seqinr)

# calculate embl version barley aa length
aa_IBSC <- read.fasta(file = "BarleyAnnos/Hordeum_vulgare.IBSC_v2.pep.all.fa", 
                      seqtype="AA")
Length_IBSC <- tibble(embl_ID = names(aa_IBSC),
                      embel_Length = getLength(aa_IBSC))

# calculate PGSB version baley AA length
aa_PGSB <- read.fasta(file = "BarleyAnnos/Barley_Morex_V2_gene_annotation_PGSB.all.fasta", 
                      seqtype="AA")
Length_PGSB <- tibble(PGSB_ID = names(aa_PGSB),
                      PGSB_Length = getLength(aa_PGSB))

# load blastp results and attach length info
IBSCmatch <- read_delim("BarleyAnnos/PGSBtoIBSC_Barly.csv", 
                      delim = "\t",
                      comment = "#",
                      col_names = FALSE) %>% 
  set_colnames(c("PGSB_ID", "embl_ID", "% identity", "alignment length", 
                 "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", 
                 "evalue", "bit score")) %>% 
  filter(mismatches == 0 & `gap opens` == 0) %>% 
  left_join(Length_IBSC, by = "embl_ID") %>% 
  left_join(Length_PGSB, by = "PGSB_ID") %>% 
  select(-mismatches, -`gap opens`, -evalue, -`bit score`) %>% 
  mutate(embl_per = `alignment length`/embel_Length,
         PGSB_per = `alignment length`/PGSB_Length)

IBSCmatch %<>% 
  filter(embl_per >= 0.80 & PGSB_per >= 0.80) # 28603

plyr::count(IBSCmatch,'embl_ID') %>% filter(freq >1) %>% arrange(-freq) %>% head()
# HORVU3Hr1G087170.1 -> 43 different matching to PGSB_ID, PGSB aa sequence are identicial