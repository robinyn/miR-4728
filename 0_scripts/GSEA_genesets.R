library(tidyverse)
library(biomaRt)
library(stringr)

setwd("~/Dev/mir-4728/5_GSEA")

targetscan = read_tsv("targetscan_predictions.txt")

mir_4728 = targetscan %>% 
  filter(`Representative miRNA`=="hsa-miR-4728-3p") %>% 
  filter(as.numeric(`Cumulative weighted context++ score`) < -0.1)

mir_21 = targetscan %>% 
  filter(`Representative miRNA`=="hsa-miR-21-5p") %>% 
  filter(as.numeric(`Cumulative weighted context++ score`) < -0.1)

mir_4728$`Transcript ID` = mir_4728$`Transcript ID` %>% 
  str_replace_all("\\.[0-9]*", "")

mir_21$`Transcript ID` = mir_21$`Transcript ID` %>% 
  str_replace_all("\\.[0-9]*", "")

mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="grch37.ensembl.org")

idMap = getBM(attributes=c('ensembl_transcript_id',
                           'ensembl_gene_id'),
              filters='ensembl_transcript_id',
              values=unique(c(mir_4728$`Transcript ID`, mir_21$`Transcript ID`)),
              mart=mart)

mir_4728 = mir_4728 %>% 
  mutate(transcriptID = `Transcript ID`) %>% 
  mutate(gs_name = `Representative miRNA`) %>% 
  dplyr::select(!c(`Transcript ID`, `Representative miRNA`)) %>% 
  left_join(idMap, by=c("transcriptID"="ensembl_transcript_id")) %>% 
  dplyr::select(!transcriptID) %>% 
  mutate(ensembl_gene = ensembl_gene_id) %>% 
  dplyr::select(!ensembl_gene_id) %>% 
  arrange(desc(`Number of 6mer sites`)) %>% 
  top_n(500) %>% 
  dplyr::select(gs_name, ensembl_gene)

mir_21 = mir_21 %>% 
  mutate(transcriptID = `Transcript ID`) %>% 
  mutate(gs_name = `Representative miRNA`) %>% 
  dplyr::select(!c(`Transcript ID`, `Representative miRNA`)) %>% 
  left_join(idMap, by=c("transcriptID"="ensembl_transcript_id")) %>% 
  dplyr::select(!transcriptID) %>% 
  mutate(ensembl_gene = ensembl_gene_id) %>% 
  dplyr::select(!ensembl_gene_id) %>% 
  arrange(desc(`Number of 6mer sites`)) %>% 
  top_n(500) %>% 
  dplyr::select(gs_name, ensembl_gene)

miRNA_targetsGS = rbind(mir_4728, mir_21)

transcription_factor_targets <- msigdbr::msigdbr() %>% 
  filter(gs_subcat == "TFT:GTRD") %>% 
  dplyr::select(c(gs_name, ensembl_gene))


