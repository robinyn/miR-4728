library(tidyverse)
library(biomaRt)
library(stringr)

setwd("~/Dev/mir-4728/5_GSEA")

# Read gene list
IRES_genes_RAW = read_tsv("~/Dev/mir-4728/5_GSEA/genesets/human_IRES_info.txt")
TOP_genes_RAW = read_tsv("~/Dev/mir-4728/5_GSEA/genesets/5p_TOP_genes.txt", col_names = FALSE)
targetscan = read_tsv("~/Dev/mir-4728/5_GSEA/genesets/targetscan_predictions.txt")

mir_4728 = targetscan %>% 
  filter(`Representative miRNA`=="hsa-miR-4728-3p") %>% 
  filter(as.numeric(`Cumulative weighted context++ score`) < -0.4)

mir_21 = targetscan %>% 
  filter(`Representative miRNA`=="hsa-miR-21-5p") %>% 
  filter(as.numeric(`Cumulative weighted context++ score`) < -0.4)

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
  #arrange(desc(`Number of 6mer sites`)) %>% 
  #top_n(500) %>% 
  dplyr::select(gs_name, ensembl_gene_id)

mir_21 = mir_21 %>% 
  mutate(transcriptID = `Transcript ID`) %>% 
  mutate(gs_name = `Representative miRNA`) %>% 
  dplyr::select(!c(`Transcript ID`, `Representative miRNA`)) %>% 
  left_join(idMap, by=c("transcriptID"="ensembl_transcript_id")) %>% 
  dplyr::select(!transcriptID) %>% 
  #arrange(desc(`Number of 6mer sites`)) %>% 
  #top_n(500) %>% 
  dplyr::select(gs_name, ensembl_gene_id)

IRES_genes = getBM(attributes=c('hgnc_symbol',
                                'ensembl_gene_id'),
                   filters='hgnc_symbol',
                   values=unique(IRES_genes_RAW$`Gene symbol`),
                   mart=mart) 

TOP_genes = getBM(attributes=c('hgnc_symbol',
                               'ensembl_gene_id'),
                  filters='hgnc_symbol',
                  values=unique(TOP_genes_RAW$X1),
                  mart=mart) 

miRNA_targetsGS = rbind(mir_4728, mir_21)

# Reformat gene sets
TOP_genes = TOP_genes %>% 
  mutate(gs_name="TOP_genes") %>% 
  dplyr::select(gs_name,ensembl_gene_id)

IRES_genes = IRES_genes %>% 
  mutate(gs_name="IRES_genes") %>% 
  dplyr::select(gs_name,ensembl_gene_id)

custom_gs = TOP_genes %>% 
  rbind(IRES_genes) %>% 
  rbind(miRNA_targetsGS)




