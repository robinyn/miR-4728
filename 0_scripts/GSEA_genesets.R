# ============================================================================================================
# Title: GSEA_genesets.R
# Author: Euisuk Robin Han
# Description: A script for the preparation of gene sets for the GSEA
# Date: 20/Apr/23
# ============================================================================================================
library(tidyverse)
library(biomaRt)
library(stringr)

setwd("~/Dev/mir-4728/5_GSEA")

# Read gene lists
IRES_genes_RAW = read_tsv("~/Dev/mir-4728/5_GSEA/genesets/IRES_helena.txt")
TOP_genes_RAW = read_tsv("~/Dev/mir-4728/5_GSEA/genesets/5p_TOP_genes.txt", col_names = FALSE)
targetscan = read_tsv("~/Dev/mir-4728/5_GSEA/genesets/targetscan_predictions.txt")
unibind = read_tsv("~/Dev/mir-4728/5_GSEA/genesets/unibind_genes.txt", col_names = FALSE)

# Filter TargetScan predictions to select miR-4728-3p/miR-21-5p targets and filter by prediction reliability
mir_4728 = targetscan %>% 
  filter(`Representative miRNA`=="hsa-miR-4728-3p") %>% 
  filter(as.numeric(`Cumulative weighted context++ score`) < -0.1) %>% 
  filter(as.numeric(`Total num nonconserved sites`) + as.numeric(`Number of 6mer sites`) > 2)

mir_21 = targetscan %>% 
  filter(`Representative miRNA`=="hsa-miR-21-5p") %>% 
  filter(as.numeric(`Cumulative weighted context++ score`) < -0.1) %>% 
  filter(as.numeric(`Total num nonconserved sites`) + as.numeric(`Number of 6mer sites`) >2)

# Remove version numbers from transcript IDs
mir_4728$`Transcript ID` = mir_4728$`Transcript ID` %>% 
  str_replace_all("\\.[0-9]*$", "")

mir_21$`Transcript ID` = mir_21$`Transcript ID` %>% 
  str_replace_all("\\.[0-9]*$", "")

# Convert transcript IDs to gene IDs
mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="grch37.ensembl.org")

idMap = getBM(attributes=c('ensembl_transcript_id',
                           'ensembl_gene_id'),
              filters='ensembl_transcript_id',
              values=unique(c(mir_4728$`Transcript ID`, mir_21$`Transcript ID`)),
              mart=mart)

# Reformat into gene set form required for GSEA
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

# Convert gene symbols in IRES/5'TOP gene sets to gene IDs 
IRES_genes = getBM(attributes=c('hgnc_symbol',
                                'ensembl_gene_id'),
                   filters='hgnc_symbol',
                   values=unique(IRES_genes_RAW$ensembl_gene_id),
                   mart=mart)

TOP_genes = getBM(attributes=c('hgnc_symbol',
                               'ensembl_gene_id'),
                  filters='hgnc_symbol',
                  values=unique(TOP_genes_RAW$X1),
                  mart=mart) 

# Combine miR-4728-3p/miR-21-5p gene sets into one
miRNA_targetsGS = rbind(mir_4728, mir_21)

# Reformat gene sets
TOP_genes = TOP_genes %>% 
  mutate(gs_name="TOP_genes") %>% 
  dplyr::select(gs_name,ensembl_gene_id)

IRES_genes = IRES_genes %>% 
  mutate(gs_name="IRES_genes") %>% 
  dplyr::select(gs_name,ensembl_gene_id)

custom_gs = TOP_genes %>% 
  rbind(IRES_genes)

colnames(unibind) = c("gs_name", "ensemble_gene_id")
unibind$ensemble_gene_id = unibind$ensemble_gene_id %>% 
  str_replace_all("\\.[0-9]*$", "")

unibind_genes = unibind %>% 
  dplyr::group_by(gs_name)


  




