# ============================================================================================================
# Title: GSEA_results_analysis.R
# Author: Euisuk Robin Han
# Description: A script for analyzing the GSEA results further. Only meant to be used internally. 
# Date: 18/May/23
# ============================================================================================================
library(tidyverse)

polysome = read_tsv("~/Dev/mir-4728/4_DE/polysome/polysome_DESeq.txt")
ribosome = read_tsv("~/Dev/mir-4728/4_DE/ribosome/ribosome_DESeq.txt")
monosome = read_tsv("~/Dev/mir-4728/4_DE/monosome/monosome_DESeq.txt")

polysome_GOBP = read_tsv("~/Dev/mir-4728/5_GSEA/results/DESeq2/polysome/tables/enrichment_GO_BP_transcription.tsv")
monosome_GOBP = read_tsv("~/Dev/mir-4728/5_GSEA/results/DESeq2/monosome/tables/enrichment_GO_BP_transcription.tsv")
ribosome_GOBP = read_tsv("~/Dev/mir-4728/5_GSEA/results/DESeq2/ribosome/tables/enrichment_GO_BP_transcription.tsv")

polysome_inflammatory_genes = polysome_GOBP$core_enrichment[polysome_GOBP$ID=="GOBP_INFLAMMATORY_RESPONSE"] %>% 
  str_split("/") %>% 
  unlist()
monosome_inflammatory_genes = monosome_GOBP$core_enrichment[monosome_GOBP$ID=="GOBP_INFLAMMATORY_RESPONSE"] %>% 
  str_split("/") %>% 
  unlist()
ribosome_inflammatory_genes = ribosome_GOBP$core_enrichment[ribosome_GOBP$ID=="GOBP_INFLAMMATORY_RESPONSE"] %>% 
  str_split("/") %>% 
  unlist()

inflammatory_genes = intersect(monosome_inflammatory_genes, polysome_inflammatory_genes) %>% 
  intersect(ribosome_inflammatory_genes)

polysome_immune_genes = polysome_GOBP$core_enrichment[polysome_GOBP$ID=="GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS"] %>% 
  str_split("/") %>% 
  unlist()
monosome_immune_genes = monosome_GOBP$core_enrichment[monosome_GOBP$ID=="GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS"] %>% 
  str_split("/") %>% 
  unlist()
ribosome_immune_genes = ribosome_GOBP$core_enrichment[ribosome_GOBP$ID=="GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS"] %>% 
  str_split("/") %>% 
  unlist()

immune_genes = intersect(monosome_immune_genes, polysome_immune_genes) %>% 
  intersect(ribosome_immune_genes)

polysome_cytokine_genes = polysome_GOBP$core_enrichment[polysome_GOBP$ID=="GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION"] %>% 
  str_split("/") %>% 
  unlist()
monosome_cytokine_genes = monosome_GOBP$core_enrichment[monosome_GOBP$ID=="GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION"] %>% 
  str_split("/") %>% 
  unlist()
ribosome_cytokine_genes = ribosome_GOBP$core_enrichment[ribosome_GOBP$ID=="GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION"] %>% 
  str_split("/") %>% 
  unlist()

cytokine_genes = intersect(monosome_cytokine_genes, polysome_cytokine_genes) %>% 
  intersect(ribosome_cytokine_genes)

all_genes = inflammatory_genes %>% intersect(immune_genes) #%>% intersect(cytokine_genes)

mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="grch37.ensembl.org")

idMap = getBM(attributes=c('hgnc_symbol',
                           'ensembl_gene_id'),
              filters='ensembl_gene_id',
              values=all_genes,
              mart=mart)

