library(stringr)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(msigdbr)

setwd("~/Dev/mir-4728/5_GSEA/results/DESeq2/ribosome/tables")

polysome = read_tsv("~/Dev/mir-4728/4_DE/polysome/polysome_DESeq.txt")
ribosome = read_tsv("~/Dev/mir-4728/4_DE/ribosome/ribosome_DESeq.txt")
monosome = read_tsv("~/Dev/mir-4728/4_DE/monosome/monosome_DESeq.txt")

plot_table = ribosome 

fcTable_transcription = plot_table %>%
  dplyr::select(geneID, totalRNA_lfc)
  
fcTable_translation = plot_table %>% 
  dplyr::select(geneID, TE_lfc)

# ========================================== Transcription =================================================
fcTable = fcTable_transcription

entrezIDs=mapIds(org.Hs.eg.db, keys=fcTable$geneID, keytype="ENSEMBL", column="ENTREZID", multiVals="first")

geneList = as.numeric(fcTable[[ncol(fcTable)]])
names(geneList) = as.character(fcTable$geneID)
geneList = sort(geneList, decreasing=TRUE)

# miRNA target sites 

enrichment_miRNA_transcription = GSEA(geneList = geneList,
                                      TERM2GENE = miRNA_targetsGS,
                                      pvalueCutoff = 1) 

# IRES / 5'TOP

enrichment_custom_transcription = GSEA(geneList=geneList,
                                       TERM2GENE=custom_gs,
                                       pvalueCutoff = 1)

# REACTOME

geneset = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_reactome_transcription = GSEA(geneList = geneList,
                                         TERM2GENE = geneset,
                                         pvalueCutoff = 1)

# GO:BP

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_BP_transcription = GSEA(geneList = geneList,
                                      TERM2GENE = geneset,
                                      pvalueCutoff = 1)

# GO:MF

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_MF_transcription = GSEA(geneList = geneList,
                                      TERM2GENE = geneset,
                                      pvalueCutoff = 1)

# GO:CC

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_CC_transcription = GSEA(geneList = geneList,
                                      TERM2GENE = geneset,
                                      pvalueCutoff = 1)

# Transcription factor targets

geneset = msigdbr::msigdbr() %>% 
  filter(gs_subcat == "TFT:GTRD") %>% 
  dplyr::select(c(gs_name, ensembl_gene))

enrichment_TFT_transcription = GSEA(geneList = geneList,
                                    TERM2GENE = geneset,
                                    pvalueCutoff = 1)

# HALLMARK

geneset = msigdbr::msigdbr() %>% 
  filter(gs_cat == "H") %>% 
  dplyr::select(c(gs_name, ensembl_gene))

enrichment_hallmark_transcription = GSEA(geneList = geneList,
                                    TERM2GENE = geneset,
                                    pvalueCutoff = 1)

# Unibind TFT

geneset = unibind_genes

enrichment_unibind_transcription = GSEA(geneList = geneList,
                                        TERM2GENE = geneset,
                                        pvalueCutoff = 1)

# KEGG

KEGG_fcTable = fcTable

KEGG_fcTable$geneID = entrezIDs[KEGG_fcTable$geneID]

KEGG_fcTable = KEGG_fcTable %>% drop_na()

geneList = KEGG_fcTable[[ncol(KEGG_fcTable)]]
names(geneList) = KEGG_fcTable[[1]]
geneList = sort(geneList, decreasing=TRUE)

enrichment_KEGG_transcription <- gseKEGG(geneList = geneList,
                                         organism = 'hsa',
                                         pvalueCutoff = 1,
                                         verbose = FALSE)

enrichment_reactome_transcription@result$Description = enrichment_reactome_transcription@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("REACTOME ")

enrichment_GO_BP_transcription@result$Description = enrichment_GO_BP_transcription@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("GOBP ")

enrichment_GO_MF_transcription@result$Description = enrichment_GO_MF_transcription@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("GOMF ")

enrichment_GO_CC_transcription@result$Description = enrichment_GO_CC_transcription@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("GOCC ")

enrichment_hallmark_transcription@result$Description = enrichment_hallmark_transcription@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("HALLMARK ")

enrichment_TFT_transcription@result$Description = enrichment_TFT_transcription@result$Description %>% 
  str_replace_all("_", " ")

# ========================================== Translation =================================================
fcTable = fcTable_translation

entrezIDs=mapIds(org.Hs.eg.db, keys=fcTable$geneID, keytype="ENSEMBL", column="ENTREZID", multiVals="first")

geneList = as.numeric(fcTable[[ncol(fcTable)]])
names(geneList) = as.character(fcTable$geneID)
geneList = sort(geneList, decreasing=TRUE)

# miRNA target sites 

enrichment_miRNA_translation = GSEA(geneList = geneList,
                                    TERM2GENE = miRNA_targetsGS,
                                    pvalueCutoff = 1) 

# IRES / 5'TOP

enrichment_custom_translation = GSEA(geneList=geneList,
                                     TERM2GENE=custom_gs,
                                     pvalueCutoff = 1)

# REACTOME

geneset = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_reactome_translation = GSEA(geneList = geneList,
                                       TERM2GENE = geneset,
                                       pvalueCutoff = 1)

# GO:BP

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_BP_translation = GSEA(geneList = geneList,
                                    TERM2GENE = geneset,
                                    pvalueCutoff = 1)

# GO:MF

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_MF_translation = GSEA(geneList = geneList,
                                    TERM2GENE = geneset,
                                    pvalueCutoff = 1)

# GO:CC

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_CC_translation = GSEA(geneList = geneList,
                                    TERM2GENE = geneset,
                                    pvalueCutoff = 1)

# Transcription factor targets

geneset = msigdbr::msigdbr() %>% 
  filter(gs_subcat == "TFT:GTRD") %>% 
  dplyr::select(c(gs_name, ensembl_gene))

enrichment_TFT_translation = GSEA(geneList = geneList,
                                  TERM2GENE = geneset,
                                  pvalueCutoff = 1)

# HALLMARK

geneset = msigdbr::msigdbr() %>% 
  filter(gs_cat == "H") %>% 
  dplyr::select(c(gs_name, ensembl_gene))

enrichment_hallmark_translation = GSEA(geneList = geneList,
                                         TERM2GENE = geneset,
                                         pvalueCutoff = 1)

# Unibind TFT

geneset = unibind_genes

enrichment_unibind_translation = GSEA(geneList = geneList,
                                        TERM2GENE = geneset,
                                        pvalueCutoff = 1)

# KEGG

KEGG_fcTable = fcTable

KEGG_fcTable$geneID = entrezIDs[KEGG_fcTable$geneID]

KEGG_fcTable = KEGG_fcTable %>% drop_na()

geneList = KEGG_fcTable[[ncol(KEGG_fcTable)]]
names(geneList) = KEGG_fcTable[[1]]
geneList = sort(geneList, decreasing=TRUE)

enrichment_KEGG_translation <- gseKEGG(geneList = geneList,
                                       organism = 'hsa',
                                       pvalueCutoff = 1,
                                       verbose = FALSE)

enrichment_reactome_translation@result$Description = enrichment_reactome_translation@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("REACTOME ")

enrichment_GO_BP_translation@result$Description = enrichment_GO_BP_translation@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("GOBP ")

enrichment_GO_MF_translation@result$Description = enrichment_GO_MF_translation@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("GOMF ")

enrichment_GO_CC_translation@result$Description = enrichment_GO_CC_translation@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("GOCC ")

enrichment_hallmark_translation@result$Description = enrichment_hallmark_translation@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("HALLMARK ")

enrichment_TFT_translation@result$Description = enrichment_TFT_translation@result$Description %>% 
  str_replace_all("_", " ")

gseaplot(enrichment_miRNA_translation, geneSetID = "hsa-miR-4728-3p")
gseaplot(enrichment_miRNA_transcription, geneSetID = "hsa-miR-4728-3p")
gseaplot(enrichment_miRNA_translation, geneSetID = "hsa-miR-21-5p")
gseaplot(enrichment_miRNA_transcription, geneSetID = "hsa-miR-21-5p")

write.table(enrichment_KEGG_transcription, sep="\t", file="enrichment_KEGG_transcription.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_reactome_transcription, sep="\t", file="enrichment_reactome_transcription.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_GO_BP_transcription, sep="\t", file="enrichment_GO_BP_transcription.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_GO_MF_transcription, sep="\t", file="enrichment_GO_MF_transcription.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_GO_CC_transcription, sep="\t", file="enrichment_GO_CC_transcription.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_TFT_transcription, sep="\t", file="enrichment_TFT_transcription.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_unibind_transcription, sep="\t", file="enrichment_unibind_transcription.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_miRNA_transcription, sep="\t", file="enrichment_miRNA_transcription.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_custom_transcription, sep="\t", file="enrichment_IRES_TOP_transcription.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_hallmark_transcription, sep="\t", file="enrichment_HALLMARK_transcription.tsv", row.names=FALSE, quote=FALSE)

write.table(enrichment_KEGG_translation, sep="\t", file="enrichment_KEGG_translation.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_reactome_translation, sep="\t", file="enrichment_reactome_translation.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_GO_BP_translation, sep="\t", file="enrichment_GO_BP_translation.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_GO_MF_translation, sep="\t", file="enrichment_GO_MF_translation.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_GO_CC_translation, sep="\t", file="enrichment_GO_CC_translation.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_TFT_translation, sep="\t", file="enrichment_TFT_translation.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_unibind_translation, sep="\t", file="enrichment_unibind_translation.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_miRNA_translation, sep="\t", file="enrichment_miRNA_translation.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_custom_translation, sep="\t", file="enrichment_IRES_TOP_translation.tsv", row.names=FALSE, quote=FALSE)
write.table(enrichment_hallmark_translation, sep="\t", file="enrichment_HALLMARK_translation.tsv", row.names=FALSE, quote=FALSE)