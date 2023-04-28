library(stringr)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(msigdbr)

# load("~/Dev/mir-4728/4_DE/polysome/polysome_anota2seq.robj")
setwd("~/Dev/mir-4728/5_GSEA")

# fcTable_transcription = anota2seqGetDeltaData(ads, "full", analysis="total mRNA", selContrast = 1) %>% 
#   as.data.frame() %>% 
#   rownames_to_column("geneID")
# 
# fcTable_transcription$geneID = fcTable_transcription$geneID %>% 
#   str_remove("\\.[0-9]*$")
# 
# fcTable_translation = anota2seqGetDeltaData(ads, "full", analysis="translation", selContrast = 1) %>% 
#   as.data.frame() %>% 
#   rownames_to_column("geneID")
# 
# fcTable_translation$geneID = fcTable_translation$geneID %>% 
#   str_remove("\\.[0-9]*$")

fcTable_transcription = plot_table %>% 
  dplyr::select(geneID, transcription_fc)

fcTable_translation = plot_table %>% 
  dplyr::select(geneID, translation_fc)

# ========================================== Transcription =================================================
fcTable = fcTable_transcription

entrezIDs=mapIds(org.Hs.eg.db, keys=fcTable$geneID, keytype="ENSEMBL", column="ENTREZID", multiVals="first")

geneList = as.numeric(fcTable[[ncol(fcTable)]])
names(geneList) = as.character(fcTable$geneID)
geneList = sort(geneList, decreasing=TRUE)

# miRNA target sites 

enrichment_miRNA_transcription = GSEA(geneList = geneList,
                                      TERM2GENE = miRNA_targetsGS,
                                      pvalueCutoff = 0.05,
                                      minGSSize = 1,
                                      maxGSSize = 100000) 

# IRES

enrichment_IRES_transcription= GSEA(geneList=geneList,
                                    TERM2GENE=IRES_genes,
                                    pvalueCutoff = 0.05)

# 5'TOP

enrichment_TOP_transcription = GSEA(geneList=geneList,
                                    TERM2GENE=TOP_genes,
                                    pvalueCutoff = 0.05)

# miRNA target / IRES / 5'TOP

enrichment_custom_transcription = GSEA(geneList=geneList,
                                       TERM2GENE=custom_gs,
                                       pvalueCutoff = 0.05)

# REACTOME

geneset = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_reactome_transcription = GSEA(geneList = geneList,
                                         TERM2GENE = geneset,
                                         pvalueCutoff = 0.05)

# GO:BP

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_BP_transcription = GSEA(geneList = geneList,
                                      TERM2GENE = geneset,
                                      pvalueCutoff = 0.05)

# GO:MF

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_MF_transcription = GSEA(geneList = geneList,
                                      TERM2GENE = geneset,
                                      pvalueCutoff = 0.05)

# GO:CC

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_CC_transcription = GSEA(geneList = geneList,
                                      TERM2GENE = geneset,
                                      pvalueCutoff = 0.05)

# Transcription factor targets

geneset = msigdbr::msigdbr() %>% 
  filter(gs_subcat == "TFT:GTRD") %>% 
  dplyr::select(c(gs_name, ensembl_gene))

enrichment_TFT_transcription = GSEA(geneList = geneList,
                                    TERM2GENE = geneset,
                                    pvalueCutoff = 0.05)

KEGG_fcTable = fcTable

KEGG_fcTable$geneID = entrezIDs[KEGG_fcTable$geneID]

KEGG_fcTable = KEGG_fcTable %>% drop_na()

geneList = KEGG_fcTable[[ncol(KEGG_fcTable)]]
names(geneList) = KEGG_fcTable[[1]]
geneList = sort(geneList, decreasing=TRUE)

enrichment_KEGG_transcription <- gseKEGG(geneList = geneList,
                                         organism = 'hsa',
                                         minGSSize = 120,
                                         pvalueCutoff = 0.05,
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
                                    pvalueCutoff = 0.05,
                                    minGSSize = 1,
                                    maxGSSize = 100000) 

# IRES

enrichment_IRES_translation = GSEA(geneList=geneList,
                                   TERM2GENE=IRES_genes,
                                   pvalueCutoff = 0.05)

# 5'TOP

enrichment_TOP_translation = GSEA(geneList=geneList,
                                  TERM2GENE=TOP_genes,
                                  pvalueCutoff = 0.05)

# miRNA target / IRES / 5'TOP

enrichment_custom_translation = GSEA(geneList=geneList,
                                     TERM2GENE=custom_gs,
                                     pvalueCutoff = 0.05)

# REACTOME

geneset = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_reactome_translation = GSEA(geneList = geneList,
                                       TERM2GENE = geneset,
                                       pvalueCutoff = 0.05)

# GO:BP

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_BP_translation = GSEA(geneList = geneList,
                                    TERM2GENE = geneset,
                                    pvalueCutoff = 0.05)

# GO:MF

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_MF_translation = GSEA(geneList = geneList,
                                    TERM2GENE = geneset,
                                    pvalueCutoff = 0.05)

# GO:CC

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_CC_translation = GSEA(geneList = geneList,
                                    TERM2GENE = geneset,
                                    pvalueCutoff = 0.05)

# Transcription factor targets

geneset = msigdbr::msigdbr() %>% 
  filter(gs_subcat == "TFT:GTRD") %>% 
  dplyr::select(c(gs_name, ensembl_gene))

enrichment_TFT_translation = GSEA(geneList = geneList,
                                  TERM2GENE = geneset,
                                  pvalueCutoff = 0.05)

KEGG_fcTable = fcTable

KEGG_fcTable$geneID = entrezIDs[KEGG_fcTable$geneID]

KEGG_fcTable = KEGG_fcTable %>% drop_na()

geneList = KEGG_fcTable[[ncol(KEGG_fcTable)]]
names(geneList) = KEGG_fcTable[[1]]
geneList = sort(geneList, decreasing=TRUE)

enrichment_KEGG_translation <- gseKEGG(geneList = geneList,
                                       organism = 'hsa',
                                       minGSSize = 120,
                                       pvalueCutoff = 0.05,
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

enrichment_TFT_translation@result$Description = enrichment_TFT_translation@result$Description %>% 
  str_replace_all("_", " ")

#write.table(enrichment_KEGG_translation, sep="\t", file="enrichment_KEGG_translation.tsv", row.names=FALSE, quote=FALSE)
#write.table(enrichment_reactome_translation, sep="\t", file="enrichment_reactome_translation.tsv", row.names=FALSE, quote=FALSE)
#write.table(enrichment_GO_BP_translation, sep="\t", file="enrichment_GO_BP_translation.tsv", row.names=FALSE, quote=FALSE)
#write.table(enrichment_GO_MF_translation, sep="\t", file="enrichment_GO_MF_translation.tsv", row.names=FALSE, quote=FALSE)
#write.table(enrichment_GO_CC_translation, sep="\t", file="enrichment_GO_CC_translation.tsv", row.names=FALSE, quote=FALSE)
#write.table(enrichment_TFT_translation, sep="\t", file="enrichment_TFT_translation.tsv", row.names=FALSE, quote=FALSE)
#write.table(enrichment_custom_translation, sep="\t", file="enrichment_IRES_TOP_miRNA_translation.tsv", row.names=FALSE, quote=FALSE)