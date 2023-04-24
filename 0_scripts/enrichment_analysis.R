library(stringr)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(msigdbr)

load("~/Dev/mir-4728/4_DE/polysome/polysome_anota2seq.robj")
setwd("~/Dev/mir-4728/5_GSEA")

fcTable = anota2seqGetDeltaData(ads, "full", analysis="translated mRNA", selContrast = 1) %>% 
  as.data.frame() %>% 
  rownames_to_column("geneID")

fcTable$geneID = fcTable$geneID %>% 
  str_remove("\\.[0-9]*$")

entrezIDs=mapIds(org.Hs.eg.db, keys=fcTable$geneID, keytype="ENSEMBL", column="ENTREZID", multiVals="first")

geneList = as.numeric(fcTable[[2]])
names(geneList) = as.character(fcTable$geneID)
geneList = sort(geneList, decreasing=TRUE)

# miRNA target sites

enrichment_miRNA = GSEA(geneList = geneList,
                           TERM2GENE = miRNA_targetsGS,
                           pvalueCutoff = 0.05,
                        minGSSize = 1,
                        maxGSSize = 100000) 

# REACTOME

geneset = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_reactome = GSEA(geneList = geneList,
                           TERM2GENE = geneset,
                           pvalueCutoff = 0.05)

# GO:BP

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_BP = GSEA(geneList = geneList,
                           TERM2GENE = geneset,
                           pvalueCutoff = 0.05)

# GO:MF

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_MF = GSEA(geneList = geneList,
                        TERM2GENE = geneset,
                        pvalueCutoff = 0.05)

# GO:CC

geneset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_GO_CC = GSEA(geneList = geneList,
                        TERM2GENE = geneset,
                        pvalueCutoff = 0.05)

# Transcription factor targets

geneset = msigdbr::msigdbr() %>% 
  filter(gs_subcat == "TFT:GTRD") %>% 
  dplyr::select(c(gs_name, ensembl_gene))

enrichment_TFT = GSEA(geneList = geneList,
                      TERM2GENE = geneset,
                      pvalueCutoff = 0.05)

# KEGG

KEGG_fcTable = fcTable

KEGG_fcTable$geneID = entrezIDs[KEGG_fcTable$geneID]

KEGG_fcTable = KEGG_fcTable %>% drop_na()

geneList = KEGG_fcTable$deltaP
names(geneList) = KEGG_fcTable$geneID
geneList = sort(geneList, decreasing=TRUE)

enrichment_KEGG <- gseKEGG(geneList = geneList,
               organism = 'hsa',
               minGSSize = 120,
               pvalueCutoff = 0.05,
               verbose = FALSE)

enrichment_reactome@result$Description = enrichment_reactome@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("REACTOME ")

enrichment_GO_BP@result$Description = enrichment_GO_BP@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("GOBP ")

enrichment_GO_MF@result$Description = enrichment_GO_MF@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("GOMF ")

enrichment_GO_CC@result$Description = enrichment_GO_CC@result$Description %>%
  str_replace_all("_", " ") %>%
  str_remove_all("GOCC ")

dotplot(enrichment_KEGG, showCategory=15, split=".sign", font.size=7, title="KEGG") +
  facet_grid(.~.sign)

dotplot(enrichment_reactome, showCategory=15, split=".sign", font.size=7, title="REACTOME") +
  facet_grid(.~.sign) + scale_y_discrete(labels=function(x) str_wrap(x, width=70))

dotplot(enrichment_GO_BP, showCategory=15, split=".sign", font.size=7, title="GO:BP") +
  facet_grid(.~.sign) + scale_y_discrete(labels=function(x) str_wrap(x, width=70))

dotplot(enrichment_GO_MF, showCategory=15, split=".sign", font.size=7, title="GO:MF") +
  facet_grid(.~.sign) + scale_y_discrete(labels=function(x) str_wrap(x, width=70))

dotplot(enrichment_GO_CC, showCategory=15, split=".sign", font.size=7, title="GO:CC") +
  facet_grid(.~.sign) + scale_y_discrete(labels=function(x) str_wrap(x, width=70))

dotplot(enrichment_TFT, showCategory=15, split=".sign", font.size=7, title="Transcription factor biding sites") +
  facet_grid(.~.sign) + scale_y_discrete(labels=function(x) str_wrap(x, width=70))


# write.table(enrichment_KEGG, sep="\t", file="enrichment_KEGG.tsv", row.names=FALSE, quote=FALSE)
# write.table(enrichment_reactome, sep="\t", file="enrichment_reactome.tsv", row.names=FALSE, quote=FALSE)