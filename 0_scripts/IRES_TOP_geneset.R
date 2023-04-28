library(tidyverse)
library(biomaRt)
library(stringr)

# Read gene list
IRES_genes_RAW = read_tsv("~/Dev/mir-4728/5_GSEA/genesets/human_IRES_info.txt")
TOP_genes_RAW = read_tsv("~/Dev/mir-4728/5_GSEA/genesets/5p_TOP_genes.txt", col_names = FALSE)

# Retrieve corresponding ensembl IDs
mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

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

