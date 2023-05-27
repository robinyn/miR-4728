# ============================================================================================================
# Title: TENT4B_PDCD4.R
# Author: Euisuk Robin Han
# Description: A script for the visualization of DGE analysis results for TENT4B and PDCD4. For internal use only
# Date: 8/May/23
# ============================================================================================================
library(stringr)
library(ggplot2)
library(tidyverse)

setwd("~/Dev/mir-4728/4_DE")

# Read DGE analysis results
polysome = read_tsv("polysome/polysome_DESeq.txt")
monosome = read_tsv("monosome/monosome_DESeq.txt")
ribosome = read_tsv("ribosome/ribosome_DESeq.txt")

# ** The script currently only plots polysome and ribosome data **

# Ensembl gene IDs for TENT4B and PDCD4
TENT4B_ensembl = "ENSG00000121274"
PDCD4_ensembl = "ENSG00000150593"

# Find DGE results for TENT4B and PDCD4
polysome_dat = polysome %>% 
  filter(geneID==TENT4B_ensembl |
           geneID==PDCD4_ensembl)
ribosome_dat = ribosome %>% 
  filter(geneID==TENT4B_ensembl |
           geneID==PDCD4_ensembl) 

polysome_lfc = polysome_dat %>% 
  dplyr::select(geneID, TE_lfc, totalRNA_lfc, polysome_lfc, mode) %>% 
  pivot_longer(cols = !c(geneID, mode),
               names_sep = "_",
               names_to = c("type", "discard"),
               values_to = "lfc") %>% 
  dplyr::select(!discard)

polysome_padj = polysome_dat %>% 
  dplyr::select(geneID, TE_padj, totalRNA_padj, polysome_padj) %>% 
  pivot_longer(cols = !c(geneID),
               names_sep = "_",
               names_to = c("type", "discard"),
               values_to = "padj") %>% 
  dplyr::select(!discard)

polysome_dat = polysome_lfc %>% 
  merge(polysome_padj, by=c("geneID", "type"))

polysome_dat$geneID[polysome_dat$geneID=="ENSG00000121274"] = "TENT4B"
polysome_dat$geneID[polysome_dat$geneID=="ENSG00000150593"] = "PDCD4"

ribosome_lfc = ribosome_dat %>% 
  dplyr::select(geneID, TE_lfc, totalRNA_lfc, polysome_lfc, mode) %>% 
  pivot_longer(cols = !c(geneID, mode),
               names_sep = "_",
               names_to = c("type", "discard"),
               values_to = "lfc") %>% 
  dplyr::select(!discard)

ribosome_padj = ribosome_dat %>% 
  dplyr::select(geneID, TE_padj, totalRNA_padj, polysome_padj) %>% 
  pivot_longer(cols = !c(geneID),
               names_sep = "_",
               names_to = c("type", "discard"),
               values_to = "padj") %>% 
  dplyr::select(!discard)

ribosome_dat = ribosome_lfc %>% 
  merge(ribosome_padj, by=c("geneID", "type"))

ribosome_dat$geneID[ribosome_dat$geneID=="ENSG00000121274"] = "TENT4B"
ribosome_dat$geneID[ribosome_dat$geneID=="ENSG00000150593"] = "PDCD4"

ribosome_dat$type[ribosome_dat$type=="polysome"] = "RPF"

polysome_anot = polysome_dat %>% 
  dplyr::select(geneID, mode) %>% 
  unique()

ribosome_anot = ribosome_dat %>% 
  dplyr::select(geneID, mode) %>% 
  unique()

# Plot data
p = ggplot() +
  geom_bar(data=polysome_dat, stat="identity", aes(x=type, y=lfc, fill=padj)) +
  scale_fill_gradient(low="red", high="grey", guide="colorbar") +
  geom_label(data=polysome_anot, aes(x=3, y=1, label=mode)) +
  facet_wrap(vars(geneID)) +
  labs(x="", y="Log2FC", title="Monosome")

print(p)

p = ggplot() +
  geom_bar(data=ribosome_dat, stat="identity", aes(x=type, y=lfc, fill=padj)) +
  scale_fill_gradient(low="red", high="grey", guide="colorbar") +
  geom_label(data=ribosome_anot, aes(x=3, y=1, label=mode)) +
  facet_wrap(vars(geneID)) +
  labs(x="", y="Log2FC", title="Ribosome")

print(p)