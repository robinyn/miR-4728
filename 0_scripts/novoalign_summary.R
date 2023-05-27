# ============================================================================================================
# Title: novoalign_summary.R
# Author: Euisuk Robin Han
# Description: A script for the visualization of alignment results from Novoalign
# Date: 17/Apr/23
# ============================================================================================================
library(tidyverse)
library(ggplot2)
library(stringr)

setwd("~/Dev/mir-4728")

# Read the parsed alignment summary
genome_align = read.csv("2_alignments/ribosome_novoalign/novoalign_ribosome_genome_results.tsv", sep="\t")
rRNA_align = read.csv("2_alignments/ribosome_novoalign/novoalign_ribosome_rRNA_results.tsv", sep="\t")

# Calculate total alignment rate for each samples
genome_align = genome_align %>% 
  mutate(alignment_rate = ((tot_seq_count - no_count)/tot_seq_count)*100)

rRNA_align = rRNA_align %>% 
  mutate(alignment_rate = ((tot_seq_count - no_count)/tot_seq_count)*100)

# Select data to plot
dat = genome_align %>% 
  dplyr::select(!alignment_rate)

# Reformat sample names to make them more legible
dat$sampleID = dat$sampleID %>% 
  str_replace_all("_", " ")

# Pivot data for plotting
subset_dat = dat %>% 
  pivot_longer(!sampleID, names_to = "type", values_to = "values") %>% 
  filter(type!="tot_seq_count")

# Reformat alignment categories to be more intuitive
subset_dat$type = subset_dat$type %>% 
  str_replace_all("multi_count", "Aligned >1 times") %>% 
  str_replace_all("uniq_count", "Aligned 1 time") %>% 
  str_replace_all("no_count", "Aligned 0 time")

# Generate plot
p = ggplot() +
  geom_bar(data=subset_dat, aes(x=sampleID, y=values, fill=type),
           stat="identity", position="stack") +
  geom_point(data=dat, aes(x=sampleID, y=tot_seq_count)) +
  labs(x="Sample", y="Number of reads") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        plot.title=element_text(hjust=0.5, size=14))

print(p)