library(tidyverse)
library(ggplot2)
library(stringr)

setwd("~/Dev/mir-4728")

genome_align = read.csv("2_alignments/ribosome_novoalign/novoalign_ribosome_genome_results.tsv", sep="\t")
rRNA_align = read.csv("2_alignments/ribosome_novoalign/novoalign_ribosome_rRNA_results.tsv", sep="\t")

genome_align = genome_align %>% 
  mutate(alignment_rate = ((tot_seq_count - no_count)/tot_seq_count)*100)

rRNA_align = rRNA_align %>% 
  mutate(alignment_rate = ((tot_seq_count - no_count)/tot_seq_count)*100)

dat = genome_align %>% 
  dplyr::select(!alignment_rate)

subset_dat = dat %>% 
  pivot_longer(!sampleID, names_to = "type", values_to = "values") %>% 
  filter(type!="tot_seq_count")

p = ggplot() +
  geom_bar(data=subset_dat, aes(x=sampleID, y=values, fill=type),
           stat="identity", position="stack") +
  geom_point(data=dat, aes(x=sampleID, y=tot_seq_count)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        plot.title=element_text(hjust=0.5))

print(p)