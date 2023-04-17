library(tidyverse)
library(ggplot2)
library(stringr)

setwd("~/Dev/mir-4728")

dat = read.csv("2_alignments/parsed_output.tsv", sep="\t")

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