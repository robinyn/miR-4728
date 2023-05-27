library(tidyverse)
library(ggplot2)
library(stringr)

setwd("~/Dev/mir-4728")

dat = read.csv("2_alignments/polysome_HISAT2_final/parsed_output.tsv", sep="\t")

dat$SampleID=dat$SampleID %>% 
  str_remove("_S[0-9]*$") %>% 
  str_replace_all("aso", "") %>% 
  str_replace_all("_", " ") 

# Overall alignment rate
p = ggplot(dat, aes(x=SampleID, y=Overall.alignment.rate)) + 
  geom_bar(stat="identity", fill="darkgreen") +
  labs(x="Sample", y = "Overall alignment rate (%)", title="Alignment rate") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        plot.title=element_text(hjust=0.5, size=14))

print(p)

# Paired/Unpaired
subset = dat %>% 
  select(SampleID, Total.pairs, Total.unpaired.reads) %>% 
  mutate(Total=Total.pairs+Total.unpaired.reads) %>% 
  mutate(Total.pairs=Total.pairs/Total*100) %>% 
  mutate(Total.unpaired.reads=Total.unpaired.reads/Total*100) %>% 
  select(SampleID, Total.pairs, Total.unpaired.reads) %>% 
  pivot_longer(cols=!SampleID, names_to = "Cat", values_to = "Values")

subset$Cat = subset$Cat %>% 
  str_replace_all("[.]", " ")

p = ggplot(subset, aes(x=SampleID, y=Values, fill=Cat)) +
  geom_bar(position="stack", stat="identity") +
  labs(x="Sample", y = "Percent of reads (%)", title="Paired/Unpaired reads") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        plot.title=element_text(hjust=0.5, size=14))

# Paired/Unpaired counts
subset = dat %>% 
  select(SampleID, Total.pairs, Total.unpaired.reads) %>% 
  pivot_longer(cols=!SampleID, names_to = "Cat", values_to = "Values")

subset$Cat = subset$Cat %>% 
  str_replace_all("[.]", " ")

p = ggplot(subset, aes(x=SampleID, y=Values, fill=Cat)) +
  geom_bar(position="stack", stat="identity") +
  labs(x="Sample", y = "Number of reads", title="Paired/Unpaired reads") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        plot.title=element_text(hjust=0.5, size=14))

print(p)

# Whole summary counts

dat = dat %>% 
  select(!c(Total.pairs, Total.unpaired.reads, Overall.alignment.rate)) %>% 
  pivot_longer(cols=!SampleID, names_to = "Cat", values_to = "Values") %>% 
  mutate(Paired = case_when(str_detect(Cat, "concordantly|discordantly")~"Paired", TRUE~"Unpaired"))

dat$Cat = dat$Cat %>% 
  str_replace_all("\\.\\.", " >") %>% 
  str_replace_all("[.]", " ") %>% 
  str_remove_all("concordantly ") %>% 
  str_remove_all(" or discordantly")

monosome = dat %>% 
  filter(str_detect(SampleID, "monosome")) %>% 
  mutate(SampleID = str_replace(SampleID, "SKBR3 monosome 4728 ", "ASO")) %>% 
  mutate(SampleID = str_replace(SampleID, "SKBR3 monosome ctrl ", "Control"))

polysome = dat %>% 
  filter(str_detect(SampleID, "polysome")) %>% 
  mutate(SampleID = str_replace(SampleID, "SKBR3 polysome 4728 ", "ASO")) %>% 
  mutate(SampleID = str_replace(SampleID, "SKBR3 polysome ctrl ", "Control"))

total = dat %>% 
  filter(str_detect(SampleID, "tot")) %>% 
  mutate(SampleID = str_replace(SampleID, "SKBR3 tot 4728 ", "ASO")) %>% 
  mutate(SampleID = str_replace(SampleID, "SKBR3 tot ctrl ", "Control"))

p = ggplot(total, aes(x=Paired, y=Values, fill=Cat)) +
  geom_bar(stat="identity", position="stack") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        plot.title=element_text(hjust=0.5, size=14)) +
  labs(x = "", y="Number of reads", title="Total RNA") +
  facet_grid(~SampleID)

print(p)

p = ggplot(monosome, aes(x=Paired, y=Values, fill=Cat)) +
  geom_bar(stat="identity", position="stack") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        plot.title=element_text(hjust=0.5, size=14)) +
  labs(x = "", y="Number of reads", title="Monosome fraction") +
  facet_grid(~SampleID)

print(p)

p = ggplot(polysome, aes(x=Paired, y=Values, fill=Cat)) +
  geom_bar(stat="identity", position="stack") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        plot.title=element_text(hjust=0.5, size=14)) +
  labs(x = "", y="Number of reads", title="Polysome fraction") +
  facet_grid(~SampleID)

print(p)