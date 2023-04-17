library(tidyverse)

setwd("~/Dev/mir-4728")

dat = read.csv("1_fastqc/polysome_summary.tsv", sep="\t")

dat$sample_name = dat$sample_name  %>%
  str_remove("_S[0-9]*") %>%
  str_replace_all("001.fastq.gz", "") %>%
  str_replace_all("_", " ")

dat = dat %>%
  mutate(fraction=case_when(grepl("monosome", sample_name) ~ "monosome",
                            grepl("polysome", sample_name) ~ "polysome",
                            grepl("tot", sample_name) ~ "total")) %>%
  mutate(fraction=as.factor(fraction))

# dat$sample_name = dat$sample_name %>%
#   str_remove(".1.fastq.gz") %>%
#   str_replace_all("_", " ")
# 
# dat = dat %>% 
#   mutate(fraction=case_when(grepl("RPF", sample_name) ~ "RPF",
#                             grepl("total", sample_name) ~ "total")) %>% 
#   mutate(fraction=as.factor(fraction))

p = ggplot(dat, aes(x=sample_name, y=total_read_count, fill=fraction)) +
  geom_bar(stat="identity") +
  labs(x="Sample", y = "Number of reads", title="Total number of reads") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        plot.title=element_text(hjust=0.5))
print(p)

p = ggplot(dat, aes(x=sample_name, y=gc_content, fill=fraction)) +
  geom_bar(stat="identity") +
  labs(x="Sample", y = "GC content (%)", title="Mean GC content") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        plot.title=element_text(hjust=0.5))
print(p)

p = ggplot(dat, aes(x=sample_name, y=mean_read_len, fill=fraction)) +
  geom_bar(stat="identity") + 
  geom_point(aes(x=sample_name, y=max_read_len), shape=6) +
  geom_point(aes(x=sample_name, y=min_read_len), shape=2) +
  labs(x="Sample", y = "Mean read length (nt)", title="Read length") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        plot.title=element_text(hjust=0.5))
print(p)