library(RColorBrewer)
library(stringr)
library(tidyverse)

setwd("~/Dev/mir-4728/5_GSEA/results/DESeq2")

transcription = read_tsv("~/Dev/mir-4728/5_GSEA/results/DESeq2/ribosome/tables/enrichment_GO_BP_transcription.tsv")
translation = read_tsv("~/Dev/mir-4728/5_GSEA/results/DESeq2/ribosome/tables/enrichment_GO_BP_translation.tsv")

# transcription = enrichment_GO_BP_transcription
# translation = enrichment_GO_BP_translation

transcription_enrichment = transcription %>% 
  dplyr::select(Description, enrichmentScore, p.adjust) %>% 
  filter(p.adjust <0.05)
translation_enrichment = translation %>% 
  dplyr::select(Description, enrichmentScore, p.adjust) %>% 
  filter(p.adjust <0.05)

# Order by padj
transcription_enrichment = transcription_enrichment[order(transcription_enrichment$p.adjust),]
translation_enrichment = translation_enrichment[order(translation_enrichment$p.adjust),]

# Order by enrichment score
# transcription_enrichment = transcription_enrichment[order(transcription_enrichment$enrichmentScore),]
# translation_enrichment = translation_enrichment[order(translation_enrichment$enrichmentScore),]

# Remove cutoff
cutoff=25

plot_data_transcription = transcription_enrichment[1:cutoff,] %>% 
  drop_na()

plot_data_translation = translation_enrichment[1:cutoff,] %>% 
  drop_na()

plot_data_transcription$Description = factor(plot_data_transcription$Description, 
                                             levels=unique(plot_data_transcription$Description))
plot_data_translation$Description = factor(plot_data_translation$Description, 
                                             levels=unique(plot_data_translation$Description))

plot_data_transcription = plot_data_transcription %>% 
  mutate(type="transcription")

plot_data_translation = plot_data_translation %>% 
  mutate(type="translation")

plot_data = plot_data_transcription %>% 
  rbind(plot_data_translation)

p = ggplot() +
  geom_bar(data=plot_data, aes(x=enrichmentScore, y=Description, fill=p.adjust), stat="identity") +
  scale_fill_gradient(low="royalblue2", high="grey") +
  scale_y_discrete(limits=rev, labels=function(x) str_wrap(x, 60)) +
  scale_x_continuous(limits=c(-1, 1)) +
  labs(x="Enrichment score", y="Pathways") +
  facet_wrap(vars(type))

print(p)

