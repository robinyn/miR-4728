library(RColorBrewer)
library(stringr)

transcription_enrichment = enrichment_KEGG_transcription@result
translation_enrichment = enrichment_KEGG_translation@result

# Order by padj
transcription_enrichment = transcription_enrichment[order(transcription_enrichment$p.adjust),]
translation_enrichment = translation_enrichment[order(translation_enrichment$p.adjust),]

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
  scale_fill_gradient(low="red", high="blue") +
  scale_y_discrete(limits=rev, labels=function(x) str_wrap(x, 60)) +
  scale_x_continuous(limits=c(-0.8, 0.8)) +
  labs(x="Enrichment score", y="Pathways") +
  facet_wrap(vars(type))

print(p)

