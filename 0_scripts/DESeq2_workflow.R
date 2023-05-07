library(DESeq2)
library(edgeR)
library(limma)
library(ggplot2)
library(tidyverse)
library(stringr)

setwd("~/Dev/mir-4728/3_counts/")

# ======================== POLYSOME ========================
# Import raw data
raw_dat = read.csv("polysome_gene_counts.txt", sep="\t")
sample_names = c("SK0001","SK0002","SK0003","SK0004","SK0005","SK0006","SK0007","SK0008","SK0009","SK0010",
                 "SK0011","SK0012","SK0013","SK0014","SK0015","SK0016","SK0017","SK0018")

colnames(raw_dat) = c("GeneID",sample_names)

raw_dat$GeneID = raw_dat$GeneID %>% 
  str_remove("\\.[0-9]*$")

dat = raw_dat %>% 
  column_to_rownames("GeneID") %>% 
  dplyr::select(1:12)

dat = dat %>%  
  filter(!grepl("__", rownames(dat))) %>% 
  as.matrix()
  
dat[dat==0] = NA
dat = dat[complete.cases(dat),]

# Create sample annotations
sample_names = sample_names[1:12]

sample_conditions = rep(c("control", "control", "control", "treatment", "treatment", "treatment"),2)
sample_type = c(rep("total", 6), rep("polysome", 6))

sample_table = data.frame(condition = as.factor(sample_conditions), type = as.factor(sample_type))
row.names(sample_table) = sample_names

# Create DESeq object
dds = DESeqDataSetFromMatrix(countData = dat,
                             colData = sample_table,
                             design = ~ condition + type + condition:type)

# Estimate size factors and normalize
dds = estimateSizeFactors(dds)

normalized_counts = counts(dds, normalized=TRUE) %>% 
  as.data.frame()

# Generate PCA for quality control
rlog_counts = rlogTransformation(dds, blind=TRUE)

print(plotPCA(rlog_counts, intgroup=c("condition", "type")))

# Perform DGE analysis
dds = DESeq(dds)

# Create a model matrix for contrasts
mod_mat = model.matrix(design(dds), colData(dds))

total_control = colMeans(mod_mat[dds$type == "total" & dds$condition == "control", ])
total_treatment = colMeans(mod_mat[dds$type == "total" & dds$condition == "treatment", ])
polysome_control = colMeans(mod_mat[dds$type == "polysome" & dds$condition == "control", ])
polysome_treatment = colMeans(mod_mat[dds$type == "polysome" & dds$condition == "treatment", ])

# total mRNA - total mRNA : transcriptional regulation
# total_control - total_treatment

# interaction : translational efficiency
# (polysome_control - polysome_treatment) - (total_control - total_treatment)

# Translational regulation 
deseq_translation = results(dds, contrast=(total_control - total_treatment)-(polysome_control - polysome_treatment))
deseq_translation = deseq_translation[order(deseq_translation$padj),] %>% 
  as.data.frame() %>% 
  rownames_to_column("geneID")

deseq_translation$geneID = deseq_translation$geneID %>% 
  str_remove("\\.[0-9]*$")

# Transcriptional regulation
deseq_transcription = results(dds, contrast=total_treatment-total_control)
deseq_transcription = deseq_transcription[order(deseq_transcription$padj),] %>% 
  as.data.frame() %>% 
  rownames_to_column("geneID")

deseq_transcription$geneID = deseq_transcription$geneID %>% 
  str_remove("\\.[0-9]*$")

# Create log2FC table for all genes for both transcriptional/translational regulation
deseq_transcription_fc = deseq_transcription %>% 
  dplyr::select(geneID, log2FoldChange)

deseq_translation_fc = deseq_translation %>% 
  dplyr::select(geneID, log2FoldChange)

colnames(deseq_transcription_fc) = c("geneID", "transcription_fc")
colnames(deseq_translation_fc) = c("geneID", "translation_fc")

# Create a plot table for log2FC vs log2FC plot
plot_table = deseq_transcription_fc %>% 
  merge(deseq_translation_fc, by="geneID")

# Filter genes to assign regulatory modes
deseq_transcription = deseq_transcription %>% 
  drop_na() %>% 
  filter(padj < 0.05)

deseq_translation = deseq_translation %>% 
  drop_na() %>% 
  filter(padj < 0.05)

plot_table = plot_table %>% 
  mutate(mode=NA)

plot_table$mode[!(plot_table$geneID %in% deseq_transcription$geneID) & 
                  (plot_table$geneID %in% deseq_translation$geneID)] = "translation"
plot_table$mode[(plot_table$geneID %in% deseq_transcription$geneID) & 
                  (plot_table$geneID %in% deseq_translation$geneID) & 
                  (sign(plot_table$transcription_fc)==sign(plot_table$translation_fc))] = "transcription"
plot_table$mode[is.na(plot_table$mode)] = "background"

plot_table$mode = factor(plot_table$mode, levels=c("background","translation", "transcription"))

# Create FC plot
p = ggplot() +
  geom_point(data=plot_table %>% arrange(mode), aes(x=transcription_fc, y=translation_fc, color=as.factor(mode))) +
  scale_x_continuous(limits=c(-6,6)) +
  scale_y_continuous(limits=c(-6,6)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_vline(xintercept=log2(1.2), linetype="solid") +
  geom_vline(xintercept=-log2(1.2), linetype="solid") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_hline(yintercept=log2(1.2), linetype="solid") +
  geom_hline(yintercept=-log2(1.2), linetype="solid") +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  geom_abline(intercept=log2(1.2), slope=1, linetype="solid") +
  geom_abline(intercept=-log2(1.2), slope=1, linetype="solid")

print(p)