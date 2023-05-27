library(DESeq2)
library(edgeR)
library(limma)
library(ggplot2)
library(tidyverse)
library(stringr)

setwd("~/Dev/mir-4728/3_counts/")

sel_Data = "ribosome"

if(sel_Data=="polysome"){
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
  sample_table$type = relevel(sample_table$type, sel_Data)
  row.names(sample_table) = sample_names 
  
}else if(sel_Data=="monosome"){
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
    dplyr::select(1:6, 13:18)
  
  dat = dat %>%  
    filter(!grepl("__", rownames(dat))) %>% 
    as.matrix()
  
  dat[dat==0] = NA
  dat = dat[complete.cases(dat),]
  
  # Create sample annotations
  sample_names = c(sample_names[1:6], sample_names[13:18])
  
  sample_conditions = rep(c("control", "control", "control", "treatment", "treatment", "treatment"),2)
  sample_type = c(rep("total", 6), rep("monosome", 6))
  
  sample_table = data.frame(condition = as.factor(sample_conditions), type = as.factor(sample_type))
  sample_table$type = relevel(sample_table$type, sel_Data)
  row.names(sample_table) = sample_names 
  
}else if(sel_Data=="ribosome"){
# ======================== RIBOSOME ========================
# Import raw data
  raw_dat = read.csv("ribosome_gene_counts_yes.txt", sep="\t")
  sample_names = colnames(raw_dat)
  sample_names = sample_names[!grepl("geneID", sample_names)]
  
  colnames(raw_dat) = c("GeneID",sample_names)
  
  raw_dat$GeneID = raw_dat$GeneID %>% 
    str_remove("\\.[0-9]*$")
  
  dat = raw_dat %>% 
    column_to_rownames("GeneID") %>% 
    dplyr::select(c(1, 3:12))
  
  dat = dat %>%  
    filter(!grepl("__", rownames(dat))) %>% 
    as.matrix()
  
  dat[dat==0] = NA
  dat = dat[complete.cases(dat),]
  
  # Create sample annotations
  sample_names = sample_names[c(1,3:12)]
  
  sample_conditions = c("treatment", "treatment","control", "control", "control",
                        "treatment", "treatment", "treatment","control", "control", "control")
  sample_type = c(rep("ribosome", 5), rep("total", 6))
  
  sample_table = data.frame(condition = as.factor(sample_conditions), type = as.factor(sample_type))
  row.names(sample_table) = sample_names
  
# ========================================================== 
}

dat_TE = dat
dat_totalRNA = dat[,rownames(sample_table)[which(sample_table$type=="total")]]
dat_polysome = dat[,rownames(sample_table)[which(sample_table$type==sel_Data)]]

# Create DESeq objects
dds_TE = DESeqDataSetFromMatrix(countData = dat_TE,
                             colData = sample_table,
                             design = ~ condition + type + condition:type)

dds_totalRNA = DESeqDataSetFromMatrix(countData = dat_totalRNA,
                                colData = sample_table[which(sample_table$type=="total"),],
                                design = ~ condition)

dds_polysome = DESeqDataSetFromMatrix(countData = dat_polysome,
                                colData = sample_table[which(sample_table$type==sel_Data),],
                                design = ~ condition)

# Preliminary QC
# Estimate size factors and normalize
dds_TE = estimateSizeFactors(dds_TE)
dds_totalRNA = estimateSizeFactors(dds_totalRNA)
dds_polysome = estimateSizeFactors(dds_polysome)

dds_TE = estimateDispersions(dds_TE)
dds_totalRNA = estimateDispersions(dds_totalRNA)
dds_polysome = estimateDispersions(dds_polysome)

normalized_counts_TE = counts(dds_TE, normalized=TRUE) %>% 
  as.data.frame()

normalized_counts_totalRNA = counts(dds_totalRNA, normalized=TRUE) %>% 
  as.data.frame()

normalized_counts_polysome = counts(dds_polysome, normalized=TRUE) %>% 
  as.data.frame()

# Generate PCA for quality control
rlog_counts_TE = rlogTransformation(dds_TE, blind=TRUE)
rlog_counts_totalRNA = rlogTransformation(dds_totalRNA, blind=TRUE)
rlog_counts_polysome = rlogTransformation(dds_polysome, blind=TRUE)

pca_obj = plotPCA(rlog_counts_TE, intgroup=c("condition", "type"), returnData=TRUE)

percentVar = round(100 * attr(pca_obj, "percentVar"))

colnames(pca_obj)[4] = "Condition"
colnames(pca_obj)[5] = "Type"

pca_plot = ggplot(pca_obj, aes(x=PC1, y=PC2, color=Condition, shape=Type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

print(pca_plot)

# print(plotPCA(rlog_counts_totalRNA, intgroup=c("condition")))
# print(plotPCA(rlog_counts_polysome, intgroup=c("condition")))

# Perform DGE analysis
dds_TE = DESeq(dds_TE)
dds_totalRNA = DESeq(dds_totalRNA)
dds_polysome = DESeq(dds_polysome)

# Create a model matrix for contrasts TE
mod_mat = model.matrix(design(dds_TE), colData(dds_TE))

total_control = colMeans(mod_mat[dds_TE$type == "total" & dds_TE$condition == "control", ])
total_treatment = colMeans(mod_mat[dds_TE$type == "total" & dds_TE$condition == "treatment", ])
polysome_control = colMeans(mod_mat[dds_TE$type == sel_Data & dds_TE$condition == "control", ])
polysome_treatment = colMeans(mod_mat[dds_TE$type == sel_Data & dds_TE$condition == "treatment", ])

# Extract TE results
deseq_TE = results(dds_TE, contrast=(total_control - total_treatment)-(polysome_control - polysome_treatment))
deseq_TE = deseq_TE[order(deseq_TE$padj),] %>% 
  as.data.frame() %>% 
  rownames_to_column("geneID")

deseq_TE$geneID = deseq_TE$geneID %>% 
  str_remove("\\.[0-9]*$")

# Create a model matrix for contrasts total RNA
mod_mat = model.matrix(design(dds_totalRNA), colData(dds_totalRNA))

total_control = colMeans(mod_mat[dds_totalRNA$type == "total" & dds_totalRNA$condition == "control", ])
total_treatment = colMeans(mod_mat[dds_totalRNA$type == "total" & dds_totalRNA$condition == "treatment", ])

# Extract total RNA results
deseq_totalRNA = results(dds_totalRNA, contrast=total_treatment-total_control)
#deseq_totalRNA = lfcShrink(dds_totalRNA, contrast=total_treatment-total_control, res=deseq_totalRNA)

deseq_totalRNA = deseq_totalRNA[order(deseq_totalRNA$padj),] %>% 
  as.data.frame() %>% 
  rownames_to_column("geneID")

deseq_totalRNA$geneID = deseq_totalRNA$geneID %>% 
  str_remove("\\.[0-9]*$")

# Create a model matrix for contrasts polysome/ribosome
mod_mat = model.matrix(design(dds_polysome), colData(dds_polysome))

polysome_control = colMeans(mod_mat[dds_polysome$type == sel_Data & dds_polysome$condition == "control", ])
polysome_treatment = colMeans(mod_mat[dds_polysome$type == sel_Data & dds_polysome$condition == "treatment", ])

# Extract polysome/ribosome results
deseq_polysome = results(dds_polysome, contrast=polysome_treatment-polysome_control)
#deseq_polysome = lfcShrink(dds_polysome, contrast=polysome_treatment-polysome_control, res=deseq_polysome)

deseq_polysome = deseq_polysome[order(deseq_polysome$padj),] %>% 
  as.data.frame() %>% 
  rownames_to_column("geneID")

deseq_polysome$geneID = deseq_polysome$geneID %>% 
  str_remove("\\.[0-9]*$")

# Reformat results for plotting
deseq_TE = deseq_TE %>% 
  dplyr::select(geneID, log2FoldChange, padj)
colnames(deseq_TE) = c("geneID", "TE_lfc", "TE_padj")

deseq_totalRNA = deseq_totalRNA %>% 
  dplyr::select(geneID, log2FoldChange, padj)
colnames(deseq_totalRNA) = c("geneID", "totalRNA_lfc", "totalRNA_padj")

deseq_polysome = deseq_polysome %>% 
  dplyr::select(geneID, log2FoldChange, padj)
colnames(deseq_polysome) = c("geneID",  "polysome_lfc", "polysome_padj")

# Merge results into one table
plot_table = deseq_TE %>% 
  merge(deseq_totalRNA, by="geneID") %>% 
  merge(deseq_polysome, by="geneID") %>% 
  mutate(mode="background")

# Remove NAs 
# plot_table = plot_table %>%
#   drop_na()

# Assign regulatory modes
plot_table$mode[which(plot_table$TE_padj > 0.05 &
                        plot_table$totalRNA_padj < 0.05 &
                        plot_table$polysome_padj < 0.05)] = "transcription"

plot_table$mode[which(plot_table$TE_padj < 0.05 &
                        plot_table$totalRNA_padj > 0.05 &
                        plot_table$polysome_padj < 0.05)] = "translation"

plot_table$mode[which(plot_table$TE_padj < 0.05 &
                        plot_table$totalRNA_padj < 0.05 &
                        plot_table$polysome_padj < 0.05 &
                        sign(plot_table$TE_lfc) != sign(plot_table$totalRNA_lfc))] = "buffered"
plot_table$mode[which(plot_table$TE_padj < 0.05 &
                        plot_table$totalRNA_padj < 0.05 &
                        plot_table$polysome_padj > 0.05)] = "buffered"

color_palette = c(background = "gray", 
                  buffered = "orange", 
                  transcription = "steelblue1", 
                  translation = "tomato2")

# Create FC plot
p = ggplot() +
  geom_point(data=plot_table %>% arrange(mode), aes(x=totalRNA_lfc, y=polysome_lfc, color=as.factor(mode))) +
  scale_x_continuous(limits=c(-6,6)) +
  scale_y_continuous(limits=c(-6,6)) +
  scale_color_manual(values=color_palette) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_vline(xintercept=log2(1.2), linetype="solid") +
  geom_vline(xintercept=-log2(1.2), linetype="solid") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_hline(yintercept=log2(1.2), linetype="solid") +
  geom_hline(yintercept=-log2(1.2), linetype="solid") +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  labs(x = "Total RNA log2FC", y = "Monosome-associated RNA log2FC", color = "Regulatory mode")

print(p)

write.table(plot_table, "~/Dev/mir-4728/4_DE/polysome/DESeq2/monosome_DESeq.txt", sep="\t", row.names=FALSE, quote=FALSE)
