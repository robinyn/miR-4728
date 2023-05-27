# ============================================================================================================
# Title: anota2seq_workflow.R
# Author: Euisuk Robin Han
# Description: A script for the DGE analysis of polysome profiling data using anota2seq
# Date: 12/Apr/23
# ============================================================================================================
library(tidyverse)
library(edgeR)
library(limma)
library(ggplot2)
library(anota2seq)
library(DESeq2)

setwd("~/Dev/mir-4728")

# ======================== POLYSOME ========================
# Import gene counts
raw_dat = read.csv("3_counts/polysome_gene_counts.txt", sep="\t")
sample_names = c("SK0001","SK0002","SK0003","SK0004","SK0005","SK0006","SK0007","SK0008","SK0009","SK0010",
                 "SK0011","SK0012","SK0013","SK0014","SK0015","SK0016","SK0017","SK0018")

colnames(raw_dat) = c("GeneID",sample_names)

dat = raw_dat %>% 
  column_to_rownames("GeneID")

dat = dat %>%
  filter(!grepl("__", rownames(dat)))

# Generate sample annotation
sample_anot = data.frame(sample = sample_names, 
                         fraction = c(rep("total", 6),
                                      rep("polysome", 6),
                                      rep("monosome", 6)),
                         condition = c(rep(c("control","control","control", "treatment", "treatment", "treatment"), 3)),
                         replicate = c(rep(c(1,2,3), 6))) %>% 
  column_to_rownames("sample")

# Quality control

# Create a PCA object 
PCA_obj = as.matrix(dat)

# Remove zeros
PCA_obj = PCA_obj[!apply(PCA_obj, 1, function(x)any(x==0)),]

# Normalize (TMM) and transform
PCA_obj_Norm = voom(calcNormFactors(DGEList(PCA_obj)))$E

# Remove transcripts with low variability across all samples 
# Standard approach is to only use transcripts in the highest quantile of SD. 

# Calculate SD 
transcripts_SD = apply(PCA_obj_Norm, 1, sd)

# Calculate quantiles
sdQuantiles = quantile(transcripts_SD)

# Select transcripts with SD in the highest quantile
PCA_obj_filt = PCA_obj_Norm[transcripts_SD>sdQuantiles[4],]

# Transpose matrix
PCA_obj_Trans = t(PCA_obj_filt)

# PCA
PCA_out = prcomp(PCA_obj_Trans)

# Visualize PCA
prop_Var = ggplot(data = NULL) +
  geom_bar(stat = "identity", aes(x=names(summary(PCA_out)$importance[2,]), 
                                  y=summary(PCA_out)$importance[2,])) + 
  geom_hline(yintercept=0.05, col="red") +
  scale_x_discrete(limits=x_names) +
  labs(x="Component", y="Proportion of variance")
  
print(prop_Var)

cumm_Var = ggplot(data = NULL) +
  geom_bar(stat = "identity", aes(x=names(summary(PCA_out)$importance[3,]),
                                  y=summary(PCA_out)$importance[3,])) +
  geom_hline(yintercept=0.9, col="red") +
  scale_x_discrete(limits=x_names) +
  labs(x="Component", y="Cumulative proportion of variance")

print(cumm_Var)

# PCA Plot object
PCA_plot = merge(PCA_out$x, sample_anot, by="row.names")

PC1vPC2 = ggplot(data=PCA_plot, aes(x=PC1, y=PC2, shape=fraction, col=condition)) +
  geom_point() +
  geom_text(aes (label=replicate,vjust=1,hjust=1))

print(PC1vPC2)

PC1vPC3 = ggplot(data=PCA_plot, aes(x=PC1, y=PC3, shape=fraction, col=condition)) +
  geom_point() +
  geom_text(aes (label=replicate,vjust=1,hjust=1))

print(PC1vPC3)

# Prep data for anota2seq
polysome_dat = dat %>% 
  dplyr::select(row.names(sample_anot[sample_anot$fraction==frac_sel,])) %>% 
  as.matrix()

total_dat = dat %>% 
  dplyr::select(row.names(sample_anot[sample_anot$fraction=="total",])) %>% 
  as.matrix()

pheno_vec = sample_anot[sample_anot$fraction==frac_sel, "condition"]

# Generate anota2seq dataset
ads = anota2seqDataSetFromMatrix(
  dataP = polysome_dat,
  dataT = total_dat,
  phenoVec = pheno_vec,
  dataType = "RNAseq",
  normalize = TRUE,
  transformation = "TMM-log2",
  filterZeroGenes = TRUE
)

# Perform QC
ads = anota2seqPerformQC(ads)

# Outlier test
ads = anota2seqResidOutlierTest(ads)

# Create contrast matrix
myContrast = matrix(nrow=length(levels (as.factor(pheno_vec))),
                     ncol=length(levels(as.factor(pheno_vec)))-1)

rownames(myContrast) = levels(as.factor (pheno_vec))

myContrast[,1] = c(-1,1)

# Run anota2seq analysis
ads = anota2seqAnalyze(ads, contrasts=myContrast)

# Select significant genes
ads = anota2seqSelSigGenes(ads,
                            maxPAdj = 0.05,
                            selDeltaPT = log2(1.2),
                            selDeltaTP = log2(1.2),
                            selDeltaP = log2(1.2),
                            selDeltaT = log2(1.2),
                            minSlopeTranslation = -1,
                            maxSlopeTranslation = 2,
                            minSlopeBuffering = -1,
                            maxSlopeBuffering = 0)

# Set regulatory modes based on analysis
ads = anota2seqRegModes(ads)

# Plot FC plot
anota2seqPlotFC(ads, selContrast = 1, plotToFile = FALSE)

# Get output
#results_table = anota2seqGetOutput(ads, output = "singleDf", selContrast = 1)

# Save ads
#save(ads, file="4_DE/polysome/polysome_anota2seq.robj")