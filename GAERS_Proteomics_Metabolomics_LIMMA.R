# This is for the single-omics analysis of both proteomics and metabolomics data using limma
# The input tables are required to have samples as columns and features as rows. 

# Setting up Working Directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Loading libraries
library(edgeR)
library(limma)
library(tidyverse)
library(ggplot2)
library(ggrepel)

#----Proteomics Single-Omics----
# Read in the file
dt_tbl <- read.csv("Prot_SCx_3wk_noOut_timpt.csv", stringsAsFactors = F, 
                   sep = ",", row.names = 1)
# Transposing
dt_tbl <- as.data.frame(t(dt_tbl))
# Making Design Table and Contrast
# Design table
samp = as.factor(c(rep("GAERS", 11), rep("NEC", 11)))
# Generating design matrix
des <- model.matrix(~0 + samp) # 0 means no intercept for linear model
colnames(des) = gsub("samp","",colnames(des))
# Making contrasts
cont.rule <- c("GAERS - NEC")
contrast = makeContrasts(contrasts = cont.rule, levels =  des)
# Fitting the linear model and performing DE analysis
# Fit the expression matrix to a linear model
fit <- lmFit(dt_tbl, des)
# Compute contrast
fit_contrast <- contrasts.fit(fit, contrast)
# Bayes statistics of differential expression
fit_contrast <- eBayes(fit_contrast, trend = T)
# We can also plot the mean-variance relationship to see if there is a trend in the expression of the features
plotSA(fit_contrast) 
title("Mean-Variance Relationship")
# Generate a table of all proteins
tab <- topTable(fit_contrast, n=Inf, adjust.method = "BH")
# Summary of number of differentially expressed genes
result <- decideTests(fit_contrast, adjust.method = "BH")
summary(result)
# Plotting Volcano Plot
# Add an additional column detailing if the DE is significantly up-regulated, 
# significantly down-regulated, or not significantly differentially expressed.
tab_mod <- tab %>%
  mutate("Sig DE" = case_when(logFC >= log2(1.5) & adj.P.Val <= 0.05 ~ "Up",
                              logFC <= -log2(1.5) & adj.P.Val <= 0.05 ~ "Down",
                              TRUE ~ "Not Sig"))
# Saving as .csv file
write.csv(tab_mod, "Prot_Tha_16wks_DEAll_FC1.5.csv", row.names = T)
# Saving the plot
pdf("Prot_Tha_16wk_DE_Plots.pdf")
# Then, we plot the volcano plot
ggplot(tab_mod, aes(logFC, -log(P.Value,10))) +
  geom_point(aes(color = `Sig DE`), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"(P-Value)")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(data = tab_mod,
                  aes(logFC, -log(P.Value,10), 
                      label = ifelse(tab_mod$`Sig DE` !="Not Sig", rownames(tab_mod), "")),
                  size = 3,
                  max.overlaps = 50) +
  ggtitle("Proteomics DE Scx 3wk GAERS v NEC, FC = 1.5")
dev.off()
# Saving DE proteins
# Subsetting
sub <- subset(tab_mod, tab_mod$`Sig DE` != "Not Sig")
# Adding the gene names
names <- read.csv("//ad.monash.edu/home/User028/dcho0009/Desktop/r/Genes and Accession.csv")
sub_names <- subset(names, names$ID %in% rownames(sub))
sub_names <- sub_names[match(rownames(sub), sub_names$ID),]
sub <- add_column(sub, `Gene Name` = sub_names$Gene.names, .before = 1)
# Saving the table
write.csv(sub,"Prot_SCx_16wks_DEOnly.csv")

#---- Metabolomics Single-Omics ----
# This is performed in a similar fashion to proteomics
# Loading Files
dt_tbl <- read.csv("//ad.monash.edu/home/User028/dcho0009/Desktop/r/Met_Tha_16wk_noOut_timpt.csv", 
                   stringsAsFactors = F, sep = ",", row.names = 1)
dt_tbl <- as.data.frame(t(dt_tbl))
# Making design table
samp = as.factor(c(rep("GAERS", 10), rep("NEC", 11))) #Change number of repetitions as required.
# Generating design matrix
des <- model.matrix(~0 + samp)
colnames(des) = gsub("samp","",colnames(des))
# Making contrasts
cont.rule <- c("GAERS - NEC")
contrast = makeContrasts(contrasts = cont.rule, levels =  des)
# Analysis
fit <- lmFit(dt_tbl, des)
fit_contrast <- contrasts.fit(fit, contrast)
fit_contrast <- eBayes(fit_contrast, trend = T)
# Checking for trend in the expression of the features
plotSA(fit_contrast) 
title("Mean-Variance Relationship")
# Generate a table of all proteins
tab <- topTable(fit_contrast, n=Inf, adjust.method = "BH")
# Summary of number of differentially expressed genes
result <- decideTests(fit_contrast, adjust.method = "BH")
summary(result)
# Plotting a volcano plot using ggplot2
# Editing data table
tab_mod <- tab %>%
  mutate("Sig DE" = case_when(logFC >= log2(1.5) & adj.P.Val <= 0.05 ~ "Up",
                              logFC <= -log2(1.5) & adj.P.Val <= 0.05 ~ "Down",
                              TRUE ~ "Not Sig"))
# Saving the plot
pdf("Met_Tha_16wk_DE_Plots.pdf")
# Plotting volcano plot
ggplot(tab_mod, aes(logFC, -log(P.Value,10))) +
  geom_point(aes(color = `Sig DE`), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"(P-Value)")) +
  #lims(x = c(-5, 5), y = c(0, 10)) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(data = tab_mod,
                  aes(logFC, -log(P.Value,10), 
                      label = ifelse(tab_mod$`Sig DE` !="Not Sig", rownames(tab_mod), "")),
                  size = 3) +
  ggtitle("Metabolomics DE Tha 16wk GAERS v NEC")
dev.off()

# Saving as .csv file
write.csv(tab_mod, "Met_Tha_16wks_DEAll_FC1.5.csv")
# Subset
sub <- subset(tab_mod, tab_mod$`Sig DE` != "Not Sig")
write.csv(sub,"Met_SCx_16wks_DEOnly.csv")