# This is for the edgeR analysis of the 16-week transcriptomic dataset

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Loading libraries
library(edgeR)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

# Reading in files
count <- read.csv("GAERS_Tha_Count.csv")
# Reorder samples
count <- count %>%
  relocate(contains("GAERS"), .after = 2)

# Making DGEList for analysis
count_list <- count[,-c(1:2)]  # Getting count dataframe
lab <- c(rep("GAERS",12), rep("NEC",12)) # Grouping information
# Sample information
samps <- colnames(count_list) 
# Setting DGEList
list <- DGEList(counts = count_list, group = lab, genes = count[,c(1:2)],
                samples = samps)
list$samples # to see library size of the samples

# Filtering 
keep <- filterByExpr(list) # Getting list of genes to remove
table(keep)
list_filt <- list[keep, , keep.lib.sizes = F] # Filtering

# Density Plots
# General information
nsamples <- ncol(list)
col.density <- colorRampPalette(brewer.pal(12, name = "Paired"))(nsamples)
# Plotting unfiltered data
lcpm_dens <- cpm(list, log = T)
plot(density(lcpm_dens[, 1], bw = 0.4), col = col.density[1], lwd = 2, ylim = c(0, 1),
     las = 2, main = "", xlab = "")
title(main = "Raw data", xlab = "Log-cpm")
for (i in 2:nsamples) {
  den <- density(lcpm_dens[, i], bw = 0.4)
  lines(den$x, den$y, col = col.density[i], lwd = 2)
}
legend("topright", samps, text.col = col.density, bty = "n", ncol = 5)
Raw_dens <- recordPlot()
# Plotting filtered data
lcpm_filt <- cpm(list_filt, log = T)
plot(density(lcpm_filt[, 1], bw = 0.4), col = col.density[1], lwd = 2, ylim = c(0, 1),
     las = 2, main = "", xlab = "")
title(main = "Filtered data", xlab = "Log-cpm")
for (i in 2:nsamples) {
  den <- density(lcpm_filt[, i], bw = 0.4)
  lines(den$x, den$y, col = col.density[i], lwd = 2)
}
legend("topright", samps, text.col = col.density, bty = "n", ncol = 5)
Filt_dens <- recordPlot()

# Normalisation
list_filt <- normLibSizes(list_filt, method = "TMM")
list_filt$samples
# Plotting boxplots
# Before normalisation
lcpm <- cpm(list, log = T)
par(mar=c(8,4,4,1))
boxplot(lcpm, las=2)
title(main = "Unormalised data", ylab = "Log-cpm")
Unnorm_boxplt <- recordPlot()
# After normalisation
lcpm2 <- cpm(list_filt, log = T)
par(mar=c(8,4,4,1))
boxplot(lcpm2, las=2)
title(main = "Normalised data", ylab = "Log-cpm")
norm_boxplt <- recordPlot()

# Check the expression level of the genes for a sample relative to other samples
plotMD(list_filt, column = 2)
abline(h=0, col="red", lty=2, lwd=2)

# Outlier removal
# Sample clustering
par(mar=c(5,4,4,2))
plotMDS(list_filt)
# Removal
outliers <- c("S86_GAERS_1","S97_GAERS_7") #c("S11_GAERS_7","S60_GAERS_1")
lt_noOut <- list_filt[ ,which(!list_filt$samples$samples %in% outliers)]
# Check
lt_noOut$samples
head(lt_noOut[["counts"]])
# New MDS
plotMDS(lt_noOut)

#----Exporting for MOFA----
# Getting Normalised dataframe
lt_noOut_cpm <- cpm(lt_noOut, log = T)

# Checking data distribution
nsamp_new <- length(colnames(lt_noOut_cpm))
plot(density(lt_noOut_cpm[, 1], bw = 0.4), col = col.density[1], lwd = 2, ylim = c(0, 1),
     las = 2, main = "", xlab = "")
title(main = "Filtered + Normalised data", xlab = "Log-cpm")
for (i in 2:nsamp_new) {
  den <- density(lt_noOut_cpm[, i], bw = 0.4)
  lines(den$x, den$y, col = col.density[i], lwd = 2)
}
legend("topright", colnames(lt_noOut_cpm), text.col = col.density, bty = "n", ncol = 5)
fin_dens <- recordPlot()

# Changing rownames so that it is gene names
rownames(lt_noOut_cpm) <- lt_noOut$genes$gene_name
# Changing column names
mofa_lab <- read.csv("GAERS Transcriptomics/GAERS_NewLabels.csv", stringsAsFactors = F, header = T)
common_col <- intersect(mofa_lab$label, colnames(lt_noOut_cpm)) # Finding common names
mofa_lab <- mofa_lab %>%
  subset(label %in% common_col) %>%
  arrange(factor(label, levels = colnames(lt_noOut_cpm)))
colnames(lt_noOut_cpm)  <- mofa_lab$new.label # Replacing names
colnames(lt_noOut_cpm) <- gsub("scx_16 weeks","16s",colnames(lt_noOut_cpm))

# Re-ordering by variances to get the top transcripts (do this before turning lt_noOut_cpm to dataframe)
library(matrixStats)
o <- order(rowVars(lt_noOut_cpm), decreasing = TRUE)
cpm_ord <- lt_noOut_cpm[o,]
cpm_ord <- as.data.frame(t(cpm_ord))
write.csv(cpm_ord, "PTE_Transcriptomics_CPMord.csv")

# Re-structuring
lt_noOut_cpm <- lt_noOut_cpm %>%
  as.data.frame(.) %>%
  select(order(colnames(.)))
lt_noOut_cpm <- as.data.frame(t(lt_noOut_cpm))
# Writing
write.csv(lt_noOut_cpm, "PTE_Transcriptomics_MOFA.csv",  row.names = T)

#----Saving Plots and Data----
pdf("JoAnne_Transcriptomics_Pre-processing Plots.pdf")
Raw_dens
Filt_dens
Unnorm_boxplt
norm_boxplt
par(mar=c(5,4,4,2))
plotMDS(list_filt)
plotMDS(lt_noOut)
fin_dens
dev.off()

# Saving as RData
save.image(file = "JoAnne_Transcriptomics_Preprocessing.RData")
# Loading
#load("GAERS Transcriptomics/GAERS_Scx_Transcriptomics_Preprocessing.RData")

#----Differential Expression----
# Since there are only two groups, a simple t-test can be performed

# Getting design matrix and contrast
lab_noOut <- as.factor(c(rep("GAERS",10), rep("NEC",12))) # Making strains factors

des <- model.matrix(~0 + lab_noOut) # model matrix
rownames(des) <- colnames(lt_noOut)
colnames(des) <- gsub("lab_noOut","", colnames(des))

cont.rule <- c("GAERS - NEC")
contrast = makeContrasts(contrasts = cont.rule, levels =  des) # making contrast

# Dispersion Estimate
lt_noOut <- estimateDisp(lt_noOut, des, robust = T)
lt_noOut$common.dispersion
plotBCV(lt_noOut)
# Fitting Model
fit <- glmQLFit(lt_noOut, des, robust = T)
plotQLDisp(fit)

# Differential Expression
qlf <- glmQLFTest(fit,contrast=contrast)
tab <- topTags(qlf, adjust.method = "BH", sort.by = "PValue", n = Inf)
is.de <- decideTestsDGE(qlf)
summary(is.de)
# Specific cut-off. Change 'lfc' as required
FC1.5 <- decideTestsDGE(qlf, adjust.method = "BH", p.value = 0.05, lfc = 0.585)
sum_FC1.5 <- as.data.frame(summary(FC1.5))
colnames(sum_FC1.5) <- c("Direction", "Pair", "Number")
write.csv(sum_FC1.5, file = "GAERS_Tha_Trans_16wks_DESummary_FC1.5.csv", row.names = F)

# Saving differentially expressed transcripts
tab_mod <- tab$table %>%
  mutate("Sig DE" = case_when(logFC >= log2(1.5) & FDR <= 0.05 ~ "Up",
                              logFC <= -log2(1.5) & FDR <= 0.05 ~ "Down",
                              TRUE ~ "Not Sig"))
# Saving as .csv file
write.csv(tab_mod, "GAERS_Tha_Trans_16wks_DEAll_FC1.5.csv", row.names = T)

# Saving as RData
save.image(file = "GAERS_Tha_Transcriptomics_DEAnalysis.RData")

#----Volcano Plot----
# Can plot all at the same time. Here we use two different FC cut-off.
pdf("GAERS_Tha_Transcriptomics_DE Plots.pdf")
plotBCV(lt_noOut)
plotQLDisp(fit)
plotMD(qlf, status = is.de, values=c(1,-1), col=c("red","blue"),
       legend = "topright", main = "GAERS - NEC")
# Fold-change of 1.5
tab_mod <- tab$table %>%
  mutate("Sig DE" = case_when(logFC >= log2(1.5) & FDR <= 0.05 ~ "Up",
                              logFC <= -log2(1.5) & FDR <= 0.05 ~ "Down",
                              TRUE ~ "Not Sig"))
top_n <- 10
top_genes <- rbind(
  tab_mod %>%
    filter(`Sig DE` == "Up") %>%
    arrange(FDR, desc(abs(logFC))) %>%
    head(top_n),
  tab_mod %>%
    filter(`Sig DE` == "Down") %>%
    arrange(FDR, desc(abs(logFC))) %>%
    head(top_n)  
)
ggplot(tab_mod, aes(logFC, -log(PValue,10))) +
  geom_point(aes(color = `Sig DE`), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"(P-Value)")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(data = top_genes,
                  aes(logFC, -log(PValue,10), 
                      label = top_genes$gene_name),
                  size = 3,
                  max.overlaps = 50) +
  ggtitle("Transcriptomics DE Tha 16wk GAERS v NEC, FC = 1.5")
# Fold change of 2
tab_mod <- tab$table %>%
  mutate("Sig DE" = case_when(logFC >= log2(2) & FDR <= 0.05 ~ "Up",
                              logFC <= -log2(2) & FDR <= 0.05 ~ "Down",
                              TRUE ~ "Not Sig"))
top_n <- 10
top_genes <- rbind(
  tab_mod %>%
    filter(`Sig DE` == "Up") %>%
    arrange(FDR, desc(abs(logFC))) %>%
    head(top_n),
  tab_mod %>%
    filter(`Sig DE` == "Down") %>%
    arrange(FDR, desc(abs(logFC))) %>%
    head(top_n)  
)
# Plot the volcano plot
ggplot(tab_mod, aes(logFC, -log(PValue,10))) +
  geom_point(aes(color = `Sig DE`), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"(P-Value)")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(data = top_genes,
                  aes(logFC, -log(PValue,10), 
                      label = top_genes$gene_name),
                  size = 3,
                  max.overlaps = 50) +
  ggtitle("Transcriptomics DE Tha 16wk GAERS v NEC, FC = 2")
dev.off()