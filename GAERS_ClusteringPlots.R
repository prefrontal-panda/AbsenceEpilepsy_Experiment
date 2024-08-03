# This is to perform several clusting plots to identify sample outliers before statistical testing
# These are: PCA, dendrogram, and heatmap

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Loading libraries
library(ggfortify)
library(plotly)
library(ggplot2)
library(tidyverse)
library(ggrepel)

# Loading files
prot <- read.csv("Prot_tha_3wk_ED1_noNA_timpt.csv", stringsAsFactors = F, sep = ",", row.names = 1)
class <- c(rep("GAERS",12), rep("NEC",12)) # Adding species class
df_prot <- add_column(prot, strain = class, .before = 1)

# Metabolomics
met = read.csv("Met_SCx_7wk_timpt.csv", stringsAsFactors = F, sep = ",", row.names = 1)
class <- c(rep("GAERS",11), rep("NEC",12))
df_met <- add_column(met, strain = class, .before = 1)

# Transcriptomics
trans <- read.csv("PTE_Transcriptomics_PCA.csv", stringsAsFactors = F, row.names = 1, sep = ",")
trans <- as.data.frame(t(trans))
class <- c(rep("Naive",10), rep("PTE",10), rep("Sham",12), rep("TBI",13))
df_trans <- add_column(trans, strain = class, .before = 1)

#----PCA----
df <- df_prot[2:4292]
pca_res <- prcomp(df, scale. = F)
autoplot(pca_res, data = df_prot) + 
  geom_point(aes(color = strain)) + 
  stat_ellipse(aes(color = strain)) +
  geom_text_repel(label = rownames(df_prot), max.overlaps = 35) +
  ggtitle("Proteomics Tha 3 wks")
# Finding the total variance explained by each PC
var_exp <- pca_res$sdev^2/sum(pca_res$sdev^2)
# Creating Scree Plot
qplot(c(1:24), var_exp) + # Change the numbers in c() to the length of var_exp
  geom_line() +
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Metabolomics Scx 7 wks") +
  ylim(0, 0.25)
# Finding exact percentage total variance explained by each principal component
print(var_exp)

# Removing Outliers
df_prot1 <- df_prot[-c(13,20,24),]
df1 <- df_prot1[2:4292]
pca_res1 <- prcomp(df1, scale. = F)
autoplot(pca_res1, data = df_prot1) + 
  geom_point(aes(color = strain)) + 
  stat_ellipse(aes(color = strain)) +
  geom_text_repel(label = rownames(df_prot1), max.overlaps = 35) +
  ggtitle("Proteomics Tha 7 wks (No NEC 1,5,9)")
# Variance plot
var_exp1 <- pca_res1$sdev^2/sum(pca_res1$sdev^2)
qplot(c(1:22), var_exp1) + # Change the numbers in c() to the length of var_exp
  geom_line() +
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Metabolomics Scx 7 wks (No Outliers)") +
  ylim(0, 0.25)
print(var_exp)

#----Dendrogram----
# Function from ggplot. This also resuslts in the plotting of a heatmap
heatmap(as.matrix(df),labCol = F)

# Can also try hclust
dist <- dist(as.matrix(df), method = "euclidian", diag = T) # Calculating the euclidian distance between samples
hc <- hclust(dist) # Performing hierarchical clustering
plot(hc) # Plotting the dendrogram

#----MDS----
# Setting as matrix
data_raw<-as.matrix(prot)
# Clustering Analysis
distance <- dist(data_raw) 
mds1<-cmdscale(distance, k=2) 
plot(mds1, type='n')
text(mds1, labels=rownames(data_raw), cex=0.6, adj=0.5)

# Using these plots, we can select the best outliers to remove and save the file file.
write.csv(df1, "Prot_Tha_16wk_noOut_timpt.csv")