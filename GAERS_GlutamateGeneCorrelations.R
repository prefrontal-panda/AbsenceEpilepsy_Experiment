# This is to perform a correlation analysis between Nit2 and other genes involves in glutamate metabolism

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Loading libraries
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggdendroplot)

# Loading dataframes
dt <- read.csv("GAERS_Tha_Transcriptomics_MOFA.csv", stringsAsFactors = F, sep = ",", row.names = 1)
ids <- read.csv("Glutamate_Met_Genes.csv", stringsAsFactors = F,  sep = ",")

# Subsetting the main dataframe
glut_gene <- dt[,names(dt) %in% ids$Gene_Symbol]

# Separating out to GAERS and NECs
#gaers <- glut_gene[c(1:10),]
#nec <- glut_gene[c(11:22),]
glut_gene <- glut_gene[c(11:22),]

#----Normality Testing----
# Shapiro test
norm <- lapply(glut_gene, shapiro.test)
sink("NEC_Tha_ShapiroTest.txt")
print(norm)
sink()
# Q-Q plot
pdf("NEC_Tha_QQPlot_Glutamate Metabolism Genes.pdf")
for (gene in names(glut_gene)) {
  plt <- ggqqplot(glut_gene[[gene]], ylab = gene, 
                  title = paste("Q-Q Plot of", gene),
                  ggtheme = theme_classic())
  print(plt)
}
dev.off()

#----Correlation between genes----
# To iterate across dataframe
# Create an empty data frame to store the results
corr_res <- data.frame(
  Column1 = character(),  # Assuming the column names are strings
  Column2 = character(),
  Correlation = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Saving to pdf file
pdf("NEC_Tha_CorrelationResults.pdf")
# Loop through pairs of columns
for (col1 in colnames(glut_gene)) {
  for (col2 in colnames(glut_gene)) {
    if (col1 != col2) {
      # Perform correlation test
      cor_test_result <- cor.test(glut_gene[[col1]], glut_gene[[col2]], method = "spearman")
      
      # Store results in the data frame
      result <- data.frame(
        Gene_1 = col1,
        Gene_2 = col2,
        Correlation = cor_test_result$estimate,
        P_Value = cor_test_result$p.value
      )
      # Append the result to the data frame
      corr_res <- rbind(corr_res, result)
      
      # Plot using ggscatter
      plt <- ggscatter(glut_gene, x = col1, y = col2,
                       color = "black", shape = 19, size = 3,
                       add.params = list(color = "blue", fill = "lightgray"),
                       add = "reg.line",
                       conf.int = TRUE, 
                       cor.coef = TRUE, 
                       cor.coeff.args = list(method = "spearman", label.sep = "\n"), #, label.x.npc = 0.8, label.y.npc = 1),
                       cor.method = "spearman",
                       xlab = col1, ylab = col2,
                       title = paste0("Spearman Correlation between ", col1, " and ", col2, " in GAERS rats"),
                       ggtheme = theme_classic())
      print(plt)
    }
  }
}
dev.off()
write.csv(corr_res, file = "NEC_Tha_CorrelationResults.csv", row.names = F)

# Getting those with cut-offs p < 0.05
df <- read.csv("Glutamate Genes Correlation/NEC_Tha_CorrelationResults.csv",
               stringsAsFactors = F, header = T)

sig_corr <- df[df$P_Value < 0.05 & (df$Correlation > 0.7 | df$Correlation < -0.7),]
write.csv(sig_corr, "NEC_Tha_SigCorr.csv", row.names = F)

#------- Correlation Heatmap ----------
# Dataset needs to be in long format (i.e. colnames = gene_1, gene_2, correlation)

# Reading in file
data <- read.csv("GAERS_Tha_CorrelationResults.csv", stringsAsFactors = F, header = T)
data_nec <- read.csv("NEC_Tha_CorrelationResults.csv", stringsAsFactors = F, header = T)
# Heatmap
pdf("GAERS_Tha_Heatmap.pdf")
ggplot(data, aes(x = Gene_1, 
                 y = factor(Gene_2, levels = rev(levels(factor(Gene_2)))), 
                 fill = Correlation)) +
  geom_tile() + 
  scale_x_discrete(position = "top") +
  ylab("") + 
  xlab("") +
  ggtitle("Glutamate Gene Correlation in GAERS Tha") +
  theme(axis.text.x.top = element_text(angle = 90))#, margin = margin(b = 10, r = 20, unit = "pt")))
ggplot(nit2_gaers, aes(x = Gene_1,
                       y = factor(Gene_2, levels = rev(levels(factor(Gene_2)))),
                       fill = Correlation)) +
  geom_tile() + 
  ylab("") + 
  xlab("") +
  ggtitle("Glutamate Gene Correlation with NIT2 in GAERS Tha") 
dev.off()

#----Jaccard Similarity Analysis----
# Sub-setting just Nit2
nit2_gaers <- data %>%
  subset(Gene_1 == "Nit2") %>%
  arrange(Gene_2)
nit2_nec <- data_nec %>%
  subset(Gene_1 == "Nit2") %>%
  arrange(Gene_2)
# Rounding to one d.p.
nit2_gaers_round <- nit2_gaers
nit2_gaers_round$Correlation <- round(nit2_gaers_round$Correlation, 2)

nit2_nec_round <- nit2_nec
nit2_nec_round$Correlation <- round(nit2_nec_round$Correlation, 2)

# Function
jaccard <- function(setA, setB) {
  intersection = length(intersect(setA, setB))
  union = length(setA) + length(setB) - intersection
  return (intersection/union)
}
# Running
jaccard(nit2_gaers_round$Correlation, nit2_nec_round$Correlation)
# Scx: 0.07936508
# Thalamus: 0.07936508

# Running Jaccard on every gene in the dataframe
# Rounding correlation values 
data$Correlation <- round(data$Correlation, 2)
data_nec$Correlation <- round(data_nec$Correlation, 2)
# Splitting into individual gene sets
list_df <- split(data, data$Gene_1)
list_nec <- split(data_nec, data_nec$Gene_1)
# Running jaccard for the individual sets
library(foreach)
jac_lst <- foreach(i = list_df, j = list_nec) %do% {
  #print(i)
  #print(j)
  jaccard(i$Correlation, j$Correlation)
}
# Changing the names
names(jac_lst) <- names(list_df)
# Saving
sink("Tha_AllGlutGenes_Jaccard.txt")
print(jac_lst)
sink()

# ---- Heatmap + dendrogram -----
# Reading in file
data <- read.csv("GAERS_Scx_CorrelationResults.csv", stringsAsFactors = F, header = T)
data_nec <- read.csv("NEC_Scx_CorrelationResults.csv", stringsAsFactors = F, header = T)
# Sub-setting just Nit2
colnames(data) <- colnames(data_nec)

nit2_gaers <- data %>%
  subset(Gene_1 == "Nit2") %>%
  arrange(Gene_2)
nit2_nec <- data_nec %>%
  subset(Gene_1 == "Nit2") %>%
  arrange(Gene_2)
# Adding columns
nit2_gaers <- nit2_gaers %>% 
  mutate(pair = paste(Gene_1, Gene_2, sep = "-"),
         strain = "GAERS",
         .before = 1)
nit2_nec <- nit2_nec %>% 
  mutate(pair = paste(Gene_1, Gene_2, sep = "-"),
         strain = "NEC",
         .before = 1)

# Making dataframe
nit2_all <- data.frame(GAERS = nit2_gaers$Correlation,
                       NEC = nit2_nec$Correlation)
row.names(nit2_all) <- nit2_gaers$pair

# perform hierarchical clustering
rowclus <- hclust(dist(nit2_all))    #cluster the rows
colclus <- hclust(dist(t(nit2_all))) #cluster the columns

# getting data ready to be used by ggplot
hm <- hmReady(nit2_all, colclus=colclus, rowclus=rowclus)

# plot the heatmap
plt_htmp <- ggplot() + 
  geom_tile(data=hm, aes(x=variable, y=y, fill=value)) + #heatmap
  scale_fill_gradientn(colors=hmGradient()) + #options for heatmap
  #geom_dendro(colclus, ylim=c(16.5, 20)) + #upper dendrogram
  geom_dendro(rowclus, xlim=c(2.5,3.5), pointing="side",  #side dendrogram
              dendrocut = 30) +             
  ylab("") + 
  xlab("") +
  theme_classic() +  #design
  ggtitle("Nit2 Pairwise Correlation in the Somatosensory Cortex")
# Save
ggsave(filename = "Nit2_Pairwise_Heatmap_Scx.tiff", plot = plt_htmp, 
       width = 24, height = 24, units = "cm", device='tiff', dpi=300)
